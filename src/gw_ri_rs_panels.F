!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Grid-panel streaming for the GW RI-RS implementations.
!>
!>        The real-space grid index l is a contraction index of every RI-RS quantity, so it can be
!>        processed in batches ("panels") of contiguous grid row-blocks without changing the
!>        result. Streaming over panels is what keeps the grid x grid intermediates
!>        (phi G phi^T, Z W Z^T, their Hadamard products) from ever being materialized for the
!>        whole system: each panel is built, consumed and discarded.
!>
!>        On top of that, the optional cutoff radii restrict a panel to the grid blocks that are
!>        geometrically near it. The planner sizes the panels against the measured free memory and
!>        the 32-bit DBCSR message-length limit, and the resulting partition is memoized in bs_env
!>        so that chi, Sigma^x and Sigma^c all use the identical one.
!>
!>        Distances are measured through bs_env%ri_rs%cell, so the same routines serve the
!>        molecular and the periodic Gamma-point driver: for a molecule cell%perd = 0 and pbc() is
!>        the identity, which reduces every test below to a plain Euclidean distance.
! **************************************************************************************************
MODULE gw_ri_rs_panels
   USE cp_dbcsr_api, ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_info, dbcsr_get_occupation, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE gw_tensor_large_cell_Gamma,      ONLY: G_occ_vir
   USE gw_utils_dbcsr,                  ONLY: hadamard_product,&
                                              release_dbcsr_topology_and_matrices,&
                                              setup_square_topology
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              max_line_length
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_panels'

   PUBLIC :: build_G_ao, build_geo_template_panel, collect_used_col_blocks, &
             contract_grid_panels, contract_grid_panels_sigma_c, create_product_matrix, &
             extract_grid_panel, extract_masked_blocks, mask_grid_blocks_near_panel, &
             panel_mem_estimate_GB, panel_template_elems, plan_grid_panels, &
             print_ri_rs_memory_estimate, reserve_blocks_within_radius, resolve_grid_panels

CONTAINS

! **************************************************************************************************
!> \brief Marks the grid blocks whose centroid lies within cutoff of the bounding box of the
!>        panel [blk0, blk1]'s chunk centroids.
!> \param centroids ...
!> \param blk0 ...
!> \param blk1 ...
!> \param cutoff ...
!> \param used ...
! **************************************************************************************************
   SUBROUTINE mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, used)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      INTEGER, INTENT(IN)                                :: blk0, blk1
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mask_grid_blocks_near_panel'

      INTEGER                                            :: c, handle, k
      REAL(KIND=dp)                                      :: cutoff2, d2, dx
      REAL(KIND=dp), DIMENSION(3)                        :: hi, lo

      CALL timeset(routineN, handle)

      cutoff2 = cutoff**2
      lo(:) = MINVAL(centroids(:, blk0:blk1), DIM=2)
      hi(:) = MAXVAL(centroids(:, blk0:blk1), DIM=2)

      ALLOCATE (used(SIZE(centroids, 2)))
      DO c = 1, SIZE(centroids, 2)
         d2 = 0.0_dp
         DO k = 1, 3
            dx = MAX(0.0_dp, lo(k) - centroids(k, c), centroids(k, c) - hi(k))
            d2 = d2 + dx*dx
         END DO
         used(c) = (d2 <= cutoff2)
      END DO

      CALL timestop(handle)

   END SUBROUTINE mask_grid_blocks_near_panel

! **************************************************************************************************
!> \brief Exact allocated-element count of the geo template of panel [blk0, blk1]: the very same
!>        per-block-pair centroid test as build_geo_template_panel, so this is the true DBCSR
!>        data size of A_pan/B_pan/C_pan (DBCSR stores whole blocks).
!> \param r_blk_sizes ...
!> \param centroids ...
!> \param used ...
!> \param blk0 ...
!> \param blk1 ...
!> \param cutoff ...
!> \param nze_tmpl ...
! **************************************************************************************************
   SUBROUTINE panel_template_elems(r_blk_sizes, centroids, used, blk0, blk1, cutoff, nze_tmpl)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: used
      INTEGER, INTENT(IN)                                :: blk0, blk1
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER(KIND=int_8), INTENT(OUT)                   :: nze_tmpl

      CHARACTER(LEN=*), PARAMETER :: routineN = 'panel_template_elems'

      INTEGER                                            :: c, handle, ib, n_used
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: used_idx
      REAL(KIND=dp)                                      :: cutoff2

      CALL timeset(routineN, handle)

      ! Compress the near mask once so the pair loop only visits candidate columns.
      n_used = COUNT(used)
      ALLOCATE (used_idx(n_used))
      n_used = 0
      DO c = 1, SIZE(used)
         IF (used(c)) THEN
            n_used = n_used + 1
            used_idx(n_used) = c
         END IF
      END DO

      cutoff2 = cutoff**2
      nze_tmpl = 0_int_8
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(blk0, blk1, n_used, used_idx, centroids, cutoff2, &
      !$OMP r_blk_sizes) PRIVATE(ib, c) REDUCTION(+:nze_tmpl)
      DO ib = blk0, blk1
         DO c = 1, n_used
            IF (SUM((centroids(:, ib) - centroids(:, used_idx(c)))**2) <= cutoff2) THEN
               nze_tmpl = nze_tmpl + INT(r_blk_sizes(ib), int_8)*INT(r_blk_sizes(used_idx(c)), int_8)
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE panel_template_elems

! **************************************************************************************************
!> \brief Per-rank peak memory (GB) of one panel step of the neighborhood-restricted
!>        contractions: three grid x grid panels of the template size (A_pan, B_pan, C_pan)
!>        plus the grid x RI intermediates (tmp2 and the accumulation operand) and the
!>        grid x AO intermediate (tmpA), whose column support is the panel's geometric
!>        neighborhood fraction f_near = width/n_grid. Shared by the panel planner and
!> \param nze_tmpl ...
!> \param pan_rows ...
!> \param width ...
!> \param n_grid_total ...
!> \param n_RI ...
!> \param n_ao ...
!> \param n_procs ...
!> \param mem_GB ...
! **************************************************************************************************
   SUBROUTINE panel_mem_estimate_GB(nze_tmpl, pan_rows, width, n_grid_total, n_RI, n_ao, &
                                    n_procs, mem_GB)

      INTEGER(KIND=int_8), INTENT(IN)                    :: nze_tmpl
      INTEGER, INTENT(IN)                                :: pan_rows, width, n_grid_total, n_RI, &
                                                            n_ao, n_procs
      REAL(KIND=dp), INTENT(OUT)                         :: mem_GB

      REAL(KIND=dp)                                      :: f_near

      f_near = REAL(width, dp)/REAL(MAX(n_grid_total, 1), dp)
      mem_GB = (3.0_dp*REAL(nze_tmpl, dp) + &
                REAL(pan_rows, dp)*f_near*(2.0_dp*REAL(n_RI, dp) + REAL(n_ao, dp)))* &
               8.0_dp/REAL(MAX(n_procs, 1), dp)*1.0E-9_dp

   END SUBROUTINE panel_mem_estimate_GB

! **************************************************************************************************
!> \brief Plans the panel boundaries for the streaming contractions. Panels grow by whole grid
!>        row-blocks towards ~panel_size rows. When the neighborhood restriction is active
!>        (centroids+cutoff), each candidate panel is additionally checked against
!>          (a) the 32-bit message bound with the panel's TRUE occupancy
!>          (b) the per-rank memory budget: panel_mem_estimate_GB <= mem_budget_GB.
!> \param bs_env ...
!> \param r_blk_sizes ...
!> \param panel_size ...
!> \param min_dim ...
!> \param pan_first ...
!> \param pan_last ...
!> \param centroids ...
!> \param cutoff ...
!> \param n_RI ...
!> \param n_ao ...
!> \param n_procs ...
!> \param mem_budget_GB ...
! **************************************************************************************************
   SUBROUTINE plan_grid_panels(bs_env, r_blk_sizes, panel_size, min_dim, pan_first, pan_last, &
                               centroids, cutoff, n_RI, n_ao, n_procs, mem_budget_GB)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes
      INTEGER, INTENT(IN)                                :: panel_size, min_dim
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pan_first, pan_last
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff
      INTEGER, INTENT(IN), OPTIONAL                      :: n_RI, n_ao, n_procs
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: mem_budget_GB

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'plan_grid_panels'

      INTEGER                                            :: blk0, blk1, handle, ib, n_grid_blocks, &
                                                            n_grid_total, n_panels, rows_acc, &
                                                            TARGET, width
      INTEGER(KIND=int_8)                                :: msg, nze_tmpl, side
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: tmp_first, tmp_last
      LOGICAL                                            :: fits, use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: used
      REAL(KIND=dp)                                      :: f_near, mem_GB

      CALL timeset(routineN, handle)

      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp
      IF (use_cutoff) THEN
         CPASSERT(PRESENT(n_RI) .AND. PRESENT(n_ao) .AND. PRESENT(n_procs))
      END IF

      n_grid_blocks = SIZE(r_blk_sizes)
      n_grid_total = SUM(r_blk_sizes)
      ALLOCATE (tmp_first(n_grid_blocks), tmp_last(n_grid_blocks))

      n_panels = 0
      blk0 = 1
      DO WHILE (blk0 <= n_grid_blocks)
         TARGET = panel_size
         DO
            rows_acc = 0
            blk1 = blk0
            DO ib = blk0, n_grid_blocks
               rows_acc = rows_acc + r_blk_sizes(ib)
               blk1 = ib
               IF (rows_acc >= TARGET) EXIT
            END DO
            IF (.NOT. use_cutoff .OR. blk1 == blk0) EXIT
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, used)
            width = SUM(r_blk_sizes, MASK=used)
            CALL panel_template_elems(r_blk_sizes, centroids, used, blk0, blk1, cutoff, nze_tmpl)
            f_near = REAL(width, dp)/REAL(MAX(n_grid_total, 1), dp)
            side = INT(REAL(rows_acc, dp)*f_near*REAL(MAX(n_RI, n_ao), dp), int_8)
            msg = MAX(nze_tmpl, side)/INT(MAX(min_dim, 1), int_8)
            fits = (msg <= bs_env%dbcsr_msg_elem_limit/4)
            IF (fits .AND. PRESENT(mem_budget_GB)) THEN
               IF (mem_budget_GB > 0.0_dp) THEN
                  CALL panel_mem_estimate_GB(nze_tmpl, rows_acc, width, n_grid_total, &
                                             n_RI, n_ao, n_procs, mem_GB)
                  fits = (mem_GB <= mem_budget_GB)
               END IF
            END IF
            ! Panel size is bounded only by the message-overflow and memory checks above; there is
            ! no neighborhood-width (f_near) cap. mp_waitall is dominated by the NUMBER of panel
            ! multiplies, so fewer/larger panels are cheaper here.
            IF (fits) EXIT
            TARGET = MAX(1, MIN(TARGET, rows_acc)/2)
         END DO
         n_panels = n_panels + 1
         tmp_first(n_panels) = blk0
         tmp_last(n_panels) = blk1
         blk0 = blk1 + 1
      END DO

      ALLOCATE (pan_first(n_panels), pan_last(n_panels))
      pan_first(:) = tmp_first(1:n_panels)
      pan_last(:) = tmp_last(1:n_panels)
      DEALLOCATE (tmp_first, tmp_last)

      CALL timestop(handle)

   END SUBROUTINE plan_grid_panels

! **************************************************************************************************
!> \brief Resolves the panel boundaries for the streaming contractions from the bs_env settings:
!> \param bs_env ...
!> \param mat_phi_mu_l ...
!> \param pan_first ...
!> \param pan_last ...
! **************************************************************************************************
   SUBROUTINE resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: pan_first, pan_last

      CHARACTER(LEN=*), PARAMETER :: routineN = 'resolve_grid_panels'

      INTEGER                                            :: handle, min_dim, n_grid_total, npcols, &
                                                            nprows, panel_size, safe_max
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      LOGICAL                                            :: use_cutoff
      REAL(KIND=dp)                                      :: mem_avail_GB, mem_budget_GB
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      IF (ALLOCATED(bs_env%ri_rs%pan_first)) THEN
         ALLOCATE (pan_first, SOURCE=bs_env%ri_rs%pan_first)
         ALLOCATE (pan_last, SOURCE=bs_env%ri_rs%pan_last)
         CALL timestop(handle)
         RETURN
      END IF

      use_cutoff = bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .AND. &
                   ALLOCATED(bs_env%ri_rs%chunk_centroids)

      ! MIN(nprows, npcols) is the divisor that bounds the worst-rank Cannon message: a
      ! P x n_grid panel is replicated into block row strips (P/nprows x n_grid) or column
      ! strips (P x n_grid/npcols) during multiply_cannon, so the largest single-rank
      ! message is ~ P*n_grid / MIN(nprows,npcols) elements.
      CALL dbcsr_get_info(mat_phi_mu_l, nfullrows_total=n_grid_total, row_blk_size=r_blk_sizes, &
                          distribution=dist)
      CALL dbcsr_distribution_get(dist, nprows=nprows, npcols=npcols)
      min_dim = MAX(MIN(nprows, npcols), 1)

      ! Panel height such that NO per-rank DBCSR message can overflow the 32-bit length field
      ! (see bs_env%dbcsr_msg_elem_limit): requiring the worst-rank message to stay under
      ! 0.5 * HUGE(int_4) gives the safe height P_safe = 0.5 * HUGE(int_4) * min_dim / n_grid.
      IF (use_cutoff) THEN
         safe_max = n_grid_total
      ELSE
         safe_max = INT(0.5_dp*REAL(bs_env%dbcsr_msg_elem_limit, dp)*REAL(min_dim, dp)/ &
                        REAL(n_grid_total, dp))
         safe_max = MAX(1, MIN(safe_max, n_grid_total))
      END IF

      ! Start from a single whole-grid panel, clamped to the overflow-safe ceiling; with the
      ! cutoff active the planner splits it further for the message/memory limits.
      panel_size = MAX(1, safe_max)

      IF (use_cutoff) THEN
         ! Half of the measured free memory as panel budget.
         CALL mp_mem_avail_per_rank_GB(bs_env%para_env, mem_avail_GB)
         mem_budget_GB = 0.5_dp*mem_avail_GB
         CALL plan_grid_panels(bs_env, r_blk_sizes, panel_size, min_dim, pan_first, pan_last, &
                               centroids=bs_env%ri_rs%chunk_centroids, &
                               cutoff=bs_env%ri_rs%cutoff_radius_v_w, &
                               n_RI=bs_env%n_RI, n_ao=bs_env%n_ao, &
                               n_procs=bs_env%para_env%num_pe, mem_budget_GB=mem_budget_GB)
      ELSE
         CALL plan_grid_panels(bs_env, r_blk_sizes, panel_size, min_dim, pan_first, pan_last)
      END IF

      ALLOCATE (bs_env%ri_rs%pan_first, SOURCE=pan_first)
      ALLOCATE (bs_env%ri_rs%pan_last, SOURCE=pan_last)

      CALL timestop(handle)

   END SUBROUTINE resolve_grid_panels

! **************************************************************************************************
!> \brief Estimates and prints per-process memory requirements for the RI-RS GW calculation.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE print_ri_rs_memory_estimate(bs_env)

!$    USE OMP_LIB, ONLY: omp_get_max_threads

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_ri_rs_memory_estimate'

      CHARACTER(LEN=max_line_length)                     :: msg
      INTEGER                                            :: handle, iatom, ipan, l, &
                                                            max_n_ao_used, max_n_local_grid, &
                                                            n_ao_used_atom, n_grid_total, &
                                                            n_local_grid, n_loc_ri_max, n_procs, &
                                                            n_RI, n_threads, natom, pan_rows, &
                                                            pan_width
      INTEGER(KIND=int_8)                                :: nze_tmpl
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      LOGICAL                                            :: use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used
      REAL(KIND=dp)                                      :: cutoff_ri, mem_avail_GB, mem_D_local_GB, &
                                                            mem_dlp_GB, mem_pan_GB, mem_panels_GB, &
                                                            mem_phi_local_GB, mem_Z_lP_GB, &
                                                            mem_Zlp_peak_GB, pos_P(3)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      n_grid_total = bs_env%ri_rs%n_grid_points
      CALL dbcsr_get_info(bs_env%ri_rs%mat_phi_mu_l, row_blk_size=r_blk_sizes)
      CPASSERT(SUM(r_blk_sizes) == n_grid_total)
      n_RI = bs_env%n_RI
      n_procs = bs_env%para_env%num_pe

      ! Z_lP upper bound: dense n_grid × n_RI, distributed evenly across all ranks.
      ! The actual sparse Z_lP is smaller due to the per-atom locality cutoff.
      mem_Z_lP_GB = REAL(n_grid_total, dp)*REAL(n_RI, dp)*8.0_dp/ &
                    REAL(n_procs, dp)*1.0E-9_dp

      ! Peak panel memory during Σ^c: two G panels (A_occ, A_vir) + one W panel plus the
      ! grid × RI / grid × AO intermediates. With the CUTOFF_RADIUS_RL_W restriction the panel
      ! matrices only allocate the geo-template blocks, so use the same nze-aware model as the
      ! panel planner (panel_mem_estimate_GB); without the cutoff, dense panel_rows × n_grid.
      ! Plus the n_RI × n_RI W_aux matrix. Distributed over n_procs ranks.
      use_cutoff = bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .AND. &
                   ALLOCATED(bs_env%ri_rs%chunk_centroids)
      CALL resolve_grid_panels(bs_env, bs_env%ri_rs%mat_phi_mu_l, pan_first, pan_last)
      mem_panels_GB = 0.0_dp
      DO ipan = 1, SIZE(pan_first)
         pan_rows = SUM(r_blk_sizes(pan_first(ipan):pan_last(ipan)))
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(bs_env%ri_rs%chunk_centroids, pan_first(ipan), &
                                             pan_last(ipan), bs_env%ri_rs%cutoff_radius_v_w, &
                                             grid_used)
            pan_width = SUM(r_blk_sizes, MASK=grid_used)
            CALL panel_template_elems(r_blk_sizes, bs_env%ri_rs%chunk_centroids, &
                                      grid_used, pan_first(ipan), pan_last(ipan), &
                                      bs_env%ri_rs%cutoff_radius_v_w, nze_tmpl)
            CALL panel_mem_estimate_GB(nze_tmpl, pan_rows, pan_width, n_grid_total, &
                                       n_RI, bs_env%n_ao, n_procs, mem_pan_GB)
         ELSE
            pan_width = n_grid_total
            mem_pan_GB = (3.0_dp*REAL(pan_rows, dp)*REAL(pan_width, dp) + &
                          2.0_dp*REAL(pan_rows, dp)*REAL(n_RI, dp))* &
                         8.0_dp/REAL(n_procs, dp)*1.0E-9_dp
         END IF
         mem_panels_GB = MAX(mem_panels_GB, mem_pan_GB)
      END DO
      mem_panels_GB = mem_panels_GB + &
                      REAL(n_RI, dp)*REAL(n_RI, dp)*8.0_dp/REAL(n_procs, dp)*1.0E-9_dp

      ! Z_lP SOLVE peak (compute_Z_lP). For the atom P with the largest integration
      ! sphere, one rank holds simultaneously:
      !   D_local  : n_local_grid x n_local_grid  (dense D'_ll', BLAS path only; O(n_local_grid^2))
      !   phi_local: n_local_grid x n_ao_used  (AOs reaching into the sphere only)
      !   d_lp     : n_local_grid x n_loc_ri, replicated once + one private copy per OMP thread
      ! n_local_grid = # grid points within cutoff_ri(P) = r_c(RI metric) + r_RI(P).
      ! This is NOT evenly distributed: n_local_grid depends on the
      ! local density of atoms/grid, so the rank owning the densest atom peaks well above the
      ! average. We report the worst-case (max over atoms) as a per-rank upper bound.
      particle_set => bs_env%ri_rs%particle_set
      CPASSERT(ASSOCIATED(particle_set))
      natom = bs_env%n_atom

      max_n_local_grid = 0
      n_loc_ri_max = 0
      max_n_ao_used = 0
      DO iatom = 1, natom
         cutoff_ri = bs_env%ri_metric%cutoff_radius + bs_env%ri_rs%radius_ri_per_atom(iatom)
         pos_P(:) = particle_set(iatom)%r(:)
         n_local_grid = 0
         DO l = 1, n_grid_total
            IF (SUM((bs_env%ri_rs%grid_points(1:3, l) - pos_P(1:3))**2) <= cutoff_ri**2) THEN
               n_local_grid = n_local_grid + 1
            END IF
         END DO
         max_n_local_grid = MAX(max_n_local_grid, n_local_grid)
         CALL get_n_ao_in_sphere(bs_env, iatom, cutoff_ri, n_ao_used_atom)
         max_n_ao_used = MAX(max_n_ao_used, n_ao_used_atom)
         n_loc_ri_max = MAX(n_loc_ri_max, &
                            bs_env%i_RI_end_from_atom(iatom) - bs_env%i_RI_start_from_atom(iatom) + 1)
      END DO

      n_threads = 1
!$    n_threads = omp_get_max_threads()

      ! D_local: dense on one rank (single-rank BLAS path). This is the upper bound; atoms whose
      ! solve does not fit are distributed automatically over a rank subgroup (ScaLAPACK).
      mem_D_local_GB = REAL(max_n_local_grid, dp)**2*8.0_dp*1.0E-9_dp
      mem_phi_local_GB = REAL(max_n_local_grid, dp)*REAL(max_n_ao_used, dp)*8.0_dp*1.0E-9_dp
      mem_dlp_GB = REAL(max_n_local_grid, dp)*REAL(n_loc_ri_max, dp)*8.0_dp* &
                   REAL(1 + n_threads, dp)*1.0E-9_dp
      mem_Zlp_peak_GB = mem_D_local_GB + mem_phi_local_GB + mem_dlp_GB

      ! Available memory per process = node MemLikelyFree / ranks-per-node, min across ranks
      ! (0 on non-Linux => warnings suppressed below).
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, mem_avail_GB)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ' '
         WRITE (bs_env%unit_nr, '(T2,A)') 'RI-RS memory estimate per MPI process:'
         WRITE (bs_env%unit_nr, '(T4,A,F37.2,A)') &
            'Available memory per process (system)', mem_avail_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F18.2,A)') &
            'Required for Z_lP (dense upper bound; actual is sparser)', mem_Z_lP_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F25.2,A)') &
            'Required for χ, W, Σ panels (peak per panel step)', mem_panels_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,F17.2,A)') &
            'Required for Z_lP solve peak (D_local+ϕ, worst-case atom)', mem_Zlp_peak_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T69,I12)') &
            'Worst-case local-grid number of grid points:', max_n_local_grid
         WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
            'Worst-case memory D_local:', mem_D_local_GB, ' GB'

      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_Z_lP_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated memory for Z_lP, ", mem_Z_lP_GB, " GB per process, exceeds the "// &
            "available ", mem_avail_GB, " GB. Z_lP (n_grid x n_RI) is distributed across all "// &
            "MPI ranks, so add nodes or use fewer MPI ranks per node."
         CPWARN(TRIM(msg))
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_panels_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated peak memory of the chi/W/Sigma panels, ", mem_panels_GB, &
            " GB per process, exceeds the available ", mem_avail_GB, " GB. Panel memory "// &
            "scales roughly as 3*panel_size*n_grid/n_procs, so add nodes, use fewer MPI ranks "// &
            "per node."
         CPWARN(TRIM(msg))
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. mem_Zlp_peak_GB > mem_avail_GB) THEN
         WRITE (msg, '(A,F0.2,A,F0.2,A)') &
            "The estimated peak memory of the Z_lP solve, ", mem_Zlp_peak_GB, &
            " GB per process, exceeds the available ", mem_avail_GB, &
            " GB. The per-atom matrix D'_ll' dominates and "// &
            "scales as n_local_grid^2, and it is not balanced across ranks: the rank owning "// &
            "the atom with the largest integration sphere peaks well above the average. "// &
            "Such atoms are distributed automatically over a rank subgroup via ScaLAPACK. "// &
            "If the run still runs out of memory, use fewer MPI ranks per node so that each "// &
            "rank has more memory for the peak atom."
         CPWARN(TRIM(msg))
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_ri_rs_memory_estimate

! **************************************************************************************************
!> \brief Creates an empty (panel_chunks x neighborhood_chunks) DBCSR matrix with zero blocks
!>        pre-allocated only where |centroid(panel_row r) - centroid(column c)| <= cutoff.
!>        Used with retain_sparsity=.TRUE. in the subsequent dbcsr_multiply so distant blocks
!>        of the grid-basis panels (ϕ G ϕ^T, Z W Z^T, ...) are never computed at all.
!> \param L_pan ...
!> \param L_full ...
!> \param centroids ...
!> \param cutoff ...
!> \param blk0 ...
!> \param A_template ...
!> \param col_map ...
! **************************************************************************************************
   SUBROUTINE build_geo_template_panel(L_pan, L_full, centroids, cutoff, blk0, A_template, col_map)
      TYPE(dbcsr_type), INTENT(IN)                       :: L_pan, L_full
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centroids
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      INTEGER, INTENT(IN)                                :: blk0
      TYPE(dbcsr_type), INTENT(OUT)                      :: A_template
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: col_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_geo_template_panel'

      INTEGER                                            :: c, cg, cs, handle, my_pcol, my_prow, &
                                                            n_grid_blks, n_pan_blks, npcols, &
                                                            nprows, r, rs
      INTEGER, DIMENSION(:), POINTER                     :: grid_blk_sizes, pan_blk_sizes
      REAL(KIND=dp)                                      :: cutoff2
      REAL(KIND=dp), ALLOCATABLE                         :: zero_blk(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      cutoff2 = cutoff**2
      CALL dbcsr_get_info(L_pan, nblkrows_total=n_pan_blks, row_blk_size=pan_blk_sizes)
      CALL dbcsr_get_info(L_full, nblkrows_total=n_grid_blks, row_blk_size=grid_blk_sizes)

      ! create_product_matrix assigns row r to process MOD(r-1,nprows) and
      ! col c to MOD(c-1,npcols), so we can determine local ownership analytically.
      CALL create_product_matrix(L_pan, L_full, 'N', 'T', A_template)
      CALL dbcsr_get_info(A_template, distribution=dist)
      CALL dbcsr_distribution_get(dist, nprows=nprows, npcols=npcols, &
                                  myprow=my_prow, mypcol=my_pcol)

      ALLOCATE (zero_blk(MAXVAL(pan_blk_sizes(1:n_pan_blks)), &
                         MAXVAL(grid_blk_sizes(1:n_grid_blks))))
      zero_blk(:, :) = 0.0_dp

      DO r = 1, n_pan_blks
         IF (MOD(r - 1, nprows) /= my_prow) CYCLE
         rs = pan_blk_sizes(r)
         DO c = 1, n_grid_blks
            IF (MOD(c - 1, npcols) /= my_pcol) CYCLE
            cg = c
            IF (PRESENT(col_map)) cg = col_map(c)
            IF ((centroids(1, blk0 + r - 1) - centroids(1, cg))**2 + &
                (centroids(2, blk0 + r - 1) - centroids(2, cg))**2 + &
                (centroids(3, blk0 + r - 1) - centroids(3, cg))**2 <= cutoff2) THEN
               cs = grid_blk_sizes(c)
               CALL dbcsr_put_block(A_template, r, c, zero_blk(1:rs, 1:cs))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(A_template)

      DEALLOCATE (zero_blk)
      CALL timestop(handle)

   END SUBROUTINE build_geo_template_panel

! **************************************************************************************************
!> \brief Slices a contiguous range of grid row-blocks [blk0, blk1] out of a (grid x n) DBCSR
!>        matrix into a new (P x n) panel matrix: iterate the source's local blocks, put the
!>        in-range ones into the panel with a remapped row-block index, then finalize. Row-block
!>        index i of the panel corresponds to source row-block blk0+i-1.
!> \param mat_full ...
!> \param blk0 ...
!> \param blk1 ...
!> \param mat_panel ...
! **************************************************************************************************
   SUBROUTINE extract_grid_panel(mat_full, blk0, blk1, mat_panel)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_full
      INTEGER, INTENT(IN)                                :: blk0, blk1
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_panel

      CHARACTER(LEN=*), PARAMETER :: routineN = 'extract_grid_panel'

      INTEGER                                            :: handle, ib, jb, npb
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_full, col_dist_full, &
                                                            row_blk_full, row_blk_pan, &
                                                            row_dist_full, row_dist_pan
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_pan
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_full, distribution=dist_full, &
                          row_blk_size=row_blk_full, col_blk_size=col_blk_full)
      CALL dbcsr_distribution_get(dist_full, row_dist=row_dist_full, col_dist=col_dist_full)

      npb = blk1 - blk0 + 1
      ALLOCATE (row_dist_pan(npb), row_blk_pan(npb))
      row_dist_pan(:) = row_dist_full(blk0:blk1)
      row_blk_pan(:) = row_blk_full(blk0:blk1)

      CALL dbcsr_distribution_new(dist_pan, template=dist_full, &
                                  row_dist=row_dist_pan, col_dist=col_dist_full)
      CALL dbcsr_create(mat_panel, name="grid_panel", dist=dist_pan, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_pan, col_blk_size=col_blk_full)

      CALL dbcsr_iterator_start(iter, mat_full)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         IF (ib < blk0 .OR. ib > blk1) CYCLE
         CALL dbcsr_put_block(mat_panel, ib - blk0 + 1, jb, blk)
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(mat_panel)

      CALL dbcsr_distribution_release(dist_pan)
      DEALLOCATE (row_dist_pan, row_blk_pan)

      CALL timestop(handle)

   END SUBROUTINE extract_grid_panel

! **************************************************************************************************
!> \brief Marks which column blocks of a DBCSR matrix carry at least one non-zero block anywhere
!>        (global union). Used to restrict the inner index of the panel multiplies to the
!>        AO/RI atoms that actually touch the panel (exact: dropped rows only meet zeros).
!> \param matrix ...
!> \param para_env ...
!> \param used ...
! **************************************************************************************************
   SUBROUTINE collect_used_col_blocks(matrix, para_env, used)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_used_col_blocks'

      INTEGER                                            :: handle, ib, jb, nblkcols
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iused
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nblkcols_total=nblkcols)
      ALLOCATE (iused(nblkcols))
      iused(:) = 0

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         iused(jb) = 1
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL para_env%sum(iused)

      ALLOCATE (used(nblkcols))
      used(:) = (iused(:) > 0)
      DEALLOCATE (iused)

      CALL timestop(handle)

   END SUBROUTINE collect_used_col_blocks

! **************************************************************************************************
!> \brief Copies the flagged block rows (compress_rows=.TRUE.) or block columns (.FALSE.) of a
!>        DBCSR matrix into a compressed matrix. The subset keeps the parent's process assignment
!>        along the compressed dimension, so every block stays on its owning rank: the extraction
!>        is purely local (zero communication), like extract_grid_panel.
!> \param mat_full ...
!> \param used ...
!> \param mat_out ...
!> \param compress_rows ...
!> \param blk_map ...
! **************************************************************************************************
   SUBROUTINE extract_masked_blocks(mat_full, used, mat_out, compress_rows, blk_map)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_full
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: used
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_out
      LOGICAL, INTENT(IN)                                :: compress_rows
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: blk_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'extract_masked_blocks'

      INTEGER                                            :: handle, ib, jb, n_blk, n_sub, r
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inv_map
      INTEGER, DIMENSION(:), POINTER                     :: blk_full, blk_sub, col_blk_full, &
                                                            col_dist_full, dist_full_1d, &
                                                            dist_sub_1d, row_blk_full, &
                                                            row_dist_full
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_distribution_type)                      :: dist_full, dist_sub
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_full, distribution=dist_full, &
                          row_blk_size=row_blk_full, col_blk_size=col_blk_full)
      CALL dbcsr_distribution_get(dist_full, row_dist=row_dist_full, col_dist=col_dist_full)

      IF (compress_rows) THEN
         blk_full => row_blk_full
         dist_full_1d => row_dist_full
      ELSE
         blk_full => col_blk_full
         dist_full_1d => col_dist_full
      END IF
      n_blk = SIZE(blk_full)

      n_sub = COUNT(used)
      CPASSERT(n_sub > 0)
      ALLOCATE (inv_map(n_blk), blk_sub(n_sub), dist_sub_1d(n_sub))
      IF (PRESENT(blk_map)) ALLOCATE (blk_map(n_sub))
      inv_map(:) = 0
      r = 0
      DO ib = 1, n_blk
         IF (used(ib)) THEN
            r = r + 1
            inv_map(ib) = r
            blk_sub(r) = blk_full(ib)
            dist_sub_1d(r) = dist_full_1d(ib)
            IF (PRESENT(blk_map)) blk_map(r) = ib
         END IF
      END DO

      IF (compress_rows) THEN
         CALL dbcsr_distribution_new(dist_sub, template=dist_full, &
                                     row_dist=dist_sub_1d, col_dist=col_dist_full)
         CALL dbcsr_create(mat_out, name="row_subset", dist=dist_sub, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=blk_sub, col_blk_size=col_blk_full)
      ELSE
         CALL dbcsr_distribution_new(dist_sub, template=dist_full, &
                                     row_dist=row_dist_full, col_dist=dist_sub_1d)
         CALL dbcsr_create(mat_out, name="col_subset", dist=dist_sub, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_full, col_blk_size=blk_sub)
      END IF

      CALL dbcsr_iterator_start(iter, mat_full)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk)
         IF (compress_rows) THEN
            IF (inv_map(ib) > 0) CALL dbcsr_put_block(mat_out, inv_map(ib), jb, blk)
         ELSE
            IF (inv_map(jb) > 0) CALL dbcsr_put_block(mat_out, ib, inv_map(jb), blk)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_finalize(mat_out)

      CALL dbcsr_distribution_release(dist_sub)
      DEALLOCATE (inv_map, blk_sub, dist_sub_1d)

      CALL timestop(handle)

   END SUBROUTINE extract_masked_blocks

! **************************************************************************************************
!> \brief Pre-seeds a square blocked DBCSR matrix with zero blocks only for block pairs whose
!>        centers lie within radius, for use with copy_fm_to_dbcsr(keep_sparsity=T) or
!>        dbcsr_multiply(retain_sparsity=T). Consumer: the RT-BSE CUTOFF_RADIUS_W0 truncation of
!>        the grid-basis W^0 (grid-blocked, centers = chunk_centroids).
!> \param matrix ...
!> \param centers block positions, one column per block row/column of matrix
!> \param radius truncation radius, same units as centers (bohr)
! **************************************************************************************************
   SUBROUTINE reserve_blocks_within_radius(matrix, centers, radius)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: centers
      REAL(KIND=dp), INTENT(IN)                          :: radius

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reserve_blocks_within_radius'

      INTEGER                                            :: handle, i, j, my_pcol, my_prow, &
                                                            nblkcols, nblkrows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      REAL(KIND=dp)                                      :: radius2
      REAL(KIND=dp), ALLOCATABLE                         :: zero_blk(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows, nblkcols_total=nblkcols, &
                          row_blk_size=row_blk, col_blk_size=col_blk, distribution=dist)
      CALL dbcsr_distribution_get(dist, row_dist=row_dist, col_dist=col_dist, &
                                  myprow=my_prow, mypcol=my_pcol)
      CPASSERT(nblkrows == SIZE(centers, 2))
      CPASSERT(nblkcols == SIZE(centers, 2))

      radius2 = radius**2
      ALLOCATE (zero_blk(MAXVAL(row_blk(1:nblkrows)), MAXVAL(col_blk(1:nblkcols))))
      zero_blk(:, :) = 0.0_dp

      DO i = 1, nblkrows
         IF (row_dist(i) /= my_prow) CYCLE
         DO j = 1, nblkcols
            IF (col_dist(j) /= my_pcol) CYCLE
            IF ((centers(1, i) - centers(1, j))**2 + (centers(2, i) - centers(2, j))**2 + &
                (centers(3, i) - centers(3, j))**2 <= radius2) THEN
               CALL dbcsr_put_block(matrix, i, j, zero_blk(1:row_blk(i), 1:col_blk(j)))
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(matrix)

      DEALLOCATE (zero_blk)
      CALL timestop(handle)

   END SUBROUTINE reserve_blocks_within_radius

! **************************************************************************************************
!> \brief Creates the (empty) result matrix of op(mat_left) * op(mat_right) with the correct block
!>        structure and a distribution on the shared process grid, ready to be filled by
!>        dbcsr_multiply. Row structure comes from op(left), column structure from op(right).
!> \param mat_left ...
!> \param mat_right ...
!> \param transa 'N' or 'T' applied to mat_left
!> \param transb 'N' or 'T' applied to mat_right
!> \param mat_out ...
! **************************************************************************************************
   SUBROUTINE create_product_matrix(mat_left, mat_right, transa, transb, mat_out)

      TYPE(dbcsr_type), INTENT(IN)                       :: mat_left, mat_right
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_out

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_product_matrix'

      INTEGER                                            :: handle, i, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_l, col_blk_r, out_col_blk, &
                                                            out_col_dist, out_row_blk, &
                                                            out_row_dist, row_blk_l, row_blk_r
      TYPE(dbcsr_distribution_type)                      :: dist_l, dist_out

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_left, distribution=dist_l, row_blk_size=row_blk_l, col_blk_size=col_blk_l)
      CALL dbcsr_get_info(mat_right, row_blk_size=row_blk_r, col_blk_size=col_blk_r)
      CALL dbcsr_distribution_get(dist_l, nprows=nprows, npcols=npcols)

      ! block SIZES follow op(left)/op(right); DISTRIBUTIONS are freshly round-robined onto the
      ! shared process grid (a transposed operand's row-dist is NOT a valid col-dist on a
      ! non-square grid). dbcsr_multiply redistributes internally, so any valid mapping works.
      IF (transa == 'N') THEN
         out_row_blk => row_blk_l
      ELSE
         out_row_blk => col_blk_l
      END IF
      IF (transb == 'N') THEN
         out_col_blk => col_blk_r
      ELSE
         out_col_blk => row_blk_r
      END IF

      ALLOCATE (out_row_dist(SIZE(out_row_blk)), out_col_dist(SIZE(out_col_blk)))
      DO i = 1, SIZE(out_row_blk)
         out_row_dist(i) = MOD(i - 1, nprows)
      END DO
      DO i = 1, SIZE(out_col_blk)
         out_col_dist(i) = MOD(i - 1, npcols)
      END DO

      CALL dbcsr_distribution_new(dist_out, template=dist_l, &
                                  row_dist=out_row_dist, col_dist=out_col_dist)
      CALL dbcsr_create(mat_out, name="panel_product", dist=dist_out, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=out_row_blk, col_blk_size=out_col_blk)
      CALL dbcsr_distribution_release(dist_out)
      DEALLOCATE (out_row_dist, out_col_dist)

      CALL timestop(handle)

   END SUBROUTINE create_product_matrix

! **************************************************************************************************
!> \brief Builds the AO-space Green's function operator G^occ/vir_µν (AO x AO DBCSR)
!> \param bs_env ...
!> \param tau ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
!> \param template ...
!> \param matrix_G_ao ...
! **************************************************************************************************
   SUBROUTINE build_G_ao(bs_env, tau, ispin, occ, vir, template, matrix_G_ao)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occ, vir
      TYPE(dbcsr_type), INTENT(INOUT)                    :: template
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix_G_ao

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_G_ao'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, dist_row_ao
      TYPE(cp_fm_type), POINTER                          :: fm_G
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao

      CALL timeset(routineN, handle)

      IF (occ) THEN
         fm_G => bs_env%fm_Gocc
      ELSE
         fm_G => bs_env%fm_Gvir
      END IF

      CALL G_occ_vir(bs_env, tau, fm_G, ispin, occ=occ, vir=vir)

      CALL setup_square_topology(template, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      CALL dbcsr_create(matrix_G_ao, name="G_ao", dist=dist_ao_ao, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=blk_ao, col_blk_size=blk_ao)

      CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(matrix_G_ao, bs_env%eps_filter)

      ! release only the topology; keep matrix_G_ao for the caller
      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao)

      CALL timestop(handle)

   END SUBROUTINE build_G_ao

! **************************************************************************************************
!> \brief Panel-streaming evaluation of  out += scale * L_out^T (A_grid ∘ B_grid) L_out,
!>        where A_grid = L_A M_A L_A^T and B_grid = L_B M_B L_B^T, WITHOUT ever forming the full
!>        grid x grid objects. The grid (row) index is processed in panels of ~panel_size rows; for
!>        each panel only P x grid slabs are built, Hadamard-multiplied, and contracted into the
!>        (small) output. Algebraically identical to L_out^T (A_grid ∘ B_grid) L_out summed over
!>        grid rows, so the result matches the non-streamed path to eps_filter.
!>
!>        Mapping (L in {phi (grid x AO), Z (grid x RI)}, M the AO/RI-space operator):
!>          chi :   L_A=L_B=phi, M_A=G_occ_ao, M_B=G_vir_ao, L_out=Z   -> RI x RI
!>          Sig :   L_A=phi (M_A=D/G), L_B=Z (M_B=V/W),       L_out=phi -> AO x AO
!> \param L_A ...
!> \param M_A ...
!> \param L_B ...
!> \param M_B ...
!> \param L_out ...
!> \param mat_out ...
!> \param scale ...
!> \param eps ...
!> \param para_env ...
!> \param pan_first ...
!> \param pan_last ...
!> \param lb_eq_la ...
!> \param lout_eq_la ...
!> \param zero_out ...
!> \param centroids ...
!> \param cutoff ...
!> \param grid_occupation ...
! **************************************************************************************************
   SUBROUTINE contract_grid_panels(L_A, M_A, L_B, M_B, L_out, mat_out, scale, eps, para_env, &
                                   pan_first, pan_last, lb_eq_la, lout_eq_la, zero_out, &
                                   centroids, cutoff, grid_occupation)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_A
      TYPE(dbcsr_type), INTENT(INOUT)                    :: M_A
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_B
      TYPE(dbcsr_type), INTENT(INOUT)                    :: M_B
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: L_out
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_out
      REAL(KIND=dp), INTENT(IN)                          :: scale, eps
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pan_first, pan_last
      LOGICAL, INTENT(IN)                                :: lb_eq_la, lout_eq_la, zero_out
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: grid_occupation

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_grid_panels'

      INTEGER                                            :: blk0, blk1, handle, ipan, n_grid_total, &
                                                            ncols_pan, nrows_pan
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap
      LOGICAL                                            :: use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, usedA, usedB
      TYPE(dbcsr_type)                                   :: A_pan, B_pan, C_pan, LA_pan, LA_panC, &
                                                            LB_pan, LB_panC, Lout_pan, MA_sub, &
                                                            MB_sub, tmp2, tmpA, tmpB
      TYPE(dbcsr_type), POINTER                          :: RB_A, RB_B, RB_out
      TYPE(dbcsr_type), TARGET                           :: LA_near, LB_near, Lout_near

      CALL timeset(routineN, handle)

      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp
      IF (PRESENT(grid_occupation)) grid_occupation = 0.0_dp

      CALL dbcsr_get_info(L_A, nfullrows_total=n_grid_total)

      IF (zero_out) CALL dbcsr_set(mat_out, 0.0_dp)

      DO ipan = 1, SIZE(pan_first)
         blk0 = pan_first(ipan)
         blk1 = pan_last(ipan)

         ! phi/Z panel slices (P x n)
         CALL extract_grid_panel(L_A, blk0, blk1, LA_pan)
         IF (.NOT. lb_eq_la) CALL extract_grid_panel(L_B, blk0, blk1, LB_pan)
         IF (.NOT. lout_eq_la) CALL extract_grid_panel(L_out, blk0, blk1, Lout_pan)

         ! Which AO/RI atoms (column blocks) actually touch this panel: the inner index of
         ! every multiply below is restricted to them, so only the matching rows of the
         ! system-wide operators M_A/M_B ever enter Cannon (exact: dropped rows meet zeros).
         CALL collect_used_col_blocks(LA_pan, para_env, usedA)
         IF (.NOT. lb_eq_la) THEN
            CALL collect_used_col_blocks(LB_pan, para_env, usedB)
         ELSE
            IF (ALLOCATED(usedB)) DEALLOCATE (usedB)
            ALLOCATE (usedB, SOURCE=usedA)
         END IF
         IF (.NOT. (ANY(usedA) .AND. ANY(usedB))) THEN
            ! empty panel slice: its Hadamard contribution is exactly zero on all ranks
            CALL dbcsr_release(LA_pan)
            IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_pan)
            IF (.NOT. lout_eq_la) CALL dbcsr_release(Lout_pan)
            CYCLE
         END IF

         ! Grid rows within reach of the panel: with the CUTOFF_RADIUS_RL_W truncation only
         ! they can appear as columns of the panel products / inner rows of the L_out multiply,
         ! so the system-wide phi/Z right operands are cut down to this neighborhood slice
         ! (local extraction, zero communication; exact w.r.t. the geo template).
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, grid_used)
            CALL extract_masked_blocks(L_A, grid_used, LA_near, compress_rows=.TRUE., blk_map=gmap)
            IF (.NOT. lb_eq_la) CALL extract_masked_blocks(L_B, grid_used, LB_near, compress_rows=.TRUE.)
            IF (.NOT. lout_eq_la) CALL extract_masked_blocks(L_out, grid_used, Lout_near, compress_rows=.TRUE.)
            RB_A => LA_near
         ELSE
            RB_A => L_A
         END IF
         IF (lb_eq_la) THEN
            RB_B => RB_A
         ELSE IF (use_cutoff) THEN
            RB_B => LB_near
         ELSE
            RB_B => L_B
         END IF
         IF (lout_eq_la) THEN
            RB_out => RB_A
         ELSE IF (use_cutoff) THEN
            RB_out => Lout_near
         ELSE
            RB_out => L_out
         END IF

         ! A_pan = LA_pan * M_A * L_A^T   (P x grid_near).
         ! When cutoff is active, A_pan is pre-seeded with only nearby blocks via
         ! build_geo_template_panel, and the multiply uses retain_sparsity to skip
         ! computing distant blocks entirely (exact: they are zero by locality).
         CALL extract_masked_blocks(LA_pan, usedA, LA_panC, compress_rows=.FALSE.)
         CALL extract_masked_blocks(M_A, usedA, MA_sub, compress_rows=.TRUE.)
         CALL create_product_matrix(LA_panC, MA_sub, 'N', 'N', tmpA)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, LA_panC, MA_sub, 0.0_dp, tmpA, filter_eps=eps)
         CALL dbcsr_release(MA_sub)
         IF (use_cutoff) THEN
            CALL build_geo_template_panel(LA_pan, LA_near, centroids, cutoff, blk0, A_pan, &
                                          col_map=gmap)
         ELSE
            CALL create_product_matrix(tmpA, RB_A, 'N', 'T', A_pan)
         END IF
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_A, 0.0_dp, A_pan, &
                             filter_eps=eps, retain_sparsity=use_cutoff)
         CALL dbcsr_release(tmpA)

         ! Grid-basis occupation of A_pan = ϕ G ϕ^T, accumulated over ALL panels into the
         ! occupation of the (never formed) full grid x grid object:
         !   sum_panels nnz(A_pan) / n_grid^2,  with nnz = occ * pan_rows * pan_cols.
         ! Panel-independent by construction -- a single-panel sample would instead report the
         ! local neighbor count of whichever region happens to land in that panel.
         IF (PRESENT(grid_occupation)) THEN
            CALL dbcsr_get_info(A_pan, nfullrows_total=nrows_pan, nfullcols_total=ncols_pan)
            grid_occupation = grid_occupation + dbcsr_get_occupation(A_pan)* &
                              REAL(ncols_pan, dp)*REAL(nrows_pan, dp)/ &
                              (REAL(n_grid_total, dp)*REAL(n_grid_total, dp))
         END IF

         ! B_pan = LB_pan * M_B * L_B^T   (P x grid_near); reuse the L_A slices when L_B == L_A.
         ! B_pan is pre-populated with A_pan's block structure so that
         ! retain_sparsity forces the final multiply to fill only those blocks (exact for ∘).
         IF (lb_eq_la) THEN
            CALL extract_masked_blocks(M_B, usedA, MB_sub, compress_rows=.TRUE.)
            CALL create_product_matrix(LA_panC, MB_sub, 'N', 'N', tmpB)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, LA_panC, MB_sub, 0.0_dp, tmpB, filter_eps=eps)
         ELSE
            CALL extract_masked_blocks(LB_pan, usedB, LB_panC, compress_rows=.FALSE.)
            CALL extract_masked_blocks(M_B, usedB, MB_sub, compress_rows=.TRUE.)
            CALL create_product_matrix(LB_panC, MB_sub, 'N', 'N', tmpB)
            CALL dbcsr_multiply('N', 'N', 1.0_dp, LB_panC, MB_sub, 0.0_dp, tmpB, filter_eps=eps)
            CALL dbcsr_release(LB_panC)
         END IF
         CALL dbcsr_release(MB_sub)
         CALL dbcsr_create(B_pan, template=A_pan)
         CALL dbcsr_copy(B_pan, A_pan)
         CALL dbcsr_set(B_pan, 0.0_dp)
         ! The output pattern is already fixed. Omitting redundant on-the-fly filtering also
         ! avoids overflowing DBCSR's single-precision screening norms for conditioned Z fits.
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_B, 0.0_dp, B_pan, &
                             retain_sparsity=.TRUE.)
         CALL dbcsr_release(tmpB)
         CALL dbcsr_release(LA_panC)

         ! C_pan = scale * (A_pan ∘ B_pan)   (P x grid_near)
         CALL dbcsr_create(C_pan, template=A_pan)
         CALL hadamard_product(A_pan, B_pan, C_pan, scale)
         CALL dbcsr_release(A_pan)
         CALL dbcsr_release(B_pan)

         ! tmp2 = C_pan * L_out   (P x n_out; inner index restricted to the neighborhood)
         CALL create_product_matrix(C_pan, RB_out, 'N', 'N', tmp2)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, RB_out, 0.0_dp, tmp2, filter_eps=eps)
         CALL dbcsr_release(C_pan)

         ! mat_out += L_out_pan^T * tmp2   (accumulate: beta = 1)
         IF (lout_eq_la) THEN
            CALL dbcsr_multiply('T', 'N', 1.0_dp, LA_pan, tmp2, 1.0_dp, mat_out, filter_eps=eps)
         ELSE
            CALL dbcsr_multiply('T', 'N', 1.0_dp, Lout_pan, tmp2, 1.0_dp, mat_out, filter_eps=eps)
            CALL dbcsr_release(Lout_pan)
         END IF
         CALL dbcsr_release(tmp2)
         IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_pan)
         CALL dbcsr_release(LA_pan)
         IF (use_cutoff) THEN
            CALL dbcsr_release(LA_near)
            IF (.NOT. lb_eq_la) CALL dbcsr_release(LB_near)
            IF (.NOT. lout_eq_la) CALL dbcsr_release(Lout_near)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE contract_grid_panels

! **************************************************************************************************
!> \brief Σ^c-specific panel loop: computes both the occupied (neg) and virtual (pos) contributions
!>        in a single pass over grid panels, forming W_pan = Z_panel × W_aux × Z^T only ONCE per
!>        panel and reusing it for both the G^occ and G^vir Hadamard contractions.
!>
!>        Computes:
!>          mat_Sigma_neg = ϕ^T ( (ϕ G^occ ϕ^T) ∘ (Z W^MIC Z^T) ) ϕ
!>          mat_Sigma_pos = ϕ^T ( (ϕ G^vir ϕ^T) ∘ (Z W^MIC Z^T) ) ϕ
!>
!> \param mat_phi ...
!> \param mat_Z ...
!> \param mat_G_occ_ao ...
!> \param mat_G_vir_ao ...
!> \param mat_W_aux ...
!> \param mat_Sigma_neg ...
!> \param mat_Sigma_pos ...
!> \param eps ...
!> \param para_env ...
!> \param pan_first ...
!> \param pan_last ...
!> \param centroids ...
!> \param cutoff ...
! **************************************************************************************************
   SUBROUTINE contract_grid_panels_sigma_c(mat_phi, mat_Z, mat_G_occ_ao, mat_G_vir_ao, &
                                           mat_W_aux, mat_Sigma_neg, mat_Sigma_pos, eps, para_env, &
                                           pan_first, pan_last, centroids, cutoff)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: mat_phi, mat_Z
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_G_occ_ao, mat_G_vir_ao, mat_W_aux, &
                                                            mat_Sigma_neg, mat_Sigma_pos
      REAL(KIND=dp), INTENT(IN)                          :: eps
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pan_first, pan_last
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: centroids
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_grid_panels_sigma_c'

      INTEGER                                            :: blk0, blk1, handle, ipan
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap
      LOGICAL                                            :: use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, used_ao, used_ri
      TYPE(dbcsr_type)                                   :: A_occ, A_vir, C_pan, G_occ_sub, &
                                                            G_vir_sub, phi_pan, phi_panC, tmp2, &
                                                            tmpA, tmpB, W_pan, W_sub, Z_pan, Z_panC
      TYPE(dbcsr_type), POINTER                          :: RB_phi, RB_Z
      TYPE(dbcsr_type), TARGET                           :: phi_near, Z_near

      CALL timeset(routineN, handle)

      use_cutoff = PRESENT(centroids) .AND. PRESENT(cutoff)
      IF (use_cutoff) use_cutoff = cutoff > 0.0_dp

      CALL dbcsr_set(mat_Sigma_neg, 0.0_dp)
      CALL dbcsr_set(mat_Sigma_pos, 0.0_dp)

      DO ipan = 1, SIZE(pan_first)
         blk0 = pan_first(ipan)
         blk1 = pan_last(ipan)

         CALL extract_grid_panel(mat_phi, blk0, blk1, phi_pan)
         CALL extract_grid_panel(mat_Z, blk0, blk1, Z_pan)

         ! AO/RI atoms touching this panel: only the matching rows of G_occ/G_vir/W ever
         ! enter the multiplies below (exact: dropped rows meet zero columns of the panel).
         CALL collect_used_col_blocks(phi_pan, para_env, used_ao)
         CALL collect_used_col_blocks(Z_pan, para_env, used_ri)
         IF (.NOT. (ANY(used_ao) .AND. ANY(used_ri))) THEN
            CALL dbcsr_release(phi_pan)
            CALL dbcsr_release(Z_pan)
            CYCLE
         END IF

         ! Neighborhood slices of phi/Z (grid rows within cutoff of the panel): they replace
         ! the system-wide right operands in every multiply (local extraction, zero comm).
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(centroids, blk0, blk1, cutoff, grid_used)
            CALL extract_masked_blocks(mat_phi, grid_used, phi_near, compress_rows=.TRUE., blk_map=gmap)
            CALL extract_masked_blocks(mat_Z, grid_used, Z_near, compress_rows=.TRUE.)
            RB_phi => phi_near
            RB_Z => Z_near
         ELSE
            RB_phi => mat_phi
            RB_Z => mat_Z
         END IF

         CALL extract_masked_blocks(phi_pan, used_ao, phi_panC, compress_rows=.FALSE.)
         CALL extract_masked_blocks(Z_pan, used_ri, Z_panC, compress_rows=.FALSE.)
         CALL extract_masked_blocks(mat_G_occ_ao, used_ao, G_occ_sub, compress_rows=.TRUE.)
         CALL extract_masked_blocks(mat_G_vir_ao, used_ao, G_vir_sub, compress_rows=.TRUE.)
         CALL extract_masked_blocks(mat_W_aux, used_ri, W_sub, compress_rows=.TRUE.)

         ! A_occ = phi_pan × G_occ × phi^T  (built first so W_pan can inherit its pattern).
         ! With cutoff active, A_occ is pre-seeded with geo-local blocks so that the
         ! phi^T multiply uses retain_sparsity and never computes distant blocks.
         CALL create_product_matrix(phi_panC, G_occ_sub, 'N', 'N', tmpA)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_occ_sub, 0.0_dp, tmpA, filter_eps=eps)
         IF (use_cutoff) THEN
            CALL build_geo_template_panel(phi_pan, phi_near, centroids, cutoff, blk0, A_occ, &
                                          col_map=gmap)
         ELSE
            CALL create_product_matrix(tmpA, RB_phi, 'N', 'T', A_occ)
         END IF
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_phi, 0.0_dp, A_occ, &
                             filter_eps=eps, retain_sparsity=use_cutoff)
         CALL dbcsr_release(tmpA)
         CALL dbcsr_release(G_occ_sub)

         ! A_vir = phi_pan × G_vir × phi^T  (same pre-screen as A_occ)
         CALL create_product_matrix(phi_panC, G_vir_sub, 'N', 'N', tmpA)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_vir_sub, 0.0_dp, tmpA, filter_eps=eps)
         IF (use_cutoff) THEN
            CALL build_geo_template_panel(phi_pan, phi_near, centroids, cutoff, blk0, A_vir, &
                                          col_map=gmap)
         ELSE
            CALL create_product_matrix(tmpA, RB_phi, 'N', 'T', A_vir)
         END IF
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_phi, 0.0_dp, A_vir, &
                             filter_eps=eps, retain_sparsity=use_cutoff)
         CALL dbcsr_release(tmpA)
         CALL dbcsr_release(G_vir_sub)

         ! W_pan = Z_pan × W_aux × Z^T  (computed once, reused for both Σ^c terms).
         ! W_pan is pre-seeded with the union of A_occ and A_vir block patterns so that
         ! retain_sparsity forces the multiply to fill only those blocks: exact since W
         ! outside G_occ∪G_vir is multiplied by zero in the Hadamard.
         CALL create_product_matrix(Z_panC, W_sub, 'N', 'N', tmpB)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, Z_panC, W_sub, 0.0_dp, tmpB, filter_eps=eps)
         CALL dbcsr_create(W_pan, template=A_occ)
         CALL dbcsr_copy(W_pan, A_occ)
         CALL dbcsr_add(W_pan, A_vir, 1.0_dp, 1.0_dp)
         CALL dbcsr_set(W_pan, 0.0_dp)
         ! The retained pattern is exact for the following Hadamard products; screening it
         ! again is redundant and may overflow DBCSR's single-precision block-norm product.
         CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpB, RB_Z, 0.0_dp, W_pan, &
                             retain_sparsity=.TRUE.)
         CALL dbcsr_release(tmpB)
         CALL dbcsr_release(W_sub)
         CALL dbcsr_release(phi_panC)
         CALL dbcsr_release(Z_panC)

         ! Σ^c_neg: ϕ^T ( A_occ ∘ W_pan ) ϕ
         CALL dbcsr_create(C_pan, template=A_occ)
         CALL hadamard_product(A_occ, W_pan, C_pan, 1.0_dp)
         CALL dbcsr_release(A_occ)
         CALL create_product_matrix(C_pan, RB_phi, 'N', 'N', tmp2)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, RB_phi, 0.0_dp, tmp2, filter_eps=eps)
         CALL dbcsr_release(C_pan)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, phi_pan, tmp2, 1.0_dp, mat_Sigma_neg, filter_eps=eps)
         CALL dbcsr_release(tmp2)

         ! Σ^c_pos: ϕ^T ( A_vir ∘ W_pan ) ϕ  — W_pan reused
         CALL dbcsr_create(C_pan, template=A_vir)
         CALL hadamard_product(A_vir, W_pan, C_pan, 1.0_dp)
         CALL dbcsr_release(A_vir)
         CALL create_product_matrix(C_pan, RB_phi, 'N', 'N', tmp2)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, RB_phi, 0.0_dp, tmp2, filter_eps=eps)
         CALL dbcsr_release(C_pan)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, phi_pan, tmp2, 1.0_dp, mat_Sigma_pos, filter_eps=eps)
         CALL dbcsr_release(tmp2)

         CALL dbcsr_release(W_pan)
         CALL dbcsr_release(Z_pan)
         CALL dbcsr_release(phi_pan)
         IF (use_cutoff) THEN
            CALL dbcsr_release(phi_near)
            CALL dbcsr_release(Z_near)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE contract_grid_panels_sigma_c

! **************************************************************************************************
   SUBROUTINE get_n_ao_in_sphere(bs_env, atom_P, cutoff_ri, n_ao_used)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom_P
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_ri
      INTEGER, INTENT(OUT)                               :: n_ao_used

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_n_ao_in_sphere'

      INTEGER                                            :: handle, ri_atom

      CALL timeset(routineN, handle)

      n_ao_used = 0
      DO ri_atom = 1, bs_env%n_atom
         IF (NORM2(bs_env%ri_rs%particle_set(ri_atom)%r(:) - &
                   bs_env%ri_rs%particle_set(atom_P)%r(:)) > &
             bs_env%ri_rs%radius_ao_per_atom(ri_atom) + cutoff_ri) CYCLE
         n_ao_used = n_ao_used + bs_env%i_ao_end_from_atom(ri_atom) - &
                     bs_env%i_ao_start_from_atom(ri_atom) + 1
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_n_ao_in_sphere

END MODULE gw_ri_rs_panels
