!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using RI-RS Approximation for molecules
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_large_cell_Gamma_ri_rs
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_add, dbcsr_copy, dbcsr_filter, dbcsr_finalize, dbcsr_iterator_stop, &
        dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale, &
        dbcsr_put_block, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE gw_large_cell_gamma,             ONLY: Fourier_transform_w_to_t,&
                                              compute_MinvVsqrt_Vsqrt,&
                                              compute_QP_energies,&
                                              delete_unnecessary_files,&
                                              compute_V_k_by_lattice_sum,&
                                              compute_cfm_W_ikp_freq_j,&
                                              compute_fm_chi_Gamma_freq,&
                                              create_fm_W_MIC_time,&
                                              fill_fm_Sigma_c_Gamma_time,&
                                              fm_write,&
                                              get_W_MIC,&
                                              multiply_fm_W_MIC_time_with_Minv_Gamma
   USE gw_non_periodic_ri_rs,           ONLY: atomic_basis_at_grid_point,&
                                              build_G_ao,&
                                              build_geo_template_panel,&
                                              collect_used_col_blocks,&
                                              compute_Sigma_c_and_QP_energies,&
                                              compute_Sigma_x,&
                                              compute_coeff_Z_lP,&
                                              create_product_matrix,&
                                              extract_grid_panel,&
                                              extract_masked_blocks,&
                                              get_mat_chi_Gamma_tau,&
                                              hadamard_product,&
                                              mask_grid_blocks_near_panel,&
                                              precompute_ri_rs_radii,&
                                              print_matrix_occupation,&
                                              print_ri_rs_memory_estimate,&
                                              resolve_grid_panels,&
                                              ri_rs_grid_assembler
   USE gw_utils,                        ONLY: de_init_bs_env
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE input_constants,                 ONLY: ri_rs_w_grid,&
                                              rtp_method_bse
   USE mathconstants,                   ONLY: twopi,&
                                              z_one,&
                                              z_zero
   USE post_scf_bandstructure_utils,    ONLY: MIC_contribution_from_ikp
   USE rpa_gw_kpoints_util,             ONLY: cp_cfm_power
   USE kinds,                           ONLY: dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_large_cell_Gamma_ri_rs'

   ! Bound on |n_i| for a nearest-image triple n = nint(h^-1 (r_l - r_l')). Atomic grid points
   ! extend beyond the cell faces, so r_l - r_l' can span more than one cell and |n_i| = 1 is
   ! NOT a safe bound; the value below is generous and only costs a (2m+1)^3 = 729-entry table.
   ! build_grid_phase_index asserts that the bound actually holds rather than clamping, because
   ! a clamped index would silently apply the wrong Bloch phase.
   INTEGER, PARAMETER, PRIVATE                        :: phase_n_max = 4

   PUBLIC :: gw_calc_large_cell_Gamma_ri_rs, &
             contract_A_B_A, &
             hadamard_product_inplace, &
             release_dbcsr_topology_and_matrices, &
             setup_square_topology

CONTAINS

! **************************************************************************************************
!> \brief GW calculation for a large periodic unit cell sampled at the Gamma point only, using
!>        the resolution of the identity in real space (RI-RS).
!>
!>        Every step except the screened interaction W is shared verbatim with the molecular
!>        driver in gw_non_periodic_ri_rs: the grid assembly, Phi, the Z_lP fit and the panel
!>        contractions for chi, Sigma^x and Sigma^c are all cell-agnostic once distances are
!>        measured through bs_env%ri_rs%cell and the periodic image sums are taken inside
!>        fill_phi_for_atom / compute_d_lp. What is genuinely periodic is get_W_MIC, which
!>        builds W on the chi/eps/W k-mesh and folds it back to the minimum image.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************

   SUBROUTINE gw_calc_large_cell_Gamma_ri_rs(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_large_cell_Gamma_ri_rs'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_grid_tau, mat_W_grid_time
      TYPE(dbcsr_type)                                   :: mat_W_pattern

      CALL timeset(routineN, handle)

      NULLIFY (mat_chi_grid_tau, mat_W_grid_time)

      ! Cache the cell so every geometric screening test measures minimum-image distances.
      ! Without this the panel/cutoff screening silently falls back to Euclidean distances,
      ! which UNDER-include under PBC (d_euclid >= d_MIC) and would drop boundary-wrapping
      ! pairs that are genuinely within the cutoff.
      CALL get_qs_env(qs_env, cell=bs_env%ri_rs%cell)

       !!========================================================================
       !! 0. Precompute AO and RI Radii
       !!    Per-atom cutoff radii from the most diffuse Gaussian primitives in
       !!    the AO and RI auxiliary basis sets. Stored in bs_env%ri_rs%
       !!    radius_ao_per_atom and radius_ri_per_atom, used for sphere-cutoff
       !!    and phi_local screening.
       !!========================================================================
      CALL precompute_ri_rs_radii(qs_env, bs_env)

       !!========================================================================
       !! 1. Grid Generation for RI-RS
       !!    (Modified Lebedev grids from Ivan Duchemin and Xavier Blase)
       !!    Generate flattened 1D array of grid points for RI-RS.
       !!    Equation: r_g(k) = R_A + r_g(A)
       !!========================================================================
      CALL ri_rs_grid_assembler(qs_env, bs_env, bs_env%ri_rs%grid_points)

       !!========================================================================
       !! 2. Atomic Basis Evaluation
       !!    Compute values of spherical atomic basis functions at grid points.
       !!    Expression: Φ_μl = Φ_μ(r_l) (mat_phi_mu_l)
       !!========================================================================
      CALL atomic_basis_at_grid_point(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)

       !!========================================================================
       !! 3. Compute RI-RS Coefficients (Z_lp)
       !!    Solve the regularized system for each atom P, where the grid domain
       !!    is restricted to r_l within a cutoff distance of atom P:
       !!    a. D_ll' = [ Σ_μ Φ_μ(r_l) Φ_μ(r_l') ]^2 (Equation 13)
       !!    b. D_lP  = Σ_{μν} Φ_μ(r_l) Φ_ν(r_l) (μν|P) (Equation 15)
       !!    c. Conditioning:
       !!       Dvec_l   = 1 / sqrt(D_ll)  (Diagonal scaling vector)
       !!       D'_ll' = Dvec_l * D_ll' * Dvec_l' + λδ_ll'
       !!       D'_lP  = Dvec_l * D_lP
       !!    d. Solve: Σ_l' D'_ll' * Z'_l'P = D'_lP (Equation 14)
       !!    e. Rescale: Z_lP = Z'_lP * Dvec_l   (Z_lP stored in mat_Z_lP)
       !!========================================================================
      CALL print_ri_rs_memory_estimate(qs_env, bs_env)

      CALL compute_coeff_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                              bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! flag the RI-RS grid as built so a subsequent RT-BSE run reuses Z_lP
      ! instead of rebuilding it
      bs_env%ri_rs%grid_built = .TRUE.

       !!========================================================================
       !! 4/5. Polarizability χ and screened interaction W^MIC.
       !!
       !!    W_TRUNCATION AUX (default): χ is contracted to the RI auxiliary basis right
       !!    away and W_PQ is truncated at atom-pair granularity. Validated -- reproduces
       !!    the analytic reference exactly -- but the auxiliary functions carry their own
       !!    extent r_RI, so the lattice sum needs a_min/2 > r_G + 2 r_RI.
       !!
       !!    W_TRUNCATION GRID: χ stays in the grid basis and W is truncated per grid point,
       !!    which converges already for a_min/2 > r_G and should therefore need a smaller
       !!    cell. Experimental; see grid_W_truncation_periodic_GW.pdf.
       !!========================================================================
      IF (bs_env%ri_rs%w_truncation == ri_rs_w_grid) THEN

         CALL get_mat_chi_grid_tau(bs_env, mat_chi_grid_tau, mat_W_pattern, &
                                   bs_env%ri_rs%mat_phi_mu_l)

         CALL compute_W_MIC_grid(bs_env, qs_env, mat_chi_grid_tau, mat_W_pattern, &
                                 mat_W_grid_time)

      ELSE

         CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                    bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

         CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

      END IF

      ! The grid coordinates are only needed for the minimum-image Bloch phases inside
      ! compute_W_MIC_grid; both branches are done with them here.
      IF (ALLOCATED(bs_env%ri_rs%grid_points)) DEALLOCATE (bs_env%ri_rs%grid_points)

       !!========================================================================
       !! 6. Compute Exact Exchange Self-Energy (Σ^x)
       !!    D_µν          = sum_n^occ C_µn C_νn
       !!    D_ll'         = sum_µν Φ_µ(r_l) D_µν Φ_ν(r_l')
       !!    V^trunc_ll'   = sum_PQ Z_lP V^trunc_PQ Z_l'Q
       !!    Σ^x_ll'       = D_ll' * V^trunc_ll'
       !!    Σ^x_λσ(k=0)   = -sum_ll' Φ_λ(r_l) Σ^x_ll' Φ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

       !!========================================================================
       !! 7. Correlation self-energy Σ^c and quasiparticle energies, iterated to
       !!    eigenvalue self-consistency when &EVGW0 is given and done in a single
       !!    pass for G0W0.
       !!
       !!    (a) W^MIC_ll'(iτ) =  sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
       !!        Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) * W^MIC_ll'(iτ), for τ < 0
       !!        Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) * W^MIC_ll'(iτ), for τ > 0
       !!        Σ^c_λσ(iτ)    =  sum_ll' Φ_λ(r_l) Σ^c_ll'(iτ) Φ_σ(r_l')
       !!    (b) Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
       !!        ϵ_nk^GW,(i) = ϵ_nk^DFT + Σ^c_nn[G^(i-1),W](ϵ_nk^GW,(i)) + Σ^x_nn - v^xc_nn
       !!========================================================================
      IF (bs_env%ri_rs%w_truncation == ri_rs_w_grid) THEN

         ! evGW0 is not available on the grid-W path yet: it would need W^MIC_ll'(iτ) to
         ! survive the eigenvalue loop, which is exactly the object we are trying to keep
         ! small. Single-shot G0W0 only; gw_utils rejects &EVGW0 + GRID at input time.
         CALL compute_Sigma_c_grid(bs_env, mat_W_grid_time, bs_env%ri_rs%mat_phi_mu_l, &
                                   fm_Sigma_c_Gamma_time)

         CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

         CALL cp_fm_release(fm_Sigma_x_Gamma)
         CALL delete_unnecessary_files(bs_env)

      ELSE

         CALL compute_Sigma_c_and_QP_energies(bs_env, qs_env, fm_W_time, fm_Sigma_x_Gamma)

      END IF

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_large_cell_Gamma_ri_rs

! **************************************************************************************************
!> \brief Evaluates atomic basis functions on a real-space grid and builds a sparse DBCSR matrix.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
! **************************************************************************************************
!> \brief Compute RI-RS Coefficients (Z_lP)
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************
!> \brief Computes the dense localized RHS d_lp(l,P) = Σ_{μν,R,S} Φ_μ(r_l)·Φ_ν(r_l)·(μν|P) for one
!>        RI atom P. OMP-threaded over (atom_j, atom_k) AO-pair blocks: per thread, sweep all
!>        (cell_R, cell_S) periodic images of (atom_j, atom_k) about atom_P at cell (0,0,0); each
!>        3c block is built by build_3c_integral_block_ctx (cached libint / sphi tables in ctx,
!>        kind-radius triangle screen → `screened` short-circuits negligible image triples), and
!>        grid-chunked pair densities are contracted into a private d_lp partial that is reduced
!>        into d_lp at the end of the parallel region.
!> \param bs_env ...
!> \param ctx shared 3c-integral context (gw_3c_ctx_create)
!> \param cell ...
!> \param phi_val Φ_μ(r_l) on the local-sphere grid (n_grid_total × n_ao)
!> \param d_lp output (n_grid_total × n_loc_ri), zeroed by the caller, accumulated here
!> \param n_grid_total number of local-sphere grid rows
!> \param n_loc_ri number of RI functions of atom_P
!> \param atom_P RI atom (pinned to cell (0,0,0))
!> \param max_ao_size ...
!> \param atom_j_mepos ...
!> \param atom_j_stride ...
! **************************************************************************************************
!> \brief Computes the χ(iτ, k=0) matrix
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************
!> \brief Computes Green's Function in grid basis
!> \param bs_env ...
!> \param tau ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
!> \param mat_phi_mu_l ...
!> \param matrix_G_grid ...
!> \param eps_filter ...
! **************************************************************************************************
!> \brief Generalized routine to compute OUT = A * B * A^T  OR  OUT = A^T * B * A using DBCSR
!> \param transA_left ...
!> \param transA_right ...
!> \param matrix_A ...
!> \param matrix_B ...
!> \param matrix_out ...
!> \param eps_filter ...
! **************************************************************************************************
!> \brief Builds the Gamma-point polarizability in the GRID basis, chi_ll'(i tau, k=0), and at the
!>        same time the block pattern that the screened interaction will be stored on.
!>
!>          chi_ll'(i tau) = sum_s g_s (phi G^occ_s phi^T)_ll' o (phi G^vir_s phi^T)_ll'
!>
!>        Unlike the auxiliary-basis path this stops before contracting with Z, because the
!>        minimum image of chi has to be taken per grid point. The matrix is still never
!>        materialized as a dense object: it is accumulated panel by panel, exactly as
!>        contract_grid_panels does, and each panel is scattered into the (sparse) accumulator
!>        with zero communication.
!>
!>        The W pattern is the union over panels, spins and time points of supp(A_occ) and
!>        supp(A_vir). That is the right bound because W enters the self-energy only through
!>        Sigma^c = G o W, so any block of W outside supp(G) is multiplied by zero. It is a much
!>        tighter bound than the a_min/2 minimum-image sphere, whose volume is 0.524 a^3 -- i.e.
!>        ~52 % dense, growing as N_grid^2 -- whereas supp(G) is set by the decay length r_G and
!>        so stays linear in system size.
!> \param bs_env ...
!> \param mat_chi_grid_tau  out, n_tau grid x grid matrices, finalized and filtered
!> \param mat_W_pattern     out, zero-valued grid x grid matrix carrying supp(G^occ) u supp(G^vir)
!> \param mat_phi_mu_l ...
! **************************************************************************************************
   SUBROUTINE get_mat_chi_grid_tau(bs_env, mat_chi_grid_tau, mat_W_pattern, mat_phi_mu_l)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_grid_tau
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_W_pattern
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: mat_phi_mu_l

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_grid_tau'

      INTEGER                                            :: blk0, blk1, handle, i, i_t, ipan, ispin, &
                                                            n_panels, npcol
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap, pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: blk_grid, dist_col_grid, dist_row_grid
      LOGICAL                                            :: have_U, use_cutoff
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, used_ao
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_grid_grid, dist_phi
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_G_occ_ao, mat_G_vir_ao
      TYPE(dbcsr_type)                                   :: A_occ, A_vir, C_pan, G_occ_sub, &
                                                            G_vir_sub, phi_pan, phi_panC, S_pan, &
                                                            tmpA, U_pan
      TYPE(dbcsr_type), POINTER                          :: RB_phi
      TYPE(dbcsr_type), TARGET                           :: phi_near

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      use_cutoff = bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .AND. &
                   ALLOCATED(bs_env%ri_rs%chunk_centroids)

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)
      n_panels = SIZE(pan_first)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,I44)') 'Number of batches for grid χ, W, Σ', n_panels
         CALL m_flush(bs_env%unit_nr)
      END IF

      ! grid x grid topology: rows inherit phi's blocking and row distribution, so panels of
      ! these matrices are co-resident with panels of phi; columns are round-robined.
      CALL dbcsr_get_info(mat_phi_mu_l, distribution=dist_phi, row_blk_size=blk_grid)
      CALL dbcsr_distribution_get(dist_phi, row_dist=dist_row_grid, npcols=npcol)
      ALLOCATE (dist_col_grid(SIZE(blk_grid)))
      DO i = 1, SIZE(blk_grid)
         dist_col_grid(i) = MOD(i - 1, npcol)
      END DO
      CALL dbcsr_distribution_new(dist_grid_grid, template=dist_phi, &
                                  row_dist=dist_row_grid, col_dist=dist_col_grid)

      NULLIFY (mat_chi_grid_tau)
      ALLOCATE (mat_chi_grid_tau(bs_env%num_time_freq_points))
      DO i_t = 1, bs_env%num_time_freq_points
         ALLOCATE (mat_chi_grid_tau(i_t)%matrix)
         CALL dbcsr_create(mat_chi_grid_tau(i_t)%matrix, "chi_grid", dist_grid_grid, &
                           dbcsr_type_no_symmetry, blk_grid, blk_grid)
      END DO
      CALL dbcsr_create(mat_W_pattern, "W_grid_pattern", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)

      ! AO-space Green's functions for every (tau, spin), built once. They are AO x AO and
      ! therefore small; precomputing them lets the panel loop be the OUTER one, which is what
      ! guarantees that each grid block of the accumulators is written exactly once.
      NULLIFY (mat_G_occ_ao, mat_G_vir_ao)
      ALLOCATE (mat_G_occ_ao(bs_env%num_time_freq_points, bs_env%n_spin))
      ALLOCATE (mat_G_vir_ao(bs_env%num_time_freq_points, bs_env%n_spin))
      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%imag_time_points(i_t)
         DO ispin = 1, bs_env%n_spin
            ALLOCATE (mat_G_occ_ao(i_t, ispin)%matrix)
            ALLOCATE (mat_G_vir_ao(i_t, ispin)%matrix)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, &
                            mat_G_occ_ao(i_t, ispin)%matrix)
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, &
                            mat_G_vir_ao(i_t, ispin)%matrix)
         END DO
      END DO

      DO ipan = 1, n_panels
         blk0 = pan_first(ipan)
         blk1 = pan_last(ipan)

         CALL extract_grid_panel(mat_phi_mu_l, blk0, blk1, phi_pan)
         CALL collect_used_col_blocks(phi_pan, bs_env%para_env, used_ao)
         IF (.NOT. ANY(used_ao)) THEN
            CALL dbcsr_release(phi_pan)
            CYCLE
         END IF

         ! Neighborhood slice of phi. gmap maps the compressed column index back to the global
         ! grid block, which the scatter into the full accumulators needs.
         IF (ALLOCATED(gmap)) DEALLOCATE (gmap)
         IF (use_cutoff) THEN
            CALL mask_grid_blocks_near_panel(bs_env%ri_rs%chunk_centroids, blk0, blk1, &
                                             bs_env%ri_rs%cutoff_radius_v_w, grid_used, &
                                             bs_env%ri_rs%cell)
            CALL extract_masked_blocks(mat_phi_mu_l, grid_used, phi_near, &
                                       compress_rows=.TRUE., blk_map=gmap)
            RB_phi => phi_near
         ELSE
            RB_phi => mat_phi_mu_l
            ALLOCATE (gmap(SIZE(blk_grid)))
            DO i = 1, SIZE(blk_grid)
               gmap(i) = i
            END DO
         END IF

         CALL extract_masked_blocks(phi_pan, used_ao, phi_panC, compress_rows=.FALSE.)

         have_U = .FALSE.

         DO i_t = 1, bs_env%num_time_freq_points

            DO ispin = 1, bs_env%n_spin

               CALL extract_masked_blocks(mat_G_occ_ao(i_t, ispin)%matrix, used_ao, G_occ_sub, &
                                          compress_rows=.TRUE.)
               CALL extract_masked_blocks(mat_G_vir_ao(i_t, ispin)%matrix, used_ao, G_vir_sub, &
                                          compress_rows=.TRUE.)

               ! A_occ = phi_pan G^occ phi^T
               CALL create_product_matrix(phi_panC, G_occ_sub, 'N', 'N', tmpA)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_occ_sub, 0.0_dp, tmpA, &
                                   filter_eps=bs_env%eps_filter)
               IF (use_cutoff) THEN
                  CALL build_geo_template_panel(phi_pan, phi_near, bs_env%ri_rs%chunk_centroids, &
                                                bs_env%ri_rs%cutoff_radius_v_w, blk0, A_occ, &
                                                col_map=gmap, cell=bs_env%ri_rs%cell)
               ELSE
                  CALL create_product_matrix(tmpA, RB_phi, 'N', 'T', A_occ)
               END IF
               CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_phi, 0.0_dp, A_occ, &
                                   filter_eps=bs_env%eps_filter, retain_sparsity=use_cutoff)
               CALL dbcsr_release(tmpA)
               CALL dbcsr_release(G_occ_sub)

               ! A_vir = phi_pan G^vir phi^T
               CALL create_product_matrix(phi_panC, G_vir_sub, 'N', 'N', tmpA)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_vir_sub, 0.0_dp, tmpA, &
                                   filter_eps=bs_env%eps_filter)
               IF (use_cutoff) THEN
                  CALL build_geo_template_panel(phi_pan, phi_near, bs_env%ri_rs%chunk_centroids, &
                                                bs_env%ri_rs%cutoff_radius_v_w, blk0, A_vir, &
                                                col_map=gmap, cell=bs_env%ri_rs%cell)
               ELSE
                  CALL create_product_matrix(tmpA, RB_phi, 'N', 'T', A_vir)
               END IF
               CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, RB_phi, 0.0_dp, A_vir, &
                                   filter_eps=bs_env%eps_filter, retain_sparsity=use_cutoff)
               CALL dbcsr_release(tmpA)
               CALL dbcsr_release(G_vir_sub)

               ! spin contribution g_s * A_occ o A_vir, summed into S_pan before scattering so
               ! that each block of the accumulator is inserted exactly once per (panel, tau)
               CALL dbcsr_create(C_pan, template=A_occ)
               CALL hadamard_product(A_occ, A_vir, C_pan, bs_env%spin_degeneracy)
               IF (ispin == 1) THEN
                  CALL dbcsr_create(S_pan, template=C_pan)
                  CALL dbcsr_copy(S_pan, C_pan)
               ELSE
                  CALL dbcsr_add(S_pan, C_pan, 1.0_dp, 1.0_dp)
               END IF
               CALL dbcsr_release(C_pan)

               ! union of supp(A_occ) and supp(A_vir) over all tau and spins of this panel
               IF (.NOT. have_U) THEN
                  CALL dbcsr_create(U_pan, template=A_occ)
                  CALL dbcsr_copy(U_pan, A_occ)
                  have_U = .TRUE.
               ELSE
                  CALL dbcsr_add(U_pan, A_occ, 1.0_dp, 1.0_dp)
               END IF
               CALL dbcsr_add(U_pan, A_vir, 1.0_dp, 1.0_dp)

               CALL dbcsr_release(A_occ)
               CALL dbcsr_release(A_vir)

            END DO ! ispin

            CALL dbcsr_filter(S_pan, bs_env%eps_filter)
            CALL accumulate_grid_panel(S_pan, blk0, mat_chi_grid_tau(i_t)%matrix, 1.0_dp, gmap)
            CALL dbcsr_release(S_pan)

         END DO ! i_t

         ! Reserve the support of G for W once per panel, with zero values.
         IF (have_U) THEN
            CALL dbcsr_set(U_pan, 0.0_dp)
            CALL reserve_grid_panel_blocks(U_pan, blk0, mat_W_pattern, gmap)
            CALL dbcsr_release(U_pan)
         END IF

         CALL dbcsr_release(phi_panC)
         CALL dbcsr_release(phi_pan)
         IF (use_cutoff) CALL dbcsr_release(phi_near)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I17,A,I3,A,F7.1,A)') &
               'Computed χ(iτ,k=0) for batch', ipan, ' /', n_panels, &
               ',    Execution time', m_walltime() - t1, ' s'
            CALL m_flush(bs_env%unit_nr)
         END IF

      END DO ! ipan

      DO i_t = 1, bs_env%num_time_freq_points
         CALL dbcsr_finalize(mat_chi_grid_tau(i_t)%matrix)
         CALL dbcsr_filter(mat_chi_grid_tau(i_t)%matrix, bs_env%eps_filter)
      END DO
      CALL dbcsr_finalize(mat_W_pattern)

      CALL print_matrix_occupation(mat_chi_grid_tau(1)%matrix, 'χ(l,l'')', &
                                   bs_env%para_env, bs_env%unit_nr)
      CALL print_matrix_occupation(mat_W_pattern, 'W(l,l'') pattern = supp(G)', &
                                   bs_env%para_env, bs_env%unit_nr)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL dbcsr_deallocate_matrix_set(mat_G_occ_ao)
      CALL dbcsr_deallocate_matrix_set(mat_G_vir_ao)
      CALL dbcsr_distribution_release(dist_grid_grid)
      DEALLOCATE (dist_col_grid)
      IF (ALLOCATED(gmap)) DEALLOCATE (gmap)

      CALL timestop(handle)

   END SUBROUTINE get_mat_chi_grid_tau

! **************************************************************************************************
   SUBROUTINE accumulate_grid_panel(mat_pan, blk0, mat_full, alpha, col_map, retain_sparsity)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_pan, mat_full
      INTEGER, INTENT(IN)                                :: blk0
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: col_map
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity

      CHARACTER(LEN=*), PARAMETER :: routineN = 'accumulate_grid_panel'

      INTEGER                                            :: col_full, handle, ib, jb, row_full
      LOGICAL                                            :: found, my_retain
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scaled
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_full, blk_pan
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      my_retain = .FALSE.
      IF (PRESENT(retain_sparsity)) my_retain = retain_sparsity

      CALL dbcsr_iterator_start(iter, mat_pan)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk_pan)

         row_full = blk0 + ib - 1
         IF (PRESENT(col_map)) THEN
            col_full = col_map(jb)
         ELSE
            col_full = jb
         END IF

         IF (my_retain) THEN
            ! Templated accumulation: blocks outside the reserved pattern are dropped. That is
            ! the intended truncation, not a lost contribution -- the template is a superset of
            ! supp(G), and anything outside it meets a zero in the Hadamard product with G.
            CALL dbcsr_get_block_p(mat_full, row_full, col_full, blk_full, found)
            IF (found) blk_full(:, :) = blk_full(:, :) + alpha*blk_pan(:, :)
         ELSE
            IF (ALLOCATED(scaled)) THEN
               IF (SIZE(scaled, 1) /= SIZE(blk_pan, 1) .OR. SIZE(scaled, 2) /= SIZE(blk_pan, 2)) &
                  DEALLOCATE (scaled)
            END IF
            IF (.NOT. ALLOCATED(scaled)) ALLOCATE (scaled(SIZE(blk_pan, 1), SIZE(blk_pan, 2)))
            scaled(:, :) = alpha*blk_pan(:, :)
            CALL dbcsr_put_block(mat_full, row_full, col_full, scaled, summation=.TRUE.)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE accumulate_grid_panel

! **************************************************************************************************
!> \brief Reserves (with zero values) the blocks of a grid panel in its full grid x grid parent:
!>          mat_full(blk0+ib-1, col_map(jb)) <- 0
!>        Unlike accumulate_grid_panel this inserts without summation, so reserving the same
!>        block repeatedly -- which happens because the support of G is collected over every
!>        time point and spin -- is idempotent rather than stacking duplicate block entries.
!>        The caller finalizes mat_full once all panels have been reserved.
!> \param mat_pan ...
!> \param blk0 ...
!> \param mat_full ...
!> \param col_map ...
! **************************************************************************************************
   SUBROUTINE reserve_grid_panel_blocks(mat_pan, blk0, mat_full, col_map)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_pan, mat_full
      INTEGER, INTENT(IN)                                :: blk0
      INTEGER, DIMENSION(:), INTENT(IN)                  :: col_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reserve_grid_panel_blocks'

      INTEGER                                            :: handle, ib, jb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: zero_blk
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_pan
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, mat_pan)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ib, jb, blk_pan)
         IF (ALLOCATED(zero_blk)) THEN
            IF (SIZE(zero_blk, 1) /= SIZE(blk_pan, 1) .OR. &
                SIZE(zero_blk, 2) /= SIZE(blk_pan, 2)) DEALLOCATE (zero_blk)
         END IF
         IF (.NOT. ALLOCATED(zero_blk)) THEN
            ALLOCATE (zero_blk(SIZE(blk_pan, 1), SIZE(blk_pan, 2)))
            zero_blk(:, :) = 0.0_dp
         END IF
         CALL dbcsr_put_block(mat_full, blk0 + ib - 1, col_map(jb), zero_blk)
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE reserve_grid_panel_blocks

! **************************************************************************************************
!> \brief Precomputes, for every element of a grid panel, a packed index encoding the nearest-image
!>        integer triple n(l,l') = nint( h^-1 (r_l - r_l') ).
!>
!>        n depends only on the geometry -- not on k, tau or omega -- yet the Bloch phase
!>        e^{i 2 pi k . n(l,l')} has to be applied on every (k, tau) and (k, omega) sweep. Computing
!>        n (a 3x3 matvec plus three NINTs) and then a COS and a SIN per element on each of those
!>        sweeps dominated the grid-W runtime. Here it is evaluated once and stored as
!>          p = (n1+m)(2m+1)^2 + (n2+m)(2m+1) + (n3+m) + 1,   m = phase_n_max
!>        so the inner loop reduces to a table lookup and two multiplies.
!>
!>        The index is carried in a DBCSR matrix sharing the panel's block structure, which keeps
!>        the lookup robust against DBCSR's block iteration order, and lets the same matrix double
!>        as the pattern template for W.
!> \param mat_pan panel whose structure is cached (the W pattern panel: a superset of chi's)
!> \param blk0 first grid row block of the panel within the parent
!> \param grid_points ...
!> \param row_off_parent cumulative row offsets of the parent grid x grid matrix
!> \param h_inv ...
!> \param mat_n_idx out, same structure as mat_pan, values = packed index p
! **************************************************************************************************
   SUBROUTINE build_grid_phase_index(mat_pan, blk0, grid_points, row_off_parent, h_inv, mat_n_idx)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_pan
      INTEGER, INTENT(IN)                                :: blk0
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_off_parent
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_n_idx

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_grid_phase_index'

      INTEGER                                            :: col, handle, i, j, l, lp, m, row, span
      INTEGER, DIMENSION(3)                              :: n
      REAL(KIND=dp), DIMENSION(3)                        :: raw
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_idx
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_create(mat_n_idx, template=mat_pan)
      CALL dbcsr_copy(mat_n_idx, mat_pan)

      m = phase_n_max
      span = 2*m + 1

      CALL dbcsr_iterator_start(iter, mat_n_idx)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk_idx)
         DO j = 1, SIZE(blk_idx, 2)
            lp = row_off_parent(col) + j
            DO i = 1, SIZE(blk_idx, 1)
               l = row_off_parent(blk0 + row - 1) + i
               raw(:) = grid_points(:, l) - grid_points(:, lp)
               n(:) = NINT(MATMUL(h_inv, raw))
               CPASSERT(MAXVAL(ABS(n)) <= m)
               blk_idx(i, j) = REAL((n(1) + m)*span*span + (n(2) + m)*span + (n(3) + m) + 1, dp)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE build_grid_phase_index

! **************************************************************************************************
!> \brief Tabulates cos/sin of the Bloch phase 2 pi k . n for every integer triple n that
!>        build_grid_phase_index can produce, indexed by the same packed index.
!> \param xkp_k fractional k-point coordinates
!> \param ctab out, cos table
!> \param stab out, sin table
! **************************************************************************************************
   SUBROUTINE build_phase_table(xkp_k, ctab, stab)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: ctab, stab

      INTEGER                                            :: i1, i2, i3, m, p, span
      REAL(KIND=dp)                                      :: arg

      m = phase_n_max
      span = 2*m + 1
      IF (.NOT. ALLOCATED(ctab)) ALLOCATE (ctab(span*span*span), stab(span*span*span))

      DO i1 = -m, m
         DO i2 = -m, m
            DO i3 = -m, m
               p = (i1 + m)*span*span + (i2 + m)*span + (i3 + m) + 1
               arg = twopi*(xkp_k(1)*i1 + xkp_k(2)*i2 + xkp_k(3)*i3)
               ctab(p) = COS(arg)
               stab(p) = SIN(arg)
            END DO
         END DO
      END DO

   END SUBROUTINE build_phase_table

! **************************************************************************************************
!> \brief Applies the k-point Bloch phase of the grid minimum image to a panel of a grid x grid
!>        matrix, splitting the result into real and imaginary parts:
!>
!>          O_ll'(k) = e^{i 2 pi k . n(l,l')} O_ll'(k=0),   n(l,l') = nint( h^-1 (r_l - r_l') )
!>
!>        n is the nearest image (the argmin) and is kept unconditionally -- there is no a_min/2
!>        cut. That is deliberate: the hard mask was found to over-truncate, discarding pairs whose
!>        nearest-image distance lies between a_min/2 and the Wigner-Seitz boundary, which the
!>        reference keeps. See grid_W_from_chi_equations.md.
!>
!>        With conj = .TRUE. the conjugate phase e^{-i 2 pi k . n} is applied and only the real
!>        part of the product is kept, which is the inverse BZ transform back to real space.
!> \param mat_re in: O(k=0) panel; out: Re O(k) (or the real part of the inverse transform)
!> \param mat_im out: Im O(k); read but not written when conj = .TRUE.
!> \param mat_n_idx packed nearest-image index from build_grid_phase_index (superset structure)
!> \param ctab cos table for this k (build_phase_table)
!> \param stab sin table for this k
!> \param conj apply the conjugate phase and keep only the real part
! **************************************************************************************************
   SUBROUTINE apply_grid_mic_phase(mat_re, mat_im, mat_n_idx, ctab, stab, conj)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_re, mat_im, mat_n_idx
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: ctab, stab
      LOGICAL, INTENT(IN)                                :: conj

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_grid_mic_phase'

      INTEGER                                            :: col, handle, i, j, p, row
      LOGICAL                                            :: found_idx, found_im
      REAL(KIND=dp)                                      :: c, s, v_im, v_re
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_idx, blk_im, blk_re
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, mat_re)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk_re)
         CALL dbcsr_get_block_p(mat_im, row, col, blk_im, found_im)
         IF (.NOT. found_im) CYCLE
         CALL dbcsr_get_block_p(mat_n_idx, row, col, blk_idx, found_idx)
         ! The index matrix is built on the union of the supports, so a block present here must
         ! be present there; a miss would mean a phase silently skipped.
         CPASSERT(found_idx)

         IF (conj) THEN
            DO j = 1, SIZE(blk_re, 2)
               DO i = 1, SIZE(blk_re, 1)
                  p = NINT(blk_idx(i, j))
                  ! Re[ e^{-i arg} (re + i im) ] = re cos + im sin
                  blk_re(i, j) = blk_re(i, j)*ctab(p) + blk_im(i, j)*stab(p)
               END DO
            END DO
         ELSE
            DO j = 1, SIZE(blk_re, 2)
               DO i = 1, SIZE(blk_re, 1)
                  p = NINT(blk_idx(i, j))
                  c = ctab(p)
                  s = stab(p)
                  v_re = blk_re(i, j)
                  v_im = blk_im(i, j)
                  blk_re(i, j) = v_re*c - v_im*s
                  blk_im(i, j) = v_re*s + v_im*c
               END DO
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE apply_grid_mic_phase

! **************************************************************************************************
!> \brief Builds the real and imaginary parts of the RI-RS fitting coefficients at k-point k:
!>          Z_lP(k) = e^{i 2 pi k . n(l,P)} Z_lP,   n(l,P) = nint( h^-1 (r_l - R_P) )
!>        The phase depends only on the grid point l and the RI atom P, so it is constant across
!>        the RI functions of P. It is required: dropping it makes the result worse and of
!>        opposite sign (23.423 vs 23.601 against a reference of 23.543), the two bracketing the
!>        reference, which confirms that the local fit coefficient of grid point l carries the
!>        image phase of atom P.
!> \param mat_Z Gamma coefficients Z_lP(k=0) (grid x RI, RI column blocks = atoms)
!> \param Zc out, Re Z_lP(k)
!> \param Zs out, Im Z_lP(k)
!> \param grid_points ...
!> \param particle_set ...
!> \param row_off_parent cumulative grid row offsets
!> \param h_inv ...
!> \param xkp_k ...
! **************************************************************************************************
   SUBROUTINE build_Z_ikp_grid(mat_Z, Zc, Zs, grid_points, particle_set, row_off_parent, &
                               h_inv, xkp_k)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_Z, Zc, Zs
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_off_parent
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_Z_ikp_grid'

      INTEGER                                            :: col, handle, i, j, l, row
      INTEGER, DIMENSION(3)                              :: n
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: arg, c, s, v
      REAL(KIND=dp), DIMENSION(3)                        :: R_P, raw
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_c, blk_s
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(Zc, mat_Z)
      CALL dbcsr_copy(Zs, mat_Z)

      CALL dbcsr_iterator_start(iter, Zc)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk_c)
         CALL dbcsr_get_block_p(Zs, row, col, blk_s, found)
         IF (.NOT. found) CYCLE

         R_P(:) = particle_set(col)%r(:)

         DO i = 1, SIZE(blk_c, 1)
            l = row_off_parent(row) + i
            raw(:) = grid_points(:, l) - R_P(:)
            n(:) = NINT(MATMUL(h_inv, raw))
            arg = twopi*(xkp_k(1)*n(1) + xkp_k(2)*n(2) + xkp_k(3)*n(3))
            c = COS(arg)
            s = SIN(arg)
            DO j = 1, SIZE(blk_c, 2)
               v = blk_c(i, j)
               blk_c(i, j) = v*c
               blk_s(i, j) = v*s
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE build_Z_ikp_grid

! **************************************************************************************************
!> \brief Multiplies W_PQ(i omega, k) by the basis-independent Brillouin-zone weight of k.
!>
!>        This is the two-grid extrapolation weight: positive on the dense mesh, negative on the
!>        coarse one, so that summing over all k realizes
!>          W^R = sum_{q in dense} v e^{iqR} W(q) - sum_{k in coarse} w e^{ikR} W(k)
!>        with v, w depending only on the two mesh sizes. It deliberately does NOT use the
!>        angular-momentum-resolved wkp / wkp_s_p / wkp_no_extra selection of the auxiliary-basis
!>        path: after the contraction to the grid there is no PQ angular momentum left, and the
!>        per-channel selection breaks the dense/coarse cancellation under the grid phase, which
!>        showed up as degeneracy-breaking noise in Sigma^c.
!> \param cfm_W ...
!> \param bs_env ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE scale_W_PQ_by_wkp(cfm_W, bs_env, ikp)

      TYPE(cp_cfm_type), INTENT(INOUT)                   :: cfm_W
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_W_PQ_by_wkp'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      cfm_W%local_data(:, :) = cfm_W%local_data(:, :)*bs_env%kpoints_chi_eps_W%wkp(ikp)

      CALL timestop(handle)

   END SUBROUTINE scale_W_PQ_by_wkp

! **************************************************************************************************
!> \brief Screened interaction of the Gamma-only periodic RI-RS scheme, truncated to the minimum
!>        image in the GRID basis rather than in the RI auxiliary basis.
!>
!>        Per k-point:
!>          chi_ll'(i tau, k) = e^{i 2 pi k . n(l,l')} chi_ll'(i tau, k=0)
!>          chi_PQ(i tau, k)  = sum_ll' Z*_lP(k) chi_ll'(i tau, k) Z_l'Q(k)
!>          chi_PQ(i w, k)    = sum_tau w_cos cos(w tau) chi_PQ(i tau, k),  projected >= 0
!>          eps, W_PQ(i w, k) = the usual dielectric construction (reused PQ routines)
!>          W~_PQ(i w, k)     = wkp(k) W_PQ(i w, k)
!>          W_ll'(i w, k)     = sum_PQ Zt_lP(k) W~_PQ(i w, k) Zt*_l'Q(k),  Zt = Z M^-1(0)
!>          W^MIC_ll'(i tau) += sum_k sum_w w_cos cos(w tau) Re[ e^{-i 2 pi k . n(l,l')} W_ll'(i w,k) ]
!>
!>        Two properties make this affordable, and both are required:
!>
!>        1. Everything grid-sized is done PANEL BY PANEL. The reference implementation kept four
!>           full grid x grid scratch matrices (chi_re, chi_im, Wgr, Wgi); here each is a panel of
!>           grid row blocks, extracted locally and discarded.
!>        2. W is stored on mat_W_pattern = supp(G^occ) u supp(G^vir), built alongside chi. This
!>           is exact because W only ever reaches the self-energy through Sigma^c = G o W, and it
!>           is what keeps the storage linear in system size; the a_min/2 minimum-image sphere
!>           would be ~52 % dense. Real and imaginary parts share that one pattern, which also
!>           removes a correctness trap: the inverse transform pairs the two through a local
!>           block lookup, so divergent patterns would silently skip blocks of the real part.
!>
!>        The nearest image n(l,l') is kept unconditionally (no a_min/2 mask): the hard mask
!>        over-truncates, discarding pairs between a_min/2 and the Wigner-Seitz boundary that the
!>        reference retains. See grid_W_from_chi_equations.md.
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_grid_tau  consumed and deallocated here
!> \param mat_W_pattern     block pattern for W; consumed and released here
!> \param mat_W_grid_time   out, W^MIC_ll'(i tau)
! **************************************************************************************************
   SUBROUTINE compute_W_MIC_grid(bs_env, qs_env, mat_chi_grid_tau, mat_W_pattern, mat_W_grid_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_grid_tau, mat_W_grid_time
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_W_pattern

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_W_MIC_grid'

      INTEGER :: blk0, blk1, handle, i_blk, i_t, ikp, ikp_batch, ikp_in_batch, ipan, j_w, &
         n_grid_blks, n_panels, n_RI, ndep
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: id_map, pan_first, pan_last, row_off
      INTEGER, DIMENSION(:), POINTER                     :: r_blk_sizes
      REAL(KIND=dp)                                      :: coeff, eps, freq_j, t1, time_i
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ctab, stab
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: n_idx_pan, Z_pan_c_k, Z_pan_s_k, &
                                                            Zt_pan_c_k, Zt_pan_s_k
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type)                                  :: cfm_chi, cfm_M_inv_V_sqrt_ikp, &
                                                            cfm_V_sqrt_ikp, cfm_W_ikp_freq_j, &
                                                            cfm_W_unscaled
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_chi_PQ_tau
      TYPE(cp_fm_type)                                   :: fm_im, fm_Minv0_work, fm_re
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_freq_all
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Minv0, fm_V_kp
      TYPE(dbcsr_type) :: chi_pan_im, chi_pan_re, M_im, M_re, mat_Minv0, mat_Z_tilde, PQ_im, &
         PQ_re, T_im, T_re, W_tmpl_pan, Wg_im, Wg_re, Wt_im, Wt_re, Z_pan_c, Z_pan_s, Zc, Zs, &
         Zt_pan_c, Zt_pan_s, Ztc, Zts
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      eps = bs_env%eps_filter
      n_RI = bs_env%n_RI

      IF (bs_env%approx_kp_extrapol .AND. bs_env%unit_nr > 0) THEN
         CPWARN("APPROX_KP_EXTRAPOL is ignored by W_TRUNCATION GRID.")
      END IF

      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat, h_inv=h_inv)

      CALL resolve_grid_panels(bs_env, bs_env%ri_rs%mat_phi_mu_l, pan_first, pan_last)
      n_panels = SIZE(pan_first)

      ! cumulative grid row offsets, and the identity column map used by the phase routine
      ! (panels taken with extract_grid_panel keep the parent's full column set)
      CALL dbcsr_get_info(mat_W_pattern, row_blk_size=r_blk_sizes)
      n_grid_blks = SIZE(r_blk_sizes)
      ALLOCATE (row_off(n_grid_blks), id_map(n_grid_blks))
      row_off(1) = 0
      DO i_blk = 2, n_grid_blks
         row_off(i_blk) = row_off(i_blk - 1) + r_blk_sizes(i_blk - 1)
      END DO
      DO i_blk = 1, n_grid_blks
         id_map(i_blk) = i_blk
      END DO

      ! W^MIC_ll'(i tau): one matrix per time point, all carrying the supp(G) pattern so that the
      ! per-(k, omega) accumulation below is a purely local block add.
      NULLIFY (mat_W_grid_time)
      ALLOCATE (mat_W_grid_time(bs_env%num_time_freq_points))
      DO i_t = 1, bs_env%num_time_freq_points
         ALLOCATE (mat_W_grid_time(i_t)%matrix)
         CALL dbcsr_create(mat_W_grid_time(i_t)%matrix, template=mat_W_pattern)
         CALL dbcsr_copy(mat_W_grid_time(i_t)%matrix, mat_W_pattern)
         CALL dbcsr_set(mat_W_grid_time(i_t)%matrix, 0.0_dp)
      END DO

      ! reusable scratch (RI-sized and grid x RI; nothing grid x grid)
      CALL dbcsr_create(Zc, template=bs_env%ri_rs%mat_Z_lP)
      CALL dbcsr_create(Zs, template=bs_env%ri_rs%mat_Z_lP)
      CALL dbcsr_create(Ztc, template=bs_env%ri_rs%mat_Z_lP)
      CALL dbcsr_create(Zts, template=bs_env%ri_rs%mat_Z_lP)
      CALL dbcsr_create(PQ_re, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_create(PQ_im, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_create(Wt_re, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_create(Wt_im, template=bs_env%mat_RI_RI%matrix)
      CALL cp_fm_create(fm_re, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_im, bs_env%fm_RI_RI%matrix_struct)

      ALLOCATE (cfm_chi_PQ_tau(bs_env%num_time_freq_points))
      DO i_t = 1, bs_env%num_time_freq_points
         CALL cp_cfm_create(cfm_chi_PQ_tau(i_t), bs_env%fm_RI_RI%matrix_struct)
      END DO

      ! RT-BSE needs W_PQ at zero frequency, which the auxiliary-basis path produces as a
      ! by-product. Accumulate the auxiliary-basis W^MIC(i omega_j) on the side (standard
      ! atom-pair MIC and standard k-weights, exactly as the AUX path does) and combine the
      ! frequency points at the end.
      IF (bs_env%rtp_method == rtp_method_bse) THEN
         ALLOCATE (fm_W_MIC_freq_all(bs_env%num_time_freq_points))
         DO j_w = 1, bs_env%num_time_freq_points
            CALL cp_fm_create(fm_W_MIC_freq_all(j_w), bs_env%fm_RI_RI%matrix_struct)
            CALL cp_fm_set_all(fm_W_MIC_freq_all(j_w), 0.0_dp)
         END DO
      END IF

      ! Outer RI-metric dressing M^-1(k=0) folded once into the back-transformation coefficients
      ! Zt = Z M^-1(0). This reproduces the reference dressing M^-1(0) W^MIC M^-1(0) while
      ! avoiding any per-k inversion of M(k); applying M^-1(k) per k instead was the source of an
      ! ~8x magnitude error in an earlier version.
      CALL cp_fm_create(fm_Minv0_work, bs_env%fm_RI_RI%matrix_struct)
      CALL RI_2c_integral_mat(qs_env, fm_Minv0, bs_env%fm_RI_RI, n_RI, bs_env%ri_metric, &
                              do_kpoints=.FALSE.)
      CALL cp_fm_power(fm_Minv0(1, 1), fm_Minv0_work, -1.0_dp, 0.0_dp, ndep)
      CALL dbcsr_create(mat_Minv0, template=bs_env%mat_RI_RI%matrix)
      CALL copy_fm_to_dbcsr(fm_Minv0(1, 1), mat_Minv0, keep_sparsity=.FALSE.)
      CALL dbcsr_create(mat_Z_tilde, template=bs_env%ri_rs%mat_Z_lP)
      CALL dbcsr_multiply("N", "N", 1.0_dp, bs_env%ri_rs%mat_Z_lP, mat_Minv0, 0.0_dp, mat_Z_tilde)
      CALL cp_fm_release(fm_Minv0)
      CALL cp_fm_release(fm_Minv0_work)
      CALL dbcsr_release(mat_Minv0)

      ! Nearest-image index per panel, evaluated once: it depends only on the geometry, not on
      ! k, tau or omega. It doubles as the pattern template for the W panels below.
      ALLOCATE (n_idx_pan(n_panels))
      ALLOCATE (Z_pan_c_k(n_panels), Z_pan_s_k(n_panels))
      ALLOCATE (Zt_pan_c_k(n_panels), Zt_pan_s_k(n_panels))
      DO ipan = 1, n_panels
         CALL extract_grid_panel(mat_W_pattern, pan_first(ipan), pan_last(ipan), W_tmpl_pan)
         CALL build_grid_phase_index(W_tmpl_pan, pan_first(ipan), bs_env%ri_rs%grid_points, &
                                     row_off, h_inv, n_idx_pan(ipan))
         CALL dbcsr_release(W_tmpl_pan)
      END DO

      DO ikp_batch = 1, bs_env%num_chi_eps_W_batches

         t1 = m_walltime()

         CALL compute_V_k_by_lattice_sum(bs_env, qs_env, fm_V_kp, ikp_batch)

         DO ikp_in_batch = 1, bs_env%nkp_chi_eps_W_batch

            ikp = (ikp_batch - 1)*bs_env%nkp_chi_eps_W_batch + ikp_in_batch
            IF (ikp > bs_env%nkp_chi_eps_W_orig_plus_extra) CYCLE

            CALL compute_MinvVsqrt_Vsqrt(bs_env, qs_env, fm_V_kp, &
                                         cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp, ikp)

            CALL bs_env%para_env%sync()
            CALL cp_fm_release(fm_V_kp(ikp, 1))
            CALL cp_fm_release(fm_V_kp(ikp, 2))

            CALL build_Z_ikp_grid(bs_env%ri_rs%mat_Z_lP, Zc, Zs, bs_env%ri_rs%grid_points, &
                                  particle_set, row_off, h_inv, &
                                  bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))
            CALL build_Z_ikp_grid(mat_Z_tilde, Ztc, Zts, bs_env%ri_rs%grid_points, &
                                  particle_set, row_off, h_inv, &
                                  bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))

            ! cos/sin of the Bloch phase for every nearest-image triple, and the Z / Ztilde
            ! panels: all depend on k but not on tau or omega, so they are hoisted out of both
            ! inner loops.
            CALL build_phase_table(bs_env%kpoints_chi_eps_W%xkp(1:3, ikp), ctab, stab)
            DO ipan = 1, n_panels
               CALL extract_grid_panel(Zc, pan_first(ipan), pan_last(ipan), Z_pan_c_k(ipan))
               CALL extract_grid_panel(Zs, pan_first(ipan), pan_last(ipan), Z_pan_s_k(ipan))
               CALL extract_grid_panel(Ztc, pan_first(ipan), pan_last(ipan), Zt_pan_c_k(ipan))
               CALL extract_grid_panel(Zts, pan_first(ipan), pan_last(ipan), Zt_pan_s_k(ipan))
            END DO

            ! ---------------------------------------------------------------------------
            ! chi_PQ(i tau, k) = Z(k)^dagger chi(k) Z(k), accumulated over grid panels
            ! ---------------------------------------------------------------------------
            DO i_t = 1, bs_env%num_time_freq_points

               CALL dbcsr_set(PQ_re, 0.0_dp)
               CALL dbcsr_set(PQ_im, 0.0_dp)

               DO ipan = 1, n_panels
                  blk0 = pan_first(ipan)
                  blk1 = pan_last(ipan)

                  ! chi panel and its Bloch phase (local slice, zero communication)
                  CALL extract_grid_panel(mat_chi_grid_tau(i_t)%matrix, blk0, blk1, chi_pan_re)
                  ! chi(k=0) is real: copy only to inherit the block pattern, then zero the
                  ! values. apply_grid_mic_phase treats (re, im) as a genuine complex pair, so a
                  ! non-zero imaginary part here would be picked up as signal.
                  CALL dbcsr_create(chi_pan_im, template=chi_pan_re)
                  CALL dbcsr_copy(chi_pan_im, chi_pan_re)
                  CALL dbcsr_set(chi_pan_im, 0.0_dp)
                  CALL apply_grid_mic_phase(chi_pan_re, chi_pan_im, n_idx_pan(ipan), &
                                            ctab, stab, .FALSE.)

                  ! M = chi_pan(k) Z(k)   [panel x RI]
                  CALL create_product_matrix(chi_pan_re, Zc, 'N', 'N', M_re)
                  CALL create_product_matrix(chi_pan_re, Zc, 'N', 'N', M_im)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, chi_pan_re, Zc, 0.0_dp, M_re, filter_eps=eps)
                  CALL dbcsr_multiply("N", "N", -1.0_dp, chi_pan_im, Zs, 1.0_dp, M_re, filter_eps=eps)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, chi_pan_re, Zs, 0.0_dp, M_im, filter_eps=eps)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, chi_pan_im, Zc, 1.0_dp, M_im, filter_eps=eps)

                  ! chi_PQ += Z_pan(k)^dagger M
                  CALL dbcsr_multiply("T", "N", 1.0_dp, Z_pan_c_k(ipan), M_re, 1.0_dp, PQ_re, &
                                      filter_eps=eps)
                  CALL dbcsr_multiply("T", "N", 1.0_dp, Z_pan_s_k(ipan), M_im, 1.0_dp, PQ_re, &
                                      filter_eps=eps)
                  CALL dbcsr_multiply("T", "N", 1.0_dp, Z_pan_c_k(ipan), M_im, 1.0_dp, PQ_im, &
                                      filter_eps=eps)
                  CALL dbcsr_multiply("T", "N", -1.0_dp, Z_pan_s_k(ipan), M_re, 1.0_dp, PQ_im, &
                                      filter_eps=eps)

                  CALL dbcsr_release(M_re)
                  CALL dbcsr_release(M_im)
                  CALL dbcsr_release(chi_pan_re)
                  CALL dbcsr_release(chi_pan_im)
               END DO ! ipan

               CALL copy_dbcsr_to_fm(PQ_re, fm_re)
               CALL copy_dbcsr_to_fm(PQ_im, fm_im)
               CALL cp_fm_to_cfm(fm_re, fm_im, cfm_chi_PQ_tau(i_t))
            END DO ! i_t

            ! ---------------------------------------------------------------------------
            ! per frequency: tau->omega, eps and W_PQ(i w, k), back to the grid, MIC, w->tau
            ! ---------------------------------------------------------------------------
            DO j_w = 1, bs_env%num_time_freq_points

               freq_j = bs_env%imag_freq_points(j_w)

               CALL cp_cfm_create(cfm_chi, bs_env%fm_RI_RI%matrix_struct)
               CALL cp_cfm_set_all(cfm_chi, z_zero)
               DO i_t = 1, bs_env%num_time_freq_points
                  time_i = bs_env%imag_time_points(i_t)
                  coeff = bs_env%weights_cos_t_to_w(j_w, i_t)*COS(time_i*freq_j)
                  CALL cp_cfm_scale_and_add(z_one, cfm_chi, &
                                            CMPLX(coeff, 0.0_dp, KIND=dp), cfm_chi_PQ_tau(i_t))
               END DO

               ! project out the negative spectrum of chi
               CALL cp_cfm_power(cfm_chi, threshold=0.0_dp, exponent=1.0_dp)

               ! eps(i w, k), W_PQ(i w, k); releases cfm_chi internally
               CALL compute_cfm_W_ikp_freq_j(bs_env, cfm_chi, cfm_V_sqrt_ikp, &
                                             cfm_M_inv_V_sqrt_ikp, cfm_W_ikp_freq_j)

               ! auxiliary-basis W(i w) on the side, for the RT-BSE kernel (unweighted here:
               ! MIC_contribution_from_ikp applies the standard k-weights itself)
               IF (bs_env%rtp_method == rtp_method_bse) THEN
                  CALL cp_cfm_create(cfm_W_unscaled, bs_env%fm_RI_RI%matrix_struct)
                  CALL cp_cfm_to_cfm(cfm_W_ikp_freq_j, cfm_W_unscaled)
                  CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_all(j_w), &
                                                 cfm_W_unscaled, ikp, &
                                                 bs_env%kpoints_chi_eps_W, "RI_AUX")
                  CALL cp_cfm_release(cfm_W_unscaled)
               END IF

               CALL scale_W_PQ_by_wkp(cfm_W_ikp_freq_j, bs_env, ikp)

               CALL cp_cfm_to_fm(cfm_W_ikp_freq_j, fm_re, fm_im)
               CALL copy_fm_to_dbcsr(fm_re, Wt_re, keep_sparsity=.FALSE.)
               CALL copy_fm_to_dbcsr(fm_im, Wt_im, keep_sparsity=.FALSE.)
               CALL cp_cfm_release(cfm_W_ikp_freq_j)

               DO ipan = 1, n_panels
                  blk0 = pan_first(ipan)
                  blk1 = pan_last(ipan)

                  ! T = Zt_pan(k) W~_PQ   [panel x RI]
                  CALL create_product_matrix(Zt_pan_c_k(ipan), Wt_re, 'N', 'N', T_re)
                  CALL create_product_matrix(Zt_pan_c_k(ipan), Wt_re, 'N', 'N', T_im)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, Zt_pan_c_k(ipan), Wt_re, 0.0_dp, T_re, &
                                      filter_eps=eps)
                  CALL dbcsr_multiply("N", "N", -1.0_dp, Zt_pan_s_k(ipan), Wt_im, 1.0_dp, T_re, &
                                      filter_eps=eps)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, Zt_pan_c_k(ipan), Wt_im, 0.0_dp, T_im, &
                                      filter_eps=eps)
                  CALL dbcsr_multiply("N", "N", 1.0_dp, Zt_pan_s_k(ipan), Wt_re, 1.0_dp, T_im, &
                                      filter_eps=eps)

                  ! W_pan(i w, k) = T Zt(k)^dagger, confined to the supp(G) pattern. The panel of
                  ! mat_W_pattern is co-resident with the accumulator by construction, so the
                  ! phase pairing and the scatter below are both purely local.
                  CALL dbcsr_create(Wg_re, template=n_idx_pan(ipan))
                  CALL dbcsr_copy(Wg_re, n_idx_pan(ipan))
                  CALL dbcsr_set(Wg_re, 0.0_dp)
                  CALL dbcsr_create(Wg_im, template=n_idx_pan(ipan))
                  CALL dbcsr_copy(Wg_im, n_idx_pan(ipan))
                  CALL dbcsr_set(Wg_im, 0.0_dp)

                  CALL dbcsr_multiply("N", "T", 1.0_dp, T_re, Ztc, 0.0_dp, Wg_re, &
                                      retain_sparsity=.TRUE.)
                  CALL dbcsr_multiply("N", "T", 1.0_dp, T_im, Zts, 1.0_dp, Wg_re, &
                                      retain_sparsity=.TRUE.)
                  CALL dbcsr_multiply("N", "T", 1.0_dp, T_im, Ztc, 0.0_dp, Wg_im, &
                                      retain_sparsity=.TRUE.)
                  CALL dbcsr_multiply("N", "T", -1.0_dp, T_re, Zts, 1.0_dp, Wg_im, &
                                      retain_sparsity=.TRUE.)

                  ! Re[ e^{-i 2 pi k . n(l,l')} W_ll'(i w, k) ], in place in Wg_re
                  CALL apply_grid_mic_phase(Wg_re, Wg_im, n_idx_pan(ipan), ctab, stab, .TRUE.)

                  DO i_t = 1, bs_env%num_time_freq_points
                     time_i = bs_env%imag_time_points(i_t)
                     coeff = bs_env%weights_cos_w_to_t(i_t, j_w)*COS(time_i*freq_j)
                     CALL accumulate_grid_panel(Wg_re, blk0, mat_W_grid_time(i_t)%matrix, &
                                                coeff, id_map, retain_sparsity=.TRUE.)
                  END DO

                  CALL dbcsr_release(Wg_re)
                  CALL dbcsr_release(Wg_im)
                  CALL dbcsr_release(T_re)
                  CALL dbcsr_release(T_im)
               END DO ! ipan

            END DO ! j_w

            DO ipan = 1, n_panels
               CALL dbcsr_release(Z_pan_c_k(ipan))
               CALL dbcsr_release(Z_pan_s_k(ipan))
               CALL dbcsr_release(Zt_pan_c_k(ipan))
               CALL dbcsr_release(Zt_pan_s_k(ipan))
            END DO

         END DO ! ikp_in_batch

         DEALLOCATE (fm_V_kp)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I8,A,I3,A,F10.1,A)') &
               'Computed grid W^MIC(iτ) for k-point batch', &
               ikp_batch, ' /', bs_env%num_chi_eps_W_batches, &
               ', Execution time', m_walltime() - t1, ' s'
            CALL m_flush(bs_env%unit_nr)
         END IF

      END DO ! ikp_batch

      CALL print_matrix_occupation(mat_W_grid_time(1)%matrix, 'W^MIC(l,l'')', &
                                   bs_env%para_env, bs_env%unit_nr)

      ! W_PQ(i omega = 0) for RT-BSE: combine the frequency points with the same weights the
      ! auxiliary-basis path uses, W(0) = sum_tau v(tau) sum_w w_cos(tau,w) cos(w tau) W(i w).
      IF (bs_env%rtp_method == rtp_method_bse) THEN
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_RI_RI%matrix_struct)
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_zero, 0.0_dp)
         DO j_w = 1, bs_env%num_time_freq_points
            freq_j = bs_env%imag_freq_points(j_w)
            coeff = 0.0_dp
            DO i_t = 1, bs_env%num_time_freq_points
               time_i = bs_env%imag_time_points(i_t)
               coeff = coeff + bs_env%imag_time_weights_freq_zero(i_t)* &
                       bs_env%weights_cos_w_to_t(i_t, j_w)*COS(time_i*freq_j)
            END DO
            CALL cp_fm_scale_and_add(1.0_dp, bs_env%fm_W_MIC_freq_zero, &
                                     coeff, fm_W_MIC_freq_all(j_w))
         END DO
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
         CALL cp_fm_release(fm_W_MIC_freq_all)
      END IF

      CALL cp_cfm_release(cfm_M_inv_V_sqrt_ikp)
      CALL cp_cfm_release(cfm_V_sqrt_ikp)
      DO i_t = 1, bs_env%num_time_freq_points
         CALL cp_cfm_release(cfm_chi_PQ_tau(i_t))
      END DO
      DEALLOCATE (cfm_chi_PQ_tau)
      CALL cp_fm_release(fm_re)
      CALL cp_fm_release(fm_im)
      CALL dbcsr_release(Zc)
      CALL dbcsr_release(Zs)
      CALL dbcsr_release(Ztc)
      CALL dbcsr_release(Zts)
      CALL dbcsr_release(PQ_re)
      CALL dbcsr_release(PQ_im)
      CALL dbcsr_release(Wt_re)
      CALL dbcsr_release(Wt_im)
      CALL dbcsr_release(mat_Z_tilde)
      CALL dbcsr_release(mat_W_pattern)
      DO ipan = 1, n_panels
         CALL dbcsr_release(n_idx_pan(ipan))
      END DO
      DEALLOCATE (n_idx_pan, Z_pan_c_k, Z_pan_s_k, Zt_pan_c_k, Zt_pan_s_k)
      IF (ALLOCATED(ctab)) DEALLOCATE (ctab, stab)
      CALL dbcsr_deallocate_matrix_set(mat_chi_grid_tau)
      DEALLOCATE (row_off, id_map)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W_MIC_grid

! **************************************************************************************************
!> \brief Correlation self-energy when W already lives in the grid basis:
!>          Sigma^c_neg = -phi^T ( (phi G^occ phi^T) o W_grid ) phi
!>          Sigma^c_pos = +phi^T ( (phi G^vir phi^T) o W_grid ) phi
!>
!>        The auxiliary-basis analogue (contract_grid_panels_sigma_c) has to build
!>        W_pan = Z_pan W_PQ Z^T with two Cannon multiplies per panel. Here W arrives already in
!>        the grid basis, so its panel is a LOCAL slice and costs no communication at all.
!>
!>        The seeding direction is deliberately inverted with respect to the auxiliary-basis
!>        routine: W's pattern defines the neighborhood, and A_occ / A_vir are created as
!>        dbcsr_create(template=W_pan) and filled with retain_sparsity. That makes the Hadamard
!>        operands co-resident and identically patterned by construction, so the local-only block
!>        lookup inside hadamard_product can never silently drop a block. It also makes this
!>        routine independent of the geometric screening entirely: the truncation is inherited
!>        from W's pattern, which was built once, correctly, in the chi stage.
!>
!>        retain_sparsity on the A_occ / A_vir multiplies is exact for the usual reason: any block
!>        of G outside supp(W) meets a zero in the Hadamard product.
!> \param bs_env ...
!> \param mat_W_grid_tau W^MIC_ll'(i tau); consumed and deallocated here
!> \param mat_phi_mu_l ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_grid(bs_env, mat_W_grid_tau, mat_phi_mu_l, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_W_grid_tau
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_grid'

      INTEGER                                            :: blk0, blk1, handle, i_t, ipan, ispin, &
                                                            n_panels
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: gmap, pan_first, pan_last
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: grid_used, used_ao
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbcsr_type) :: A_occ, A_vir, C_pan, G_occ_sub, G_vir_sub, matrix_G_occ_ao, &
         matrix_G_vir_ao, phi_pan, phi_panC, tmp2, tmpA, W_pan, W_pan_full
      TYPE(dbcsr_type), TARGET                           :: phi_near

      CALL timeset(routineN, handle)

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)
      n_panels = SIZE(pan_first)

      NULLIFY (mat_Sigma_neg_tau, mat_Sigma_pos_tau)
      ALLOCATE (mat_Sigma_neg_tau(bs_env%num_time_freq_points, bs_env%n_spin))
      ALLOCATE (mat_Sigma_pos_tau(bs_env%num_time_freq_points, bs_env%n_spin))
      DO i_t = 1, bs_env%num_time_freq_points
         DO ispin = 1, bs_env%n_spin
            ALLOCATE (mat_Sigma_neg_tau(i_t, ispin)%matrix)
            ALLOCATE (mat_Sigma_pos_tau(i_t, ispin)%matrix)
            CALL dbcsr_create(mat_Sigma_neg_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
            CALL dbcsr_create(mat_Sigma_pos_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
         END DO
      END DO

      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%imag_time_points(i_t)

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_ao)
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_ao)

            CALL dbcsr_set(mat_Sigma_neg_tau(i_t, ispin)%matrix, 0.0_dp)
            CALL dbcsr_set(mat_Sigma_pos_tau(i_t, ispin)%matrix, 0.0_dp)

            DO ipan = 1, n_panels
               blk0 = pan_first(ipan)
               blk1 = pan_last(ipan)

               CALL extract_grid_panel(mat_phi_mu_l, blk0, blk1, phi_pan)
               CALL collect_used_col_blocks(phi_pan, bs_env%para_env, used_ao)
               IF (.NOT. ANY(used_ao)) THEN
                  CALL dbcsr_release(phi_pan)
                  CYCLE
               END IF

               ! W panel: local slice, then compress to the columns it actually occupies. That
               ! column set IS the neighborhood, so no centroid screening is needed here.
               CALL extract_grid_panel(mat_W_grid_tau(i_t)%matrix, blk0, blk1, W_pan_full)
               CALL collect_used_col_blocks(W_pan_full, bs_env%para_env, grid_used)
               IF (.NOT. ANY(grid_used)) THEN
                  CALL dbcsr_release(W_pan_full)
                  CALL dbcsr_release(phi_pan)
                  CYCLE
               END IF
               CALL extract_masked_blocks(W_pan_full, grid_used, W_pan, &
                                          compress_rows=.FALSE., blk_map=gmap)
               CALL dbcsr_release(W_pan_full)
               CALL extract_masked_blocks(mat_phi_mu_l, grid_used, phi_near, compress_rows=.TRUE.)

               CALL extract_masked_blocks(phi_pan, used_ao, phi_panC, compress_rows=.FALSE.)
               CALL extract_masked_blocks(matrix_G_occ_ao, used_ao, G_occ_sub, compress_rows=.TRUE.)
               CALL extract_masked_blocks(matrix_G_vir_ao, used_ao, G_vir_sub, compress_rows=.TRUE.)

               ! A_occ inherits W_pan's pattern AND distribution
               CALL dbcsr_create(A_occ, template=W_pan)
               CALL dbcsr_copy(A_occ, W_pan)
               CALL dbcsr_set(A_occ, 0.0_dp)
               CALL create_product_matrix(phi_panC, G_occ_sub, 'N', 'N', tmpA)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_occ_sub, 0.0_dp, tmpA, &
                                   filter_eps=bs_env%eps_filter)
               CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, phi_near, 0.0_dp, A_occ, &
                                   filter_eps=bs_env%eps_filter, retain_sparsity=.TRUE.)
               CALL dbcsr_release(tmpA)
               CALL dbcsr_release(G_occ_sub)

               CALL dbcsr_create(A_vir, template=W_pan)
               CALL dbcsr_copy(A_vir, W_pan)
               CALL dbcsr_set(A_vir, 0.0_dp)
               CALL create_product_matrix(phi_panC, G_vir_sub, 'N', 'N', tmpA)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, phi_panC, G_vir_sub, 0.0_dp, tmpA, &
                                   filter_eps=bs_env%eps_filter)
               CALL dbcsr_multiply('N', 'T', 1.0_dp, tmpA, phi_near, 0.0_dp, A_vir, &
                                   filter_eps=bs_env%eps_filter, retain_sparsity=.TRUE.)
               CALL dbcsr_release(tmpA)
               CALL dbcsr_release(G_vir_sub)

               ! Sigma^c_neg += phi_pan^T ( A_occ o W_pan ) phi
               CALL dbcsr_create(C_pan, template=A_occ)
               CALL hadamard_product(A_occ, W_pan, C_pan, 1.0_dp)
               CALL dbcsr_release(A_occ)
               CALL create_product_matrix(C_pan, phi_near, 'N', 'N', tmp2)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, phi_near, 0.0_dp, tmp2, &
                                   filter_eps=bs_env%eps_filter)
               CALL dbcsr_release(C_pan)
               CALL dbcsr_multiply('T', 'N', 1.0_dp, phi_pan, tmp2, 1.0_dp, &
                                   mat_Sigma_neg_tau(i_t, ispin)%matrix, filter_eps=bs_env%eps_filter)
               CALL dbcsr_release(tmp2)

               ! Sigma^c_pos += phi_pan^T ( A_vir o W_pan ) phi
               CALL dbcsr_create(C_pan, template=A_vir)
               CALL hadamard_product(A_vir, W_pan, C_pan, 1.0_dp)
               CALL dbcsr_release(A_vir)
               CALL create_product_matrix(C_pan, phi_near, 'N', 'N', tmp2)
               CALL dbcsr_multiply('N', 'N', 1.0_dp, C_pan, phi_near, 0.0_dp, tmp2, &
                                   filter_eps=bs_env%eps_filter)
               CALL dbcsr_release(C_pan)
               CALL dbcsr_multiply('T', 'N', 1.0_dp, phi_pan, tmp2, 1.0_dp, &
                                   mat_Sigma_pos_tau(i_t, ispin)%matrix, filter_eps=bs_env%eps_filter)
               CALL dbcsr_release(tmp2)

               CALL dbcsr_release(W_pan)
               CALL dbcsr_release(phi_near)
               CALL dbcsr_release(phi_panC)
               CALL dbcsr_release(phi_pan)
            END DO ! ipan

            CALL dbcsr_scale(mat_Sigma_neg_tau(i_t, ispin)%matrix, -1.0_dp)

            CALL dbcsr_release(matrix_G_occ_ao)
            CALL dbcsr_release(matrix_G_vir_ao)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I15,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ) for time point', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
               CALL m_flush(bs_env%unit_nr)
            END IF

         END DO ! ispin
      END DO ! i_t

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)
      ! W in the grid basis is the largest object this path holds; release it here, where it
      ! stops being needed. (Unlike the auxiliary-basis Sigma^c this routine is entered once --
      ! evGW0 is not supported on the grid path -- so there is nothing to preserve for a
      ! subsequent cycle.)
      CALL dbcsr_deallocate_matrix_set(mat_W_grid_tau)
      IF (ALLOCATED(gmap)) DEALLOCATE (gmap)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_grid

! **************************************************************************************************

   SUBROUTINE contract_A_B_A(transA_left, transA_right, matrix_A, matrix_B, matrix_out, eps_filter)

      CHARACTER(LEN=1), INTENT(IN)                       :: transA_left, transA_right
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_out
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_A_B_A'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix_tmp, template=matrix_A)

      IF (transA_left == "N" .AND. transA_right == "T") THEN
         ! Path 1: Out = A * B * A^T
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_A, matrix_B, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, matrix_A, &
                             0.0_dp, matrix_out, filter_eps=eps_filter)

      ELSE IF (transA_left == "T" .AND. transA_right == "N") THEN
         ! Path 2: Out = A^T * B * A
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_B, matrix_A, &
                             0.0_dp, matrix_tmp, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_A, matrix_tmp, &
                             0.0_dp, matrix_out, filter_eps=eps_filter)
      ELSE
         CPABORT("Unsupported transposition pair in contract_A_B_A")
      END IF

      CALL dbcsr_release(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE contract_A_B_A

! **************************************************************************************************
!> \brief Computes C = A ◦ B (Element-wise Hadamard product) for sparse DBCSR matrices.
!> \param matrix_A ...
!> \param matrix_B ...
!> \param matrix_C ...
!> \param fac (Scaling factor applied to the product)
! **************************************************************************************************
!> \brief In-place Hadamard A <- fac * (A ◦ B). Value mutation only (no block insert/delete),
!>        so iterating A while writing through the block pointer is safe.
!> \param matrix_A in/out factor (overwritten by the product)
!> \param matrix_B second factor (looked up; blocks absent in B zero the A block)
!> \param fac (Scaling factor applied to the product)
! **************************************************************************************************
   SUBROUTINE hadamard_product_inplace(matrix_A, matrix_B, fac)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B
      REAL(KIND=dp), INTENT(IN)                          :: fac

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hadamard_product_inplace'

      INTEGER                                            :: col, handle, row
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk_A, blk_B
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, matrix_A)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk_A)

         CALL dbcsr_get_block_p(matrix_B, row, col, blk_B, found)

         IF (found) THEN
            blk_A(:, :) = fac*blk_A(:, :)*blk_B(:, :)
         ELSE
            ! If B is sparse here, the product is zero
            blk_A(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE hadamard_product_inplace

! **************************************************************************************************
!> \brief Compute screened Coulomb interaction matrix
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_time ...
! **************************************************************************************************

   SUBROUTINE compute_W(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W'

      INTEGER                                            :: handle, i_t, j_w
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type)                                   :: fm_M_inv_V_sqrt, fm_V, fm_V_sqrt

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL create_fm_W_MIC_time(bs_env, fm_W_time)

      ! 1. Allocate V and M matrices
      CALL cp_fm_create(fm_V, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_M_inv_V_sqrt, bs_env%fm_RI_RI%matrix_struct)

      ! Compute V and M^-1 * V^0.5
      CALL compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_M_inv_V_sqrt)

      ! 2. Loop over frequencies
      DO j_w = 1, bs_env%num_time_freq_points
         ! Fourier transformation of χ_PQ(iτ) to χ_PQ(iω_j)
         CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)

         ! ε(iω_j) = Id - V^0.5*M^-1*χ(iω_j)*M^-1*V^0.5
         ! W(iω_j) = V^0.5*(ε^-1(iω_j)-Id)*V^0.5
         CALL compute_fm_W_freq(bs_env, bs_env%fm_chi_Gamma_freq, fm_V_sqrt, &
                                fm_M_inv_V_sqrt, bs_env%fm_W_MIC_freq)

         ! Fourier transform from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
         CALL Fourier_transform_w_to_t(bs_env, fm_W_time, bs_env%fm_W_MIC_freq, j_w)
      END DO

      ! M^-1*W^MIC(iτ)*M^-1
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_W_time)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T55,A,F10.1,A)') &
            'Computed W(iτ),', ' Execution time', m_walltime() - t1, ' s'
      END IF

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      ! Cleanup
      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_V_sqrt)
      CALL cp_fm_release(fm_M_inv_V_sqrt)

      ! Marek : Fourier transform W^MIC(itau) back to get it at a specific im.frequency point - iomega = 0
      IF (bs_env%rtp_method == rtp_method_bse) THEN
         t1 = m_walltime()
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
         ! Set to zero
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_zero, 0.0_dp)
         ! Sum over all times
         DO i_t = 1, bs_env%num_time_freq_points
            ! Add the relevant structure with correct weight
            CALL cp_fm_scale_and_add(1.0_dp, bs_env%fm_W_MIC_freq_zero, &
                                     bs_env%imag_time_weights_freq_zero(i_t), fm_W_time(i_t))
         END DO
         ! Done, save to file
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
         ! Report calculation
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T55,A,F10.1,A)') &
               'Computed W(0),', ' Execution time', m_walltime() - t1, ' s'
         END IF
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W

! **************************************************************************************************
!> \brief Computes V, V^0.5, and M^-1 * V^0.5
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
! **************************************************************************************************

   SUBROUTINE compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_V, fm_V_sqrt, fm_Minv_Vsqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_MinvVsqrt'

      INTEGER                                            :: handle, info, n_RI, ndep
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      CALL cp_fm_create(fm_work, fm_V%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. Build Coulomb Matrix V(k=0) using the kp-routine but only for ikp=1
      ! -----------------------------------------------------------------------
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      ALLOCATE (mat_V_kp(1:1, 1:2))
      NULLIFY (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
      ALLOCATE (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)

      CALL dbcsr_create(mat_V_kp(1, 1)%matrix, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 1)%matrix)
      CALL dbcsr_set(mat_V_kp(1, 1)%matrix, 0.0_dp)

      ! Dummy imaginary part just to satisfy the routine
      CALL dbcsr_create(mat_V_kp(1, 2)%matrix, template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 2)%matrix)
      CALL dbcsr_set(mat_V_kp(1, 2)%matrix, 0.0_dp)

      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, "RI_AUX", cell, &
                                      particle_set, qs_kind_set, atomic_kind_set, &
                                      bs_env%size_lattice_sum_V, operator_coulomb, 1, 1)

      ! Copy real part to fm_V
      CALL copy_dbcsr_to_fm(mat_V_kp(1, 1)%matrix, fm_V)

      CALL dbcsr_deallocate_matrix(mat_V_kp(1, 1)%matrix)
      CALL dbcsr_deallocate_matrix(mat_V_kp(1, 2)%matrix)
      DEALLOCATE (mat_V_kp)

      ! -----------------------------------------------------------------------
      ! 2. Get RI-Metric Matrix M(k=0)
      ! -----------------------------------------------------------------------
      CALL RI_2c_integral_mat(qs_env, fm_M, fm_V, n_RI, bs_env%ri_metric, &
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)

      ! -----------------------------------------------------------------------
      ! 3. M -> M^-1
      ! -----------------------------------------------------------------------
      CALL cp_fm_cholesky_decompose(fm_M(1, 1), info_out=info)
      IF (info == 0) THEN
         CALL cp_fm_cholesky_invert(fm_M(1, 1))
         CALL cp_fm_uplo_to_full(fm_M(1, 1), fm_work)
      ELSE
         ! Fallback if Cholesky fails due to conditioning
         CALL cp_fm_power(fm_M(1, 1), fm_work, -1.0_dp, bs_env%eps_eigval_mat_RI, ndep)
         CALL cp_fm_to_fm(fm_work, fm_M(1, 1))
      END IF

      ! -----------------------------------------------------------------------
      ! 4. V -> V^0.5
      ! -----------------------------------------------------------------------
      CALL cp_fm_to_fm(fm_V, fm_V_sqrt)
      CALL cp_fm_cholesky_decompose(fm_V_sqrt, info_out=info)
      IF (info == 0) THEN
         CALL clean_lower_part(fm_V_sqrt)
      ELSE
         CALL cp_fm_power(fm_V, fm_V_sqrt, 0.5_dp, bs_env%eps_eigval_mat_RI, ndep)
      END IF

      ! -----------------------------------------------------------------------
      ! 5. M^-1 * V^0.5
      ! -----------------------------------------------------------------------
      CALL parallel_gemm("N", "T", n_RI, n_RI, n_RI, 1.0_dp, fm_M(1, 1), fm_V_sqrt, &
                         0.0_dp, fm_Minv_Vsqrt)

      CALL cp_fm_release(fm_M)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE compute_V_MinvVsqrt

! **************************************************************************************************
!> \brief Computes W(iω) from χ_PQ(iω_j)
!> \param bs_env ...
!> \param fm_chi_freq_j ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
!> \param fm_W_freq_j ...
! **************************************************************************************************

   SUBROUTINE compute_fm_W_freq(bs_env, fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt, fm_W_freq_j)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_W_freq_j

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_fm_W_freq'

      INTEGER                                            :: handle, info, n_RI, ndep
      TYPE(cp_fm_type)                                   :: fm_eps_freq_j, fm_work

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      CALL cp_fm_create(fm_eps_freq_j, fm_chi_freq_j%matrix_struct)
      CALL cp_fm_create(fm_work, fm_chi_freq_j%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. ε(iω_j) = Id - (M^-1 * V^0.5)^T * χ(iω_j) * (M^-1 * V^0.5)
      ! -----------------------------------------------------------------------
      ! work = χ(iω_j) * (M^-1 * V^0.5)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_chi_freq_j, fm_Minv_Vsqrt, 0.0_dp, fm_work)

      ! eps_work = (M^-1 * V^0.5)^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_Minv_Vsqrt, fm_work, 0.0_dp, fm_eps_freq_j)

      ! ε(iω_j) = Id - eps_work  -->  -eps_work + Id
      CALL fm_add_on_diag(fm_eps_freq_j, 1.0_dp)

      ! Force perfect symmetry before Cholesky to avoid info != 0 due to GEMM noise
      CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)

      ! -----------------------------------------------------------------------
      ! 2. W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      ! -----------------------------------------------------------------------

      ! a) Cholesky decomposition of ε(iω_j)
      CALL cp_fm_cholesky_decompose(fm_eps_freq_j, info_out=info)

      ! b) Inversion
      IF (info == 0) THEN
         CALL cp_fm_cholesky_invert(fm_eps_freq_j)
         CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)
      ELSE
         ! Fallback to expensive diagonalization if Cholesky fails
         CALL cp_fm_power(fm_eps_freq_j, fm_work, -1.0_dp, bs_env%eps_eigval_mat_RI, ndep)
         CALL cp_fm_to_fm(fm_work, fm_eps_freq_j)
      END IF

      ! c) ε^-1(iω_j) - Id
      CALL fm_add_on_diag(fm_eps_freq_j, -1.0_dp)

      ! d) work = (ε^-1(iω_j) - Id) * V^0.5
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_eps_freq_j, fm_V_sqrt, &
                         0.0_dp, fm_work)

      ! e) W(iw) = V^0.5^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_V_sqrt, fm_work, &
                         0.0_dp, fm_W_freq_j)

      ! Cleanup
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_eps_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_freq

! **************************************************************************************************
!> \brief Adds a real scalar value to the diagonal of a real full matrix (fm)
!> \param fm ...
!> \param alpha ...
! **************************************************************************************************

   SUBROUTINE fm_add_on_diag(fm, alpha)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm
      REAL(KIND=dp), INTENT(IN)                          :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_add_on_diag'

      INTEGER                                            :: handle, i_global, i_row, j_col, &
                                                            j_global, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            i_global = row_indices(i_row)
            IF (j_global == i_global) THEN
               fm%local_data(i_row, j_col) = fm%local_data(i_row, j_col) + alpha
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_add_on_diag

! **************************************************************************************************
!> \brief Zeroes out the strictly lower triangular part of a real matrix
!> \param fm_mat ...
! **************************************************************************************************
   SUBROUTINE clean_lower_part(fm_mat)
      TYPE(cp_fm_type)                                   :: fm_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'clean_lower_part'

      INTEGER                                            :: handle, i_row, j_col, j_global, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            IF (j_global < row_indices(i_row)) fm_mat%local_data(i_row, j_col) = 0.0_dp
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE clean_lower_part

! **************************************************************************************************
!> \brief Computes the exact exchange part of the GW self-energy
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy
!> \param bs_env ...
!> \param fm_W_time ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
!> \brief DBCSR Topology Generation
!> \param matrix_template ...
!> \param dim_type ...
!> \param square_dist ...
!> \param blk_sizes ...
!> \param mapped_dist ...
! **************************************************************************************************

   SUBROUTINE setup_square_topology(matrix_template, dim_type, square_dist, blk_sizes, mapped_dist)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      CHARACTER(LEN=*), INTENT(IN)                       :: dim_type
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: square_dist
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: blk_sizes, mapped_dist

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_square_topology'

      INTEGER                                            :: handle, i, np, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      TYPE(dbcsr_distribution_type)                      :: dist_template

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_template, distribution=dist_template, &
                          row_blk_size=row_blk, col_blk_size=col_blk)
      CALL dbcsr_distribution_get(dist_template, row_dist=row_dist, col_dist=col_dist, &
                                  nprows=nprows, npcols=npcols)

      IF (TRIM(dim_type) == 'ROW') THEN
         ! Creates ROW x ROW (e.g., Grid x Grid from mat_phi_mu_l)
         blk_sizes => row_blk
         np = npcols
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=row_dist, col_dist=mapped_dist)

      ELSE IF (TRIM(dim_type) == 'COL') THEN
         ! Creates COL x COL (e.g., Aux x Aux from mat_Z_lP)
         blk_sizes => col_blk
         np = nprows
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=mapped_dist, col_dist=col_dist)
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_square_topology

! **************************************************************************************************
!> \brief DBCSR matrices deallocation
!> \param dist ...
!> \param mapped_dist    ...
!> \param m1 ...
!> \param m2 ...
!> \param m3 ...
!> \param m4 ...
! **************************************************************************************************

   SUBROUTINE release_dbcsr_topology_and_matrices(dist, mapped_dist, m1, m2, m3, m4)

      TYPE(dbcsr_distribution_type), INTENT(INOUT), &
         OPTIONAL                                        :: dist
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL, &
         POINTER                                         :: mapped_dist
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: m1, m2, m3, m4

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_dbcsr_topology_and_matrices'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (PRESENT(dist)) CALL dbcsr_distribution_release(dist)
      IF (PRESENT(mapped_dist)) THEN
         IF (ASSOCIATED(mapped_dist)) THEN
            DEALLOCATE (mapped_dist)
            NULLIFY (mapped_dist)
         END IF
      END IF
      IF (PRESENT(m1)) CALL dbcsr_release(m1)
      IF (PRESENT(m2)) CALL dbcsr_release(m2)
      IF (PRESENT(m3)) CALL dbcsr_release(m3)
      IF (PRESENT(m4)) CALL dbcsr_release(m4)

      CALL timestop(handle)

   END SUBROUTINE release_dbcsr_topology_and_matrices

END MODULE gw_large_cell_Gamma_ri_rs
