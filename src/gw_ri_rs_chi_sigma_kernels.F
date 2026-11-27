!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared RI-RS polarizability and self-energy assembly.
! **************************************************************************************************
MODULE gw_ri_rs_chi_sigma_kernels
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_distribution_type,&
                                              dbcsr_filter,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_scale,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE gw_tensor_large_cell_Gamma,      ONLY: delete_unnecessary_files
   USE gw_utils,                        ONLY: G_occ_vir
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: cfm_ikp_from_fm_Gamma
   USE qs_environment_types,            ONLY: qs_environment_type
   USE ri_rs_matrices,                  ONLY: build_ri_rs_bare_aux_interaction,&
                                              contract_grid_panels,&
                                              print_matrix_occupation,&
                                              release_dbcsr_topology_and_matrices,&
                                              resolve_grid_panels,&
                                              setup_square_topology
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_chi_sigma_kernels'
   PUBLIC :: get_mat_chi_Gamma_tau, compute_Sigma_x, compute_Sigma_c, &
             sigma_mo_diag_type, sigma_mo_diag_at_k, sigma_mo_diag_release

! **************************************************************************************************
!> \brief MO diagonals of Σ^c at the band/DOS k-points (Γ-only RI-RS GW). Each dense
!>        Σ^c_µν(±i|τ_j|) is reduced to these diagonals and released as soon as it is built:
!>          Σ^c_nn(k, ±iτ) = Σ_µ C*_µn(k) [Σ^c(k, ±iτ) C(k)]_µn,  Σ^c_µν(k) = Σ_R e^{ik·R} Σ^c_µν(R)
!>        with the MIC phase of cfm_ikp_from_fm_Gamma; repeated k-points are evaluated once.
! **************************************************************************************************
   TYPE sigma_mo_diag_type
      TYPE(qs_environment_type), POINTER                 :: qs_env => NULL()
      INTEGER                                            :: n_kp_distinct = 0
      ! ikp_distinct(i_kd): first occurrence of distinct k-point i_kd in bs_env%kpoints_DOS
      ! i_kd_of_ikp(ikp): the distinct k-point that ikp copies
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ikp_distinct, i_kd_of_ikp
      ! C(k), solutions of H(k) C(k) = S(k) C(k) ϵ(k), per distinct k-point and spin
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:, :)    :: cfm_mos
      ! Σ^c_nn(k, ±i|τ_j|): (n_ao, n_t, positive_time/negative_time, n_kp_distinct, n_spin)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :, :) :: Sigma_c_n_time
      ! 3D Γ-point correction ΔΣ^c_n(±i|τ_j|): (n_mo, n_t, 2, n_spin), added in the AO basis as
      ! (S C) diag(ΔΣ_n) (S C)ᵀ before the k-transform; not allocated without the correction
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dSigma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_SC
      ! work: Σ(k) and Σ(k) C(k) (complex AO x AO), (S C) diag(ΔΣ) (real AO x AO)
      TYPE(cp_cfm_type)                                  :: cfm_ikp, cfm_work
      TYPE(cp_fm_type)                                   :: fm_C_scaled
   END TYPE sigma_mo_diag_type
CONTAINS

! **************************************************************************************************
!> \brief Computes χ_PQ(iτ) from the occupied and virtual Green's functions.
!> \param bs_env ...
!> \param mat_chi_Gamma_tau Response matrices χ_PQ(iτ) to be computed
!> \param cutoff grid-centroid distance cutoff; zero disables geometric truncation
! **************************************************************************************************
   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, mat_chi_Gamma_tau, cutoff)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER                                            :: handle, i_t, ispin, n_panels
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      REAL(KIND=dp)                                      :: first_grid_occ, t1, tau
      TYPE(dbcsr_type)                                   :: matrix_G_occ_ao(1), matrix_G_vir_ao

      CALL timeset(routineN, handle)

      ! Panel boundaries for the grid-streaming contraction.
      ! The panels are identical for χ, Σ^x and Σ^c, so the count is reported once here for all three stages.

      CALL resolve_grid_panels(bs_env%ri_rs, bs_env%para_env, bs_env%dbcsr_msg_elem_limit, &
                               pan_first, pan_last, cutoff, unit_nr=bs_env%unit_nr)

      n_panels = SIZE(pan_first)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T74,I9)') &
            'Number of batches for χ, Σ matrices', n_panels
         WRITE (bs_env%unit_nr, '(A)') ' '
         CALL m_flush(bs_env%unit_nr)
      END IF

      ! =========================================================================
      ! IMAGINARY TIME LOOP
      ! χ_PQ(iτ) = Σ_s g_s · Z^T ( (ϕ G^occ_s ϕ^T) ∘ (ϕ G^vir_s ϕ^T) ) Z
      ! (g_s = spin degeneracy)
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         t1 = m_walltime()
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)

         DO ispin = 1, bs_env%n_spin

            ! AO-space Green's functions G^occ_µν, G^vir_µν (dense AO x AO, small)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_occ_ao(1))
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_vir_ao)

            ! χ_PQ += g_s · Z^T ( (ϕ G^occ ϕ^T) ∘ (ϕ G^vir ϕ^T) ) Z
            IF (i_t == 1) THEN
               CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_G_occ_ao, &
                                         L_B=bs_env%ri_rs%mat_phi_mu_l, M_B=matrix_G_vir_ao, &
                                         L_out=bs_env%ri_rs%mat_Z_lP, mat_out=mat_chi_Gamma_tau(i_t:i_t), &
                                         scale=bs_env%spin_degeneracy, eps=bs_env%eps_filter, &
                                         para_env=bs_env%para_env, &
                                         pan_first=pan_first, pan_last=pan_last, &
                                         lb_eq_la=.TRUE., lout_eq_la=.FALSE., &
                                         zero_out=(ispin == 1), &
                                         cell=bs_env%ri_rs%cell, centroids=bs_env%ri_rs%chunk_centroids, &
                                         cutoff=cutoff, grid_occupation=first_grid_occ)
            ELSE
               CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_G_occ_ao, &
                                         L_B=bs_env%ri_rs%mat_phi_mu_l, M_B=matrix_G_vir_ao, &
                                         L_out=bs_env%ri_rs%mat_Z_lP, mat_out=mat_chi_Gamma_tau(i_t:i_t), &
                                         scale=bs_env%spin_degeneracy, eps=bs_env%eps_filter, &
                                         para_env=bs_env%para_env, &
                                         pan_first=pan_first, pan_last=pan_last, &
                                         lb_eq_la=.TRUE., lout_eq_la=.FALSE., &
                                         zero_out=(ispin == 1), &
                                         cell=bs_env%ri_rs%cell, centroids=bs_env%ri_rs%chunk_centroids, &
                                         cutoff=cutoff)
            END IF

            CALL dbcsr_release(matrix_G_occ_ao(1))
            CALL dbcsr_release(matrix_G_vir_ao)

         END DO ! ispin

         ! Sparsity reports
         IF (i_t == 1) THEN
            CALL print_matrix_occupation(bs_env%ri_rs%mat_Z_lP, 'Z_lP', bs_env%para_env, bs_env%unit_nr, &
                                         bs_env%dbcsr_msg_elem_limit)
            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,T73,F7.2,A)') &
                  'Percentage of non-zero matrix elements in G_ll'', χ_ll'', W_ll''', &
                  first_grid_occ*100.0_dp, ' %'
               CALL m_flush(bs_env%unit_nr)
            END IF
            CALL print_matrix_occupation(mat_chi_Gamma_tau(i_t)%matrix, 'χ_PQ', &
                                         bs_env%para_env, bs_env%unit_nr, bs_env%dbcsr_msg_elem_limit, &
                                         suffix=' for time point 1')
            IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
         END IF

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I13,A,I3,A,F7.1,A)') &
               'Computed χ(iτ,k=0) for time point', i_t, ' /', bs_env%num_time_freq_points, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! i_t

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE get_mat_chi_Gamma_tau

! **************************************************************************************************
!> \brief Computes the exact-exchange part of the GW self-energy:
!>          D_μν        = Σ_n^occ C_μn C_νn                     (density matrix = G^occ at τ=0)
!>          V^tr_PQ     = M^-1 (P|Q)_trunc M^-1                 (truncated Coulomb, RI basis)
!>          Σ^x_λσ(k=0) = -Σ_ll' ϕ_λ(r_l) [ (ϕ D ϕ^T)_ll' ∘ (Z V^tr Z^T)_ll' ] ϕ_σ(r_l')
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param cutoff grid-centroid distance cutoff; zero disables geometric truncation
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma, cutoff)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, ispin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      REAL(KIND=dp)                                      :: t1
      TYPE(dbcsr_p_type)                                 :: mat_Sigma_x_Gamma(1)
      TYPE(dbcsr_type)                                   :: matrix_D_ao(1), matrix_V_aux
      TYPE(dbcsr_type), TARGET                           :: matrix_Sigma_x

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      mat_Sigma_x_Gamma(1)%matrix => matrix_Sigma_x
      CALL dbcsr_create(matrix_Sigma_x, template=bs_env%mat_ao_ao%matrix)

      CALL resolve_grid_panels(bs_env%ri_rs, bs_env%para_env, bs_env%dbcsr_msg_elem_limit, &
                               pan_first, pan_last, cutoff, unit_nr=bs_env%unit_nr)

      ! =========================================================================
      ! 1. COMPUTE V^tr_PQ (RI x RI)
      ! =========================================================================
      CALL build_ri_rs_bare_aux_interaction(qs_env, bs_env%fm_Minv_Gamma, bs_env%trunc_coulomb, &
                                            bs_env%auto_ri, bs_env%ri_rs%mat_Z_lP, matrix_V_aux)
      CALL dbcsr_filter(matrix_V_aux, bs_env%eps_filter)

      ! =========================================================================
      ! 2. SPIN LOOP FOR EXACT EXCHANGE
      ! Σ^x_λσ = -Σ_ll' ϕ_λ(r_l) ( D_ll' V^tr_ll' ) ϕ_σ(r_l')
      !        = -ϕ^T ( (ϕ D ϕ^T) ∘ (Z V^tr Z^T) ) ϕ
      ! =========================================================================
      DO ispin = 1, bs_env%n_spin

         ! AO-space density matrix D_µν = G^occ at τ = 0
         CALL build_G_ao(bs_env, 0.0_dp, ispin, .TRUE., .FALSE., &
                         bs_env%ri_rs%mat_phi_mu_l, matrix_D_ao(1))

         CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_D_ao, &
                                   L_B=bs_env%ri_rs%mat_Z_lP, M_B=matrix_V_aux, &
                                   L_out=bs_env%ri_rs%mat_phi_mu_l, mat_out=mat_Sigma_x_Gamma, &
                                   scale=1.0_dp, eps=bs_env%eps_filter, &
                                   para_env=bs_env%para_env, &
                                   pan_first=pan_first, pan_last=pan_last, &
                                   lb_eq_la=.FALSE., lout_eq_la=.TRUE., zero_out=.TRUE., &
                                   cell=bs_env%ri_rs%cell, centroids=bs_env%ri_rs%chunk_centroids, &
                                   cutoff=cutoff)
         CALL dbcsr_scale(matrix_Sigma_x, -1.0_dp)

         CALL dbcsr_release(matrix_D_ao(1))

         ! Data I/O and Export to CP2K Full Matrices
         CALL copy_dbcsr_to_fm(matrix_Sigma_x, fm_Sigma_x_Gamma(ispin))

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      ! =========================================================================
      ! 3. CLEANUP
      ! =========================================================================
      CALL dbcsr_release(matrix_V_aux)
      CALL dbcsr_release(matrix_Sigma_x)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy on the imaginary-time grid:
!>          Σ^c_λσ(iτ<0) = -Σ_ll' ϕ_λ(r_l) [ (ϕ G^occ ϕᵀ)_ll' ∘ (Z W^MIC Zᵀ)_ll' ] ϕ_σ(r_l')
!>          Σ^c_λσ(iτ>0) = +Σ_ll' ϕ_λ(r_l) [ (ϕ G^vir ϕᵀ)_ll' ∘ (Z W^MIC Zᵀ)_ll' ] ϕ_σ(r_l')
!>        with G^occ/vir_µν(i|τ|) = Σ_n^occ/vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_νn, l, l' grid points and
!>        W^MIC = M⁻¹ W M⁻¹ (M⁺ W M⁺ in the Γ-only driver).
!> \param bs_env band-structure environment (grid, Z_lP, ϕ_µ(r_l), MOs, time grid)
!> \param fm_W_time W^MIC_PQ(iτ_j), one dense RI x RI matrix per time point
!> \param fm_Sigma_c_Gamma_time Σ^c_λσ(±i|τ_j|) at Γ, (n_t, positive/negative time, n_spin);
!>        allocated here unless mo_diag is present
!> \param cutoff grid-centroid distance cutoff; zero disables geometric truncation
!> \param keep_w_and_files keep W and scratch files for another Sigma_c evaluation
!> \param release_w_per_time release W^MIC(iτ_j) right after it was copied to DBCSR (a single
!>        Σ^c evaluation, G0W0); scratch files are still governed by keep_w_and_files
!> \param mo_diag if present, each Σ^c_µν(±i|τ_j|) is reduced to its MO diagonals at the band/DOS
!>        k-points (sigma_mo_diag_type) and released; fm_Sigma_c_Gamma_time is then not allocated
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c(bs_env, fm_W_time, fm_Sigma_c_Gamma_time, cutoff, keep_w_and_files, &
                              release_w_per_time, mo_diag)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_w_and_files, release_w_per_time
      TYPE(sigma_mo_diag_type), INTENT(INOUT), OPTIONAL  :: mo_diag

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      CHARACTER(LEN=32)                                  :: sigma_label
      INTEGER                                            :: handle, i_sector, i_t, ispin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: blk_aux, dist_row_aux
      LOGICAL                                            :: release_w, retain_work
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(cp_fm_type), DIMENSION(2)                     :: fm_Sigma_pair
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux
      TYPE(dbcsr_p_type)                                 :: mat_Sigma_tau(2)
      TYPE(dbcsr_type)                                   :: matrix_G_ao(2), matrix_W_aux
      TYPE(dbcsr_type), TARGET                           :: matrix_Sigma_tau(2)

      CALL timeset(routineN, handle)

      retain_work = .FALSE.
      IF (PRESENT(keep_w_and_files)) retain_work = keep_w_and_files
      release_w = .FALSE.
      IF (PRESENT(release_w_per_time)) release_w = release_w_per_time

      ! =========================================================================
      ! 1. SETUP AUXILIARY TOPOLOGY AND PRE-ALLOCATE OUTPUT ARRAYS
      ! =========================================================================
      CALL setup_square_topology(bs_env%ri_rs%mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL resolve_grid_panels(bs_env%ri_rs, bs_env%para_env, bs_env%dbcsr_msg_elem_limit, &
                               pan_first, pan_last, cutoff, unit_nr=bs_env%unit_nr)

      ! Only the final FM time array persists; DBCSR outputs are streamed one time/spin at a time.
      ! With mo_diag, only the (+, −) pair of one time point and spin is held as dense matrices.
      IF (.NOT. PRESENT(mo_diag)) THEN
         ALLOCATE (fm_Sigma_c_Gamma_time(bs_env%num_time_freq_points, 2, bs_env%n_spin))
      END IF
      DO i_sector = 1, 2
         mat_Sigma_tau(i_sector)%matrix => matrix_Sigma_tau(i_sector)
      END DO

      ! =========================================================================
      ! 2. IMAGINARY TIME LOOP
      ! Σ^c_neg_λσ(iτ) = -ϕ^T ( (ϕ G^occ ϕ^T) ∘ (Z W^MIC Z^T) ) ϕ
      ! Σ^c_pos_λσ(iτ) =  ϕ^T ( (ϕ G^vir ϕ^T) ∘ (Z W^MIC Z^T) ) ϕ
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)

         CALL dbcsr_create(matrix_W_aux, "W_aux", dist_aux_aux, dbcsr_type_no_symmetry, &
                           blk_aux, blk_aux)
         CALL copy_fm_to_dbcsr(fm_W_time(i_t), matrix_W_aux, keep_sparsity=.FALSE.)
         CALL dbcsr_filter(matrix_W_aux, bs_env%eps_filter)

         ! The dense W^MIC_PQ(iτ_j) is not needed again once it lives in matrix_W_aux.
         IF (release_w) CALL cp_fm_release(fm_W_time(i_t))

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            ! AO-space Green's functions G^occ_µν, G^vir_µν (dense AO x AO, small)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_ao(1))
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_ao(2))

            DO i_sector = 1, 2
               CALL dbcsr_create(matrix_Sigma_tau(i_sector), template=bs_env%mat_ao_ao%matrix)
            END DO

            ! Σ^c_neg and Σ^c_pos in a single panel loop: W_pan = Z_panel × W × Z^T built once
            CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_G_ao, &
                                      L_B=bs_env%ri_rs%mat_Z_lP, M_B=matrix_W_aux, &
                                      L_out=bs_env%ri_rs%mat_phi_mu_l, &
                                      mat_out=mat_Sigma_tau, &
                                      scale=1.0_dp, eps=bs_env%eps_filter, para_env=bs_env%para_env, &
                                      pan_first=pan_first, pan_last=pan_last, &
                                      lb_eq_la=.FALSE., lout_eq_la=.TRUE., zero_out=.TRUE., &
                                      cell=bs_env%ri_rs%cell, centroids=bs_env%ri_rs%chunk_centroids, cutoff=cutoff)
            CALL dbcsr_scale(matrix_Sigma_tau(1), -1.0_dp)

            ! FM time sector 1 is positive, sector 2 is negative.
            DO i_sector = 1, 2
               CALL dbcsr_release(matrix_G_ao(i_sector))
               IF (PRESENT(mo_diag)) THEN
                  CALL cp_fm_create(fm_Sigma_pair(3 - i_sector), bs_env%fm_s_Gamma%matrix_struct)
                  CALL copy_dbcsr_to_fm(matrix_Sigma_tau(i_sector), fm_Sigma_pair(3 - i_sector))
               ELSE
                  CALL cp_fm_create(fm_Sigma_c_Gamma_time(i_t, 3 - i_sector, ispin), &
                                    bs_env%fm_s_Gamma%matrix_struct)
                  CALL copy_dbcsr_to_fm(matrix_Sigma_tau(i_sector), &
                                        fm_Sigma_c_Gamma_time(i_t, 3 - i_sector, ispin))
               END IF
            END DO

            ! Report the first-time outputs before releasing the current DBCSR pair.
            IF (i_t == 1) THEN
               sigma_label = 'Σ^c'
               IF (bs_env%n_spin == 2 .AND. ispin == 1) sigma_label = 'Σ^c_α'
               IF (bs_env%n_spin == 2 .AND. ispin == 2) sigma_label = 'Σ^c_β'
               CALL print_matrix_occupation(matrix_Sigma_tau(2), &
                                            TRIM(sigma_label)//'(iτ)', bs_env%para_env, bs_env%unit_nr, &
                                            bs_env%dbcsr_msg_elem_limit, suffix=' for time point 1')
               CALL print_matrix_occupation(matrix_Sigma_tau(1), &
                                            TRIM(sigma_label)//'(-iτ)', bs_env%para_env, bs_env%unit_nr, &
                                            bs_env%dbcsr_msg_elem_limit, suffix=' for time point 1')
               IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
            END IF
            DO i_sector = 1, 2
               CALL dbcsr_release(matrix_Sigma_tau(i_sector))
            END DO

            ! Σ^c_µν(±i|τ_j|) (+ 3D ΔΣ^c) -> Σ^c_nn(k, ±i|τ_j|) at every distinct k-point
            IF (PRESENT(mo_diag)) THEN
               CALL sigma_mo_diag_add(bs_env, mo_diag, i_t, ispin, fm_Sigma_pair)
               DO i_sector = 1, 2
                  CALL cp_fm_release(fm_Sigma_pair(i_sector))
               END DO
            END IF

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I15,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ) for time point', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

         CALL dbcsr_release(matrix_W_aux)

      END DO ! i_t

      ! The DBCSR topology is local to this contraction.
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      ! The non-periodic evGW0 path keeps these across Sigma_c evaluations.
      IF (.NOT. retain_work) THEN
         CALL cp_fm_release(fm_W_time)
         CALL delete_unnecessary_files(bs_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief Adds the 3D Γ-point correction to one Σ^c_µν(±i|τ_j|) pair of one spin and stores its MO
!>        diagonals at every distinct band/DOS k-point:
!>          Σ^c_µν(±i|τ_j|) += Σ_n (S C)_µn ΔΣ^c_n(±i|τ_j|) (S C)_νn     (if mo_diag%dSigma)
!>          Σ^c_nn(k, ±i|τ_j|) = Σ_µ C*_µn(k) [Σ^c(k, ±i|τ_j|) C(k)]_µn
!> \param bs_env band-structure environment (AO and MO matrix layouts)
!> \param mo_diag MO-diagonal storage, C(k) and work matrices
!> \param i_t time point
!> \param ispin spin
!> \param fm_Sigma_pair Σ^c_µν(+i|τ_j|), Σ^c_µν(−i|τ_j|) (positive_time, negative_time);
!>        overwritten by the corrected matrices
! **************************************************************************************************
   SUBROUTINE sigma_mo_diag_add(bs_env, mo_diag, i_t, ispin, fm_Sigma_pair)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(sigma_mo_diag_type), INTENT(INOUT)            :: mo_diag
      INTEGER, INTENT(IN)                                :: i_t, ispin
      TYPE(cp_fm_type), DIMENSION(2), INTENT(INOUT)      :: fm_Sigma_pair

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sigma_mo_diag_add'

      INTEGER                                            :: handle, i_kd, i_sector, n_ao, n_mo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: scaling

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      ! work matrices of this step only: Σ(k), Σ(k) C(k) (complex) and, for the 3D shift,
      ! (S C) diag(ΔΣ) (real)
      CALL cp_cfm_create(mo_diag%cfm_ikp, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(mo_diag%cfm_work, bs_env%fm_s_Gamma%matrix_struct)
      IF (ALLOCATED(mo_diag%dSigma)) THEN
         CALL cp_fm_create(mo_diag%fm_C_scaled, bs_env%fm_mo_coeff_Gamma(1)%matrix_struct)
      END IF
      DO i_sector = 1, 2
         ! ΔΣ^c_µν(±iτ) = (S C) diag(ΔΣ_n) (S C)ᵀ, folded onto the k-points like Σ^c (Cᵀ S C = 1)
         IF (ALLOCATED(mo_diag%dSigma)) THEN
            n_mo = SIZE(mo_diag%dSigma, 1)
            ALLOCATE (scaling(n_ao))
            scaling(:) = 0.0_dp
            scaling(1:n_mo) = mo_diag%dSigma(:, i_t, i_sector, ispin)
            CALL cp_fm_to_fm(mo_diag%fm_SC(ispin), mo_diag%fm_C_scaled)
            CALL cp_fm_column_scale(mo_diag%fm_C_scaled, scaling)
            CALL parallel_gemm('N', 'T', n_ao, n_ao, n_mo, 1.0_dp, mo_diag%fm_C_scaled, &
                               mo_diag%fm_SC(ispin), 1.0_dp, fm_Sigma_pair(i_sector))
            DEALLOCATE (scaling)
         END IF
         DO i_kd = 1, mo_diag%n_kp_distinct
            CALL sigma_mo_diag_at_k(bs_env, mo_diag, fm_Sigma_pair(i_sector), i_kd, ispin, &
                                    mo_diag%Sigma_c_n_time(:, i_t, i_sector, i_kd, ispin))
         END DO
      END DO

      CALL cp_cfm_release(mo_diag%cfm_ikp)
      CALL cp_cfm_release(mo_diag%cfm_work)
      CALL cp_fm_release(mo_diag%fm_C_scaled)

      CALL timestop(handle)

   END SUBROUTINE sigma_mo_diag_add

! **************************************************************************************************
!> \brief Real MO diagonal of a real Γ AO operator at distinct k-point i_kd, with one complex
!>        n_AO³ product; the real part selects the Hermitian part of the operator:
!>          O_µν(k) = Σ_R e^{ik·R} O_µν(R)   (MIC phase, cfm_ikp_from_fm_Gamma)
!>          X = O(k) C(k),   o_n = Re Σ_µ C*_µn(k) X_µn
!> \param bs_env band-structure environment (k-point list kpoints_DOS)
!> \param mo_diag C(k), k-point map and work matrices
!> \param fm_Gamma O_µν at Γ (the sum over images of the MIC)
!> \param i_kd distinct k-point
!> \param ispin spin of C(k)
!> \param diag o_n, n = 1 .. n_ao
! **************************************************************************************************
   SUBROUTINE sigma_mo_diag_at_k(bs_env, mo_diag, fm_Gamma, i_kd, ispin, diag)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(sigma_mo_diag_type), INTENT(INOUT)            :: mo_diag
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Gamma
      INTEGER, INTENT(IN)                                :: i_kd, ispin
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: diag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sigma_mo_diag_at_k'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: diag_loc
      INTEGER                                            :: handle, i_col, i_row, n, n_ao, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      CALL cfm_ikp_from_fm_Gamma(mo_diag%cfm_ikp, fm_Gamma, mo_diag%ikp_distinct(i_kd), &
                                 mo_diag%qs_env, bs_env%kpoints_DOS, "ORB")
      ! X = O(k) C(k)
      CALL parallel_gemm('N', 'N', n_ao, n_ao, n_ao, z_one, mo_diag%cfm_ikp, &
                         mo_diag%cfm_mos(i_kd, ispin), z_zero, mo_diag%cfm_work)
      ! o_n = Σ_µ C*_µn X_µn on the local elements (C and X share the layout), summed over ranks
      CALL cp_cfm_get_info(mo_diag%cfm_work, nrow_local=nrow_local, ncol_local=ncol_local, &
                           col_indices=col_indices)
      ALLOCATE (diag_loc(n_ao))
      diag_loc(:) = (0.0_dp, 0.0_dp)
      DO i_col = 1, ncol_local
         n = col_indices(i_col)
         DO i_row = 1, nrow_local
            diag_loc(n) = diag_loc(n) + CONJG(mo_diag%cfm_mos(i_kd, ispin)%local_data(i_row, i_col))* &
                          mo_diag%cfm_work%local_data(i_row, i_col)
         END DO
      END DO
      CALL mo_diag%cfm_work%matrix_struct%para_env%sum(diag_loc)
      diag(1:n_ao) = REAL(diag_loc, KIND=dp)
      DEALLOCATE (diag_loc)

      CALL timestop(handle)

   END SUBROUTINE sigma_mo_diag_at_k

! **************************************************************************************************
!> \brief Releases the matrices and arrays of a sigma_mo_diag_type.
!> \param mo_diag MO-diagonal storage to release
! **************************************************************************************************
   SUBROUTINE sigma_mo_diag_release(mo_diag)
      TYPE(sigma_mo_diag_type), INTENT(INOUT)            :: mo_diag

      INTEGER                                            :: i_kd, ispin

      IF (ALLOCATED(mo_diag%cfm_mos)) THEN
         DO ispin = 1, SIZE(mo_diag%cfm_mos, 2)
            DO i_kd = 1, SIZE(mo_diag%cfm_mos, 1)
               CALL cp_cfm_release(mo_diag%cfm_mos(i_kd, ispin))
            END DO
         END DO
         DEALLOCATE (mo_diag%cfm_mos)
      END IF
      IF (ALLOCATED(mo_diag%fm_SC)) CALL cp_fm_release(mo_diag%fm_SC)
      CALL cp_cfm_release(mo_diag%cfm_ikp)
      CALL cp_cfm_release(mo_diag%cfm_work)
      CALL cp_fm_release(mo_diag%fm_C_scaled)
      IF (ALLOCATED(mo_diag%ikp_distinct)) DEALLOCATE (mo_diag%ikp_distinct)
      IF (ALLOCATED(mo_diag%i_kd_of_ikp)) DEALLOCATE (mo_diag%i_kd_of_ikp)
      IF (ALLOCATED(mo_diag%Sigma_c_n_time)) DEALLOCATE (mo_diag%Sigma_c_n_time)
      IF (ALLOCATED(mo_diag%dSigma)) DEALLOCATE (mo_diag%dSigma)
      NULLIFY (mo_diag%qs_env)
      mo_diag%n_kp_distinct = 0

   END SUBROUTINE sigma_mo_diag_release

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
      INTEGER, CONTIGUOUS, DIMENSION(:), POINTER         :: blk_ao, dist_row_ao
      TYPE(cp_fm_type)                                   :: fm_G
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao

      CALL timeset(routineN, handle)

      ! Dense G^occ/vir_µν(i|τ|), local to this routine, only for the conversion to DBCSR
      CALL cp_fm_create(fm_G, bs_env%fm_s_Gamma%matrix_struct)
      CALL G_occ_vir(bs_env, tau, fm_G, ispin, occ=occ, vir=vir)

      CALL setup_square_topology(template, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      CALL dbcsr_create(matrix_G_ao, name="G_ao", dist=dist_ao_ao, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=blk_ao, col_blk_size=blk_ao)

      CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.FALSE.)
      CALL cp_fm_release(fm_G)
      CALL dbcsr_filter(matrix_G_ao, bs_env%eps_filter)

      ! release only the topology; keep matrix_G_ao for the caller
      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao)

      CALL timestop(handle)

   END SUBROUTINE build_G_ao

END MODULE gw_ri_rs_chi_sigma_kernels
