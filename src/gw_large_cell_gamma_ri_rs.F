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
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_set, dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
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
                                              compute_fm_chi_Gamma_freq,&
                                              create_fm_W_MIC_time,&
                                              fm_write,&
                                              get_W_MIC,&
                                              multiply_fm_W_MIC_time_with_Minv_Gamma
   USE gw_non_periodic_ri_rs,           ONLY: atomic_basis_at_grid_point,&
                                              compute_Sigma_c_and_QP_energies,&
                                              compute_Sigma_x,&
                                              compute_coeff_Z_lP,&
                                              get_mat_chi_Gamma_tau,&
                                              precompute_ri_rs_radii,&
                                              print_ri_rs_memory_estimate,&
                                              ri_rs_grid_assembler
   USE gw_utils,                        ONLY: de_init_bs_env
   USE input_constants,                 ONLY: rtp_method_bse
   USE kinds,                           ONLY: dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE machine,                         ONLY: m_walltime
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

      CALL timeset(routineN, handle)

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
       !! 4. Compute Independent-Particle Polarizability (χ)
       !!    G^occ_µλ(i|τ|)  = sum_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^vir_µλ(i|τ|)  = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^occ_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^occ_µν Φ_ν(r_l')
       !!    G^vir_ll'(i|τ|) = sum_µν Φ_µ(r_l) G^vir_µν Φ_ν(r_l')
       !!    χ_ll'(iτ)       = G^occ_ll'(i|τ|) * G^vir_ll'(i|τ|)
       !!    χ_PQ(iτ)        = sum_ll' Z_lP χ_ll'(iτ) Z_l'Q
       !!========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

       !!========================================================================
       !! 5. Compute Screened Interaction (W^MIC)
       !!    χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W^MIC_PQ(iτ)
       !!========================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

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
      CALL compute_Sigma_c_and_QP_energies(bs_env, qs_env, fm_W_time, fm_Sigma_x_Gamma)

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
