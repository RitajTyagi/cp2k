!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using the RI-RS approximation for a large periodic unit cell sampled at the Gamma
!>        point only.
!>
!>        Every stage except the screened interaction is shared with the molecular driver in
!>        gw_ri_rs_non_periodic: the grid, phi, the Z_lP fit and the panel-streamed chi, Sigma^x
!>        and Sigma^c are all cell-agnostic once distances are measured through the minimum image
!>        and the periodic image sums are taken inside evaluate_ao_basis_on_points and
!>        compute_d_lp. For a molecule cell%perd = 0 and those sums collapse to a single term, so
!>        the same code serves both. What is genuinely periodic is get_W_MIC, which builds W on
!>        the chi/eps/W k-mesh and folds it back to the minimum image.
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_ri_rs_large_cell_Gamma
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE gw_ri_rs_grid_setup_main,        ONLY: setup_ri_rs_grid
   USE gw_ri_rs_compute_Z_lP,           ONLY: compute_Z_lP
   USE gw_ri_rs_non_periodic,           ONLY: atomic_basis_at_grid_point,&
                                              compute_Sigma_c,&
                                              compute_Sigma_x,&
                                              get_mat_chi_Gamma_tau,&
                                              init_evGW0,&
                                              print_evGW0_band_edges,&
                                              reset_and_clean_bs_env,&
                                              update_eigenvalues_G
   USE gw_ri_rs_panels,                 ONLY: print_ri_rs_memory_estimate
   USE gw_ri_rs_utils,                  ONLY: precompute_ri_rs_radii
   USE gw_tensor_large_cell_Gamma,      ONLY: compute_QP_energies,&
                                              delete_unnecessary_files,&
                                              get_W_MIC
   USE gw_utils,                        ONLY: de_init_bs_env
   USE input_constants,                 ONLY: G0W0
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_print_mem_per_rank
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_large_cell_Gamma'

   PUBLIC :: gw_calc_ri_rs_large_cell_Gamma

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using the RI-RS formalism for a large periodic unit cell at Gamma.
!> \param qs_env ...
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************
   SUBROUTINE gw_calc_ri_rs_large_cell_Gamma(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_large_cell_Gamma'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time

      CALL timeset(routineN, handle)

      ! ==============================================================================
      ! 0. Precompute AO and RI radii
      !    Per-atom cutoff radii from the most diffuse Gaussian primitives of
      !    the AO ("ORB") and RI auxiliary ("RI_AUX") basis sets:
      !      α_min,ao = min { ζ_ao | ζ_ao > 10⁻³ },  α_min,ri analogous
      !      r_ao = sqrt( -ln(ε) / α_min,ao )   (radius_ao_per_atom)
      !      r_ri = sqrt( -ln(ε) / α_min,ri )   (radius_ri_per_atom)
      ! ==============================================================================
      CALL precompute_ri_rs_radii(bs_env)

      ! ==============================================================================
      ! 1. Grid generation for RI-RS
      !    Modified Lebedev atomic grids (Duchemin & Blase), one per atom,
      !    concatenated into a flat global list:  r_l = R_A + r_l^(A)
      ! ==============================================================================
      CALL setup_ri_rs_grid(bs_env, bs_env%ri_rs%grid_points)

      ! ==============================================================================
      ! 2a. Atomic basis evaluation on the grid (grid x AO matrix)
      !     ϕ_μl = Σ_R ϕ_μ^R(r_l)   (mat_phi_mu_l), Bloch-summed over the images of μ
      ! ==============================================================================
      CALL atomic_basis_at_grid_point(bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)

      ! ==============================================================================
      ! 2b. Print the memory estimate for the RI-RS calculation
      ! ==============================================================================
      CALL print_ri_rs_memory_estimate(bs_env)

      ! ==============================================================================
      ! 3. RI-RS fitting coefficients Z_lP (grid x RI matrix)
      !    Per-atom regularized solve, restricted to grid points r_l within a
      !    cutoff distance of atom P:
      !      a. D_ll' = [ Σ_μ ϕ_μ(r_l) ϕ_μ(r_l') ]²
      !      b. D_lP  = Σ_μν Σ_RS ϕ_μ(r_l) ϕ_ν(r_l) (μ^R ν^S|P)
      !      c. Jacobi conditioning with d_l = 1/sqrt(D_ll):
      !         D'_ll' = d_l D_ll' d_l' + λδ_ll' ,  D'_lP = d_l D_lP
      !      d. Solve Σ_l' D'_ll' Z'_l'P = D'_lP
      !      e. Rescale Z_lP = d_l Z'_lP                        (mat_Z_lP)
      ! ==============================================================================
      CALL compute_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                        bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! flag the RI-RS grid as built so a subsequent RT-BSE run reuses Z_lP
      ! instead of rebuilding it
      bs_env%ri_rs%grid_built = .TRUE.

      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, &
                                 "after the RI-RS grid setup")

      ! ==============================================================================
      ! 4. Polarizability matrix χ on the imaginary-time grid, streamed over grid panels
      !      G^occ_µλ(i|τ|)  = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^occ_ll'(i|τ|) = Σ_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')  (G^vir analogous)
      !      χ_ll'(iτ)       = G^occ_ll'(i|τ|) ∘ G^vir_ll'(i|τ|)
      !      χ_PQ(iτ)        = Σ_ll' Z_lP χ_ll'(iτ) Z_l'Q
      ! ==============================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! ==============================================================================
      ! 5. Screened Coulomb interaction W, truncated to the minimum image in the RI
      !    auxiliary basis. This is the only genuinely periodic stage:
      !      χ_PQ(iτ,k=0) -> χ_PQ(iω,k) -> ε_PQ(iω,k) -> W_PQ(iω,k) -> W^MIC_PQ(iτ)
      ! ==============================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

      ! ==============================================================================
      ! 6. Exact-exchange self-energy Σ^x
      !      D_µν        = Σ_n^occ C_µn C_νn                (density matrix)
      !      V^tr_ll'    = Σ_PQ Z_lP V^tr_PQ Z_l'Q          (truncated Coulomb)
      !      Σ^x_λσ(k=0) = -Σ_ll' ϕ_λ(r_l) [ D_ll' ∘ V^tr_ll' ] ϕ_σ(r_l')
      ! ==============================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

      ! ==============================================================================
      ! 7. Correlation self-energy Σ^c and quasiparticle energies, iterated until
      !    eigenvalue self-consistency if &EVGW0 is given, single pass for G0W0.
      ! ==============================================================================
      CALL compute_Sigma_c_and_QP_energies(bs_env, qs_env, fm_W_time, fm_Sigma_x_Gamma)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_large_cell_Gamma

! **************************************************************************************************
!> \brief Correlation self-energy and quasiparticle energies, iterated to eigenvalue
!>        self-consistency in G (evGW0).
!>
!>        The loop body is identical to the molecular one in gw_ri_rs_non_periodic; only the
!>        quasiparticle step differs, because the periodic one has to take the Gamma matrices to
!>        the band and DOS k-points by the minimum image convention and diagonalize there. That
!>        is why this is a separate routine rather than a shared one: the two
!>        compute_QP_energies cannot both be in scope under the same name.
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_time ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_and_QP_energies(bs_env, qs_env, fm_W_time, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time, fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_and_QP_energies'

      INTEGER                                            :: handle, i_iter, n_iter
      LOGICAL                                            :: converged
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2)                        :: e_fermi_dft
      REAL(KIND=dp), DIMENSION(3, 2)                     :: band_prev
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

      converged = .FALSE.

      CALL init_evGW0(bs_env, n_iter, band_prev, eigenval_scf_Gamma_dft, e_fermi_dft)

      ! evGW0 self-consistency loop; for G0W0 the loop ends after one iteration
      DO i_iter = 1, n_iter

         bs_env%ri_rs%evgw0_i_iter = i_iter

         ! W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         ! Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
         ! Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
         ! Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
         CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                              bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time)

         ! Σ^c_λσ(iτ) -> Σ^c_nn(ϵ) at every band/DOS k-point
         ! ϵ_nk^GW,(i) = ϵ_nk^DFT + Σ^c_nn[G^(i-1),W](ϵ_nk^GW,(i)) + Σ^x_nn - v^xc_nn
         CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

         IF (bs_env%gw_flavour == G0W0) EXIT

         CALL print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

         IF (i_iter == n_iter .OR. converged) EXIT

         ! eigenvalues ϵ_n to be updated in G
         CALL update_eigenvalues_G(bs_env)

      END DO

      CALL cp_fm_release(fm_W_time)

      CALL reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      CALL delete_unnecessary_files(bs_env)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_and_QP_energies

END MODULE gw_ri_rs_large_cell_Gamma
