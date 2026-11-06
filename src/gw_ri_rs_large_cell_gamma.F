!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using the RI-RS approximation for large periodic cells at Gamma
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_ri_rs_large_cell_Gamma
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE gw_compute_V_and_W,              ONLY: get_W_MIC
   USE gw_gamma_correction,             ONLY: gamma_corr_add_to_sigma,&
                                              gamma_corr_build_coulomb_maps,&
                                              gamma_corr_response,&
                                              gamma_corr_tabulate
   USE gw_ri_rs_chi_sigma_kernels,      ONLY: compute_Sigma_c,&
                                              compute_Sigma_x,&
                                              get_mat_chi_Gamma_tau
   USE gw_ri_rs_non_periodic,           ONLY: init_evGW0,&
                                              print_evGW0_band_edges,&
                                              print_ri_rs_memory_estimate,&
                                              reset_and_clean_bs_env,&
                                              update_eigenvalues_G
   USE gw_tensor_large_cell_Gamma,      ONLY: compute_QP_energies,&
                                              delete_unnecessary_files
   USE gw_utils,                        ONLY: de_init_bs_env
   USE input_constants,                 ONLY: G0W0
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_print_mem_per_rank
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE ri_rs_setup,                     ONLY: ensure_ri_rs_representation
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_large_cell_Gamma'

   PUBLIC :: gw_calc_ri_rs_large_cell_Gamma

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for large periodic cells at Gamma
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

      ! Build the common grid, phi and Z representation with the system's boundary conditions.
      CALL ensure_ri_rs_representation(qs_env, bs_env)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
      CALL print_ri_rs_memory_estimate(bs_env)

      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, &
                                 label='Memory per MPI process after computing Z_lP:')

       !!========================================================================
       !! 4. Compute Independent-Particle Polarizability (χ)
       !!    G^occ_µλ(i|τ|)  = sum_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^vir_µλ(i|τ|)  = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^occ_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')
       !!    G^vir_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^vir_µν ϕ_ν(r_l')
       !!    χ_ll'(iτ)       = G^occ_ll'(i|τ|) * G^vir_ll'(i|τ|)
       !!    χ_PQ(iτ)        = sum_ll' Z_lP χ_ll'(iτ) Z_l'Q
       !!========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, cutoff=0.0_dp)

       !!========================================================================
       !! 4b. Gamma-point head-and-wing correction of W (2D slabs and 3D bulk cells whose k-mesh
       !!     for χ, ε, W contains q = 0, where the Coulomb operator between charged auxiliary
       !!     functions is singular).
       !!     2D: c_P, Π = Id - c cᵀ, V_b = Π V(0) Π       -- separate the charge channel
       !!     3D: V(0) = Ewald sum without its G = 0 term  -- finite, no projector
       !!     Y_αβ(iτ), X_lα(iτ)                           -- optical head and Cartesian wings
       !!     B_PQ = δ_PQ + F_densityᵀ χ_PQ(Γ) F_density   -- the regular body, reused as is
       !!     S_αβ = Y - Eᵀ B⁻¹ E,  T_Pα                   -- block inversion in closed form
       !!     K_h(R), K_w(R)                               -- Brillouin-zone cell integrals
       !!     ΔW^c_PQ(iτ) is then added to W^MIC_PQ(iτ) in RI space inside get_W_MIC (atom-pair
       !!     minimum-image cell integrals, gamma_corr_add_to_W), before M⁻¹ W M⁻¹.
       !!     It has to run here, after χ_PQ(Γ,iτ) exists and before get_W_MIC consumes it.
       !!========================================================================
      IF (bs_env%ri_rs%gamma_corr%enabled) THEN
         ! c, Pi, V_b^(+-1/2) and the maps F_density = M^-1 Pi V_b^(1/2),
         ! G_extern = V(0) Pi V_b^(-1/2) (gwp42 Eqs. 63-67; Pi = 1 in 3D), replicated or
         ! distributed depending on the memory estimate.
         CALL gamma_corr_build_coulomb_maps(qs_env, bs_env)
         ! Head S_ab(i omega) and RI wing vector T_Pa(i omega) from the velocity-gauge dipoles
         ! and chi_PQ(Gamma, i tau) (gwp42 Eqs. 48-49, 71-74).
         CALL gamma_corr_response(bs_env, qs_env, bs_env%mat_chi_Gamma_tau)
         ! 2D: cell integrals K_h(R, i omega), K_w(R, i omega) on a 128 x 128 mesh of in-plane
         ! minimum-image displacements, interpolated at atom-pair separations in get_W_MIC.
         ! 3D: per-direction coefficients of the analytic cell integrals (n_tab unused).
         CALL gamma_corr_tabulate(bs_env, n_tab=128)
         ! From here on get_W_MIC may add DeltaW to W^MIC (gamma_corr_add_to_W asserts this).
         bs_env%ri_rs%gamma_corr%ready = .TRUE.
      END IF

       !!========================================================================
       !! 5. Compute Screened Interaction (W^MIC)
       !!    χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W^MIC_PQ(iτ)
       !!========================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

       !!========================================================================
       !! 6. Compute Exact Exchange Self-Energy (Σ^x)
       !!    D_µν          = sum_n^occ C_µn C_νn
       !!    D_ll'         = sum_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
       !!    V^trunc_ll'   = sum_PQ Z_lP V^trunc_PQ Z_l'Q
       !!    Σ^x_ll'       = D_ll' * V^trunc_ll'
       !!    Σ^x_λσ(k=0)   = -sum_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma, cutoff=0.0_dp)

       !!========================================================================
       !! 7. Compute Correlation Self-Energy (Σ^c) and 8. Quasiparticle Energies,
       !!    iterated until eigenvalue self-consistency if &EVGW0 is given and
       !!    done in a single pass for G0W0.
       !!    W^MIC_ll'(iτ) =  sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
       !!    Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) * W^MIC_ll'(iτ), for τ < 0
       !!    Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) * W^MIC_ll'(iτ), for τ > 0
       !!    Σ^c_λσ(iτ)    =  sum_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
       !!    Σ^c_λσ(iτ)   -> Σ^c_nn(ϵ)
       !!    ϵ_nk^GW       = ϵ_nk^DFT + Σ^c_nn(ϵ) + Σ^x_nn - v^xc_nn
       !!========================================================================
      CALL compute_Sigma_c_and_QP_energies(bs_env, qs_env, fm_W_time, fm_Sigma_x_Gamma)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_large_cell_Gamma

! **************************************************************************************************
!> \brief Correlation self-energy and quasiparticle energies, iterated to eigenvalue
!>        self-consistency in G (evGW0).
!>
!>        The loop body is the molecular one of gw_ri_rs_non_periodic; only the quasiparticle
!>        step differs, because the periodic one has to carry the Gamma matrices to the band and
!>        DOS k-points by the minimum image convention and diagonalize there. That is why this
!>        is a separate routine rather than a shared one: the two compute_QP_energies cannot
!>        both be in scope under the same name.
!>
!>        For G0W0 this runs once with G^(0) built from the DFT eigenvalues. For evGW0 the
!>        Green's function is rebuilt from the quasiparticle energies of all states until the
!>        quasiparticle HOMO, LUMO and HOMO-LUMO gap change by less than EPS_ITER between two
!>        cycles, or MAX_ITER cycles are spent. W^MIC, Σ^x and everything computed before this
!>        routine stay frozen, which is what keep_w_and_files=.TRUE. preserves across cycles;
!>        the cleanup it suppresses is done once here instead.
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

      ! evGW0 self-consistency loop; for G0W0 the loop is terminated after one iteration
      DO i_iter = 1, n_iter

         bs_env%evgw0_i_iter = i_iter

         ! W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         ! Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
         ! Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
         ! Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
         CALL compute_Sigma_c(bs_env, fm_W_time, fm_Sigma_c_Gamma_time, cutoff=0.0_dp, &
                              keep_w_and_files=.TRUE.)

         ! 3D Gamma-point correction: head, head-dipole and wings are added to Σ^c in the AO
         ! basis through the orbital-pair kernel Ŵ_nm (gamma_corr_add_to_sigma); the 2D
         ! correction was already added to W^MIC inside get_W_MIC.
         IF (bs_env%ri_rs%gamma_corr%enabled .AND. .NOT. bs_env%ri_rs%gamma_corr%slab) &
            CALL gamma_corr_add_to_sigma(bs_env, qs_env, fm_Sigma_c_Gamma_time)

         ! Σ^c_λσ(iτ) -> Σ^c_nn(ϵ) at every band and DOS k-point
         ! ϵ_nk^GW,(i) = ϵ_nk^DFT + Σ^c_nn[G^(i-1),W](ϵ_nk^GW,(i)) + Σ^x_nn - v^xc_nn
         CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

         IF (bs_env%gw_flavour == G0W0) EXIT

         CALL print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

         IF (i_iter == n_iter .OR. converged) EXIT

         ! eigenvalues ϵ_n to be updated in G:
         ! G^occ_µλ(i|τ|) = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         ! G^vir_µλ(i|τ|) = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         CALL update_eigenvalues_G(bs_env)

      END DO

      CALL cp_fm_release(fm_W_time)

      CALL reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      CALL delete_unnecessary_files(bs_env)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_and_QP_energies

END MODULE gw_ri_rs_large_cell_Gamma
