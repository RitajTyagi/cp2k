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

MODULE gw_ri_rs_non_periodic
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_data_size, dbcsr_get_info, dbcsr_get_occupation, dbcsr_p_type, dbcsr_put_block, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set,&
                                              max_elements_per_block
   USE cp_fm_basic_linalg,              ONLY: cp_fm_add_on_diag,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_uplo_to_full
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE gw_ri_rs_compute_Z_lP,           ONLY: compute_Z_lP
   USE gw_ri_rs_grid_setup_main,        ONLY: setup_ri_rs_grid
   USE gw_ri_rs_panels,                 ONLY: build_G_ao,&
                                              contract_grid_panels,&
                                              contract_grid_panels_sigma_c,&
                                              print_ri_rs_memory_estimate,&
                                              reserve_blocks_within_radius,&
                                              resolve_grid_panels
   USE gw_ri_rs_utils,                  ONLY: evaluate_ao_basis_on_points,&
                                              precompute_ri_rs_radii
   USE gw_tensor_large_cell_Gamma,      ONLY: Fourier_transform_w_to_t,&
                                              compute_fm_chi_Gamma_freq,&
                                              create_fm_W_MIC_time,&
                                              delete_unnecessary_files,&
                                              fill_fm_Sigma_c_Gamma_time,&
                                              fm_write
   USE gw_utils,                        ONLY: analyt_conti_and_print,&
                                              de_init_bs_env
   USE gw_utils_dbcsr,                  ONLY: release_dbcsr_topology_and_matrices,&
                                              setup_square_topology
   USE gw_utils_fm,                     ONLY: fm_contract_ABA,&
                                              fm_invert,&
                                              fm_sqrt
   USE input_constants,                 ONLY: G0W0,&
                                              evGW0,&
                                              rtp_method_bse
   USE kinds,                           ONLY: dp,&
                                              max_line_length
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_print_mem_per_rank
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: get_all_VBM_CBM_bandgaps
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE time_frequency_grids,            ONLY: negative_time,&
                                              positive_time,&
                                              transform_time_to_frequency
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_non_periodic'

   PUBLIC :: gw_calc_ri_rs_non_periodic, atomic_basis_at_grid_point, &
             reserve_blocks_within_radius

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for molecules
!> \param qs_env ...
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************
   SUBROUTINE gw_calc_ri_rs_non_periodic(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_non_periodic'

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
      !     ϕ_μl = ϕ_μ(r_l)   (mat_phi_mu_l)
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
      !      b. D_lP  = Σ_μν ϕ_μ(r_l) ϕ_ν(r_l) (μν|P)
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
                                 label='Memory per MPI process after computing Z_lP:')

      ! ==============================================================================
      ! 4. Polarizability matrix χ on the imaginary-time grid
      !      G^occ_µλ(i|τ|)  = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^vir_µλ(i|τ|)  = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^occ_ll'(i|τ|) = Σ_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')  (G^vir analogous)
      !      χ_ll'(iτ)       = G^occ_ll'(i|τ|) ∘ G^vir_ll'(i|τ|) (element-wise)
      !      χ_PQ(iτ)        = Σ_ll' Z_lP χ_ll'(iτ) Z_l'Q
      ! ==============================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      ! ==============================================================================
      ! 5. Screened Coulomb interaction W (RI basis)
      !      χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W_PQ(iτ)
      ! ==============================================================================
      CALL compute_W(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

      ! ==============================================================================
      ! 6. Exact-exchange self-energy Σ^x
      !      D_µν          = Σ_n^occ C_µn C_νn                (density matrix)
      !      D_ll'         = Σ_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
      !      V^tr_ll'      = Σ_PQ Z_lP V^tr_PQ Z_l'Q          (truncated Coulomb)
      !      Σ^x_ll'       = D_ll' ∘ V^tr_ll'
      !      Σ^x_λσ(k=0)   = -Σ_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
      ! ==============================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

      ! ==============================================================================
      ! 7. Correlation self-energy Σ^c and quasiparticle energies, iterated
      !    until eigenvalue self-consistency if the &EVGW0 section is given and
      !    done in a single pass for G0W0.
      !
      !    (a) W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
      !        Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
      !        Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
      !        Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
      !    (b) Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
      !        ϵ_n^GW,(i) = ϵ_n^DFT + Σ^c_nn[G^(i-1),W](ϵ_n^GW,(i)) + Σ^x_nn - v^xc_nn
      ! ==============================================================================
      CALL compute_Sigma_c_and_QP_energies(bs_env, fm_W_time, fm_Sigma_x_Gamma)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_non_periodic

! **************************************************************************************************
!> \brief Correlation self-energy and quasiparticle energies, iterated to eigenvalue
!>        self-consistency in G (evGW0).
!>
!>        For G0W0 this runs once with G^(0) built from the DFT eigenvalues. For evGW0 the
!>        Green's function is rebuilt from the quasiparticle energies of all states until the
!>        quasiparticle HOMO, LUMO and HOMO-LUMO gap change by less than EPS_ITER between two
!>        cycles, or MAX_ITER cycles are spent. W, Σ^x and everything computed before this
!>        routine stay frozen: they either do not depend on the eigenvalues at all (Σ^x is G
!>        at τ = 0, a pure density matrix) or are held fixed by construction in GW0. The
!>        imaginary-time grid is fixed as well, since W(iτ) lives on it.
!>
!> \param bs_env ...
!> \param fm_W_time ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_and_QP_energies(bs_env, fm_W_time, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
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

      ! evGW0 self-consistency loop; for G0W0, loop is terminated after one iteration
      DO i_iter = 1, n_iter

         bs_env%ri_rs%evgw0_i_iter = i_iter

         ! W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         ! Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
         ! Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
         ! Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
         CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                              bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time)

         ! Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
         ! ϵ_n^GW,(i) = ϵ_n^DFT + Σ^c_nn[G^(i-1),W](ϵ_n^GW,(i)) + Σ^x_nn - v^xc_nn
         CALL compute_QP_energies(bs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

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

! **************************************************************************************************
!> \brief Sets up the evGW0 eigenvalue self-consistency loop: the cycle count, the eigenvalues
!>        the first Green's function is built from, and the DFT reference that the loop
!>        overwrites. A G0W0 run reduces to a single cycle and needs none of it.
!> \param bs_env ...
!> \param n_iter ...
!> \param band_prev ...
!> \param eigenval_scf_Gamma_dft ...
!> \param e_fermi_dft ...
! **************************************************************************************************
   SUBROUTINE init_evGW0(bs_env, n_iter, band_prev, eigenval_scf_Gamma_dft, e_fermi_dft)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(OUT)                               :: n_iter
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(OUT)        :: band_prev
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: e_fermi_dft

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_evGW0'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      n_iter = 1
      band_prev(:, :) = 0.0_dp
      e_fermi_dft(:) = 0.0_dp

      IF (bs_env%gw_flavour == evGW0) THEN
         n_iter = bs_env%ri_rs%evgw0_iter

         ! eigenvalues currently in G; the first cycle starts is a plain G0W0 step
         bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_scf(:, :, :)
         ! the loop overwrites these; post-GW printing expects the DFT values back
         ALLOCATE (eigenval_scf_Gamma_dft, SOURCE=bs_env%eigenval_scf_Gamma)
         e_fermi_dft(:) = bs_env%e_fermi(:)
      END IF

      CALL timestop(handle)

   END SUBROUTINE init_evGW0

! **************************************************************************************************
!> \brief Feeds the quasiparticle energies of the current evGW0 cycle back into the Green's
!>        function used by the next one, and re-centres the Fermi level between the new
!>        quasiparticle HOMO and LUMO.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE update_eigenvalues_G(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'update_eigenvalues_G'

      INTEGER                                            :: handle, i_mo, ispin

      CALL timeset(routineN, handle)

      ! record this cycle's evGW0 result; it is also what the next Green's function is built from
      bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_GW(:, :, :)

      DO ispin = 1, bs_env%n_spin
         ! Update all physical states; exclude linear-dependency placeholders.
         DO i_mo = 1, bs_env%n_mo_retained
            bs_env%eigenval_scf_Gamma(i_mo, ispin) = bs_env%eigenval_GW(i_mo, 1, ispin)
         END DO
         bs_env%e_fermi(ispin) = &
            0.5_dp*(bs_env%eigenval_GW(bs_env%n_occ(ispin), 1, ispin) + &
                    bs_env%eigenval_GW(bs_env%n_occ(ispin) + 1, 1, ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE update_eigenvalues_G

! **************************************************************************************************
!> \brief Restores the DFT reference that the evGW0 loop overwrote and cleanup
!> \param bs_env ...
!> \param eigenval_scf_Gamma_dft ...
!> \param e_fermi_dft ...
!> \param fm_Sigma_x_Gamma ..
! **************************************************************************************************
   SUBROUTINE reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: e_fermi_dft
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reset_and_clean_bs_env'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (bs_env%gw_flavour == evGW0) THEN
         ! Update array with final evGW0 result
         bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_GW(:, :, :)
         ! put the DFT reference back for the post-GW DOS/band-edge printing
         bs_env%eigenval_scf_Gamma(:, :) = eigenval_scf_Gamma_dft(:, :)
         bs_env%e_fermi(:) = e_fermi_dft(:)
         DEALLOCATE (eigenval_scf_Gamma_dft)
         CALL cp_fm_release(fm_Sigma_x_Gamma)
      END IF

      CALL timestop(handle)

   END SUBROUTINE reset_and_clean_bs_env

! **************************************************************************************************
!> \brief Print the quasiparticle HOMO, LUMO and HOMO-LUMO gap of the current evGW0 cycle
!> \param bs_env ...
!> \param band_prev ...
!> \param i_iter ...
!> \param n_iter ...
!> \param converged ...
! **************************************************************************************************
   SUBROUTINE print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(INOUT)      :: band_prev
      INTEGER, INTENT(IN)                                :: i_iter, n_iter
      LOGICAL, INTENT(OUT)                               :: converged

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_evGW0_band_edges'

      INTEGER                                            :: handle, homo, ispin, u
      REAL(KIND=dp)                                      :: max_delta
      REAL(KIND=dp), DIMENSION(3)                        :: band

      CALL timeset(routineN, handle)

      u = bs_env%unit_nr

      converged = (i_iter > 1)
      max_delta = 0.0_dp

      IF (u > 0) THEN
         WRITE (u, '(A)') ' '
         WRITE (u, '(T2,A)') REPEAT('-', 79)
         WRITE (u, '(T2,A,I4,A,I4)') 'evGW0 cycle', i_iter, ' /', n_iter
         WRITE (u, '(T2,A)') REPEAT('-', 79)
      END IF

      DO ispin = 1, bs_env%n_spin

         homo = bs_env%n_occ(ispin)
         band(1) = bs_env%eigenval_GW(homo, 1, ispin)
         band(2) = bs_env%eigenval_GW(homo + 1, 1, ispin)
         band(3) = band(2) - band(1)

         IF (i_iter > 1) THEN
            max_delta = MAX(max_delta, MAXVAL(ABS(band(:) - band_prev(:, ispin))))
         END IF

         IF (u > 0) THEN
            IF (bs_env%n_spin == 2) WRITE (u, '(T2,A,I0)') 'Spin ', ispin
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 HOMO (eV)', band(1)*evolt
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 LUMO (eV)', band(2)*evolt
            WRITE (u, '(T2,A,T61,F20.3)') 'evGW0 HOMO-LUMO gap (eV)', band(3)*evolt
         END IF

         band_prev(:, ispin) = band(:)

      END DO

      IF (i_iter > 1) THEN
         converged = (max_delta < bs_env%ri_rs%evgw0_eps_iter)
         IF (u > 0) WRITE (u, '(T2,A,T61,F20.6)') 'Max. change to previous cycle (eV)', &
            max_delta*evolt
      END IF

      IF (u > 0) WRITE (u, '(T2,A)') REPEAT('-', 79)

      IF (converged) THEN
         IF (u > 0) THEN
            WRITE (u, '(A)') ' '
            WRITE (u, '(T2,A,I4,A)') &
               'evGW0 eigenvalue self-consistency reached in', i_iter, ' cycles.'
            WRITE (u, '(A)') ' '
         END IF
      ELSE IF (i_iter == n_iter) THEN
         CALL cp_warn(__LOCATION__, &
                      "The evGW0 eigenvalue self-consistency cycle did not converge "// &
                      "within MAX_ITER cycles. The reported quasiparticle energies are "// &
                      "those of the last cycle.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_evGW0_band_edges

! **************************************************************************************************
!> \brief Evaluates the AO basis on the RI-RS grid and stores it as the sparse DBCSR matrix
!>        ϕ_μl = ϕ_μ(r_l)  (rows = grid points in atom-aligned blocks of at most
!>        max_elements_per_block points, columns = one block per atom's full AO set).
!>        Grid points outside the reach of an atom's most
!>        diffuse Gaussian (or the CUTOFF_RADIUS_RL_AO) are skipped, and only blocks
!>        with at least one element > eps_filter are stored. This locality is the source of
!>        ALL grid-dimension sparsity used downstream. Also caches the per-chunk centroids
!>        needed by the optional CUTOFF_RADIUS_RL_W operator truncation.
!> \param bs_env ...
!> \param ri_rs_grid_points ...
!> \param mat_phi_mu_l ...
! **************************************************************************************************
   SUBROUTINE atomic_basis_at_grid_point(bs_env, ri_rs_grid_points, mat_phi_mu_l)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: ri_rs_grid_points(:, :)
      TYPE(dbcsr_type), INTENT(OUT)                      :: mat_phi_mu_l

      CHARACTER(LEN=*), PARAMETER :: routineN = 'atomic_basis_at_grid_point'

      INTEGER :: bs_eff, c_size, handle, i, i_blk, ia, iatom, ikind, natom, npcol, nprow, &
         num_grid_chunks, r_end, r_start, remaining, run, safe_max
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_row_start
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, r_blk_sizes, row_dist, sizes_AO
      REAL(KIND=dp)                                      :: r2_threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_col_buffer
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist_ks, dist
      TYPE(gto_basis_set_type), POINTER                  :: ao_basis
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      natom = bs_env%n_atom
      cell => bs_env%ri_rs%cell
      para_env => bs_env%para_env
      particle_set => bs_env%ri_rs%particle_set
      sizes_AO => bs_env%sizes_AO
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(SIZE(sizes_AO) == natom)
      CPASSERT(SIZE(ri_rs_grid_points, 2) == bs_env%ri_rs%n_grid_points)

      ! =========================================================================
      ! 1. SETUP DBCSR MATRIX TOPOLOGY
      ! =========================================================================

      ! B. Define Row Block Sizes: atom-aligned blocks (a block never spans two atoms' grid
      ! runs), each atom's run subdivided into blocks of at most bs_eff points.

      ! Fetch CP2K's default process grid configuration
      CALL dbcsr_get_info(bs_env%mat_ao_ao%matrix, distribution=dbcsr_dist_ks)
      CALL dbcsr_distribution_get(dbcsr_dist_ks, nprows=nprow, npcols=npcol)

      ! Overflow-safe upper bound on the block size (see bs_env%dbcsr_msg_elem_limit).
      safe_max = INT(0.5_dp*REAL(bs_env%dbcsr_msg_elem_limit, dp)* &
                     REAL(MAX(MIN(nprow, npcol), 1), dp)/ &
                     REAL(bs_env%ri_rs%n_grid_points, dp))
      safe_max = MAX(1, safe_max)
      ! Block size = CP2K's global max_elements_per_block (GLOBAL/DBCSR input; default 32),
      ! overflow-capped.
      bs_eff = MAX(1, MIN(max_elements_per_block, safe_max))

      ! Count the atom-aligned blocks, then fill r_blk_sizes and each block's starting grid row.
      num_grid_chunks = 0
      DO ia = 1, natom
         run = bs_env%ri_rs%grid_atom_boundaries(ia + 1) - bs_env%ri_rs%grid_atom_boundaries(ia)
         IF (run > 0) num_grid_chunks = num_grid_chunks + (run + bs_eff - 1)/bs_eff
      END DO
      ALLOCATE (r_blk_sizes(num_grid_chunks), blk_row_start(num_grid_chunks))
      i_blk = 0
      r_start = 1
      DO ia = 1, natom
         remaining = bs_env%ri_rs%grid_atom_boundaries(ia + 1) - bs_env%ri_rs%grid_atom_boundaries(ia)
         DO WHILE (remaining > 0)
            i_blk = i_blk + 1
            r_blk_sizes(i_blk) = MIN(bs_eff, remaining)
            blk_row_start(i_blk) = r_start
            r_start = r_start + r_blk_sizes(i_blk)
            remaining = remaining - r_blk_sizes(i_blk)
         END DO
      END DO

      IF (bs_env%unit_nr > 0) THEN
         ! T71 compensates for the two two-byte Greek characters in the label.
         WRITE (bs_env%unit_nr, '(T2,A,T71,I12)') &
            'RI-RS grid row-blocks of ϕ_μ(r_l)', num_grid_chunks
         WRITE (bs_env%unit_nr, '(T2,A,T69,I12)') 'RI-RS grid points per block (max)', bs_eff
      END IF

      ! Cache per-chunk centroids for the optional CUTOFF_RADIUS_RL_W / CUTOFF_RADIUS_W0 block
      ! truncations. Left unallocated otherwise, so PRESENT(centroids) stays .FALSE. at the
      ! contract_grid_panels* calls.
      IF (bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp .OR. &
          bs_env%ri_rs%cutoff_radius_w0 > 0.0_dp) THEN
         IF (ALLOCATED(bs_env%ri_rs%chunk_centroids)) DEALLOCATE (bs_env%ri_rs%chunk_centroids)
         ALLOCATE (bs_env%ri_rs%chunk_centroids(3, num_grid_chunks))
         DO i_blk = 1, num_grid_chunks
            r_start = blk_row_start(i_blk)
            r_end = r_start + r_blk_sizes(i_blk) - 1
            bs_env%ri_rs%chunk_centroids(1, i_blk) = &
               SUM(ri_rs_grid_points(1, r_start:r_end))/REAL(r_blk_sizes(i_blk), dp)
            bs_env%ri_rs%chunk_centroids(2, i_blk) = &
               SUM(ri_rs_grid_points(2, r_start:r_end))/REAL(r_blk_sizes(i_blk), dp)
            bs_env%ri_rs%chunk_centroids(3, i_blk) = &
               SUM(ri_rs_grid_points(3, r_start:r_end))/REAL(r_blk_sizes(i_blk), dp)
         END DO
      END IF

      ! C. Build Custom Mappings using Round-Robin across the 2D process grid

      ALLOCATE (row_dist(num_grid_chunks))
      DO i = 1, num_grid_chunks
         row_dist(i) = MOD(i - 1, nprow)
      END DO

      ALLOCATE (col_dist(natom))
      DO i = 1, natom
         col_dist(i) = MOD(i - 1, npcol)
      END DO

      ! E. Create the DBCSR Distribution and Initialize the Matrix
      CALL dbcsr_distribution_new(dist, template=dbcsr_dist_ks, &
                                  row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(mat_phi_mu_l, name="phi_val_sparse", dist=dist, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=r_blk_sizes, col_blk_size=sizes_AO)

      ! =========================================================================
      ! 2. STREAM DATA DIRECTLY INTO SPARSE MATRIX
      ! =========================================================================
      ! Iterate over the atoms assigned to this specific MPI rank
      DO iatom = para_env%mepos + 1, natom, para_env%num_pe

         c_size = sizes_AO(iatom)

         ! Allocate a temporary dense buffer just for this specific atom
         ALLOCATE (atom_col_buffer(bs_env%ri_rs%n_grid_points, c_size))
         atom_col_buffer = 0.0_dp

         ! Evaluate the basis functions on the grid. Skip grid points outside
         ! the spatial extent of the most diffuse AO Gaussian on iatom; beyond
         ! that radius the contribution is guaranteed below eps_filter. A positive
         ! CUTOFF_RADIUS_RL_AO overrides this with a user-defined hard cutoff.
         IF (bs_env%ri_rs%cutoff_radius_ri_ao > 0.0_dp) THEN
            r2_threshold = bs_env%ri_rs%cutoff_radius_ri_ao**2
         ELSE
            r2_threshold = bs_env%ri_rs%radius_ao_per_atom(iatom)**2
         END IF
         ikind = particle_set(iatom)%atomic_kind%kind_number
         ao_basis => bs_env%basis_set_AO(ikind)%gto_basis_set
         CALL evaluate_ao_basis_on_points(atom_col_buffer, ri_rs_grid_points, ao_basis, &
                                          particle_set(iatom)%r, cell, cutoff_squared=r2_threshold)

         ! Slice the dense column into the atom-aligned grid row-blocks and insert into DBCSR
         DO i_blk = 1, num_grid_chunks
            r_start = blk_row_start(i_blk)
            r_end = r_start + r_blk_sizes(i_blk) - 1

            ! Apply dynamic sparsity filtering: Only store blocks with physical significance
            IF (MAXVAL(ABS(atom_col_buffer(r_start:r_end, 1:c_size))) > bs_env%eps_filter) THEN
               CALL dbcsr_put_block(mat_phi_mu_l, row=i_blk, col=iatom, &
                                    block=atom_col_buffer(r_start:r_end, 1:c_size))
            END IF
         END DO

         DEALLOCATE (atom_col_buffer)

      END DO

      CALL dbcsr_finalize(mat_phi_mu_l)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
      CALL print_matrix_occupation(mat_phi_mu_l, 'ϕ_μ(r_l)', bs_env)

      ! -------------------------------------------------------------------------
      ! CLEANUP
      ! -------------------------------------------------------------------------
      DEALLOCATE (r_blk_sizes, row_dist, col_dist, blk_row_start)
      CALL dbcsr_distribution_release(dist)

      CALL timestop(handle)

   END SUBROUTINE atomic_basis_at_grid_point

! **************************************************************************************************
!> \brief Computes χ_PQ(iτ) from the occupied and virtual Green's functions.
!> \param bs_env ...
!> \param mat_chi_Gamma_tau Response matrices χ_PQ(iτ) to be computed
!> \param mat_phi_mu_l AO values ϕ_μ(r_l) on the RI-RS grid
!> \param mat_Z_lP RI-RS fitting matrix Z_lP
! **************************************************************************************************
   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, mat_chi_Gamma_tau, mat_phi_mu_l, mat_Z_lP)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER                                            :: handle, i_t, ispin, n_panels
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      REAL(KIND=dp)                                      :: grid_occ, t1, tau
      TYPE(dbcsr_type)                                   :: matrix_G_occ_ao, matrix_G_vir_ao

      CALL timeset(routineN, handle)

      ! Panel boundaries for the grid-streaming contraction.
      ! The panels are identical for χ, Σ^x and Σ^c, so the count is reported once here for all three stages.

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

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
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_ao)
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_ao)

            ! χ_PQ += g_s · Z^T ( (ϕ G^occ ϕ^T) ∘ (ϕ G^vir ϕ^T) ) Z
            CALL contract_grid_panels(L_A=mat_phi_mu_l, M_A=matrix_G_occ_ao, &
                                      L_B=mat_phi_mu_l, M_B=matrix_G_vir_ao, &
                                      L_out=mat_Z_lP, mat_out=mat_chi_Gamma_tau(i_t)%matrix, &
                                      scale=bs_env%spin_degeneracy, eps=bs_env%eps_filter, &
                                      para_env=bs_env%para_env, &
                                      pan_first=pan_first, pan_last=pan_last, &
                                      lb_eq_la=.TRUE., lout_eq_la=.FALSE., &
                                      zero_out=(ispin == 1), &
                                      centroids=bs_env%ri_rs%chunk_centroids, &
                                      cutoff=bs_env%ri_rs%cutoff_radius_v_w, &
                                      grid_occupation=grid_occ)

            CALL dbcsr_release(matrix_G_occ_ao)
            CALL dbcsr_release(matrix_G_vir_ao)

         END DO ! ispin

         ! Sparsity reports
         IF (i_t == 1) THEN
            CALL print_matrix_occupation(mat_Z_lP, 'Z_lP', bs_env)
            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,T73,F7.2,A)') &
                  'Percentage of non-zero matrix elements in G_ll'', χ_ll'', W_ll''', &
                  grid_occ*100.0_dp, ' %'
               CALL m_flush(bs_env%unit_nr)
            END IF
            CALL print_matrix_occupation(mat_chi_Gamma_tau(i_t)%matrix, 'χ_PQ', bs_env, &
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
!> \brief Computes the screened Coulomb interaction on the imaginary-time grid, entirely in the
!>        RI auxiliary (PQ) basis:
!>          χ_PQ(iω) = Σ_τ w(ω,τ) cos(ωτ) χ_PQ(iτ)                  (cosine transform)
!>          ε(iω)    = Id - V^0.5 M^-1 χ(iω) M^-1 V^0.5             (dielectric function)
!>          W(iω)    = V^0.5 ( ε^-1(iω) - Id ) V^0.5                (correlation part only)
!>          W(iτ)    = Σ_ω w̃(τ,ω) cos(ωτ) W(iω)                    (back transform)
!>          W(iτ)   <- M^-1 W(iτ) M^-1                              (fold in the RI metric)
!>        where V is the bare Coulomb matrix and M the RI metric.
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

      ! M^-1(k=0) W^MIC(iτ) M^-1(k=0) -> fm_W_time
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_W_time)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
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
                                     bs_env%time_frequency_grid%time_weights_at_zero_frequency(i_t), fm_W_time(i_t))
         END DO
         ! Done, save to file
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
         ! Report calculation
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed W(0),', ' Execution time', m_walltime() - t1, ' s'
         END IF
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W

! **************************************************************************************************
!> \brief Computes V, V^0.5, and M^-1 V^0.5 for the RI-RS dielectric function.
!>        The Coulomb matrix V is constructed by the RI-RS k-point path. The inverse metric
!>        M^-1(k=0) is precomputed once in gw_utils and stored in bs_env.
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V Coulomb matrix V(k=0)
!> \param fm_V_sqrt symmetric factor V^0.5
!> \param fm_Minv_Vsqrt product M^-1(k=0) V^0.5
! **************************************************************************************************
   SUBROUTINE compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_V, fm_V_sqrt, fm_Minv_Vsqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_MinvVsqrt'

      INTEGER                                            :: handle
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (bs_env%auto_ri%enabled) THEN
         ! -------------------------------------------------------------------
         ! 1a. The optimized AB functions span two atoms, so their previously transformed
         !     Coulomb matrix cannot be rebuilt by the atom-local k-point integral routine.
         ! -------------------------------------------------------------------
         CALL cp_fm_to_fm(bs_env%auto_ri%V_pq, fm_V)
      ELSE
         ! -------------------------------------------------------------------
         ! 1b. Build Coulomb Matrix V(k=0) using the kp-routine but only for ikp=1.
         ! -------------------------------------------------------------------
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                         qs_kind_set=qs_kind_set)
         atomic_kind_set => bs_env%ri_rs%atomic_kind_set

         ALLOCATE (mat_V_kp(1:1, 1:2))
         NULLIFY (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         ALLOCATE (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         CALL dbcsr_create(mat_V_kp(1, 1)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_set(mat_V_kp(1, 1)%matrix, 0.0_dp)
         CALL dbcsr_create(mat_V_kp(1, 2)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 2)%matrix)
         ! The dummy imaginary part is required only by the k-point routine interface.
         CALL dbcsr_set(mat_V_kp(1, 2)%matrix, 0.0_dp)

         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig
         CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, "RI_AUX", cell, &
                                         particle_set, qs_kind_set, atomic_kind_set, &
                                         bs_env%size_lattice_sum_V, operator_coulomb, 1, 1)
         CALL copy_dbcsr_to_fm(mat_V_kp(1, 1)%matrix, fm_V)

         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 2)%matrix)
         DEALLOCATE (mat_V_kp)
      END IF

      ! -----------------------------------------------------------------------
      ! 2. V -> V^0.5.
      ! -----------------------------------------------------------------------
      CALL fm_sqrt(fm_V, fm_V_sqrt, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! -----------------------------------------------------------------------
      ! 3. M^-1(k=0) V^0.5.
      ! -----------------------------------------------------------------------
      CALL parallel_gemm("N", "T", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI, 1.0_dp, &
                         bs_env%fm_Minv_Gamma, fm_V_sqrt, 0.0_dp, fm_Minv_Vsqrt)

      CALL timestop(handle)

   END SUBROUTINE compute_V_MinvVsqrt

! **************************************************************************************************
!> \brief Computes the screened interaction at one imaginary frequency:
!>          ε(iω_j)  = Id - (M^-1 V^0.5)^T χ(iω_j) (M^-1 V^0.5)
!>          W(iω_j)  = V^0.5^T ( ε^-1(iω_j) - Id ) V^0.5
!>        ε is inverted via Cholesky; if that fails due to conditioning, via
!>        eigendecomposition (cp_fm_power) with eigenvalue filtering.
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

      INTEGER                                            :: handle, n_RI
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
      CALL cp_fm_add_on_diag(fm_eps_freq_j, 1.0_dp)

      ! Force perfect symmetry before Cholesky to avoid info != 0 due to GEMM noise
      CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)

      ! -----------------------------------------------------------------------
      ! 2. W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      ! -----------------------------------------------------------------------

      ! a) Invert ε by Cholesky decomposition or, if that fails, by diagonalization.
      CALL fm_invert(fm_eps_freq_j, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! b) ε^-1(iω_j) - Id
      CALL cp_fm_add_on_diag(fm_eps_freq_j, -1.0_dp)

      ! c) W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      CALL fm_contract_ABA(fm_V_sqrt, fm_eps_freq_j, fm_W_freq_j)

      ! Cleanup
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_eps_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_freq

! **************************************************************************************************
!> \brief Computes the exact-exchange part of the GW self-energy:
!>          D_μν        = Σ_n^occ C_μn C_νn                     (density matrix = G^occ at τ=0)
!>          V^tr_PQ     = M^-1 (P|Q)_trunc M^-1                 (truncated Coulomb, RI basis)
!>          Σ^x_λσ(k=0) = -Σ_ll' ϕ_λ(r_l) [ (ϕ D ϕ^T)_ll' ∘ (Z V^tr Z^T)_ll' ] ϕ_σ(r_l')
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(bs_env, qs_env, mat_phi_mu_l, mat_Z_lP, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, ispin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux
      TYPE(dbcsr_type)                                   :: mat_Sigma_x_Gamma, matrix_D_ao, &
                                                            matrix_V_aux

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      CALL dbcsr_create(mat_Sigma_x_Gamma, template=bs_env%mat_ao_ao%matrix)

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

      ! =========================================================================
      ! 1. COMPUTE V^tr_PQ (RI x RI)
      ! =========================================================================
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      IF (bs_env%auto_ri%enabled) THEN
         ALLOCATE (fm_Vtr_Gamma(1, 1))
         CALL cp_fm_create(fm_Vtr_Gamma(1, 1), bs_env%fm_RI_RI%matrix_struct)
         CALL cp_fm_to_fm(bs_env%auto_ri%V_pq, fm_Vtr_Gamma(1, 1))
      ELSE
         CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI%matrix_struct, bs_env%n_RI, &
                                 bs_env%trunc_coulomb)
      END IF

      ! M^-1(k=0) V^tr(τ) M^-1(k=0) -> fm_Vtr_Gamma
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(matrix_V_aux, bs_env%eps_filter)

      ! =========================================================================
      ! 2. SPIN LOOP FOR EXACT EXCHANGE
      ! Σ^x_λσ = -Σ_ll' ϕ_λ(r_l) ( D_ll' V^tr_ll' ) ϕ_σ(r_l')
      !        = -ϕ^T ( (ϕ D ϕ^T) ∘ (Z V^tr Z^T) ) ϕ
      ! =========================================================================
      DO ispin = 1, bs_env%n_spin

         ! AO-space density matrix D_µν = G^occ at τ = 0
         CALL build_G_ao(bs_env, 0.0_dp, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_D_ao)

         CALL contract_grid_panels(L_A=mat_phi_mu_l, M_A=matrix_D_ao, &
                                   L_B=mat_Z_lP, M_B=matrix_V_aux, &
                                   L_out=mat_phi_mu_l, mat_out=mat_Sigma_x_Gamma, &
                                   scale=1.0_dp, eps=bs_env%eps_filter, &
                                   para_env=bs_env%para_env, &
                                   pan_first=pan_first, pan_last=pan_last, &
                                   lb_eq_la=.FALSE., lout_eq_la=.TRUE., zero_out=.TRUE., &
                                   centroids=bs_env%ri_rs%chunk_centroids, &
                                   cutoff=bs_env%ri_rs%cutoff_radius_v_w)
         CALL dbcsr_scale(mat_Sigma_x_Gamma, -1.0_dp)

         CALL dbcsr_release(matrix_D_ao)

         ! Data I/O and Export to CP2K Full Matrices
         CALL copy_dbcsr_to_fm(mat_Sigma_x_Gamma, fm_Sigma_x_Gamma(ispin))

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
      CALL dbcsr_release(mat_Sigma_x_Gamma)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy on the imaginary-time grid:
!>          Σ^c_λσ(iτ<0) = -Σ_ll' ϕ_λ(r_l) [ (ϕ G^occ ϕ^T)_ll' ∘ (Z W^MIC Z^T)_ll' ] ϕ_σ(r_l')
!>          Σ^c_λσ(iτ>0) = +Σ_ll' ϕ_λ(r_l) [ (ϕ G^vir ϕ^T)_ll' ∘ (Z W^MIC Z^T)_ll' ] ϕ_σ(r_l')
!> \param bs_env ...
!> \param fm_W_time ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c(bs_env, fm_W_time, mat_phi_mu_l, mat_Z_lP, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      INTEGER                                            :: handle, i_t, ispin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, dist_row_aux
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbcsr_type)                                   :: matrix_G_occ_ao, matrix_G_vir_ao, &
                                                            matrix_W_aux

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP AUXILIARY TOPOLOGY AND PRE-ALLOCATE OUTPUT ARRAYS
      ! =========================================================================
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL resolve_grid_panels(bs_env, mat_phi_mu_l, pan_first, pan_last)

      ! Pre-allocate local DBCSR matrices to act as targets for final output
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

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            ! AO-space Green's functions G^occ_µν, G^vir_µν (dense AO x AO, small)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_ao)
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_ao)

            ! Σ^c_neg and Σ^c_pos in a single panel loop: W_pan = Z_panel × W × Z^T built once
            CALL contract_grid_panels_sigma_c(mat_phi=mat_phi_mu_l, mat_Z=mat_Z_lP, &
                                              mat_G_occ_ao=matrix_G_occ_ao, &
                                              mat_G_vir_ao=matrix_G_vir_ao, &
                                              mat_W_aux=matrix_W_aux, &
                                              mat_Sigma_neg=mat_Sigma_neg_tau(i_t, ispin)%matrix, &
                                              mat_Sigma_pos=mat_Sigma_pos_tau(i_t, ispin)%matrix, &
                                              eps=bs_env%eps_filter, &
                                              para_env=bs_env%para_env, &
                                              pan_first=pan_first, pan_last=pan_last, &
                                              centroids=bs_env%ri_rs%chunk_centroids, &
                                              cutoff=bs_env%ri_rs%cutoff_radius_v_w)
            CALL dbcsr_scale(mat_Sigma_neg_tau(i_t, ispin)%matrix, -1.0_dp)

            CALL dbcsr_release(matrix_G_occ_ao)
            CALL dbcsr_release(matrix_G_vir_ao)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I15,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ) for time point', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

         CALL dbcsr_release(matrix_W_aux)

      END DO ! i_t

      ! -------------------------------------------------------------------------
      ! 3. FINALIZE AND CLEANUP
      ! -------------------------------------------------------------------------
      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      ! fm_W_time and the scratch files are released by the caller: in an evGW0 cycle this
      ! routine is entered once per iteration and both have to survive until it is done.
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief Builds the DBCSR distribution.
!> \param matrix_template ...
!> \param square_dist ...
!> \param blk_sizes ...
!> \param mapped_dist ...
! **************************************************************************************************
!> \brief Releases a distribution created by setup_square_topology.
!> \param dist ...
!> \param mapped_dist ...
! **************************************************************************************************
!> \brief Σ^c_λσ(iτ) -> Σ^c_nn(ϵ) and the quasi-particle levels of the non-periodic RI-RS path,
!>        ϵ_n^GW = ϵ_n^DFT + Σ^c_nn(ϵ_n^GW) + Σ^x_nn - v^xc_nn.
!> \param bs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(bs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: Sigma_c_n_freq
      INTEGER                                            :: handle, ispin, j_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_n, V_xc_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_n_time
      TYPE(cp_fm_type)                                   :: fm_ks, fm_mos, fm_s, fm_work

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_ks, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_s, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_mos, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_work, bs_env%fm_s_Gamma%matrix_struct)

      ALLOCATE (V_xc_n(bs_env%n_ao), Sigma_x_n(bs_env%n_ao))
      ALLOCATE (Sigma_c_n_time(bs_env%n_ao, bs_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_n_freq(bs_env%n_ao, bs_env%num_freq_points_fit))

      DO ispin = 1, bs_env%n_spin

         ! 1. Roothaan-Hall H^KS_µν C_νn = S_µν C_νn ϵ_n
         CALL cp_fm_to_fm(bs_env%fm_ks_Gamma(ispin), fm_ks)
         CALL cp_fm_to_fm(bs_env%fm_s_Gamma, fm_s)
         CALL cp_fm_geeig(fm_ks, fm_s, fm_mos, bs_env%eigenval_scf(:, 1, ispin), fm_work)

         ! 2. v^xc_µν -> v^xc_nn  and  Σ^x_µν -> Σ^x_nn
         CALL to_Gamma_and_mo_real(V_xc_n, bs_env%fm_V_xc_Gamma(ispin), fm_mos)
         CALL to_Gamma_and_mo_real(Sigma_x_n, fm_Sigma_x_Gamma(ispin), fm_mos)

         ! 3. Σ^c_µν(+/-i|τ_j|) -> Σ^c_nn(+/-i|τ_j|)
         DO j_t = 1, bs_env%num_time_freq_points
            CALL to_Gamma_and_mo_real(Sigma_c_n_time(:, j_t, positive_time), &
                                      fm_Sigma_c_Gamma_time(j_t, positive_time, ispin), fm_mos)
            CALL to_Gamma_and_mo_real(Sigma_c_n_time(:, j_t, negative_time), &
                                      fm_Sigma_c_Gamma_time(j_t, negative_time, ispin), fm_mos)
         END DO

         ! 4. Σ^c_nn(iτ) -> Σ^c_nn(iω)
         CALL transform_time_to_frequency(bs_env%time_frequency_grid, Sigma_c_n_time, Sigma_c_n_freq)
         Sigma_c_n_freq(1:bs_env%n_occ(ispin), :) = CONJG(Sigma_c_n_freq(1:bs_env%n_occ(ispin), :))

         ! 5. Analytic continuation Σ^c_nn(iω) -> Σ^c_nn(ϵ) and the QP levels
         CALL analyt_conti_and_print(bs_env, Sigma_c_n_freq, Sigma_x_n, V_xc_n, &
                                     bs_env%eigenval_scf(:, 1, ispin), 1, ispin)

      END DO ! ispin

      CALL get_all_VBM_CBM_bandgaps(bs_env)

      IF (bs_env%gw_flavour == G0W0) CALL cp_fm_release(fm_Sigma_x_Gamma)
      CALL cp_fm_release(fm_Sigma_c_Gamma_time)

      CALL cp_fm_release(fm_ks)
      CALL cp_fm_release(fm_s)
      CALL cp_fm_release(fm_mos)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

! **************************************************************************************************
!> \brief AO -> MO transform of a Γ-point matrix
!> \param array_n ...
!> \param fm_Gamma ...
!> \param fm_mos ...
! **************************************************************************************************
   SUBROUTINE to_Gamma_and_mo_real(array_n, fm_Gamma, fm_mos)

      REAL(KIND=dp), DIMENSION(:)                        :: array_n
      TYPE(cp_fm_type)                                   :: fm_Gamma, fm_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'to_Gamma_and_mo_real'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: fm_mo

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mo, fm_Gamma%matrix_struct)

      ! A_nn' = Σ_μν C_μn A_μν C_νn'
      CALL fm_contract_ABA(fm_mos, fm_Gamma, fm_mo)

      CALL cp_fm_get_diag(fm_mo, array_n)

      CALL cp_fm_release(fm_mo)

      CALL timestop(handle)

   END SUBROUTINE to_Gamma_and_mo_real

! **************************************************************************************************
!> \brief Counts the AO functions whose radial support intersects an RI fitting sphere.
!> \param bs_env ...
!> \param atom_P Atom at the center of the RI fitting sphere
!> \param cutoff_ri Radius of the RI fitting sphere
!> \param n_ao_used Number of intersecting AO functions
! **************************************************************************************************
!> \brief Prints the percentage of non-zero elements in a distributed RI-RS matrix.
!> \param matrix Distributed matrix whose occupation is reported
!> \param label Mathematical matrix label used in the output
!> \param bs_env ...
!> \param suffix Optional text appended to the matrix label
! **************************************************************************************************
   SUBROUTINE print_matrix_occupation(matrix, label, bs_env, suffix)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: suffix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_matrix_occupation'

      CHARACTER(LEN=32)                                  :: output_format
      CHARACTER(LEN=max_line_length)                     :: msg, output_label
      INTEGER                                            :: handle, i, unicode_shift
      REAL(KIND=dp)                                      :: frac_2p31, max_loc, occ

      CALL timeset(routineN, handle)

      occ = dbcsr_get_occupation(matrix)
      max_loc = REAL(dbcsr_get_data_size(matrix), dp)
      CALL bs_env%para_env%max(max_loc)

      IF (bs_env%unit_nr > 0) THEN
         frac_2p31 = max_loc/REAL(bs_env%dbcsr_msg_elem_limit, dp)
         output_label = 'Percentage of non-zero matrix elements in '//TRIM(label)
         IF (PRESENT(suffix)) output_label = TRIM(output_label)//TRIM(suffix)
         ! Fortran counts UTF-8 bytes, whereas the terminal displays each Greek letter in one
         ! column. Shift the absolute output tab once for every continuation byte.
         unicode_shift = 0
         DO i = 1, LEN_TRIM(output_label)
            IF (IAND(IACHAR(output_label(i:i)), 192) == 128) unicode_shift = unicode_shift + 1
         END DO
         WRITE (output_format, '(A,I0,A)') '(T2,A,T', 72 + unicode_shift, ',F7.2,A)'
         WRITE (bs_env%unit_nr, output_format) TRIM(output_label), occ*100.0_dp, ' %'
         IF (frac_2p31 > 0.5_dp) THEN
            WRITE (msg, '(3A,F0.2,A)') &
               "The largest per-rank message of ", TRIM(label), " reaches ", frac_2p31, &
               " of the 32-bit limit that DBCSR uses for its message length. Beyond it the "// &
               "length overflows and multiply_cannon fails. Reduce the per-rank block size, "// &
               "for instance with more MPI ranks."
            CPWARN(TRIM(msg))
         END IF
         CALL m_flush(bs_env%unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_matrix_occupation

END MODULE gw_ri_rs_non_periodic
