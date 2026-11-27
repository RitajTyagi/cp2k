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
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_release,&
                                              cp_fm_struct_retain,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output
   USE gw_compute_V_and_W,              ONLY: get_W_MIC
   USE gw_gamma_correction,             ONLY: gamma_corr_build_coulomb_maps,&
                                              gamma_corr_classify_atoms,&
                                              gamma_corr_n_tab_rem,&
                                              gamma_corr_response,&
                                              gamma_corr_sigma_shifts,&
                                              gamma_corr_tabulate,&
                                              gamma_corr_tail_operator
   USE gw_nevanlinna,                   ONLY: gw_qp_summary_type
   USE gw_ri_rs_chi_sigma_kernels,      ONLY: compute_Sigma_c,&
                                              compute_Sigma_x,&
                                              get_mat_chi_Gamma_tau,&
                                              sigma_mo_diag_at_k,&
                                              sigma_mo_diag_release,&
                                              sigma_mo_diag_type
   USE gw_ri_rs_non_periodic,           ONLY: init_evGW0,&
                                              print_evGW0_band_edges,&
                                              print_ri_rs_memory_estimate,&
                                              reset_and_clean_bs_env,&
                                              update_eigenvalues_G
   USE gw_tensor_large_cell_Gamma,      ONLY: delete_unnecessary_files
   USE gw_utils,                        ONLY: analyt_conti_and_print,&
                                              de_init_bs_env
   USE input_constants,                 ONLY: G0W0
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_mem_used_per_rank_GB,&
                                              mp_print_mem_per_rank
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: cfm_ikp_from_fm_Gamma,&
                                              get_all_VBM_CBM_bandgaps
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE ri_rs_grid_initialization,       ONLY: rs_ao_ratio_of_kind
   USE ri_rs_setup,                     ONLY: ensure_ri_rs_representation
   USE time_frequency_grids,            ONLY: transform_time_to_frequency
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_large_cell_Gamma'

   PUBLIC :: gw_calc_ri_rs_large_cell_Gamma

CONTAINS

! **************************************************************************************************
!> \brief G0W0 / evGW0 with the RI-RS approximation for large periodic cells at the Γ point, with
!>        the Γ-point head-and-wing correction (2D slabs, 3D bulk); at the band/DOS k-points
!>          ϵ_nk^GW = ϵ_nk^DFT + Re Σ^c_nn(k, ϵ_nk^GW) + Σ^x_nn(k) − v^xc_nn(k).
!>        M⁺ is the relative-cut pseudo-inverse of the RI metric at Γ.
!> \param qs_env QS environment (cell, basis sets, integrals, input)
!> \param bs_env band-structure environment containing the GW parameters, MOs and grids
! **************************************************************************************************

   SUBROUTINE gw_calc_ri_rs_large_cell_Gamma(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_large_cell_Gamma'
      INTEGER, PARAMETER                                 :: n_tab_slab = 128

      INTEGER                                            :: handle, i_work
      REAL(KIND=dp)                                      :: dense_GB, panels_GB
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ao
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time

      CALL timeset(routineN, handle)

      fm_struct_ao => bs_env%fm_work_mo(1)%matrix_struct
      CALL cp_fm_struct_retain(fm_struct_ao)
      DO i_work = 2, SIZE(bs_env%fm_work_mo)
         CALL cp_fm_release(bs_env%fm_work_mo(i_work))
      END DO
      CALL cp_cfm_release(bs_env%cfm_work_mo)
      CALL cp_cfm_release(bs_env%cfm_work_mo_2)
      CALL cp_fm_release(bs_env%fm_Gocc)
      CALL cp_fm_release(bs_env%fm_Gvir)
      CALL cp_fm_release(bs_env%fm_RI_RI)

      CALL gamma_only_memory_estimate(qs_env, bs_env, dense_GB)

      ! ==============================================================================
      ! 1.-3. RI-RS grid r_l, AOs on the grid and fit coefficients Z_lP (periodic cell)
      !      ϕ_µ(r_l)  = Σ_R φ_µ(r_l − A_µ − R)                   (Γ Bloch sums)
      !      D_ll'     = [Σ_µ ϕ_µ(r_l) ϕ_µ(r_l')]²
      !      d_lP      = Σ_µν ϕ_µ(r_l) ϕ_ν(r_l) (µν|P)_m           (RI metric m)
      !      Σ_l' D_ll' Z_l'P = d_lP                               (Tikhonov-regularized)
      !    and the Γ-point correction set-up (gamma_corr_init)
      ! ==============================================================================
      CALL ensure_ri_rs_representation(qs_env, bs_env)
      CALL print_ri_rs_memory_estimate(bs_env, panels_GB=panels_GB)

      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, &
                                 label='Memory per MPI process after computing Z_lP:')

      ! Peak of the whole run: the dense terms of the governing phase + the memory in use after
      ! the Z_lP fit (grids, ϕ, Z, setup matrices) + the χ/W/Σ panel term.
      CALL gamma_only_peak_estimate(bs_env, dense_GB, panels_GB)

      ! ==============================================================================
      ! 4. Polarizability matrix χ on the imaginary-time grid
      !      G^occ_µλ(i|τ|)  = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^vir_µλ(i|τ|)  = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
      !      G^occ_ll'(i|τ|) = Σ_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')  (G^vir analogous)
      !      χ_ll'(iτ)       = G^occ_ll'(i|τ|) ∘ G^vir_ll'(i|τ|) (element-wise)
      !      χ_PQ(iτ)        = Σ_ll' Z_lP χ_ll'(iτ) Z_l'Q
      ! ==============================================================================
      ! CUTOFF_RADIUS_RL_W > 0: only grid-block pairs with |c_b − c_b'|_mic ≤ r_c enter the
      ! point-wise products of χ, Σ^x and Σ^c (not the dense W_PQ(iτ)); otherwise all pairs.
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, cutoff=bs_env%ri_rs%cutoff_radius_v_w)

      ! ==============================================================================
      ! 5. Γ-point head-and-wing correction
      !      V(0)      = Ewald sum without its G = 0 term (2D Ewald sum for slabs)
      !      ĉ = c/|c|, c_P = ∫ϕ_P;  Π = 1 − ĉ ĉᵀ (2D), Π = 1 (3D);  V_b = Π V(0) Π
      !      F_density = M⁺_c Π V_b^(1/2),   G_extern = V(0) Π V_b^(-1/2)
      !      µ^α_ai    = ⟨a|[H,r^α]|i⟩/(ϵ_a − ϵ_i)          (velocity-gauge dipoles)
      !      Y_αβ(iτ)  = (2/A) Σ_ia µ^α_ai µ^β_ai e^(-(ϵ_a-ϵ_i)τ)
      !      X_lα(iτ)  = (2/√A) Σ_ia ψ_i(r_l) ψ_a(r_l) µ^α_ai e^(-(ϵ_a-ϵ_i)τ)
      !      B(iω)     = 1 + F_densityᵀ χ(Γ,iω) F_density,   E(iω) = F_densityᵀ Zᵀ X(iω)
      !      S_αβ(iω)  = Y_αβ − [Eᵀ B⁻¹ E]_αβ,               T_Pα(iω) = [G_extern B⁻¹ E]_Pα
      !      2D: K_h, K_w,α, K_h1, K_h2, K_w1 (R, iω), integrals over the Wigner-Seitz cell
      !          of the reciprocal supercell lattice, tabulated on in-plane displacements R
      !      3D: head constant w(iω), A_αβ(iω), B_αβ(iω), 1/ε0(iω) − 1 = −4πā/(1 + 4πā),
      !          ā = tr S/3;  tail operator K = V^tr − V^E − ξ c cᵀ
      !    (A: cell area in 2D, volume Ω in 3D; α, β Cartesian; i occupied, a empty.)
      ! ==============================================================================
      ! check_gw_input admits PERIODIC XY or XYZ only and setup_kpoints_chi_eps_W sets the
      ! single-point mesh, so gamma_corr_init has switched the correction on.
      CPASSERT(bs_env%ri_rs%gamma_corr%enabled)
      ! V(0), c, Π, V_b^(±1/2) and the maps F_density = M⁺_c Π V_b^(1/2),
      ! G_extern = V(0) Π V_b^(-1/2) (Π = 1 in 3D), all as distributed n_RI x n_RI matrices.
      CALL gamma_corr_build_coulomb_maps(qs_env, bs_env)
      ! Head S_αβ(iω) = Y_αβ − [Eᵀ B⁻¹ E]_αβ, RI wing vector T_Pα(iω) = [G_extern B⁻¹ E]_Pα and
      ! the finite-q pair weights that give the head its |q| dependence.
      CALL gamma_corr_response(bs_env, qs_env, bs_env%mat_chi_Gamma_tau)
      ! 2D: cell integrals K_h(R,iω), K_w,α(R,iω) on an n_tab_slab x n_tab_slab mesh of in-plane
      ! minimum-image displacements R; 3D: head constant and the A, B tensors by quadrature.
      CALL gamma_corr_tabulate(bs_env, n_tab=n_tab_slab, qs_env=qs_env)
      ! From here on get_W_MIC may add ΔW to W^MIC (gamma_corr_add_folded_to_W asserts this).
      bs_env%ri_rs%gamma_corr%ready = .TRUE.
      ! 3D sphere-tail operator K = V^tr − V^E − ξ c cᵀ (built while gc%fm_V_zero exists), used
      ! in get_W_MIC. A W read from file lacks the tail while the B:dd term of Σ^c is dropped.
      IF (.NOT. bs_env%ri_rs%gamma_corr%slab .AND. bs_env%ri_rs%gamma_corr%tail_3d) THEN
         IF (bs_env%all_W_exist) THEN
            CPABORT('TAIL_3D: W^MIC would be read from file; delete the W files and rerun')
         END IF
         CALL gamma_corr_tail_operator(qs_env, bs_env, bs_env%ri_rs%gamma_corr%fm_K_tail)
      END IF

      ! ==============================================================================
      ! 6. Screened Coulomb interaction, real at k = 0 and folded with M⁺ (W_f = W^MIC_PQ)
      !      χ_PQ(iω_k) = Σ_j γ^cos_kj χ_PQ(iτ_j)
      !      ε(iω)      = 1 + V(0)^(1/2) M⁺ χ(iω) M⁺ V(0)^(1/2)
      !      W^c(iω)    = V(0)^(1/2) (ε⁻¹(iω) − 1) V(0)^(1/2)
      !      W_f(iτ_j)  = M⁺ [Σ_k γ̃^cos_jk W^c(iω_k)] M⁺
      !    with the Γ-point correction added in the same folded representation:
      !      2D: W_f(iτ) += M⁺_c ΔW(iτ) M⁺_c = F K(iτ) Fᵀ, the head and wings
      !            ΔW_PQ = (1/A) c_P c_Q K_h(R_PQ) + (1/√A) Σ_α [T_Pα c_Q − c_P T_Qα] K_w,α(R_PQ)
      !                    + the out-of-plane moment terms with K_h1, K_h2, K_w1
      !      3D: W_f(iτ) += [Σ_k γ̃^cos(τ, ω_k) (1/ε0(iω_k) − 1)] M⁺ K M⁺   (sphere tail);
      !          head and wings enter Σ^c in step 8
      ! ==============================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

      ! ==============================================================================
      ! 7. Exact-exchange self-energy Σ^x (truncated Coulomb operator V^tr, radius R_c)
      !      D_µν          = Σ_n^occ C_µn C_νn                (density matrix)
      !      D_ll'         = Σ_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
      !      V^tr_ll'      = Σ_PQ Z_lP (M⁺ V^tr M⁺)_PQ Z_l'Q
      !      Σ^x_ll'       = D_ll' ∘ V^tr_ll'
      !      Σ^x_λσ(k=0)   = -Σ_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
      ! ==============================================================================
      CALL compute_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma, cutoff=bs_env%ri_rs%cutoff_radius_v_w)

      ! ==============================================================================
      ! 8. Correlation self-energy Σ^c and quasiparticle energies, in a single pass for G0W0
      !    and iterated to eigenvalue self-consistency for evGW0 (&EVGW0 section)
      !
      !    (a) W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
      !        Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
      !        Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
      !        Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
      !        3D: Σ^c_λσ(±i|τ|) += Σ_n (SC)_λn ΔΣ^c_n(±i|τ|) (SC)_σn  (head and wings)
      !    (b) Σ^c_nn(k, iτ) = Σ_µ C*_µn(k) [Σ^c(k, iτ) C(k)]_µn  at the band/DOS k-points
      !        Σ^c_nn(k, iτ) -> Σ^c_nn(k, iω) -> Σ^c_nn(k, ϵ)   (ANALYTIC_CONTINUATION)
      !        ϵ_nk^GW,(i) = ϵ_nk^DFT + Σ^c_nn[G^(i-1),W](ϵ_nk^GW,(i)) + Σ^x_nn(k) - v^xc_nn(k)
      ! ==============================================================================
      CALL compute_Sigma_c_and_QP_energies(bs_env, qs_env, fm_W_time, fm_Sigma_x_Gamma)

      ! recreate the AO and RI work matrices released at the start
      DO i_work = 2, SIZE(bs_env%fm_work_mo)
         CALL cp_fm_create(bs_env%fm_work_mo(i_work), fm_struct_ao)
      END DO
      CALL cp_cfm_create(bs_env%cfm_work_mo, fm_struct_ao)
      CALL cp_cfm_create(bs_env%cfm_work_mo_2, fm_struct_ao)
      CALL cp_fm_create(bs_env%fm_Gocc, fm_struct_ao)
      CALL cp_fm_create(bs_env%fm_Gvir, fm_struct_ao)
      CALL cp_fm_struct_release(fm_struct_ao)
      CALL cp_fm_create(bs_env%fm_RI_RI, bs_env%fm_struct_RI_RI)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_large_cell_Gamma

! **************************************************************************************************
!> \brief Memory estimate of the Γ-only RI-RS GW run per MPI process, printed before the RI-RS
!>        representation is built: peak = max(W step, Σ^c and QP step, Γ correction) + replicated
!>        arrays, with dense matrices in units of R = 8 n_RI²/n_proc and A = 8 n_AO²/n_proc bytes;
!>        aborts when it cannot fit.
!> \param qs_env QS environment, for the atoms (form-factor classes) and the input
!> \param bs_env band-structure environment
!> \param dense_GB the dense and replicated peak of the governing phase per process (GB)
! **************************************************************************************************
   SUBROUTINE gamma_only_memory_estimate(qs_env, bs_env, dense_GB)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: dense_GB

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_only_memory_estimate'
      INTEGER, PARAMETER :: n_ao_dos_warn = 10000, n_tab_rem = gamma_corr_n_tab_rem, n_w_extra = 9
      REAL(KIND=dp), PARAMETER :: bytes_dp = 8.0_dp, eps_kp_duplicate = 1.0E-10_dp, &
         mem_fraction = 0.9_dp, to_GB = 1.0E-9_dp, W_restart_warn_GB = 100.0_dp

      CHARACTER(LEN=1024)                                :: msg
      INTEGER                                            :: handle, ikp, jkp, n_class, n_k, n_mo, &
                                                            n_occ, n_proc, n_t, n_t_max, n_vir
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_class, rep_atom
      LOGICAL                                            :: dos_default, duplicate, slab, write_W
      REAL(KIND=dp) :: A, corr_GB, d_grid, d_occ, d_vir, distributed_GB, mem_avail_GB, peak_GB, &
         qp_GB, qp_work, R, rem_GB, repl_GB, sig_GB, w_GB, W_restart_GB, w_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: d_class
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      n_proc = bs_env%para_env%num_pe
      n_t = bs_env%num_time_freq_points
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = bs_env%n_ao
      n_occ = MAXVAL(bs_env%n_occ(1:bs_env%n_spin))
      n_vir = n_mo - MINVAL(bs_env%n_occ(1:bs_env%n_spin))
      ! N_r = α N_AO grid points; α = RS_AO_RATIO is the same for every kind outside full-k GW
      d_grid = rs_ao_ratio_of_kind(bs_env, 1)*REAL(bs_env%n_ao, dp)
      d_occ = REAL(n_occ, dp)
      d_vir = REAL(n_vir, dp)
      ! the Γ-point correction is always on in this driver: slab (2D) or bulk (3D) version
      slab = ALL(bs_env%periodic == [1, 1, 0])

      ! R and A in GB per process
      R = bytes_dp*REAL(bs_env%n_RI, dp)**2/REAL(n_proc, dp)*to_GB
      A = bytes_dp*REAL(bs_env%n_ao, dp)**2/REAL(n_proc, dp)*to_GB

      ! distinct k-points of the quasiparticle step and its work relative to the W step
      n_k = 0
      dos_default = .FALSE.
      IF (ASSOCIATED(bs_env%kpoints_DOS)) THEN
         DO ikp = 1, bs_env%kpoints_DOS%nkp
            duplicate = .FALSE.
            DO jkp = 1, ikp - 1
               IF (ALL(ABS(bs_env%kpoints_DOS%xkp(:, ikp) - bs_env%kpoints_DOS%xkp(:, jkp)) &
                       < eps_kp_duplicate)) duplicate = .TRUE.
            END DO
            IF (.NOT. duplicate) n_k = n_k + 1
         END DO
         dos_default = ANY(bs_env%nkp_grid_DOS_input < 0)
      END IF
      ! distributed dense matrices per phase
      w_GB = REAL(n_t + n_w_extra, dp)*R + REAL(2 + 3*bs_env%n_spin, dp)*A
      qp_GB = REAL(n_t + 4, dp)*R + REAL(8 + 2*bs_env%n_spin + 2*n_k*bs_env%n_spin, dp)*A
      IF (.NOT. slab) qp_GB = qp_GB + REAL(bs_env%n_spin + 1, dp)*A
      corr_GB = 11.0_dp*R + bytes_dp*(2.0_dp*d_grid*d_occ + REAL(bs_env%n_ao, dp)*d_occ + &
                                      3.0_dp*d_occ*d_vir)/REAL(n_proc, dp)*to_GB
      distributed_GB = MAX(w_GB, qp_GB, corr_GB)

      ! replicated arrays of the correction: X(iτ), X(iω) and the slab remainder tables
      repl_GB = 2.0_dp*3.0_dp*d_grid*REAL(n_t, dp)*bytes_dp*to_GB
      rem_GB = 0.0_dp
      n_class = 0
      IF (slab) THEN
         CALL gamma_corr_classify_atoms(qs_env, n_class, atom_class, d_class, rep_atom)
         DEALLOCATE (atom_class, d_class, rep_atom)
         rem_GB = REAL(n_tab_rem**2, dp)*REAL(n_t, dp)* &
                  (REAL(n_class*(n_class + 1)/2, dp) + 2.0_dp*REAL(n_class, dp))*bytes_dp*to_GB
      END IF
      ! Σ^c_nn(k, ±iτ_j) at the distinct k-points, replicated
      sig_GB = 2.0_dp*REAL(n_t, dp)*REAL(bs_env%n_ao, dp)*REAL(n_k*bs_env%n_spin, dp)*bytes_dp*to_GB
      repl_GB = repl_GB + sig_GB
      peak_GB = distributed_GB + repl_GB + rem_GB
      dense_GB = peak_GB

      ! Dense memory held during the panel contractions (the larger of the χ and Σ^c loops, plus
      ! the replicated arrays), reserved before resolve_grid_panels sets the panel budget.
      bs_env%ri_rs%panel_reserve_GB = MAX(REAL(n_t, dp)*R + REAL(2 + 3*bs_env%n_spin, dp)*A, qp_GB) + &
                                      repl_GB + rem_GB

      qp_work = REAL(n_k, dp)*REAL(2*n_t + 2, dp)*8.0_dp*REAL(bs_env%n_ao, dp)**3
      w_work = 8.0_dp*REAL(n_t, dp)*REAL(bs_env%n_RI, dp)**3

      ! size of the W(iτ) restart files, if &PRINT%RESTART is active
      logger => cp_get_default_logger()
      CALL get_qs_env(qs_env, input=input)
      write_W = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART'), &
                      cp_p_file)
      W_restart_GB = REAL(n_t, dp)*bytes_dp*REAL(bs_env%n_RI, dp)**2*to_GB

      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, mem_avail_GB)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') 'Γ-only RI-RS GW memory estimate per MPI process:'
         WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
            'Available memory per process (system)', mem_avail_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
            'W step, dense n_RI x n_RI and n_AO x n_AO matrices', w_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
            'Self-energy and quasiparticle step (dense)', qp_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T70,F9.2,A)') &
            'Γ-point correction: maps and response (dense)', corr_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T71,F9.2,A)') &
            'Replicated wing grid vectors X(iτ), X(iω)', repl_GB - sig_GB, ' GB'
         IF (slab) THEN
            WRITE (bs_env%unit_nr, '(T4,A,T71,I10)') 'Form-factor classes of the slab correction', &
               n_class
            WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
               'Replicated zone-remainder tables', rem_GB, ' GB'
         END IF
         WRITE (bs_env%unit_nr, '(T4,A,T72,F9.2,A)') &
            'Replicated MO diagonals Σ^c_nn(k, ±iτ)', sig_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T71,I10)') &
            'Quasiparticle step: distinct band and DOS k-points', n_k
         ! candidate peaks, each with all replicated arrays; the larger one governs
         WRITE (bs_env%unit_nr, '(T4,A,T69,F9.2,A)') &
            'Estimated peak in the W step', w_GB + repl_GB + rem_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,I0,A,T70,F9.2,A)') &
            'Estimated peak in the Σ^c and QP step (', n_k, ' k-points)', &
            qp_GB + repl_GB + rem_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T70,F9.2,A)') &
            'Estimated peak in the Γ-point correction', corr_GB + repl_GB + rem_GB, ' GB'
         IF (w_GB >= MAX(qp_GB, corr_GB)) THEN
            WRITE (bs_env%unit_nr, '(T4,A,T75,A)') 'Phase that sets the peak', 'W step'
         ELSE IF (qp_GB >= corr_GB) THEN
            WRITE (bs_env%unit_nr, '(T4,A,T71,A)') 'Phase that sets the peak', 'Σ^c and QP'
         ELSE
            WRITE (bs_env%unit_nr, '(T4,A,T74,A)') 'Phase that sets the peak', 'Γ corr.'
         END IF
         WRITE (bs_env%unit_nr, '(T4,A,T71,F8.1,A)') &
            'Quasiparticle step work relative to the W step', 100.0_dp*qp_work/MAX(w_work, 1.0_dp), ' %'
         IF (write_W) WRITE (bs_env%unit_nr, '(T4,A,T70,F9.2,A)') &
            'W(iτ) restart files, all written by one process', W_restart_GB, ' GB'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      IF (write_W .AND. W_restart_GB > W_restart_warn_GB) THEN
         WRITE (msg, '(A,F0.1,A)') 'The W(iτ) restart files of this run take ', W_restart_GB, &
            ' GB and are written through one MPI process. Switch them off with '// &
            '&PROPERTIES&BANDSTRUCTURE&GW&PRINT&RESTART OFF if no restart is needed.'
         CPWARN(TRIM(msg))
      END IF

      IF (dos_default .AND. bs_env%n_ao > n_ao_dos_warn .AND. n_k > 1) THEN
         WRITE (msg, '(A,I0,A)') 'DOS%KPOINTS is not set, so the default DOS mesh adds its k-points '// &
            'to the quasiparticle step (', n_k, ' distinct k-points in total), each a full pass '// &
            'over all Σ^c(iτ). For Γ-point results set &DOS KPOINTS 1 1 1 in &BANDSTRUCTURE.'
         CPWARN(TRIM(msg))
      END IF

      ! Abort on the certain part only: the n_RI² and n_AO² terms scale as 1/n_proc, the
      ! replicated arrays do not.
      IF (mem_avail_GB > 0.0_dp .AND. peak_GB > mem_fraction*mem_avail_GB) THEN
         ! W step: (n_t + n_w_extra) R + replicated ≤ mem_fraction × available
         n_t_max = INT((mem_fraction*mem_avail_GB - repl_GB - rem_GB)/R) - n_w_extra
         IF (mem_fraction*mem_avail_GB > repl_GB + rem_GB) THEN
            WRITE (msg, '(A,F0.1,A,F0.1,A,I0,A,I0,A)') &
               'The Γ-only RI-RS GW run needs about ', peak_GB, ' GB per MPI process (dense matrices '// &
               'and replicated arrays), but only ', mem_avail_GB, ' GB are available. The dense '// &
               'matrices are distributed: about ', &
               CEILING(REAL(n_proc, dp)*distributed_GB/(mem_fraction*mem_avail_GB - repl_GB - rem_GB)), &
               ' MPI processes with the same memory per process would fit (now ', n_proc, &
               '). Alternatively use fewer time points or fewer MPI processes per node.'
         ELSE
            WRITE (msg, '(A,F0.1,A,F0.1,A)') &
               'The Γ-only RI-RS GW run needs about ', peak_GB, ' GB per MPI process, but only ', &
               mem_avail_GB, ' GB are available, and the arrays replicated on every process alone '// &
               'do not fit: use fewer MPI processes per node or fewer time points.'
         END IF
         IF (n_t_max < n_t .AND. w_GB >= MAX(qp_GB, corr_GB)) THEN
            WRITE (msg, '(2A,I0,A)') TRIM(msg), ' The W step alone allows at most ', &
               MAX(n_t_max, 0), ' time points on this allocation.'
         END IF
         CPABORT(TRIM(msg))
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_only_memory_estimate

! **************************************************************************************************
!> \brief Estimated peak memory of the Γ-only RI-RS GW run per MPI process, printed after the
!>        Z_lP fit: peak ≈ max(dense + in use after the fit + χ/W/Σ panels, in use after the fit
!>        + Z_lP fit bound). Aborts if it exceeds the available memory by more than 20 %.
!> \param bs_env band-structure environment
!> \param dense_GB dense and replicated peak of the governing phase (GB)
!> \param panels_GB panel term (GB)
! **************************************************************************************************
   SUBROUTINE gamma_only_peak_estimate(bs_env, dense_GB, panels_GB)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: dense_GB, panels_GB

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_only_peak_estimate'
      REAL(KIND=dp), PARAMETER                           :: mem_abort_factor = 1.2_dp, &
                                                            to_GB = 1.0E-9_dp

      CHARACTER(LEN=1024)                                :: msg
      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: mem_avail_GB, mem_used_GB, peak_GB

      CALL timeset(routineN, handle)

      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, mem_avail_GB)
      CALL mp_mem_used_per_rank_GB(bs_env%para_env, mem_used_GB)
      ! the larger of the panel phases and the Z_lP fit (sphere-local ϕ, Gram matrix and RHS on
      ! every rank; 0 if Z_lP was read from file)
      peak_GB = MAX(dense_GB + mem_used_GB + panels_GB, &
                    mem_used_GB + bs_env%ri_rs%fit_peak_bytes*to_GB)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T69,F9.2,A)') &
            'Estimated peak of the run per MPI process', peak_GB, ' GB'
         WRITE (bs_env%unit_nr, '(T4,A,T45,3F9.2)') '(dense, in use after Z_lP, panels)', &
            dense_GB, mem_used_GB, panels_GB
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      IF (mem_avail_GB > 0.0_dp .AND. peak_GB > mem_abort_factor*mem_avail_GB) THEN
         WRITE (msg, '(A,F0.1,A,F0.1,A,F0.1,A)') 'The Γ-only RI-RS GW run needs an estimated ', peak_GB, &
            ' GB per MPI process at its peak (dense matrices, memory in use after the Z_lP fit and '// &
            'the χ/W/Σ panels), more than ', mem_abort_factor, ' times the ', mem_avail_GB, &
            ' GB available. Use more MPI processes (nodes), fewer time points or fewer MPI processes per node.'
         CPABORT(TRIM(msg))
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_only_peak_estimate

! **************************************************************************************************
!> \brief Correlation self-energy and quasiparticle energies at the band/DOS k-points: one pass
!>        for G0W0, iterated to eigenvalue self-consistency in G for evGW0 (W^MIC, Σ^x frozen):
!>          ϵ_nk^GW,(i) = ϵ_nk^DFT + Re Σ^c_nn[G^(i-1),W](ϵ_nk^GW,(i)) + Σ^x_nn(k) − v^xc_nn(k)
!> \param bs_env band-structure environment (MOs, energies, time grid, Γ correction)
!> \param qs_env QS environment, for the k-point phases and the dipoles of the Γ correction
!> \param fm_W_time W^MIC_PQ(iτ_j) = M⁺ W(iτ_j) M⁺, one dense RI x RI matrix per time point;
!>        released here
!> \param fm_Sigma_x_Gamma Σ^x_λσ at Γ per spin; released after its MO diagonals are taken
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_and_QP_energies(bs_env, qs_env, fm_W_time, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time, fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_and_QP_energies'

      INTEGER                                            :: handle, i_iter, n_iter
      LOGICAL                                            :: converged
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenval_scf_Gamma_dft
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_x_n, V_xc_n
      REAL(KIND=dp), DIMENSION(2)                        :: e_fermi_dft
      REAL(KIND=dp), DIMENSION(3, 2)                     :: band_prev
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time
      TYPE(sigma_mo_diag_type)                           :: mo_diag

      CALL timeset(routineN, handle)

      converged = .FALSE.

      CALL init_evGW0(bs_env, n_iter, band_prev, eigenval_scf_Gamma_dft, e_fermi_dft)

      ! H(k) C(k) = S(k) C(k) ϵ(k), Σ^x_nn(k), v^xc_nn(k) at the distinct band/DOS k-points, once:
      ! evGW0 changes only the energies in G, not H^KS, C(k), Σ^x or v^xc
      CALL sigma_mo_diag_setup(bs_env, qs_env, fm_Sigma_x_Gamma, mo_diag, Sigma_x_n, V_xc_n)

      ! evGW0 self-consistency loop; for G0W0 the loop is terminated after one iteration
      DO i_iter = 1, n_iter

         ! 3D Γ-point correction through the orbital-pair kernel Ŵ_nm, per iteration (it depends on
         ! the energies in G): ΔΣ^c_n(∓i|τ|) = ∓ Σ_m^occ/vir e^(-|τ||ϵ_m-ϵ_F|) Ŵ_nm(iτ).
         ! The 2D correction is already part of W^MIC.
         IF (.NOT. bs_env%ri_rs%gamma_corr%slab) THEN
            CALL gamma_corr_sigma_shifts(bs_env, qs_env, mo_diag%dSigma)
         END IF

         ! W_ll'(iτ)   =  Σ_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         ! Σ^c_ll'(iτ) = -G^occ_ll'(i|τ|) ∘ W_ll'(iτ),  τ < 0
         ! Σ^c_ll'(iτ) =  G^vir_ll'(i|τ|) ∘ W_ll'(iτ),  τ > 0
         ! Σ^c_λσ(iτ)  =  Σ_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l') (+ (SC) diag(ΔΣ^c_n) (SC)ᵀ in 3D)
         ! Σ^c_nn(k, iτ) = Σ_µ C*_µn(k) [Σ^c(k, iτ) C(k)]_µn, per time point; the AO matrices
         ! are not kept
         CALL compute_Sigma_c(bs_env, fm_W_time, fm_Sigma_c_Gamma_time, cutoff=bs_env%ri_rs%cutoff_radius_v_w, &
                              keep_w_and_files=.TRUE., &
                              release_w_per_time=(bs_env%gw_flavour == G0W0), mo_diag=mo_diag)

         ! Σ^c_nn(k, iτ) -> Σ^c_nn(k, iω) -> Σ^c_nn(k, ϵ) at every band and DOS k-point
         ! ϵ_nk^GW,(i) = ϵ_nk^DFT + Σ^c_nn[G^(i-1),W](ϵ_nk^GW,(i)) + Σ^x_nn - v^xc_nn
         CALL qp_energies_from_mo_diag(bs_env, mo_diag, Sigma_x_n, V_xc_n, fm_Sigma_x_Gamma, i_iter)

         IF (bs_env%gw_flavour == G0W0) EXIT

         CALL print_evGW0_band_edges(bs_env, band_prev, i_iter, n_iter, converged)

         IF (i_iter == n_iter .OR. converged) EXIT

         ! eigenvalues ϵ_n to be updated in G:
         ! G^occ_µλ(i|τ|) = Σ_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         ! G^vir_µλ(i|τ|) = Σ_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
         CALL update_eigenvalues_G(bs_env)

      END DO

      CALL cp_fm_release(fm_W_time)
      CALL sigma_mo_diag_release(mo_diag)
      DEALLOCATE (Sigma_x_n, V_xc_n)

      CALL reset_and_clean_bs_env(bs_env, eigenval_scf_Gamma_dft, e_fermi_dft, fm_Sigma_x_Gamma)

      CALL delete_unnecessary_files(bs_env)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_and_QP_energies

! **************************************************************************************************
!> \brief Everything the streamed MO-diagonal Σ^c needs before the Σ^c loop, at the distinct
!>        band/DOS k-points: H(k) C(k) = S(k) C(k) ϵ(k) with H(k), S(k) from the Γ matrices by the
!>        MIC phase, Σ^x_nn(k) = Re[C^H(k) Σ^x(k) C(k)]_nn, v^xc_nn(k) likewise, and S C in 3D.
!> \param bs_env band-structure environment (H^KS, S, v^xc and MOs at Γ, k-point list)
!> \param qs_env QS environment, for the MIC phases of cfm_ikp_from_fm_Gamma
!> \param fm_Sigma_x_Gamma Σ^x_µν at Γ per spin; released here after its MO diagonals are taken
!> \param mo_diag set up here; filled by compute_Sigma_c
!> \param Sigma_x_n Σ^x_nn(k), (n_ao, n_kp_distinct, n_spin)
!> \param V_xc_n v^xc_nn(k), (n_ao, n_kp_distinct, n_spin)
! **************************************************************************************************
   SUBROUTINE sigma_mo_diag_setup(bs_env, qs_env, fm_Sigma_x_Gamma, mo_diag, Sigma_x_n, V_xc_n)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      TYPE(sigma_mo_diag_type), INTENT(INOUT)            :: mo_diag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: Sigma_x_n, V_xc_n

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sigma_mo_diag_setup'
      REAL(KIND=dp), PARAMETER                           :: eps_k = 1.0E-10_dp

      INTEGER                                            :: handle, i_kd, ikp, ispin, n_ao, n_kd, nkp
      TYPE(cp_cfm_type)                                  :: cfm_ks_ikp, cfm_s_ikp

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      nkp = bs_env%nkp_bs_and_DOS

      ! Distinct k-points: ikp copies the first earlier point with the same coordinates.
      ALLOCATE (mo_diag%i_kd_of_ikp(nkp), mo_diag%ikp_distinct(nkp))
      n_kd = 0
      DO ikp = 1, nkp
         mo_diag%i_kd_of_ikp(ikp) = 0
         DO i_kd = 1, n_kd
            IF (MAXVAL(ABS(bs_env%kpoints_DOS%xkp(:, ikp) - &
                           bs_env%kpoints_DOS%xkp(:, mo_diag%ikp_distinct(i_kd)))) < eps_k) THEN
               mo_diag%i_kd_of_ikp(ikp) = i_kd
               EXIT
            END IF
         END DO
         IF (mo_diag%i_kd_of_ikp(ikp) == 0) THEN
            n_kd = n_kd + 1
            mo_diag%ikp_distinct(n_kd) = ikp
            mo_diag%i_kd_of_ikp(ikp) = n_kd
         END IF
      END DO
      mo_diag%n_kp_distinct = n_kd
      mo_diag%qs_env => qs_env

      CALL cp_cfm_create(mo_diag%cfm_ikp, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(mo_diag%cfm_work, bs_env%fm_s_Gamma%matrix_struct)
      ALLOCATE (mo_diag%cfm_mos(n_kd, bs_env%n_spin))
      ALLOCATE (Sigma_x_n(n_ao, n_kd, bs_env%n_spin), V_xc_n(n_ao, n_kd, bs_env%n_spin))
      ALLOCATE (mo_diag%Sigma_c_n_time(n_ao, bs_env%num_time_freq_points, 2, n_kd, bs_env%n_spin))
      mo_diag%Sigma_c_n_time(:, :, :, :, :) = 0.0_dp

      DO ispin = 1, bs_env%n_spin
         DO i_kd = 1, n_kd
            ikp = mo_diag%ikp_distinct(i_kd)

            ! H^KS_µν(k), S_µν(k) from the Γ matrices
            CALL cfm_ikp_from_fm_Gamma(cfm_ks_ikp, bs_env%fm_ks_Gamma(ispin), &
                                       ikp, qs_env, bs_env%kpoints_DOS, "ORB")
            CALL cfm_ikp_from_fm_Gamma(cfm_s_ikp, bs_env%fm_s_Gamma, &
                                       ikp, qs_env, bs_env%kpoints_DOS, "ORB")

            ! Roothaan-Hall: H(k) C(k) = S(k) C(k) ϵ(k)
            CALL cp_cfm_create(mo_diag%cfm_mos(i_kd, ispin), bs_env%fm_s_Gamma%matrix_struct)
            CALL cp_cfm_geeig(cfm_ks_ikp, cfm_s_ikp, mo_diag%cfm_mos(i_kd, ispin), &
                              bs_env%eigenval_scf(:, ikp, ispin), mo_diag%cfm_work)

            ! v^xc_nn(k), Σ^x_nn(k)
            CALL sigma_mo_diag_at_k(bs_env, mo_diag, bs_env%fm_V_xc_Gamma(ispin), i_kd, ispin, &
                                    V_xc_n(:, i_kd, ispin))
            CALL sigma_mo_diag_at_k(bs_env, mo_diag, fm_Sigma_x_Gamma(ispin), i_kd, ispin, &
                                    Sigma_x_n(:, i_kd, ispin))
         END DO
         DO ikp = 1, nkp
            bs_env%eigenval_scf(:, ikp, ispin) = &
               bs_env%eigenval_scf(:, mo_diag%ikp_distinct(mo_diag%i_kd_of_ikp(ikp)), ispin)
         END DO
      END DO
      CALL cp_cfm_release(cfm_ks_ikp)
      CALL cp_cfm_release(cfm_s_ikp)
      ! the work matrices are created again inside every diagonal step (sigma_mo_diag_add)
      CALL cp_cfm_release(mo_diag%cfm_ikp)
      CALL cp_cfm_release(mo_diag%cfm_work)
      ! v^xc and Σ^x enter only through their MO diagonals taken above: release the dense AO
      ! matrices (bs_env%fm_V_xc_Gamma is used by no later step of the Γ-only driver)
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_release(bs_env%fm_V_xc_Gamma(ispin))
      END DO
      CALL cp_fm_release(fm_Sigma_x_Gamma)

      ! 3D: S C for Σ^c_µν += (S C) diag(ΔΣ^c_n) (S C)ᵀ per time point; ΔΣ^c_n itself is
      ! computed per evGW0 iteration (gamma_corr_sigma_shifts)
      IF (.NOT. bs_env%ri_rs%gamma_corr%slab) THEN
         ALLOCATE (mo_diag%fm_SC(bs_env%n_spin))
         DO ispin = 1, bs_env%n_spin
            CALL cp_fm_create(mo_diag%fm_SC(ispin), bs_env%fm_mo_coeff_Gamma(1)%matrix_struct)
            CALL parallel_gemm('N', 'N', n_ao, n_ao, n_ao, 1.0_dp, bs_env%fm_s_Gamma, &
                               bs_env%fm_mo_coeff_Gamma(ispin), 0.0_dp, mo_diag%fm_SC(ispin))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE sigma_mo_diag_setup

! **************************************************************************************************
!> \brief Quasiparticle energies from the stored MO diagonals at every band/DOS k-point:
!>          Σ^c_nn(k, iτ) -> Σ^c_nn(k, iω) -> Σ^c_nn(k, ϵ)   (analytic continuation)
!>          ϵ_nk^GW = ϵ_nk^DFT + Re Σ^c_nn(k, ϵ_nk^GW) + Σ^x_nn(k) − v^xc_nn(k)
!> \param bs_env band-structure environment (time-frequency grid, KS energies, k-point list)
!> \param mo_diag filled by compute_Sigma_c
!> \param Sigma_x_n Σ^x_nn(k)
!> \param V_xc_n v^xc_nn(k)
!> \param fm_Sigma_x_Gamma Σ^x_µν at Γ; released here for G0W0 if still allocated
!> \param i_iter evGW0 iteration
! **************************************************************************************************
   SUBROUTINE qp_energies_from_mo_diag(bs_env, mo_diag, Sigma_x_n, V_xc_n, fm_Sigma_x_Gamma, i_iter)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(sigma_mo_diag_type), INTENT(INOUT)            :: mo_diag
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Sigma_x_n, V_xc_n
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      INTEGER, INTENT(IN)                                :: i_iter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qp_energies_from_mo_diag'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: Sigma_c_ikp_n_freq
      INTEGER                                            :: handle, i_kd, ikp, ispin
      TYPE(gw_qp_summary_type)                           :: qp_summary

      CALL timeset(routineN, handle)

      ALLOCATE (Sigma_c_ikp_n_freq(bs_env%n_ao, bs_env%num_freq_points_fit))

      qp_summary%n_kpoints = bs_env%nkp_bs_and_DOS
      DO ispin = 1, bs_env%n_spin
         DO ikp = 1, bs_env%nkp_bs_and_DOS
            i_kd = mo_diag%i_kd_of_ikp(ikp)

            ! Σ^c_nn(k_i,iτ) -> Σ^c_nn(k_i,iω); occupied states continue from the lower half plane
            CALL transform_time_to_frequency(bs_env%time_frequency_grid_sigma, &
                                             mo_diag%Sigma_c_n_time(:, :, :, i_kd, ispin), &
                                             Sigma_c_ikp_n_freq)
            Sigma_c_ikp_n_freq(1:bs_env%n_occ(ispin), :) = &
               CONJG(Sigma_c_ikp_n_freq(1:bs_env%n_occ(ispin), :))

            ! Σ^c_nn(k_i,iω) -> Σ^c_nn(k_i,ϵ) and
            ! ϵ_nk_i^GW = ϵ_nk_i^DFT + Σ^c_nn(k_i,ϵ) + Σ^x_nn(k_i) - v^xc_nn(k_i)
            CALL analyt_conti_and_print(bs_env, Sigma_c_ikp_n_freq, Sigma_x_n(:, i_kd, ispin), &
                                        V_xc_n(:, i_kd, ispin), bs_env%eigenval_scf(:, ikp, ispin), &
                                        ikp, ispin, qp_summary, i_iter)
         END DO
      END DO

      CALL get_all_VBM_CBM_bandgaps(bs_env)

      ! Σ^x_µν was already released in sigma_mo_diag_setup; no-op for an unallocated array
      IF (bs_env%gw_flavour == G0W0) CALL cp_fm_release(fm_Sigma_x_Gamma)

      DEALLOCATE (Sigma_c_ikp_n_freq)

      CALL timestop(handle)

   END SUBROUTINE qp_energies_from_mo_diag

END MODULE gw_ri_rs_large_cell_Gamma
