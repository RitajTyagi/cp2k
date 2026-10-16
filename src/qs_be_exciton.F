!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Becke exciton (BE) model for open-shell singlet excited states on top of a
!>        state-specific (MOM/IMOM) UKS triplet.
!>
!>        For a single excitation i -> f from a closed-shell reference the four determinants
!>        {aa, ab, ba, bb} give  E_S = E^ab + K_if  and  E_T = E^ab - K_if = E^aa, so the
!>        zeroth-order singlet is E_S = E_T + 2 K_if. The virial-theorem argument of Becke adds
!>        the nonlocal singlet correlation E_nlCorr^S = -K_if, which gives the working equation
!>
!>                 E_S = E_T + K_if,     K_if = (if|if) = int int phi_i(1) phi_f(1) phi_i(2) phi_f(2) / r12
!>
!>        phi_i and phi_f are the two singly occupied (canonical alpha) orbitals of the M_S=+1
!>        triplet of the current SCF.
!> \par History
!>      10.2026 created
! **************************************************************************************************
MODULE qs_be_exciton
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Becke2018,&
                                              Roy2026,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals,&
                                              init_coulomb_local
   USE hartree_local_types,             ONLY: hartree_local_create,&
                                              hartree_local_release,&
                                              hartree_local_type
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE hfx_types,                       ONLY: hfx_create,&
                                              hfx_release,&
                                              hfx_type
   USE input_constants,                 ONLY: be_kif_analytic,&
                                              be_kif_grid
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                              local_rho_set_release,&
                                              local_rho_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace,&
                                              rho0_s_grid_create
   USE qs_rho0_methods,                 ONLY: init_rho0
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals,&
                                              calculate_rho_atom_coeff
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_be_exciton'

   PUBLIC :: be_exciton_main

CONTAINS

! **************************************************************************************************
!> \brief Driver of &PROPERTIES&BE_EXCITON: E_S = E_T + K_if for every requested orbital pair
!> \param qs_env ...
!> \param be_section the BE_EXCITON input section
! **************************************************************************************************
   SUBROUTINE be_exciton_main(qs_env, be_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: be_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'be_exciton_main'

      INTEGER                                            :: handle, i_mo, f_mo, idir, ipair, iw, &
                                                            method, nao, nmo_a, npairs, nocc_a, nocc_b
      INTEGER, DIMENSION(:), POINTER                     :: pair
      LOGICAL                                            :: do_e_gs
      REAL(KIND=dp)                                      :: d2, e_gs, e_s, e_t, kif, osc
      REAL(KIND=dp), DIMENSION(3)                        :: dip
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eig_a, occ_a, occ_b
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_vec
      TYPE(cp_fm_type)                                   :: c_f, c_i, work_vec
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_a
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipmat, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_energy_type), POINTER                      :: energy

      CALL timeset(routineN, handle)

      CALL cite_reference(Roy2026)
      CALL cite_reference(Becke2018)

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, be_section, "PRINT", extension=".beLog")

      CALL get_qs_env(qs_env, dft_control=dft_control, mos=mos, energy=energy, matrix_s=matrix_s)

      ! the BE model sits on top of the M_S=+1 UKS triplet
      IF (dft_control%nspins /= 2) &
         CPABORT("BE_EXCITON needs a spin-unrestricted (UKS) triplet reference")
      IF (dft_control%nimages /= 1) &
         CPABORT("BE_EXCITON is implemented for Gamma-point calculations only")

      CALL section_vals_val_get(be_section, "METHOD", i_val=method)
      CALL section_vals_val_get(be_section, "ORBITAL_PAIR", n_rep_val=npairs)
      IF (npairs == 0) CPABORT("BE_EXCITON: give at least one ORBITAL_PAIR")
      CALL section_vals_val_get(be_section, "E_GROUND_STATE", explicit=do_e_gs)
      e_gs = 0.0_dp
      IF (do_e_gs) CALL section_vals_val_get(be_section, "E_GROUND_STATE", r_val=e_gs)

      CALL get_mo_set(mos(1), mo_coeff=mo_coeff_a, eigenvalues=eig_a, occupation_numbers=occ_a, &
                      nao=nao, nmo=nmo_a)
      CALL get_mo_set(mos(2), occupation_numbers=occ_b)
      ! after MOM the occupations are non-aufbau: count, do not use homo
      nocc_a = COUNT(occ_a > 0.5_dp)
      nocc_b = COUNT(occ_b > 0.5_dp)

      ! E_T: total energy of the converged state-specific triplet SCF
      e_t = energy%total

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') REPEAT("-", 79)
         WRITE (iw, '(T2,A)') "BE exciton model (Roy & Ghosal, JCTC 2026):  E_S = E_T + K_if"
         WRITE (iw, '(T2,A)') REPEAT("-", 79)
         WRITE (iw, '(T2,A,T61,2I10)') "Occupied alpha / beta orbitals", nocc_a, nocc_b
         IF (method == be_kif_grid) THEN
            WRITE (iw, '(T2,A,T71,A10)') "K_if evaluated with", "GRID"
         ELSE
            WRITE (iw, '(T2,A,T71,A10)') "K_if evaluated with", "ANALYTIC"
         END IF
      END IF
      IF (nocc_a - nocc_b /= 2) &
         CPWARN("BE_EXCITON: SCF state is not an M_S=+1 triplet (N_a - N_b /= 2)")

      ! column vectors c_i and c_f in the distribution of the alpha MO coefficients
      NULLIFY (fm_struct_vec)
      CALL cp_fm_struct_create(fm_struct_vec, nrow_global=nao, ncol_global=1, &
                               template_fmstruct=mo_coeff_a%matrix_struct)
      CALL cp_fm_create(c_i, fm_struct_vec, name="BE c_i")
      CALL cp_fm_create(c_f, fm_struct_vec, name="BE c_f")
      CALL cp_fm_create(work_vec, fm_struct_vec, name="BE work")
      CALL cp_fm_struct_release(fm_struct_vec)

      ! AO dipole matrices X_g,mu,nu = <mu| g |nu>, g = x, y, z (about the origin: <f|i> = 0, so the
      ! transition dipole does not depend on the reference point)
      NULLIFY (dipmat)
      CALL dbcsr_allocate_matrix_set(dipmat, 3)
      DO idir = 1, 3
         ALLOCATE (dipmat(idir)%matrix)
         CALL dbcsr_copy(dipmat(idir)%matrix, matrix_s(1)%matrix, "BE dipole")
         CALL dbcsr_set(dipmat(idir)%matrix, 0.0_dp)
      END DO
      CALL build_local_moment_matrix(qs_env, dipmat, 1, ref_point=[0.0_dp, 0.0_dp, 0.0_dp])

      DO ipair = 1, npairs
         CALL section_vals_val_get(be_section, "ORBITAL_PAIR", i_rep_val=ipair, i_vals=pair)
         i_mo = pair(1)
         f_mo = pair(2)
         IF (i_mo < 1 .OR. i_mo > nmo_a .OR. f_mo < 1 .OR. f_mo > nmo_a) &
            CPABORT("BE_EXCITON: ORBITAL_PAIR outside the alpha MO range (ADDED_MOS?)")
         IF (i_mo == f_mo) CPABORT("BE_EXCITON: ORBITAL_PAIR needs two different orbitals")
         ! in the M_S=+1 triplet both singly occupied orbitals are alpha-occupied
         IF (occ_a(i_mo) < 0.5_dp .OR. occ_a(f_mo) < 0.5_dp) &
            CPABORT("BE_EXCITON: ORBITAL_PAIR orbitals must both be alpha-occupied")

         CALL cp_fm_to_fm(mo_coeff_a, c_i, 1, i_mo, 1)
         CALL cp_fm_to_fm(mo_coeff_a, c_f, 1, f_mo, 1)

         SELECT CASE (method)
         CASE (be_kif_grid)
            CALL be_kif_grid_eval(qs_env, c_i, c_f, kif, iw)
         CASE (be_kif_analytic)
            CALL be_kif_analytic_eval(qs_env, be_section, c_i, c_f, kif)
         CASE DEFAULT
            CPABORT("BE_EXCITON: unknown METHOD")
         END SELECT

         ! working equation of the BE model: E_S = E_T + 2 K_if + E_nlCorr^S = E_T + K_if
         e_s = e_t + kif

         ! frozen-orbital transition dipole d_g = <phi_f| g |phi_i> = c_f^T X_g c_i and the
         ! oscillator strength of the singlet  f = 2/3 omega_S sum_g |d_g|^2,  omega_S = E_S - E_GS
         DO idir = 1, 3
            CALL cp_dbcsr_sm_fm_multiply(dipmat(idir)%matrix, c_i, work_vec, ncol=1)
            CALL cp_fm_trace(c_f, work_vec, dip(idir))
         END DO
         d2 = SUM(dip**2)
         osc = 0.0_dp
         IF (do_e_gs) osc = 2.0_dp/3.0_dp*(e_s - e_gs)*d2

         IF (iw > 0) THEN
            WRITE (iw, '(/,T2,A,I6,A,I6)') "Orbital pair (alpha)  i =", i_mo, "   f =", f_mo
            WRITE (iw, '(T4,A,T51,F14.6,A,F12.6,A)') "eps_i", eig_a(i_mo), " Ha", eig_a(i_mo)*evolt, " eV"
            WRITE (iw, '(T4,A,T51,F14.6,A,F12.6,A)') "eps_f", eig_a(f_mo), " Ha", eig_a(f_mo)*evolt, " eV"
            WRITE (iw, '(T4,A,T51,F14.8,A,F12.6,A)') "BE| K_if = (if|if)", kif, " Ha", kif*evolt, " eV"
            WRITE (iw, '(T4,A,T44,F21.10,A)') "BE| E_T  (triplet SCF)", e_t, " Ha"
            WRITE (iw, '(T4,A,T44,F21.10,A)') "BE| E_S = E_T + K_if", e_s, " Ha"
            IF (do_e_gs) THEN
               WRITE (iw, '(T4,A,T61,F14.6,A)') "BE| Triplet excitation energy", (e_t - e_gs)*evolt, " eV"
               WRITE (iw, '(T4,A,T61,F14.6,A)') "BE| Singlet excitation energy", (e_s - e_gs)*evolt, " eV"
            END IF
            WRITE (iw, '(T4,A,T44,3F12.6)') "BE| <f|r|i> x, y, z [bohr]", dip
            IF (do_e_gs) WRITE (iw, '(T4,A,T61,ES14.6)') "BE| Oscillator strength (singlet)", osc
         END IF
      END DO

      CALL cp_fm_release(c_i)
      CALL cp_fm_release(c_f)
      CALL cp_fm_release(work_vec)
      CALL dbcsr_deallocate_matrix_set(dipmat)

      IF (iw > 0) WRITE (iw, '(T2,A,/)') REPEAT("-", 79)
      CALL cp_print_key_finished_output(iw, logger, be_section, "PRINT")

      CALL timestop(handle)

   END SUBROUTINE be_exciton_main

! **************************************************************************************************
!> \brief K_if on the real-space grid, GPW and GAPW.
!>
!>        D_s = 1/2 (c_i c_f^T + c_f c_i^T)   (symmetric, same density rho_if = phi_i phi_f)
!>        rho_if(r) = sum_{mu nu} D_s,mu,nu chi_mu(r) chi_nu(r)
!>                  [GAPW: = rho~_if + sum_A (rho1_A - rho1~_A), with the rho0 compensation]
!>        V_if = Poisson(rho_if),   J[D_s]_mu,nu = <mu| V_if |nu>   (+ one-centre terms in GAPW)
!>        K_if = c_i^T J[D_s] c_f = sum D_s,mu,nu (mu nu|la si) D_s,la,si = (if|if)
!>
!>        int rho_if = <i|f> = 0, so the G=0 term plays no role under PBC.
!> \param qs_env ...
!> \param c_i alpha MO i (nao x 1)
!> \param c_f alpha MO f (nao x 1)
!> \param kif K_if in Hartree
!> \param iw output unit (debug line for GPW)
! **************************************************************************************************
   SUBROUTINE be_kif_grid_eval(qs_env, c_i, c_f, kif, iw)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: c_i, c_f
      REAL(KIND=dp), INTENT(OUT)                         :: kif
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'be_kif_grid_eval'

      INTEGER                                            :: handle, ispin, natom, nspins
      LOGICAL                                            :: gapw
      REAL(KIND=dp)                                      :: e1c, ehartree
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type)                                   :: work_vec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_p, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mpga
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hartree_local_type), POINTER                  :: hartree_local
      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(pw_c1d_gs_type)                               :: rho_g, v_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: rho_r, v_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, matrix_h, matrix_p, matrix_s, mpga, dft_control, hartree_local, &
               local_rho_set, para_env, sab_orb, oce, pw_env, poisson_env, auxbas_pw_pool, &
               qs_kind_set, ks_env)

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control, matrix_s=matrix_s, &
                      para_env=para_env, pw_env=pw_env, natom=natom, &
                      atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
      gapw = dft_control%qs_control%gapw
      nspins = dft_control%nspins

      ! spin-resolved AO matrices: P = [D_s, 0], H = [J[D_s], .], sparsity of S (sab_orb)
      CALL dbcsr_allocate_matrix_set(matrix_p, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_h, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_p(ispin)%matrix, matrix_h(ispin)%matrix)
         CALL dbcsr_create(matrix_p(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(matrix_p(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(matrix_p(ispin)%matrix, 0.0_dp)
         CALL dbcsr_create(matrix_h(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(matrix_h(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(matrix_h(ispin)%matrix, 0.0_dp)
      END DO
      ! D_s = 1/2 (c_i c_f^T + c_f c_i^T)
      CALL cp_dbcsr_plus_fm_fm_t(matrix_p(1)%matrix, matrix_v=c_i, matrix_g=c_f, ncol=1, &
                                 alpha=1.0_dp, keep_sparsity=.TRUE., symmetry_mode=1)

      ! GAPW: one-centre expansion of D_s (rho1_A, rho1~_A) and its rho0 compensation density
      IF (gapw) THEN
         CALL get_qs_env(qs_env, oce=oce, sab_orb=sab_orb)
         CALL local_rho_set_create(local_rho_set)
         CALL allocate_rho_atom_internals(local_rho_set%rho_atom_set, atomic_kind_set, &
                                          qs_kind_set, dft_control, para_env)
         CALL init_rho0(local_rho_set, qs_env, dft_control%qs_control%gapw_control, zcore=0.0_dp)
         CALL rho0_s_grid_create(pw_env, local_rho_set%rho0_mpole)
         CALL hartree_local_create(hartree_local)
         CALL init_coulomb_local(hartree_local, natom)
         mpga(1:nspins, 1:1) => matrix_p(1:nspins)
         CALL calculate_rho_atom_coeff(qs_env, mpga, local_rho_set%rho_atom_set, &
                                       qs_kind_set, oce, sab_orb, para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set, do_rho0=.TRUE.)
      END IF

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)
      CALL auxbas_pw_pool%create_pw(rho_r)
      CALL auxbas_pw_pool%create_pw(rho_g)
      CALL auxbas_pw_pool%create_pw(v_r)
      CALL auxbas_pw_pool%create_pw(v_g)

      ! rho_if on the grid (soft part in GAPW)
      CALL calculate_rho_elec(ks_env=ks_env, matrix_p=matrix_p(1)%matrix, rho=rho_r, &
                              rho_gspace=rho_g, soft_valid=gapw)
      IF (gapw) THEN
         CALL pw_axpy(local_rho_set%rho0_mpole%rho0_s_gs, rho_g)
         IF (ASSOCIATED(local_rho_set%rho0_mpole%rhoz_cneo_s_gs)) &
            CALL pw_axpy(local_rho_set%rho0_mpole%rhoz_cneo_s_gs, rho_g)
      END IF

      ! V_if = 4 pi rho_if(G)/G^2;  ehartree = 1/2 int rho_if V_if  (= K_if/2 in GPW)
      CALL pw_poisson_solve(poisson_env, rho_g, ehartree, v_g)
      CALL pw_transfer(v_g, v_r)
      CALL pw_scale(v_r, v_r%pw_grid%dvol)

      ! J[D_s]_mu,nu = <mu|V_if|nu> on the grid
      CALL integrate_v_rspace(v_rspace=v_r, pmat=matrix_p(1), hmat=matrix_h(1), qs_env=qs_env, &
                              calculate_forces=.FALSE., gapw=gapw)
      IF (gapw) THEN
         ! one-centre Hartree potentials of rho1_A - rho1~_A, the rho0 (g0) term on the grid,
         ! and their AO matrix elements
         CALL Vh_1c_gg_integrals(qs_env, e1c, hartree_local%ecoul_1c, local_rho_set, para_env, &
                                 tddft=.TRUE., core_2nd=.TRUE.)
         CALL integrate_vhg0_rspace(qs_env, v_r, para_env, calculate_forces=.FALSE., &
                                    local_rho_set=local_rho_set)
         CALL update_ks_atom(qs_env, matrix_h, matrix_p, forces=.FALSE., tddft=.TRUE., &
                             rho_atom_external=local_rho_set%rho_atom_set)
      END IF

      ! K_if = c_i^T J[D_s] c_f
      CALL cp_fm_create(work_vec, c_f%matrix_struct, name="BE J c_f")
      CALL cp_dbcsr_sm_fm_multiply(matrix_h(1)%matrix, c_f, work_vec, ncol=1)
      CALL cp_fm_trace(c_i, work_vec, kif)
      CALL cp_fm_release(work_vec)

      IF (iw > 0 .AND. .NOT. gapw) &
         WRITE (iw, '(T4,A,T51,F14.8,A)') "BE| check: 2 E_H[rho_if] (GPW)", 2.0_dp*ehartree, " Ha"

      CALL auxbas_pw_pool%give_back_pw(rho_r)
      CALL auxbas_pw_pool%give_back_pw(rho_g)
      CALL auxbas_pw_pool%give_back_pw(v_r)
      CALL auxbas_pw_pool%give_back_pw(v_g)
      IF (gapw) THEN
         CALL local_rho_set_release(local_rho_set)
         CALL hartree_local_release(hartree_local)
      END IF
      CALL dbcsr_deallocate_matrix_set(matrix_p)
      CALL dbcsr_deallocate_matrix_set(matrix_h)

      CALL timestop(handle)

   END SUBROUTINE be_kif_grid_eval

! **************************************************************************************************
!> \brief K_if from four-centre integrals of the HFX module (reference path).
!>
!>        With the alpha "density" P_a = c_i c_i^T, P_b = 0, the UKS exchange matrix is
!>        H_a = -x K[P_a],  K[P]_mu,nu = sum_{la si} (mu la|nu si) P_la,si,  x = HF FRACTION, so
!>        c_f^T K[c_i c_i^T] c_f = (fi|fi) = (if|if) = K_if   =>   K_if = -c_f^T H_a c_f / x.
!>        The operator (COULOMB, TRUNCATED for PBC, ...) is set in BE_EXCITON&HF.
!> \param qs_env ...
!> \param be_section ...
!> \param c_i alpha MO i (nao x 1)
!> \param c_f alpha MO f (nao x 1)
!> \param kif K_if in Hartree
! **************************************************************************************************
   SUBROUTINE be_kif_analytic_eval(qs_env, be_section, c_i, c_f, kif)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: be_section
      TYPE(cp_fm_type), INTENT(IN)                       :: c_i, c_f
      REAL(KIND=dp), INTENT(OUT)                         :: kif

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'be_kif_analytic_eval'

      INTEGER                                            :: handle, ispin, n_rep_hf, nspins
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: fraction
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: work_vec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: hfx_section

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, cell, matrix_h, matrix_p, matrix_s, dft_control, x_data, &
               para_env, particle_set, qs_kind_set, hfx_section)

      hfx_section => section_vals_get_subs_vals(be_section, "HF")
      CALL section_vals_get(hfx_section, explicit=explicit, n_repetition=n_rep_hf)
      IF (.NOT. explicit) &
         CPABORT("BE_EXCITON METHOD ANALYTIC needs a BE_EXCITON&HF section")
      IF (n_rep_hf /= 1) CPABORT("BE_EXCITON&HF must be given exactly once")
      CALL section_vals_val_get(hfx_section, "FRACTION", r_val=fraction)
      IF (ABS(fraction) < EPSILON(1.0_dp)) CPABORT("BE_EXCITON&HF FRACTION must not be zero")

      CALL get_qs_env(qs_env, dft_control=dft_control, matrix_s=matrix_s, para_env=para_env, &
                      atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell)
      nspins = dft_control%nspins

      CALL hfx_create(x_data, para_env, hfx_section, atomic_kind_set, qs_kind_set, &
                      particle_set, dft_control, cell, orb_basis="ORB")

      CALL dbcsr_allocate_matrix_set(matrix_p, nspins)
      CALL dbcsr_allocate_matrix_set(matrix_h, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_p(ispin)%matrix, matrix_h(ispin)%matrix)
         CALL dbcsr_create(matrix_p(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(matrix_p(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(matrix_p(ispin)%matrix, 0.0_dp)
         CALL dbcsr_create(matrix_h(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(matrix_h(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(matrix_h(ispin)%matrix, 0.0_dp)
      END DO
      ! P_a = c_i c_i^T
      CALL cp_dbcsr_plus_fm_fm_t(matrix_p(1)%matrix, matrix_v=c_i, ncol=1, alpha=1.0_dp, &
                                 keep_sparsity=.TRUE.)

      ! H_a = -x K[P_a]
      CALL tddft_hfx_matrix(matrix_h, matrix_p, qs_env, update_energy=.FALSE., &
                            recalc_integrals=.TRUE., external_hfx_sections=hfx_section, &
                            external_x_data=x_data)

      ! K_if = -c_f^T H_a c_f / x
      CALL cp_fm_create(work_vec, c_f%matrix_struct, name="BE K c_f")
      CALL cp_dbcsr_sm_fm_multiply(matrix_h(1)%matrix, c_f, work_vec, ncol=1)
      CALL cp_fm_trace(c_f, work_vec, kif)
      kif = -kif/fraction
      CALL cp_fm_release(work_vec)

      CALL dbcsr_deallocate_matrix_set(matrix_p)
      CALL dbcsr_deallocate_matrix_set(matrix_h)
      CALL hfx_release(x_data)

      CALL timestop(handle)

   END SUBROUTINE be_kif_analytic_eval

END MODULE qs_be_exciton
