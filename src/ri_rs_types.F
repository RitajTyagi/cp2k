!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared real-space resolution-of-identity grid and matrix data.
! **************************************************************************************************
MODULE ri_rs_types
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_release,&
                                              dbcsr_type
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: ri_rs_env, ri_rs_grid_opt_type, rirs_grid_type, ri_rs_release, gamma_corr_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_types'

   ! Atom-relative RI-RS grids
   TYPE rirs_grid_type
      INTEGER                    :: npts = 0
      REAL(KIND=dp), ALLOCATABLE :: raw_points(:, :)
   END TYPE rirs_grid_type

   ! State of the Gamma-point head-and-wing correction of W for a Gamma-only 2D slab or 3D bulk
   ! cell.
   !
   ! At q = 0 the Coulomb operator between two CHARGED auxiliary functions is singular, so the
   ! charge channel is split off and carried analytically:
   !
   !   c_P    = integral dr phi_P0(r)                              (auxiliary charge)
   !   2D:  Pi_PQ = delta_PQ - c_P c_Q / sum_S |c_S|^2             (neutral projector)
   !        V_b   = Pi V(0) Pi                                     (neutral body, converges as R^-3)
   !        phi_0 = 1/sqrt(A), v(q) = 2 pi / q                     (the separated charge mode)
   !   3D:  V(0)  = Ewald sum with the G = 0 Fourier coefficient of 1/r omitted (finite in the
   !        charge channel, no projector), phi_0 = 1/sqrt(Omega), v(q) = 4 pi / q^2
   !
   ! and the head/wing of eps in that extended basis is integrated over the whole small
   ! Brillouin zone of the supercell instead of being sampled. See GAMMA_W_CORRECTION_PLAN.md
   ! and large_cell_gamma/GAMMA_CORRECTION_3D_PLAN.md.
   TYPE gamma_corr_type
      ! enabled: set by gamma_corr_init (PERIODIC XY or XYZ with KPOINTS_W 1 1 1); every call
      !          site of the correction tests it.
      ! slab:    PERIODIC XY (2D head 2 pi/q, tabulated in-plane kernels) when .TRUE.,
      !          PERIODIC XYZ (3D head 4 pi/q^2, kernels per atom pair) when .FALSE.
      ! ready:   set once step 4b of the Gamma-only driver has built the maps, T_RI
      !          and the tables; the Sigma^c loop asserts it before adding DeltaW.
      ! reported_projection: the V(0) - D_box c c^T summary has been printed (print once).
      LOGICAL                                          :: enabled = .FALSE.
      LOGICAL                                          :: slab = .TRUE.
      LOGICAL                                          :: ready = .FALSE.
      ! distributed: the dense n_RI x n_RI algebra (maps, dielectric body) runs on cp_fm /
      !          ScaLAPACK matrices instead of replicated LAPACK arrays (gamma_corr_init decides
      !          from the memory estimate).
      LOGICAL                                          :: distributed = .FALSE.
      LOGICAL                                          :: reported_projection = .FALSE.
      ! measure: the layer area A = |a_1 x a_2| (bohr^2, 2D) or the cell volume Omega (bohr^3,
      ! 3D); every 1/A of the slab formulas is 1/Omega in the bulk ones.
      ! b_recip(1:3, i) is the reciprocal vector b_i in Cartesian components (1/bohr): the two
      ! in-plane ones for a slab (third column zero), all three for a bulk cell.
      REAL(KIND=dp)                                    :: measure = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                   :: b_recip = 0.0_dp
      ! ewald_omega: splitting parameter of the 3D Ewald sum that builds V(0) (1/bohr); chosen
      ! from the image box of the real-space erfc sum by gamma_corr_coulomb_at_zero.
      REAL(KIND=dp)                                    :: ewald_omega = 0.0_dp
      ! Zero-wavevector Fourier value of the RI metric kernel: pi/omega^2 (erfc) or
      ! 2 pi r_c^2 (truncated Coulomb). The charge-channel right-hand side is scaled by it so
      ! that the fixed Tikhonov parameter regularizes it on the same footing as the Gaussian
      ! columns, and the fitted coefficients are divided by it again.
      REAL(KIND=dp)                                    :: metric_head = 0.0_dp
      ! V_PQ(q=0), the real-space lattice sum of the bare Coulomb operator, distributed. Built
      ! once by gamma_corr_coulomb_at_zero and handed to compute_V_k_by_lattice_sum, which
      ! would otherwise repeat the same lattice sum for the regular W path; released there.
      LOGICAL                                          :: V_zero_cached = .FALSE.
      TYPE(cp_fm_type)                                 :: fm_V_zero = cp_fm_type()
      ! Normalized auxiliary charge vector c_P / ||c||, n_RI entries, replicated, and ||c||.
      REAL(KIND=dp), ALLOCATABLE                       :: charge(:)
      REAL(KIND=dp)                                    :: charge_norm = 0.0_dp
      ! F_density = M^-1(0) Pi V_b^(1/2)   maps a fitted grid density to the response basis
      ! G_extern  = V(0) Pi V_b^(-1/2)      maps the body response to the external potential, in
      !             the RI projection representation that W^MIC_PQ is in before M^-1 W M^-1
      !             (it is M times the F_extern of the equations)
      ! Both n_RI x n_RI, replicated on every rank (replicated dense-algebra path).
      REAL(KIND=dp), ALLOCATABLE                       :: F_density(:, :), G_extern(:, :)
      ! The same two maps as distributed cp_fm matrices on the distributed path.
      TYPE(cp_fm_type)                                 :: fm_F_density = cp_fm_type()
      TYPE(cp_fm_type)                                 :: fm_G_extern = cp_fm_type()
      ! T_RI(P,a,i_w): the wing vector in the RI projection representation,
      ! T_Pa = sum_Q G_extern(P,Q) [B^-1 E]_Qa. Shape (n_RI, 3, n_freq), replicated.
      REAL(KIND=dp), ALLOCATABLE                       :: T_RI(:, :, :)
      ! S_head(a,b,i_w): the Cartesian dielectric head with the body Schur complement already
      ! removed, S = Y - E^T B^-1 E. The whole directional dependence of the correction enters
      ! through the scalar a(q_hat,i_w) = q_hat^T S q_hat.
      ! Shape (3, 3, n_freq), replicated.
      REAL(KIND=dp), ALLOCATABLE                       :: S_head(:, :, :)
      ! Cell-integral tables over the in-plane minimum-image displacement, see
      ! gamma_corr_tabulate: K_head(i_R, i_w) and K_wing(i_R, a, i_w) on a regular R mesh.
      ! The last index is the imaginary FREQUENCY; gamma_corr_add_to_W interpolates them at the
      ! atom-pair separations and transforms the result to imaginary time.
      ! n_tab: table points along a_1 and a_2; i_R = 1 + i_u + n_tab(1) i_v over [-1/2, 1/2]^2.
      INTEGER, DIMENSION(2)                            :: n_tab = 0
      REAL(KIND=dp), ALLOCATABLE                       :: K_head(:, :), K_wing(:, :, :)
      ! 3D bulk cells (slab = .FALSE.): the screening denominator 1 + 4 pi a(q_hat, i_w) is
      ! independent of |q|, so the radial integrals are analytic and no table is needed. The
      ! six-face cubature of the BZ cell (gamma_corr_bz_faces) gives the directions q_hat
      ! (dir_hat(3, n_dir)) and the radius to the cell boundary (dir_radius(n_dir)); per
      ! direction and frequency, head_coef(n_dir, n_w) and wing_coef(n_dir, n_w) are the
      ! angular weight times the head/wing prefactor (gamma_corr_direction_coefficients).
      ! gamma_corr_add_to_W_bulk contracts them with sin(L t)/t and (sin(L t) - L t cos(L t))/t^2,
      ! t = q_hat . R, at the minimum-image separation of every local atom pair.
      REAL(KIND=dp), ALLOCATABLE                       :: dir_hat(:, :), dir_radius(:)
      REAL(KIND=dp), ALLOCATABLE                       :: head_coef(:, :), wing_coef(:, :)
      ! head_pot_time(i_t): the head potential w(i tau) = K_h(0, i omega)/Omega transformed to
      ! imaginary time; gamma_corr_add_to_sigma adds -+ w S G S to Sigma^c(-+ i tau) (3D only).
      REAL(KIND=dp), ALLOCATABLE                       :: head_pot_time(:)
      ! head_pot_freq(i_w): the same head potential in imaginary frequency; A_wing(a,b,i_w) and
      ! B_hd(a,b,i_w): the cell-averaged wing and head-dipole tensors (1/Omega), see
      ! gamma_corr_add_to_sigma (3D only).
      REAL(KIND=dp), ALLOCATABLE                       :: head_pot_freq(:)
      ! q_max: largest radius of the BZ cell (1/bohr); dipoles are clipped at 1/q_max.
      REAL(KIND=dp)                                    :: q_max = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE                       :: A_wing(:, :, :), B_hd(:, :, :)
   END TYPE gamma_corr_type

   ! User input for the in-memory RI-RS grid optimization.
   TYPE ri_rs_grid_opt_type
      LOGICAL                                          :: enabled = .FALSE.
      REAL(KIND=dp)                                    :: cutoff_atomic_cluster = -1.0_dp
      REAL(KIND=dp)                                    :: rs_ao_ratio = 7.0_dp
      INTEGER                                          :: max_iter = 20
   END TYPE ri_rs_grid_opt_type

   TYPE ri_rs_env

      ! Input parameters for RI-RS
      INTEGER                                          :: grid_select = 1
      REAL(KIND=dp)                                    :: tikhonov = 1.0E-08_dp
      REAL(KIND=dp)                                    :: cutoff_radius_ri_ao = -1.0_dp
      REAL(KIND=dp)                                    :: cutoff_radius_v_w = -1.0_dp
      TYPE(ri_rs_grid_opt_type)                        :: grid_opt = ri_rs_grid_opt_type()

      ! Geometry references used by RI-RS grid setup and optimization.
      ! They remain owned by qs_env and must never be deallocated here.
      TYPE(cell_type), POINTER                         :: cell => NULL()
      TYPE(particle_type), DIMENSION(:), POINTER       :: particle_set => NULL()
      TYPE(atomic_kind_type), DIMENSION(:), POINTER    :: atomic_kind_set => NULL()

      ! Data types for cutoffs based DBCSR matrices
      REAL(KIND=dp), ALLOCATABLE                       :: chunk_centroids(:, :)
      INTEGER, ALLOCATABLE                             :: grid_atom_boundaries(:)
      INTEGER, ALLOCATABLE                             :: pan_first(:), pan_last(:)
      ! Panel plans share a fixed grid, process layout and message limit until release.
      REAL(KIND=dp)                                    :: panel_plan_cutoff = 0.0_dp

      ! Atom-specific grids used directly by the grid assembler
      TYPE(rirs_grid_type), ALLOCATABLE                :: atomic_grids(:)

      ! Data types for storing RI-RS matrices
      INTEGER                                          :: n_grid_points = 0
      TYPE(dbcsr_type)                                 :: mat_phi_mu_l = dbcsr_type()
      TYPE(dbcsr_type)                                 :: mat_Z_lP = dbcsr_type()
      REAL(KIND=dp), ALLOCATABLE                       :: grid_points(:, :)
      LOGICAL                                          :: Z_lP_exists = .FALSE.

      ! Dense Gram + phi + threaded RHS memory model, collected during fitting setup.
      INTEGER                                          :: fit_max_grid = 0
      REAL(KIND=dp)                                    :: fit_peak_bytes = 0.0_dp

      ! Spatial extent of most diffuse Gaussian primitive of AO/RI basis for each atom
      REAL(KIND=dp), ALLOCATABLE                       :: radius_ao_per_atom(:)
      REAL(KIND=dp), ALLOCATABLE                       :: radius_ri_per_atom(:)

      ! Precomputed grid-basis kernels for RT-BSE
      ! mat_V_aux_rtbse(P, Q) = truncated-Coulomb V_PQ with M^-1 sandwich (RI x RI);
      ! the Hartree grid kernel Z V Z^T is applied factorized, never materialized
      ! mat_W0_grid_rtbse(l, l') = sum_PQ  Z_lP W^{w=0}_PQ Z_l'Q
      TYPE(dbcsr_type)                                 :: mat_V_aux_rtbse = dbcsr_type()
      TYPE(dbcsr_type)                                 :: mat_W0_grid_rtbse = dbcsr_type()

      ! RTBSE%EPS_FILTER_RHO: cuts Delta-rho_AO in compute_sigma_ri_rs, also cuts SEX and Hartree
      REAL(KIND=dp)                                    :: eps_filter_rho = -1.0_dp

      ! RTBSE%CUTOFF_RADIUS_W0 (bohr): drops grid block pairs of W0 beyond this distance.
      ! BSE-only: the GW quasiparticle energies use their own &GW CUTOFF_RADIUS_* family.
      REAL(KIND=dp)                                    :: cutoff_radius_w0 = -1.0_dp

      ! Set TRUE once mat_phi_mu_l and mat_Z_lP have been populated in memory (by GW RI-RS
      ! or on-demand by the RT-BSE RI-RS kernel path).
      LOGICAL                                          :: grid_built = .FALSE.

      ! Independent toggles for auxiliary Coulomb / grid W0 availability
      LOGICAL                                          :: V_grid_built = .FALSE.
      LOGICAL                                          :: W0_grid_built = .FALSE.

      ! Gamma-point head-and-wing correction of W (Gamma-only 2D slabs and 3D bulk cells); see
      ! gamma_corr_type.
      TYPE(gamma_corr_type)                            :: gamma_corr = gamma_corr_type()

   END TYPE ri_rs_env

CONTAINS

! **************************************************************************************************
!> \brief Release owned RI-RS matrices and grid data, then clear construction state.
!>        Input settings and restart-file availability are retained. Geometry references are
!>        borrowed from qs_env and are only nullified. Repeated calls are safe.
!> \param ri_rs RI-RS input, owned grid/matrix data and borrowed geometry references
! **************************************************************************************************
   SUBROUTINE ri_rs_release(ri_rs)
      TYPE(ri_rs_env), INTENT(INOUT)                     :: ri_rs

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ri_rs_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ri_rs%W0_grid_built) CALL dbcsr_release(ri_rs%mat_W0_grid_rtbse)
      IF (ri_rs%V_grid_built) CALL dbcsr_release(ri_rs%mat_V_aux_rtbse)
      IF (ri_rs%grid_built) THEN
         CALL dbcsr_release(ri_rs%mat_phi_mu_l)
         CALL dbcsr_release(ri_rs%mat_Z_lP)
      END IF
      IF (ALLOCATED(ri_rs%grid_points)) DEALLOCATE (ri_rs%grid_points)
      IF (ALLOCATED(ri_rs%atomic_grids)) DEALLOCATE (ri_rs%atomic_grids)
      IF (ALLOCATED(ri_rs%radius_ao_per_atom)) DEALLOCATE (ri_rs%radius_ao_per_atom)
      IF (ALLOCATED(ri_rs%radius_ri_per_atom)) DEALLOCATE (ri_rs%radius_ri_per_atom)
      IF (ALLOCATED(ri_rs%chunk_centroids)) DEALLOCATE (ri_rs%chunk_centroids)
      IF (ALLOCATED(ri_rs%grid_atom_boundaries)) DEALLOCATE (ri_rs%grid_atom_boundaries)
      IF (ALLOCATED(ri_rs%pan_first)) DEALLOCATE (ri_rs%pan_first)
      IF (ALLOCATED(ri_rs%pan_last)) DEALLOCATE (ri_rs%pan_last)
      ri_rs%panel_plan_cutoff = 0.0_dp
      NULLIFY (ri_rs%atomic_kind_set, ri_rs%cell, ri_rs%particle_set)
      ri_rs%fit_max_grid = 0
      ri_rs%fit_peak_bytes = 0.0_dp
      ri_rs%n_grid_points = 0
      ri_rs%grid_built = .FALSE.
      ri_rs%V_grid_built = .FALSE.
      ri_rs%W0_grid_built = .FALSE.

      ! Gamma-point correction: release the cached V(0), free every allocatable of
      ! gamma_corr_type and clear ready, so a later GW run rebuilds the correction. enabled is
      ! left alone; gamma_corr_init resets it.
      IF (ri_rs%gamma_corr%V_zero_cached) THEN
         CALL cp_fm_release(ri_rs%gamma_corr%fm_V_zero)
         ri_rs%gamma_corr%V_zero_cached = .FALSE.
      END IF
      IF (ALLOCATED(ri_rs%gamma_corr%charge)) DEALLOCATE (ri_rs%gamma_corr%charge)
      IF (ALLOCATED(ri_rs%gamma_corr%F_density)) DEALLOCATE (ri_rs%gamma_corr%F_density)
      IF (ALLOCATED(ri_rs%gamma_corr%G_extern)) DEALLOCATE (ri_rs%gamma_corr%G_extern)
      CALL cp_fm_release(ri_rs%gamma_corr%fm_F_density)
      CALL cp_fm_release(ri_rs%gamma_corr%fm_G_extern)
      IF (ALLOCATED(ri_rs%gamma_corr%T_RI)) DEALLOCATE (ri_rs%gamma_corr%T_RI)
      IF (ALLOCATED(ri_rs%gamma_corr%S_head)) DEALLOCATE (ri_rs%gamma_corr%S_head)
      IF (ALLOCATED(ri_rs%gamma_corr%K_head)) DEALLOCATE (ri_rs%gamma_corr%K_head)
      IF (ALLOCATED(ri_rs%gamma_corr%K_wing)) DEALLOCATE (ri_rs%gamma_corr%K_wing)
      IF (ALLOCATED(ri_rs%gamma_corr%dir_hat)) DEALLOCATE (ri_rs%gamma_corr%dir_hat)
      IF (ALLOCATED(ri_rs%gamma_corr%dir_radius)) DEALLOCATE (ri_rs%gamma_corr%dir_radius)
      IF (ALLOCATED(ri_rs%gamma_corr%head_coef)) DEALLOCATE (ri_rs%gamma_corr%head_coef)
      IF (ALLOCATED(ri_rs%gamma_corr%wing_coef)) DEALLOCATE (ri_rs%gamma_corr%wing_coef)
      IF (ALLOCATED(ri_rs%gamma_corr%head_pot_time)) DEALLOCATE (ri_rs%gamma_corr%head_pot_time)
      IF (ALLOCATED(ri_rs%gamma_corr%head_pot_freq)) DEALLOCATE (ri_rs%gamma_corr%head_pot_freq)
      IF (ALLOCATED(ri_rs%gamma_corr%A_wing)) DEALLOCATE (ri_rs%gamma_corr%A_wing)
      IF (ALLOCATED(ri_rs%gamma_corr%B_hd)) DEALLOCATE (ri_rs%gamma_corr%B_hd)
      ri_rs%gamma_corr%ready = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE ri_rs_release

END MODULE ri_rs_types
