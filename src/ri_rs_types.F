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
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
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

   ! State of the Gamma-point head-and-wing correction of W for a Gamma-only 2D slab.
   !
   ! At q = 0 the real-space Coulomb lattice sum over two CHARGED auxiliary functions diverges
   ! with the summation range, so the charge channel is split off and carried analytically:
   !
   !   c_P    = integral dr phi_P0(r)                              (auxiliary charge)
   !   Pi_PQ  = delta_PQ - c_P c_Q / sum_S |c_S|^2                 (neutral projector)
   !   V_b    = Pi V(0) Pi                                         (neutral body, converges as R^-3)
   !   phi_0  = 1/sqrt(A), v(q) = 2 pi / q                         (the separated charge mode)
   !
   ! and the head/wing of eps in that extended basis is integrated over the whole small
   ! Brillouin zone of the supercell instead of being sampled. See GAMMA_W_CORRECTION_PLAN.md.
   TYPE gamma_corr_type
      LOGICAL                                          :: enabled = .FALSE.
      LOGICAL                                          :: ready = .FALSE.
      LOGICAL                                          :: reported_projection = .FALSE.
      ! Layer area A = |a_1 x a_2| and the two in-plane reciprocal vectors b_1, b_2.
      REAL(KIND=dp)                                    :: area = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 2)                   :: b_recip = 0.0_dp
      ! Zero-wavevector Fourier value of the RI metric kernel: pi/omega^2 (erfc) or
      ! 2 pi r_c^2 (truncated Coulomb). The charge-channel right-hand side is scaled by it so
      ! that the fixed Tikhonov parameter regularizes it on the same footing as the Gaussian
      ! columns, and the fitted coefficients are divided by it again.
      REAL(KIND=dp)                                    :: metric_head = 0.0_dp
      ! Normalized auxiliary charge vector c_P / ||c||, n_RI entries, replicated, and ||c||.
      REAL(KIND=dp), ALLOCATABLE                       :: charge(:)
      REAL(KIND=dp)                                    :: charge_norm = 0.0_dp
      ! F_density = M^-1(0) Pi V_b^(1/2)   maps a fitted grid density to the response basis
      ! F_extern  = M^-1(0) V(0) Pi V_b^(-1/2)  maps it to the external potential
      ! Both n_RI x n_RI, replicated on every rank.
      REAL(KIND=dp), ALLOCATABLE                       :: F_density(:, :), F_extern(:, :)
      ! Z_l0: grid coefficients of the normalized constant function, all grid rows.
      REAL(KIND=dp), ALLOCATABLE                       :: Z_l0(:)
      ! S_head(a,b,i_w): the Cartesian dielectric head with the body Schur complement already
      ! removed, S = Y - E^T B^-1 E. The whole directional dependence of the correction enters
      ! through the scalar a(q_hat,i_w) = q_hat^T S q_hat.
      REAL(KIND=dp), ALLOCATABLE                       :: S_head(:, :, :)
      ! T_wing(l,a,i_w): the grid wing map, T^w_la = sum_Q Z_lQ [F_extern B^-1 E]_Qa. The wing
      ! of W is linear in the direction, t_l(q_hat) = sum_a q_hat_a T^w_la.
      REAL(KIND=dp), ALLOCATABLE                       :: T_wing(:, :, :)
      ! Cell-integral tables over the in-plane minimum-image displacement, see
      ! gamma_corr_tabulate: K_head(i_R, i_t) and K_wing(i_R, a, i_t) on a regular R mesh.
      INTEGER, DIMENSION(2)                            :: n_tab = 0
      REAL(KIND=dp), ALLOCATABLE                       :: K_head(:, :), K_wing(:, :, :)
      ! Hot-loop layout of the same information, built once at the end of gamma_corr_tabulate.
      ! Sigma^c touches every grid PAIR, so the kernel evaluation has to cost a handful of
      ! contiguous loads rather than a frequency transform and a pbc() call per pair:
      !   K_time(1:4, i_R, i_t)   head and the three wing components already transformed to
      !                           imaginary time, interleaved so one bilinear corner is one
      !                           cache line
      !   U_grid(1:4, l, i_t)     (Z_l0, T_wing(l,1:3)) packed per grid point
      !   frac(1:2, l)            in-plane fractional coordinates of every grid point, so the
      !                           minimum image is a subtraction and an ANINT
      REAL(KIND=dp), ALLOCATABLE                       :: K_time(:, :, :), U_grid(:, :, :), &
                                                          frac(:, :)
   END TYPE gamma_corr_type

   ! User input for the in-memory RI-RS grid optimization.
   TYPE ri_rs_grid_opt_type
      LOGICAL                                          :: enabled = .FALSE.
      REAL(KIND=dp)                                    :: cutoff_atomic_cluster = 3.0_dp/angstrom
      REAL(KIND=dp)                                    :: rs_ao_ratio = 5.0_dp
      INTEGER                                          :: max_iter = 300
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

      ! Gamma-point head-and-wing correction of W (Gamma-only 2D slabs); see gamma_corr_type.
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

      IF (ALLOCATED(ri_rs%gamma_corr%charge)) DEALLOCATE (ri_rs%gamma_corr%charge)
      IF (ALLOCATED(ri_rs%gamma_corr%F_density)) DEALLOCATE (ri_rs%gamma_corr%F_density)
      IF (ALLOCATED(ri_rs%gamma_corr%F_extern)) DEALLOCATE (ri_rs%gamma_corr%F_extern)
      IF (ALLOCATED(ri_rs%gamma_corr%Z_l0)) DEALLOCATE (ri_rs%gamma_corr%Z_l0)
      IF (ALLOCATED(ri_rs%gamma_corr%S_head)) DEALLOCATE (ri_rs%gamma_corr%S_head)
      IF (ALLOCATED(ri_rs%gamma_corr%T_wing)) DEALLOCATE (ri_rs%gamma_corr%T_wing)
      IF (ALLOCATED(ri_rs%gamma_corr%K_head)) DEALLOCATE (ri_rs%gamma_corr%K_head)
      IF (ALLOCATED(ri_rs%gamma_corr%K_wing)) DEALLOCATE (ri_rs%gamma_corr%K_wing)
      IF (ALLOCATED(ri_rs%gamma_corr%K_time)) DEALLOCATE (ri_rs%gamma_corr%K_time)
      IF (ALLOCATED(ri_rs%gamma_corr%U_grid)) DEALLOCATE (ri_rs%gamma_corr%U_grid)
      IF (ALLOCATED(ri_rs%gamma_corr%frac)) DEALLOCATE (ri_rs%gamma_corr%frac)
      ri_rs%gamma_corr%ready = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE ri_rs_release

END MODULE ri_rs_types
