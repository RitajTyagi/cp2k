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

   TYPE gamma_corr_type
      LOGICAL                                          :: enabled = .FALSE.
      LOGICAL                                          :: slab = .TRUE.
      LOGICAL                                          :: ready = .FALSE.
      REAL(KIND=dp)                                    :: metric_cut = 1.0E-4_dp
      REAL(KIND=dp)                                    :: wall_time = 0.0_dp
      REAL(KIND=dp)                                    :: measure = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                   :: b_recip = 0.0_dp
      REAL(KIND=dp)                                    :: ewald_omega = 0.0_dp
      LOGICAL                                          :: V_zero_cached = .FALSE.
      TYPE(cp_fm_type)                                 :: fm_V_zero = cp_fm_type()
      REAL(KIND=dp), ALLOCATABLE                       :: charge(:)
      REAL(KIND=dp)                                    :: charge_norm = 0.0_dp
      TYPE(cp_fm_type)                                 :: fm_F_density = cp_fm_type()
      TYPE(cp_fm_type)                                 :: fm_G_extern = cp_fm_type()
      REAL(KIND=dp), ALLOCATABLE                       :: T_RI(:, :, :)
      REAL(KIND=dp), ALLOCATABLE                       :: S_head(:, :, :)
      INTEGER, DIMENSION(2)                            :: n_tab = 0
      REAL(KIND=dp), ALLOCATABLE                       :: K_head(:, :), K_wing(:, :, :)
      REAL(KIND=dp), ALLOCATABLE                       :: K_head1(:, :), K_head2(:, :), &
                                                          K_wing1(:, :, :)
      INTEGER                                          :: q_form = 0
      REAL(KIND=dp), ALLOCATABLE                       :: q_par(:, :)
      INTEGER                                          :: n_class = 0, n_tab_rem = 0
      INTEGER, ALLOCATABLE                             :: atom_class(:)
      REAL(KIND=dp), ALLOCATABLE                       :: K_rem_head(:, :, :), K_rem_wing(:, :, :, :)
      REAL(KIND=dp), ALLOCATABLE                       :: dir_hat(:, :), dir_radius(:)
      REAL(KIND=dp), ALLOCATABLE                       :: head_pot_freq(:)
      REAL(KIND=dp)                                    :: q_max = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE                       :: A_wing(:, :, :), B_hd(:, :, :)
      REAL(KIND=dp)                                    :: xi_madelung = 0.0_dp
      LOGICAL                                          :: tail_3d = .TRUE.
      REAL(KIND=dp), ALLOCATABLE                       :: tail_factor(:)
      TYPE(cp_fm_type)                                 :: fm_K_tail = cp_fm_type()
      INTEGER                                          :: n_low = 0
      TYPE(cp_fm_type)                                 :: fm_U_low = cp_fm_type()
      REAL(KIND=dp), ALLOCATABLE                       :: m_low(:)
      REAL(KIND=dp), ALLOCATABLE                       :: eigenval_dipole(:, :)
      REAL(KIND=dp), ALLOCATABLE                       :: head_gamma(:)
   END TYPE gamma_corr_type

   ! User input for the in-memory RI-RS grid optimization.
   TYPE ri_rs_grid_opt_type
      LOGICAL                                          :: enabled = .FALSE.
      REAL(KIND=dp)                                    :: cutoff_atomic_cluster = -1.0_dp
      REAL(KIND=dp)                                    :: rs_ao_ratio = -1.0_dp
      INTEGER                                          :: max_iter = 20
      LOGICAL                                          :: reproducible_ties = .TRUE.
      INTEGER, ALLOCATABLE                             :: grid_rep(:), grid_centers(:)
   END TYPE ri_rs_grid_opt_type

   TYPE ri_rs_env

      ! Input parameters for RI-RS
      INTEGER                                          :: grid_select = 1
      INTEGER                                          :: d_lP_computation = 0
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
      INTEGER, ALLOCATABLE                             :: grid_atom_order(:)
      INTEGER, ALLOCATABLE                             :: pan_first(:), pan_last(:)
      REAL(KIND=dp)                                    :: panel_plan_cutoff = 0.0_dp
      REAL(KIND=dp)                                    :: panel_reserve_GB = 0.0_dp

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
      IF (ALLOCATED(ri_rs%grid_opt%grid_rep)) DEALLOCATE (ri_rs%grid_opt%grid_rep)
      IF (ALLOCATED(ri_rs%grid_opt%grid_centers)) DEALLOCATE (ri_rs%grid_opt%grid_centers)
      IF (ALLOCATED(ri_rs%radius_ao_per_atom)) DEALLOCATE (ri_rs%radius_ao_per_atom)
      IF (ALLOCATED(ri_rs%radius_ri_per_atom)) DEALLOCATE (ri_rs%radius_ri_per_atom)
      IF (ALLOCATED(ri_rs%chunk_centroids)) DEALLOCATE (ri_rs%chunk_centroids)
      IF (ALLOCATED(ri_rs%grid_atom_boundaries)) DEALLOCATE (ri_rs%grid_atom_boundaries)
      IF (ALLOCATED(ri_rs%grid_atom_order)) DEALLOCATE (ri_rs%grid_atom_order)
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

      IF (ri_rs%gamma_corr%V_zero_cached) THEN
         CALL cp_fm_release(ri_rs%gamma_corr%fm_V_zero)
         ri_rs%gamma_corr%V_zero_cached = .FALSE.
      END IF
      IF (ALLOCATED(ri_rs%gamma_corr%charge)) DEALLOCATE (ri_rs%gamma_corr%charge)
      CALL cp_fm_release(ri_rs%gamma_corr%fm_F_density)
      CALL cp_fm_release(ri_rs%gamma_corr%fm_G_extern)
      IF (ALLOCATED(ri_rs%gamma_corr%T_RI)) DEALLOCATE (ri_rs%gamma_corr%T_RI)
      IF (ALLOCATED(ri_rs%gamma_corr%S_head)) DEALLOCATE (ri_rs%gamma_corr%S_head)
      IF (ALLOCATED(ri_rs%gamma_corr%K_head)) DEALLOCATE (ri_rs%gamma_corr%K_head)
      IF (ALLOCATED(ri_rs%gamma_corr%K_wing)) DEALLOCATE (ri_rs%gamma_corr%K_wing)
      IF (ALLOCATED(ri_rs%gamma_corr%K_head1)) DEALLOCATE (ri_rs%gamma_corr%K_head1)
      IF (ALLOCATED(ri_rs%gamma_corr%K_head2)) DEALLOCATE (ri_rs%gamma_corr%K_head2)
      IF (ALLOCATED(ri_rs%gamma_corr%K_wing1)) DEALLOCATE (ri_rs%gamma_corr%K_wing1)
      IF (ALLOCATED(ri_rs%gamma_corr%q_par)) DEALLOCATE (ri_rs%gamma_corr%q_par)
      IF (ALLOCATED(ri_rs%gamma_corr%atom_class)) DEALLOCATE (ri_rs%gamma_corr%atom_class)
      IF (ALLOCATED(ri_rs%gamma_corr%K_rem_head)) DEALLOCATE (ri_rs%gamma_corr%K_rem_head)
      IF (ALLOCATED(ri_rs%gamma_corr%K_rem_wing)) DEALLOCATE (ri_rs%gamma_corr%K_rem_wing)
      IF (ALLOCATED(ri_rs%gamma_corr%dir_hat)) DEALLOCATE (ri_rs%gamma_corr%dir_hat)
      IF (ALLOCATED(ri_rs%gamma_corr%dir_radius)) DEALLOCATE (ri_rs%gamma_corr%dir_radius)
      IF (ALLOCATED(ri_rs%gamma_corr%head_pot_freq)) DEALLOCATE (ri_rs%gamma_corr%head_pot_freq)
      IF (ALLOCATED(ri_rs%gamma_corr%A_wing)) DEALLOCATE (ri_rs%gamma_corr%A_wing)
      IF (ALLOCATED(ri_rs%gamma_corr%B_hd)) DEALLOCATE (ri_rs%gamma_corr%B_hd)
      IF (ALLOCATED(ri_rs%gamma_corr%tail_factor)) DEALLOCATE (ri_rs%gamma_corr%tail_factor)
      CALL cp_fm_release(ri_rs%gamma_corr%fm_K_tail)
      CALL cp_fm_release(ri_rs%gamma_corr%fm_U_low)
      IF (ALLOCATED(ri_rs%gamma_corr%m_low)) DEALLOCATE (ri_rs%gamma_corr%m_low)
      ri_rs%gamma_corr%n_low = 0
      IF (ALLOCATED(ri_rs%gamma_corr%eigenval_dipole)) DEALLOCATE (ri_rs%gamma_corr%eigenval_dipole)
      IF (ALLOCATED(ri_rs%gamma_corr%head_gamma)) DEALLOCATE (ri_rs%gamma_corr%head_gamma)
      ri_rs%gamma_corr%ready = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE ri_rs_release

END MODULE ri_rs_types
