!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initialize atom-owned molecular or periodic RI-RS grids by Cholesky selection.
! **************************************************************************************************
MODULE ri_rs_grid_initialization
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE input_constants,                 ONLY: do_gapw_log,&
                                              ri_rs_small_cell_full_kp
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE message_passing,                 ONLY: mp_comm_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              deallocate_grid_atom,&
                                              grid_atom_type
   USE ri_rs_types,                     ONLY: rirs_grid_type
   USE ri_rs_utils,                     ONLY: evaluate_ao_basis_on_points,&
                                              filter_grid_to_voronoi,&
                                              get_rirs_cluster_atoms
   USE string_utilities,                ONLY: uppercase
   USE util,                            ONLY: locate
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_grid_initialization'

   PUBLIC :: initialize_ri_rs_grid, find_translation_orbits, rs_ao_ratio_of_kind, &
             translation_atom_images

CONTAINS

! **************************************************************************************************
!> \brief Construct the initial RI-RS grid of every atom.
!>
!> For every atom B, the Lebedev grid points are
!>
!> (1)  r_l = R_B + ρ_s Ω_a,
!>
!> where ρ_s is a radial quadrature point and Ω_a is a Lebedev direction. Lebedev grids are
!> tabulated according to the angular degree L up to which they integrate exactly. To integrate a
!> three-centre integral (μν|P), the smallest available Lebedev grid satisfying
!>
!> (2)  L >= 2 l_AO,max + l_RI,max + ΔL
!>
!> is used. Here, l_AO,max is the maximum angular momentum of the atomic AO basis functions ϕ_μ,
!> l_RI,max is the maximum angular momentum of the auxiliary basis functions φ_P, and the internal
!> angular buffer is ΔL. Enough radial points ρ_s are used that the Lebedev grid of every atom B
!> contains at least
!>
!> (3)  N_initial^B >= α_initial N_AO^B,
!>
!> points, where N_AO^B is the number of atomic orbitals on atom B and
!> α_initial = max(30, 2 RS_AO_RATIO). The factor two provides candidates that can be discarded by
!> the atomic Voronoi filter. The requested initial RI-RS grid of atom A contains
!>
!> (4)  N_R^A = ceil(α * N_AO^A)
!>
!> points, with α = RS_AO_RATIO of the atomic kind (rs_ao_ratio_of_kind). Around every atom A, a
!> cluster is defined as
!>
!> (5)  C_A contains all atoms and, for periodic systems, their images within the cluster radius.
!>
!> The Cholesky selection uses all Lebedev grid points of every atom image B in C_A. The
!> Voronoi cell of atom A is
!>
!> (6)  V_A = {r_l : |r_l - R_A| <= |r_l - R_B| for every atom image B}.
!>
!> For cluster C_A, Cholesky selection is performed on
!>
!> (7)  D_ll' = [Σ_(μ in C_A) ϕ_μ(r_l) ϕ_μ(r_l')]^2.
!>
!> The first grid point is the point with the largest diagonal element,
!>
!> (8)  d_l^(0) = D_ll,                 q_1 = arg max_l d_l^(0).
!>
!> For every selected point q_k, the Cholesky column and diagonal are updated according to
!>
!> (9)  L_lk = [D_lq_k - Σ_(j<k) L_lj L_q_kj]/sqrt(d_q_k^(k-1)),
!>
!> (10) d_l^(k) = max(0, d_l^(k-1) - L_lk^2),
!>
!> and the next point q_(k+1) is the point with the largest d_l^(k). Points outside V_A take part
!> in Eqs. (7)-(10), but only selected points inside V_A are placed into the RI-RS grid of atom A.
!> If D_ll' reaches its numerical rank before N_R^A points have been retained, let S_A contain the
!> retained points. Every unused point r_l in V_A is assigned the distance
!>
!> (11) ρ_l = min_(q in S_A) |r_l - q|,
!>
!> and the point
!>
!> (12) q_new = arg max_(r_l in V_A and r_l not in S_A) ρ_l
!>
!> is appended to S_A. Equations (11)-(12) are repeated until Eq. (4) is satisfied. This maximin
!> distance is only a geometric selection criterion, not the three-centre-integral fitting error.
!> Clusters C_A are assigned to MPI groups; Cholesky rows are distributed within each group.
!> The completed atom grids are broadcast to all ranks.
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************
   SUBROUTINE initialize_ri_rs_grid(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_ri_rs_grid'
      INTEGER, PARAMETER :: l_additional = 2, radial_quadrature = do_gapw_log
      REAL(KIND=dp), PARAMETER :: minimum_candidate_ratio = 30.0_dp

      INTEGER                                            :: handle, ikind
      REAL(KIND=dp)                                      :: candidate_ratio
      TYPE(rirs_grid_type), ALLOCATABLE                  :: radial_lebedev_grids(:)

      CALL timeset(routineN, handle)

      ! Build the Lebedev grids of Eqs. (1)-(3) for all elements in the calculation
      ALLOCATE (radial_lebedev_grids(SIZE(bs_env%basis_set_AO)))
      DO ikind = 1, SIZE(bs_env%basis_set_AO)
         candidate_ratio = MAX(minimum_candidate_ratio, 2.0_dp*rs_ao_ratio_of_kind(bs_env, ikind))
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T69,F12.1)') 'RS_AO_RATIO of '// &
               TRIM(bs_env%ri_rs%atomic_kind_set(ikind)%element_symbol)//' ('// &
               TRIM(bs_env%basis_set_AO(ikind)%gto_basis_set%name)//')', rs_ao_ratio_of_kind(bs_env, ikind)
         END IF
         CALL build_lebedev_grid(bs_env%basis_set_AO(ikind)%gto_basis_set, &
                                 bs_env%basis_set_RI(ikind)%gto_basis_set, &
                                 candidate_ratio, l_additional, radial_quadrature, &
                                 radial_lebedev_grids(ikind))
      END DO

      ! Apply Eqs. (5)-(12) for every atom and keep points inside Voronoi volume, Eq. (6)
      CALL cholesky_selection_and_voronoi_filtering(bs_env, radial_lebedev_grids)

      CALL broadcast_ri_rs_grids(bs_env)

      CALL timestop(handle)

   END SUBROUTINE initialize_ri_rs_grid

! **************************************************************************************************
!> \brief Build one fixed-orientation Lebedev grid according to Eqs. (1)-(3).
!> \param ao ...
!> \param ri ...
!> \param ratio α_initial in Eq. (3).
!> \param l_additional ΔL in Eq. (2).
!> \param radial_quadrature ...
!> \param grid ...
! **************************************************************************************************
   SUBROUTINE build_lebedev_grid(ao, ri, ratio, l_additional, radial_quadrature, grid)
      TYPE(gto_basis_set_type), POINTER                  :: ao, ri
      REAL(KIND=dp), INTENT(IN)                          :: ratio
      INTEGER, INTENT(IN)                                :: l_additional, radial_quadrature
      TYPE(rirs_grid_type), INTENT(OUT)                  :: grid

      INTEGER                                            :: degree, ir, nang, nrad, offset, rule
      TYPE(grid_atom_type), POINTER                      :: radial_grid

      CPASSERT(ASSOCIATED(ao) .AND. ASSOCIATED(ri))
      degree = 2*MAXVAL(ao%lmax) + MAXVAL(ri%lmax) + l_additional
      rule = get_number_of_lebedev_grid(l=degree)
      nang = lebedev_grid(rule)%n
      nrad = MAX(2, CEILING(ratio*REAL(ao%nsgf, dp)/REAL(nang, dp)))

      NULLIFY (radial_grid)
      CALL allocate_grid_atom(radial_grid)
      CALL create_grid_atom(radial_grid, nrad, nang, 0, rule, radial_quadrature)
      grid%npts = nrad*nang
      ALLOCATE (grid%raw_points(3, grid%npts))
      DO ir = 1, nrad
         offset = (ir - 1)*nang
         grid%raw_points(:, offset + 1:offset + nang) = &
            radial_grid%rad(ir)*lebedev_grid(rule)%r
      END DO
      CALL deallocate_grid_atom(radial_grid)
   END SUBROUTINE build_lebedev_grid

! **************************************************************************************************
!> \brief Apply Eqs. (5)-(12) independently for every atom and keep points inside Eq. (6).
!> \param bs_env ...
!> \param radial_lebedev_grids ...
! **************************************************************************************************
   SUBROUTINE cholesky_selection_and_voronoi_filtering(bs_env, radial_lebedev_grids)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(rirs_grid_type), INTENT(IN)                   :: radial_lebedev_grids(:)

      INTEGER                                            :: iatom, icenter, ikind, my_group, &
                                                            ncenters, ngroups, npoints
      TYPE(mp_comm_type)                                 :: atom_env

      CPASSERT(.NOT. ALLOCATED(bs_env%ri_rs%atomic_grids))
      ! Translation orbits: only the representative atoms (grid_centers) get their own
      ! Cholesky selection; images receive a copy in broadcast_ri_rs_grids.
      CALL find_translation_orbits(bs_env)
      ncenters = SIZE(bs_env%ri_rs%grid_opt%grid_centers)
      ALLOCATE (bs_env%ri_rs%atomic_grids(bs_env%n_atom))
      ngroups = MIN(ncenters, bs_env%para_env%num_pe)
      my_group = MOD(bs_env%para_env%mepos, ngroups)
      CALL atom_env%from_split(bs_env%para_env, my_group)
      DO iatom = 1, bs_env%n_atom
         ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
         npoints = CEILING(rs_ao_ratio_of_kind(bs_env, ikind)* &
                           bs_env%basis_set_AO(ikind)%gto_basis_set%nsgf)
         bs_env%ri_rs%atomic_grids(iatom)%npts = npoints
         ALLOCATE (bs_env%ri_rs%atomic_grids(iatom)%raw_points(3, npoints))
      END DO
      DO icenter = 1, ncenters
         IF (MOD(icenter - 1, ngroups) == my_group) THEN
            CALL select_cholesky_grid_iatom(bs_env, radial_lebedev_grids, &
                                            bs_env%ri_rs%grid_opt%grid_centers(icenter), atom_env)
         END IF
      END DO
      CALL atom_env%free()
   END SUBROUTINE cholesky_selection_and_voronoi_filtering

! **************************************************************************************************
!> \brief Select one atom grid with Eqs. (5)-(12) and the atomic Voronoi cell in Eq. (6).
!> \param bs_env ...
!> \param radial_lebedev_grids ...
!> \param iatom ...
!> \param atom_env ...
! **************************************************************************************************
   SUBROUTINE select_cholesky_grid_iatom(bs_env, radial_lebedev_grids, iatom, atom_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(rirs_grid_type), INTENT(IN)                   :: radial_lebedev_grids(:)
      INTEGER, INTENT(IN)                                :: iatom
      TYPE(mp_comm_type), INTENT(IN)                     :: atom_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'select_cholesky_grid_iatom'

      CHARACTER(LEN=3*default_string_length)             :: failure
      INTEGER                                            :: handle, igrid_point, npoints, nretained
      INTEGER, ALLOCATABLE :: ao_point_indices(:), atom_n_ao(:), atom_point_offsets(:), &
         atom_value_offsets(:), cluster_atoms(:), cluster_cells(:, :), selected_indices(:)
      LOGICAL, ALLOCATABLE                               :: inside_voronoi(:)
      REAL(KIND=dp)                                      :: center(3)
      REAL(KIND=dp), ALLOCATABLE                         :: ao_values(:), cluster_centers(:, :), &
                                                            cluster_points(:, :)

      CALL timeset(routineN, handle)

      CALL get_rirs_cluster_atoms(bs_env%ri_rs%particle_set, bs_env%ri_rs%cell, iatom, &
                                  bs_env%ri_rs%grid_opt%cutoff_atomic_cluster, cluster_atoms, &
                                  cluster_cells, cluster_centers)
      center = cluster_centers(:, 1)
      CALL collect_cluster_lebedev_points(bs_env, radial_lebedev_grids, cluster_atoms, &
                                          cluster_centers, center, cluster_points)
      npoints = SIZE(cluster_points, 2)
      ALLOCATE (inside_voronoi(npoints))
      CALL filter_grid_to_voronoi(cluster_points, iatom, bs_env%ri_rs%particle_set, &
                                  bs_env%ri_rs%cell, mask=inside_voronoi)
      IF (COUNT(inside_voronoi) < bs_env%ri_rs%atomic_grids(iatom)%npts) THEN
         WRITE (failure, '(A,I0,A,I0,A,I0,A)') 'Atom ', iatom, ': only ', COUNT(inside_voronoi), &
            ' Voronoi grid points for ', bs_env%ri_rs%atomic_grids(iatom)%npts, &
            ' points; reduce RS_AO_RATIO.'
         CPABORT(TRIM(failure))
      END IF

      CALL evaluate_cluster_ao_values(bs_env, cluster_atoms, cluster_centers, center, cluster_points, &
                                      ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                      atom_n_ao, ao_values)
      ALLOCATE (selected_indices(bs_env%ri_rs%atomic_grids(iatom)%npts))
      CALL select_cholesky_grid_points(ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                       atom_n_ao, ao_values, npoints, cluster_points, &
                                       inside_voronoi, selected_indices, atom_env, &
                                       bs_env%ri_rs%grid_opt%reproducible_ties)
      nretained = COUNT(selected_indices > 0)
      IF (nretained == 0) THEN
         WRITE (failure, '(A,I0,A)') 'Atom ', iatom, &
            ': Cholesky selection found no numerically independent point.'
         CPABORT(TRIM(failure))
      ELSE IF (nretained < SIZE(selected_indices)) THEN
         CALL complete_grid_by_maximin_distance(cluster_points, inside_voronoi, selected_indices, &
                                                bs_env%ri_rs%grid_opt%reproducible_ties)
      END IF
      DO igrid_point = 1, SIZE(selected_indices)
         bs_env%ri_rs%atomic_grids(iatom)%raw_points(:, igrid_point) = &
            cluster_points(:, selected_indices(igrid_point))
      END DO

      CALL timestop(handle)

   END SUBROUTINE select_cholesky_grid_iatom

! **************************************************************************************************
!> \brief Complete a rank-saturated atom grid with the geometric maximin rule in Eqs. (11)-(12).
!> \param cluster_points Coordinates r_l of all existing cluster-grid points.
!> \param inside_voronoi True for points inside the atomic Voronoi cell V_A.
!> \param selected_indices Cholesky indices on entry and the completed indices on exit.
!> \param reproducible_ties Lexicographic tie rule (else the first point within 1e-8 of the maximum).
! **************************************************************************************************
   SUBROUTINE complete_grid_by_maximin_distance(cluster_points, inside_voronoi, selected_indices, &
                                                reproducible_ties)
      REAL(KIND=dp), INTENT(IN)                          :: cluster_points(:, :)
      LOGICAL, INTENT(IN)                                :: inside_voronoi(:)
      INTEGER, INTENT(INOUT)                             :: selected_indices(:)
      LOGICAL, INTENT(IN)                                :: reproducible_ties

      REAL(KIND=dp), PARAMETER :: tie_tolerance = 1.0E-8_dp, &
         tie_tolerance_reproducible = 1.0E-10_dp

      INTEGER                                            :: igrid_point, inew_point, iselected, &
                                                            nretained
      LOGICAL, ALLOCATABLE                               :: unselected_inside(:)
      REAL(KIND=dp)                                      :: distance_sq, largest_distance_sq
      REAL(KIND=dp), ALLOCATABLE                         :: nearest_distance_sq(:)

      CPASSERT(SIZE(cluster_points, 2) == SIZE(inside_voronoi))
      CPASSERT(COUNT(inside_voronoi) >= SIZE(selected_indices))

      nretained = COUNT(selected_indices > 0)
      CPASSERT(nretained > 0)
      ALLOCATE (unselected_inside(SIZE(inside_voronoi)), &
                nearest_distance_sq(SIZE(inside_voronoi)))
      unselected_inside(:) = inside_voronoi
      nearest_distance_sq(:) = HUGE(1.0_dp)

      ! Initialize ρ_l^2 in Eq. (11) from all retained Cholesky points.
      DO iselected = 1, nretained
         unselected_inside(selected_indices(iselected)) = .FALSE.
         DO igrid_point = 1, SIZE(inside_voronoi)
            IF (.NOT. unselected_inside(igrid_point)) CYCLE
            distance_sq = SUM((cluster_points(:, igrid_point) - &
                               cluster_points(:, selected_indices(iselected)))**2)
            nearest_distance_sq(igrid_point) = &
               MIN(nearest_distance_sq(igrid_point), distance_sq)
         END DO
      END DO

      DO WHILE (nretained < SIZE(selected_indices))
         largest_distance_sq = MAXVAL(nearest_distance_sq, MASK=unselected_inside)
         IF (reproducible_ties) THEN
            inew_point = pick_tied_point(nearest_distance_sq, unselected_inside, cluster_points, &
                                         tie_tolerance_reproducible*largest_distance_sq)
         ELSE
            ! Ties up to rounding go to the first point, independent of the compiler.
            inew_point = 0
            DO igrid_point = 1, SIZE(inside_voronoi)
               IF (.NOT. unselected_inside(igrid_point)) CYCLE
               IF (nearest_distance_sq(igrid_point) >= (1.0_dp - tie_tolerance)*largest_distance_sq) THEN
                  inew_point = igrid_point
                  EXIT
               END IF
            END DO
         END IF
         CPASSERT(inew_point > 0)

         nretained = nretained + 1
         selected_indices(nretained) = inew_point
         unselected_inside(inew_point) = .FALSE.
         DO igrid_point = 1, SIZE(inside_voronoi)
            IF (.NOT. unselected_inside(igrid_point)) CYCLE
            distance_sq = SUM((cluster_points(:, igrid_point) - &
                               cluster_points(:, inew_point))**2)
            nearest_distance_sq(igrid_point) = &
               MIN(nearest_distance_sq(igrid_point), distance_sq)
         END DO
      END DO
   END SUBROUTINE complete_grid_by_maximin_distance

! **************************************************************************************************
!> \brief Collect all cluster Lebedev points used in Eqs. (5) and (7).
!> \param bs_env ...
!> \param radial_lebedev_grids ...
!> \param cluster_atoms ...
!> \param cluster_centers ...
!> \param center Position of the central atom, defining atom-relative grid coordinates.
!> \param cluster_points ...
! **************************************************************************************************
   SUBROUTINE collect_cluster_lebedev_points(bs_env, radial_lebedev_grids, cluster_atoms, &
                                             cluster_centers, center, cluster_points)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(rirs_grid_type), INTENT(IN)                   :: radial_lebedev_grids(:)
      INTEGER, INTENT(IN)                                :: cluster_atoms(:)
      REAL(KIND=dp), INTENT(IN)                          :: cluster_centers(:, :), center(3)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: cluster_points(:, :)

      INTEGER                                            :: cluster_iatom, i_cluster_atom, ikind, n, &
                                                            npoints, offset

      npoints = 0
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         ikind = bs_env%ri_rs%particle_set(cluster_atoms(i_cluster_atom))%atomic_kind%kind_number
         npoints = npoints + radial_lebedev_grids(ikind)%npts
      END DO
      ALLOCATE (cluster_points(3, npoints))

      offset = 0
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         cluster_iatom = cluster_atoms(i_cluster_atom)
         ikind = bs_env%ri_rs%particle_set(cluster_iatom)%atomic_kind%kind_number
         n = radial_lebedev_grids(ikind)%npts
         cluster_points(:, offset + 1:offset + n) = radial_lebedev_grids(ikind)%raw_points + &
                                                    SPREAD(cluster_centers(:, i_cluster_atom) - &
                                                           center, 2, n)
         offset = offset + n
      END DO
   END SUBROUTINE collect_cluster_lebedev_points

! **************************************************************************************************
!> \brief Evaluate ϕ_μ(r_l) needed for D_ll' in Eq. (7), using AO locality. Atom b uses the
!>        point rows atom_point_offsets(b):atom_point_offsets(b+1)-1. Its column-major block
!>        ϕ_μ(r_l) has shape (number of point rows, atom_n_ao(b)) and starts at
!>        atom_value_offsets(b) in ao_values.
!> \param bs_env ...
!> \param cluster_atoms ...
!> \param cluster_centers ...
!> \param center Position of the central atom, defining atom-relative grid coordinates.
!> \param cluster_points ...
!> \param ao_point_indices Cluster-point index for every stored point row.
!> \param atom_point_offsets First stored point row for each atom, followed by the final bound.
!> \param atom_value_offsets First AO value for each atom, followed by the final bound.
!> \param atom_n_ao Number of AO functions for each atom.
!> \param ao_values Compact atom-local values ϕ_μ(r_l).
! **************************************************************************************************
   SUBROUTINE evaluate_cluster_ao_values(bs_env, cluster_atoms, cluster_centers, center, cluster_points, &
                                         ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                         atom_n_ao, ao_values)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cluster_atoms(:)
      REAL(KIND=dp), INTENT(IN)                          :: cluster_centers(:, :), center(3), &
                                                            cluster_points(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: ao_point_indices(:), &
                                                            atom_point_offsets(:), &
                                                            atom_value_offsets(:), atom_n_ao(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: ao_values(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_cluster_ao_values'

      INTEGER :: cluster_iatom, first_point, first_value, handle, i_cluster_atom, igrid_point, &
         ikind, last_point, last_value, nactive, npoints
      INTEGER, ALLOCATABLE                               :: point_indices(:)
      LOGICAL, ALLOCATABLE                               :: active(:)
      REAL(KIND=dp)                                      :: source_position(3)
      REAL(KIND=dp), ALLOCATABLE                         :: atom_points(:, :), atom_values(:, :), &
                                                            relative_points(:, :)
      TYPE(gto_basis_set_type), POINTER                  :: ao

      CALL timeset(routineN, handle)

      npoints = SIZE(cluster_points, 2)
      ALLOCATE (relative_points(3, npoints))
      ALLOCATE (point_indices(npoints))
      point_indices(:) = [(igrid_point, igrid_point=1, npoints)]
      ALLOCATE (active(npoints), atom_n_ao(SIZE(cluster_atoms)), &
                atom_point_offsets(SIZE(cluster_atoms) + 1), &
                atom_value_offsets(SIZE(cluster_atoms) + 1))

      atom_point_offsets(1) = 1
      atom_value_offsets(1) = 1
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         cluster_iatom = cluster_atoms(i_cluster_atom)
         ikind = bs_env%ri_rs%particle_set(cluster_iatom)%atomic_kind%kind_number
         ao => bs_env%basis_set_AO(ikind)%gto_basis_set
         CPASSERT(ao%kind_radius > 0.0_dp)
         relative_points(:, :) = cluster_points + &
                                 SPREAD(center - cluster_centers(:, i_cluster_atom), 2, npoints)

         active(:) = SUM(relative_points**2, DIM=1) &
                     <= ao%kind_radius**2
         nactive = COUNT(active)
         atom_n_ao(i_cluster_atom) = ao%nsgf
         atom_point_offsets(i_cluster_atom + 1) = atom_point_offsets(i_cluster_atom) + nactive
         atom_value_offsets(i_cluster_atom + 1) = atom_value_offsets(i_cluster_atom) + &
                                                  nactive*atom_n_ao(i_cluster_atom)
      END DO

      ALLOCATE (ao_point_indices(atom_point_offsets(SIZE(cluster_atoms) + 1) - 1), &
                ao_values(atom_value_offsets(SIZE(cluster_atoms) + 1) - 1))
      DO i_cluster_atom = 1, SIZE(cluster_atoms)
         cluster_iatom = cluster_atoms(i_cluster_atom)
         ikind = bs_env%ri_rs%particle_set(cluster_iatom)%atomic_kind%kind_number
         ao => bs_env%basis_set_AO(ikind)%gto_basis_set
         relative_points(:, :) = cluster_points + &
                                 SPREAD(center - cluster_centers(:, i_cluster_atom), 2, npoints)

         active(:) = SUM(relative_points**2, DIM=1) &
                     <= ao%kind_radius**2
         first_point = atom_point_offsets(i_cluster_atom)
         last_point = atom_point_offsets(i_cluster_atom + 1) - 1
         first_value = atom_value_offsets(i_cluster_atom)
         last_value = atom_value_offsets(i_cluster_atom + 1) - 1
         nactive = last_point - first_point + 1
         ao_point_indices(first_point:last_point) = PACK(point_indices, active)
         ALLOCATE (atom_points(3, nactive), atom_values(nactive, ao%nsgf))
         atom_points(:, :) = RESHAPE(PACK(relative_points, SPREAD(active, 1, 3)), [3, nactive])
         atom_values(:, :) = 0.0_dp
         source_position = 0.0_dp
         CALL evaluate_ao_basis_on_points(atom_values, atom_points, ao, source_position, &
                                          bs_env%ri_rs%cell, minimum_image=.FALSE.)
         ao_values(first_value:last_value) = RESHAPE(atom_values, [nactive*ao%nsgf])
         DEALLOCATE (atom_points, atom_values)
      END DO

      CALL timestop(handle)

   END SUBROUTINE evaluate_cluster_ao_values

! **************************************************************************************************
!> \brief Apply the Cholesky selection of Eqs. (7)-(10) and retain selected points in V_A.
!> \param ao_point_indices Cluster-point index for every stored point row.
!> \param atom_point_offsets First stored point row for each atom, followed by the final bound.
!> \param atom_value_offsets First AO value for each atom, followed by the final bound.
!> \param atom_n_ao Number of AO functions for each atom.
!> \param ao_values Compact atom-local values ϕ_μ(r_l).
!> \param npoints Number of points in the cluster Lebedev grids.
!> \param cluster_points Coordinates r_l of the cluster Lebedev points.
!> \param inside_voronoi True for points inside the atomic Voronoi cell V_A.
!> \param selected_indices Selected indices in V_A; zero denotes rank exhaustion.
!> \param atom_env ...
!> \param reproducible_ties Lexicographic tie rule (else the first point within 1e-8 of the maximum).
! **************************************************************************************************
   SUBROUTINE select_cholesky_grid_points(ao_point_indices, atom_point_offsets, atom_value_offsets, &
                                          atom_n_ao, ao_values, npoints, cluster_points, &
                                          inside_voronoi, selected_indices, atom_env, &
                                          reproducible_ties)
      INTEGER, INTENT(IN)                                :: ao_point_indices(:), &
                                                            atom_point_offsets(:), &
                                                            atom_value_offsets(:), atom_n_ao(:)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET      :: ao_values(:)
      INTEGER, INTENT(IN)                                :: npoints
      REAL(KIND=dp), INTENT(IN)                          :: cluster_points(3, npoints)
      LOGICAL, INTENT(IN)                                :: inside_voronoi(npoints)
      INTEGER, INTENT(OUT)                               :: selected_indices(:)
      TYPE(mp_comm_type), INTENT(IN)                     :: atom_env
      LOGICAL, INTENT(IN)                                :: reproducible_ties

      CHARACTER(LEN=*), PARAMETER :: routineN = 'select_cholesky_grid_points'
      REAL(KIND=dp), PARAMETER :: rank_eps_factor = 64.0_dp, tie_tolerance = 1.0E-8_dp, &
         tie_tolerance_reproducible = 1.0E-12_dp

      INTEGER :: capacity, first_local, first_point, first_row, first_value, handle, iblock, &
         igrid_point, irow, last_local, last_point, last_row, last_value, n_ao, new_capacity, &
         nretained, nrows, nselected, selected_point
      INTEGER, ALLOCATABLE                               :: first_row_local(:), last_row_local(:)
      REAL(KIND=dp)                                      :: initial_max, rank_tolerance, &
                                                            selected_diagonal, tie_tol
      REAL(KIND=dp), ALLOCATABLE :: diagonal(:), factor(:, :), grown(:, :), initial_diagonal(:), &
         selected_ao_values(:), selected_factor_values(:), selection_matrix_column(:)
      REAL(KIND=dp), POINTER                             :: atom_ao_values(:, :)

      CALL timeset(routineN, handle)

      ALLOCATE (diagonal(npoints), initial_diagonal(npoints), selection_matrix_column(npoints))

      ! d_l^(0) = D_ll = [Σ_μ ϕ_μ(r_l)^2]^2.
      diagonal(:) = 0.0_dp
      DO iblock = 1, SIZE(atom_n_ao)
         first_point = atom_point_offsets(iblock)
         last_point = atom_point_offsets(iblock + 1) - 1
         first_value = atom_value_offsets(iblock)
         last_value = atom_value_offsets(iblock + 1) - 1
         nrows = last_point - first_point + 1
         atom_ao_values(1:nrows, 1:atom_n_ao(iblock)) => ao_values(first_value:last_value)
         DO irow = 1, nrows
            igrid_point = ao_point_indices(first_point + irow - 1)
            diagonal(igrid_point) = diagonal(igrid_point) + SUM(atom_ao_values(irow, :)**2)
         END DO
      END DO
      diagonal(:) = diagonal**2
      IF (reproducible_ties) THEN
         ! Rank cut per point against d_l^(0) and tie tolerance against the initial maximum: the
         ! round-off of the residual d_l^(k) stays at about eps * max_l d_l^(0) as d_l^(k) shrinks.
         initial_diagonal(:) = diagonal
         initial_max = MAXVAL(initial_diagonal)
         tie_tol = tie_tolerance_reproducible*initial_max
      ELSE
         ! Residuals below N*eps*max(D_ll) are rounding noise (rank tolerance of LAPACK DPSTRF);
         ! accepting them would let the compiler's rounding decide the grid.
         rank_tolerance = REAL(npoints, dp)*EPSILON(1.0_dp)*MAXVAL(diagonal)
      END IF

      capacity = MIN(32, npoints)
      ALLOCATE (factor(npoints, capacity), selected_ao_values(MAXVAL(atom_n_ao)), &
                selected_factor_values(npoints))
      selected_indices(:) = 0
      first_local = 1 + INT((INT(npoints, int_8)*atom_env%mepos)/atom_env%num_pe)
      last_local = INT((INT(npoints, int_8)*(atom_env%mepos + 1))/atom_env%num_pe)
      ALLOCATE (first_row_local(SIZE(atom_n_ao)), last_row_local(SIZE(atom_n_ao)))
      DO iblock = 1, SIZE(atom_n_ao)
         first_point = atom_point_offsets(iblock)
         last_point = atom_point_offsets(iblock + 1) - 1
         first_row_local(iblock) = 1 + COUNT(ao_point_indices(first_point:last_point) < first_local)
         last_row_local(iblock) = COUNT(ao_point_indices(first_point:last_point) <= last_local)
      END DO
      nselected = 0
      nretained = 0
      DO WHILE (nretained < SIZE(selected_indices) .AND. nselected < npoints)
         ! Pivot q_k = argmax_l d_l^(k-1) (Eq. (8)). Lebedev candidates on a symmetric site are
         ! degenerate over a symmetry orbit up to round-off, which must not decide the grid.
         selected_diagonal = MAXVAL(diagonal)
         IF (selected_diagonal <= 0.0_dp) EXIT
         IF (reproducible_ties) THEN
            ! T = {l : d_l > 0, d_l >= max - tol}; its lexicographically smallest relative r_l,
            ! independent of absolute coordinates and of the atom order of the supercell
            selected_point = pick_tied_point(diagonal, diagonal > 0.0_dp, cluster_points, tie_tol)
         ELSE
            ! the first point within the relative tolerance of the maximum
            selected_point = FINDLOC(diagonal >= (1.0_dp - tie_tolerance)*selected_diagonal, .TRUE., DIM=1)
         END IF
         selected_diagonal = diagonal(selected_point)

         IF (nselected == capacity) THEN
            new_capacity = MIN(npoints, capacity + MAX(32, capacity/2))
            ALLOCATE (grown(npoints, new_capacity))
            grown(:, :capacity) = factor(:, :capacity)
            CALL MOVE_ALLOC(grown, factor)
            capacity = new_capacity
         END IF

         ! Form D_lq from atom-local AO products without storing the dense D_ll' matrix.
         selection_matrix_column(:) = 0.0_dp
         DO iblock = 1, SIZE(atom_n_ao)
            first_point = atom_point_offsets(iblock)
            last_point = atom_point_offsets(iblock + 1) - 1
            nrows = last_point - first_point + 1
            IF (nrows == 0) CYCLE
            irow = locate(ao_point_indices(first_point:last_point), selected_point)
            IF (irow == 0) CYCLE
            first_value = atom_value_offsets(iblock)
            last_value = atom_value_offsets(iblock + 1) - 1
            n_ao = atom_n_ao(iblock)
            atom_ao_values(1:nrows, 1:n_ao) => ao_values(first_value:last_value)
            selected_ao_values(:n_ao) = atom_ao_values(irow, :)
            first_row = first_row_local(iblock)
            last_row = last_row_local(iblock)
            ! Sum each row in the same order independently of the MPI row partition.
            DO irow = first_row, last_row
               igrid_point = ao_point_indices(first_point + irow - 1)
               selection_matrix_column(igrid_point) = &
                  selection_matrix_column(igrid_point) + &
                  DOT_PRODUCT(atom_ao_values(irow, :), selected_ao_values(:n_ao))
            END DO
         END DO
         CALL atom_env%sum(selection_matrix_column)
         selection_matrix_column(:) = selection_matrix_column**2

         ! L_lk = [D_lq_k - sum_(j<k) L_lj L_q_kj]/sqrt(d_q_k^(k-1)).
         IF (nselected > 0) THEN
            selected_factor_values(:nselected) = factor(selected_point, :nselected)
            IF (atom_env%num_pe > 1) selection_matrix_column(:first_local - 1) = 0.0_dp
            IF (atom_env%num_pe > 1) selection_matrix_column(last_local + 1:) = 0.0_dp
            DO igrid_point = first_local, last_local
               selection_matrix_column(igrid_point) = &
                  selection_matrix_column(igrid_point) - &
                  DOT_PRODUCT(factor(igrid_point, :nselected), selected_factor_values(:nselected))
            END DO
            CALL atom_env%sum(selection_matrix_column)
         END IF
         nselected = nselected + 1
         factor(:, nselected) = selection_matrix_column/SQRT(selected_diagonal)
         diagonal(:) = MAX(0.0_dp, diagonal - factor(:, nselected)**2)
         IF (reproducible_ties) THEN
            WHERE (diagonal <= rank_eps_factor*EPSILON(1.0_dp)*initial_diagonal) diagonal = 0.0_dp
         ELSE
            WHERE (diagonal <= rank_tolerance) diagonal = 0.0_dp
         END IF
         diagonal(selected_point) = 0.0_dp

         IF (inside_voronoi(selected_point)) THEN
            nretained = nretained + 1
            selected_indices(nretained) = selected_point
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE select_cholesky_grid_points

! **************************************************************************************************
!> \brief Reproducible argmax: among T = {l : eligible_l, values_l >= max_l' values_l' - tol},
!>        return the point with the lexicographically smallest coordinate relative to the cluster
!>        centre (x, then y, then z; compared as integers in units of 1e-8 bohr).
!> \param values Quantity to maximize.
!> \param eligible Candidate mask.
!> \param points Relative coordinates r_l.
!> \param tol Absolute tie tolerance.
!> \return Selected index, 0 if no point is eligible.
! **************************************************************************************************
   FUNCTION pick_tied_point(values, eligible, points, tol) RESULT(ipick)
      REAL(KIND=dp), INTENT(IN)                          :: values(:)
      LOGICAL, INTENT(IN)                                :: eligible(:)
      REAL(KIND=dp), INTENT(IN)                          :: points(:, :), tol
      INTEGER                                            :: ipick

      REAL(KIND=dp), PARAMETER                           :: quantum = 1.0E-8_dp

      INTEGER                                            :: idim, l
      INTEGER(KIND=int_8)                                :: qbest(3), ql(3)
      LOGICAL                                            :: smaller
      REAL(KIND=dp)                                      :: vmax

      ipick = 0
      IF (.NOT. ANY(eligible)) RETURN
      vmax = MAXVAL(values, MASK=eligible)
      qbest(:) = 0
      DO l = 1, SIZE(values)
         IF (.NOT. eligible(l)) CYCLE
         IF (values(l) < vmax - tol) CYCLE
         ql(:) = NINT(points(:, l)/quantum, KIND=int_8)
         IF (ipick == 0) THEN
            smaller = .TRUE.
         ELSE
            smaller = .FALSE.
            DO idim = 1, 3
               IF (ql(idim) /= qbest(idim)) THEN
                  smaller = ql(idim) < qbest(idim)
                  EXIT
               END IF
            END DO
         END IF
         IF (smaller) THEN
            ipick = l
            qbest(:) = ql
         END IF
      END DO
   END FUNCTION pick_tied_point

! **************************************************************************************************
!> \brief Communicate the independently constructed atom grids after Eqs. (5)-(12).
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE broadcast_ri_rs_grids(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: atom_mpi_rank, iatom, icenter, ncenters, &
                                                            ngroups

      ncenters = SIZE(bs_env%ri_rs%grid_opt%grid_centers)
      ngroups = MIN(ncenters, bs_env%para_env%num_pe)
      DO icenter = 1, ncenters
         iatom = bs_env%ri_rs%grid_opt%grid_centers(icenter)
         atom_mpi_rank = MOD(icenter - 1, ngroups)
         CALL bs_env%para_env%bcast(bs_env%ri_rs%atomic_grids(iatom)%raw_points, atom_mpi_rank)
      END DO
      ! Atom-relative grids are translation invariant: an image A = rep(A) + t gets the grid
      ! of its representative, r_l^A - R_A = r_l^rep(A) - R_rep(A).
      DO iatom = 1, bs_env%n_atom
         IF (bs_env%ri_rs%grid_opt%grid_rep(iatom) /= iatom) THEN
            bs_env%ri_rs%atomic_grids(iatom)%raw_points(:, :) = &
               bs_env%ri_rs%atomic_grids(bs_env%ri_rs%grid_opt%grid_rep(iatom))%raw_points
         END IF
      END DO
   END SUBROUTINE broadcast_ri_rs_grids

! **************************************************************************************************
!> \brief RS_AO_RATIO α of one atomic kind in Eq. (4).
!>
!> A positive RS_AO_RATIO applies to every kind. The default (a negative value) is 7 points per AO
!> function, except for full-k periodic GW with the small augmented MOLOPT basis sets: 14 for
!> aug-SZV-MOLOPT-GTH-tier-2, 9 for aug-DZVP-MOLOPT-GTH-tier-1 and 8 for aug-DZVP-MOLOPT-GTH-tier-2.
!> \param bs_env ...
!> \param ikind ...
!> \return α of the kind
! **************************************************************************************************
   FUNCTION rs_ao_ratio_of_kind(bs_env, ikind) RESULT(ratio)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: ikind
      REAL(KIND=dp)                                      :: ratio

      CHARACTER(LEN=default_string_length)               :: basis_name

      ratio = bs_env%ri_rs%grid_opt%rs_ao_ratio
      IF (ratio > 0.0_dp) RETURN
      ratio = 7.0_dp
      IF (bs_env%gw_implementation /= ri_rs_small_cell_full_kp) RETURN
      basis_name = bs_env%basis_set_AO(ikind)%gto_basis_set%name
      CALL uppercase(basis_name)
      SELECT CASE (TRIM(basis_name))
      CASE ("AUG-SZV-MOLOPT-GTH-TIER-2")
         ratio = 14.0_dp
      CASE ("AUG-DZVP-MOLOPT-GTH-TIER-1")
         ratio = 9.0_dp
      CASE ("AUG-DZVP-MOLOPT-GTH-TIER-2")
         ratio = 8.0_dp
      END SELECT
   END FUNCTION rs_ao_ratio_of_kind

! **************************************************************************************************
!> \brief Translation orbits of the atoms (grids of a primitive-cell atom copied to its images).
!>
!> With fractional coordinates s_A = h^-1 R_A, a vector t is a lattice translation of the
!> structure if every atom A is mapped onto an atom B(A,t) of the same kind,
!>
!> (T1) | h (s_A + t - s_B) | < eps_trans,   s_A + t - s_B reduced modulo 1 along periodic
!>                                            directions only; t has no non-periodic component.
!>
!> The candidates are t = s_B - s_1 for every atom B of the kind of atom 1; the valid ones form
!> the translation group T (including t = 0). The orbit of A is {B(A,t) : t in T} and its
!> representative is the smallest atom index in it,
!>
!> (T2) grid_rep(A) = min_{t in T} B(A,t),
!>
!> so grid_rep(A) <= A. A Gamma-only n1 x n2 x n3 supercell of a primitive cell has
!> |T| = n1 n2 n3 and the representatives are the primitive-cell atoms.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE find_translation_orbits(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER                                            :: iatom, natom, ncenters, ntrans
      INTEGER, ALLOCATABLE                               :: images(:, :), rep(:)

      natom = bs_env%n_atom
      CALL translation_atom_images(bs_env, images)
      ntrans = SIZE(images, 2)
      ! Eq. (T2)
      ALLOCATE (rep(natom))
      DO iatom = 1, natom
         rep(iatom) = MINVAL(images(iatom, :))
      END DO
      DEALLOCATE (images)

      ncenters = COUNT([(rep(iatom) == iatom, iatom=1, natom)])
      IF (ALLOCATED(bs_env%ri_rs%grid_opt%grid_rep)) DEALLOCATE (bs_env%ri_rs%grid_opt%grid_rep)
      IF (ALLOCATED(bs_env%ri_rs%grid_opt%grid_centers)) DEALLOCATE (bs_env%ri_rs%grid_opt%grid_centers)
      ALLOCATE (bs_env%ri_rs%grid_opt%grid_centers(ncenters))
      bs_env%ri_rs%grid_opt%grid_centers(:) = PACK([(iatom, iatom=1, natom)], &
                                                   [(rep(iatom) == iatom, iatom=1, natom)])
      CALL MOVE_ALLOC(rep, bs_env%ri_rs%grid_opt%grid_rep)
      ! The translations act freely, so every orbit has |T| atoms.
      CPASSERT(ncenters*ntrans == natom)

      IF (bs_env%unit_nr > 0 .AND. ntrans > 1) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T71,I10)') &
            'RI-RS grid: primitive cells in the supercell (translations)', ntrans
      END IF

   END SUBROUTINE find_translation_orbits

! **************************************************************************************************
!> \brief Atom maps B(A,t) of the lattice translations t in T of Eq. (T1).
!>
!> Atoms are binned on an m^3 cell list in fractional coordinates (m ~ N^(1/3)), so each
!> candidate t = s_B - s_1 is checked in O(N_atom) and the whole search costs O(N_atom^2).
!> \param bs_env ...
!> \param images images(A, t) = B(A, t); column 1 is the identity t = 0
! **************************************************************************************************
   SUBROUTINE translation_atom_images(bs_env, images)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: images

      REAL(KIND=dp), PARAMETER                           :: eps_trans = 1.0E-5_dp

      INTEGER                                            :: d1, d2, d3, iatom, icand, jatom, natom, &
                                                            nbin, ntrans
      INTEGER, ALLOCATABLE                               :: head(:), image(:), kinds(:), next(:), &
                                                            trans_images(:, :)
      INTEGER, DIMENSION(3)                              :: ib, jb
      LOGICAL                                            :: found, valid
      LOGICAL, DIMENSION(3)                              :: periodic
      REAL(KIND=dp), ALLOCATABLE                         :: s(:, :)
      REAL(KIND=dp), DIMENSION(3)                        :: d, q, t
      TYPE(cell_type), POINTER                           :: cell

      natom = bs_env%n_atom
      ntrans = 1
      cell => bs_env%ri_rs%cell
      periodic(:) = cell%perd(:) == 1
      IF (natom > 1 .AND. ANY(periodic)) THEN
         ALLOCATE (s(3, natom), kinds(natom), next(natom), image(natom))
         DO iatom = 1, natom
            s(:, iatom) = MATMUL(cell%h_inv, bs_env%ri_rs%particle_set(iatom)%r)
            kinds(iatom) = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
         END DO
         ! each translation maps atom 1 onto a different atom of its kind
         ALLOCATE (trans_images(natom, COUNT(kinds == kinds(1))))
         trans_images(:, 1) = [(iatom, iatom=1, natom)]
         nbin = MAX(1, FLOOR(REAL(natom, dp)**(1.0_dp/3.0_dp)))
         ALLOCATE (head(nbin**3))
         head(:) = 0
         DO iatom = 1, natom
            ib(:) = MIN(nbin - 1, FLOOR(MODULO(s(:, iatom), 1.0_dp)*nbin))
            next(iatom) = head(1 + ib(1) + nbin*(ib(2) + nbin*ib(3)))
            head(1 + ib(1) + nbin*(ib(2) + nbin*ib(3))) = iatom
         END DO
         DO icand = 2, natom
            IF (kinds(icand) /= kinds(1)) CYCLE
            t(:) = s(:, icand) - s(:, 1)
            ! no translations along non-periodic directions
            IF (ANY(.NOT. periodic .AND. ABS(t)*NORM2(cell%hmat, DIM=1) > eps_trans)) CYCLE
            valid = .TRUE.
            DO iatom = 1, natom
               q(:) = s(:, iatom) + t
               ib(:) = MIN(nbin - 1, FLOOR(MODULO(q, 1.0_dp)*nbin))
               found = .FALSE.
               ! search the 27 neighbouring bins of q for B with Eq. (T1)
               DO d3 = -1, 1
                  DO d2 = -1, 1
                     DO d1 = -1, 1
                        jb(:) = MODULO(ib + [d1, d2, d3], nbin)
                        jatom = head(1 + jb(1) + nbin*(jb(2) + nbin*jb(3)))
                        DO WHILE (jatom > 0 .AND. .NOT. found)
                           IF (kinds(jatom) == kinds(iatom)) THEN
                              d(:) = q - s(:, jatom)
                              WHERE (periodic) d = d - ANINT(d)
                              IF (NORM2(MATMUL(cell%hmat, d)) < eps_trans) THEN
                                 found = .TRUE.
                                 image(iatom) = jatom
                              END IF
                           END IF
                           jatom = next(jatom)
                        END DO
                     END DO
                  END DO
               END DO
               IF (.NOT. found) THEN
                  valid = .FALSE.
                  EXIT
               END IF
            END DO
            IF (valid) THEN
               ntrans = ntrans + 1
               trans_images(:, ntrans) = image
            END IF
         END DO
         DEALLOCATE (s, kinds, next, image, head)
         ALLOCATE (images(natom, ntrans), SOURCE=trans_images(:, 1:ntrans))
         DEALLOCATE (trans_images)
      ELSE
         ALLOCATE (images(natom, 1))
         images(:, 1) = [(iatom, iatom=1, natom)]
      END IF

   END SUBROUTINE translation_atom_images

END MODULE ri_rs_grid_initialization
