!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to compute the Coulomb integral V_(alpha beta)(k) for a k-point k using lattice
!>        summation in real space. These integrals are e.g. needed in periodic RI for RPA, GW
!> \par History
!>       2018.05 created [Jan Wilhelm]
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE kpoint_coulomb_2c
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc,&
                                              get_cell,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release_p, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE generic_shg_integrals,           ONLY: int_operators_r12_ab_shg
   USE generic_shg_integrals_init,      ONLY: contraction_matrix_shg
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: fourpi,&
                                              gaussi,&
                                              pi,&
                                              twopi,&
                                              z_one
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_coulomb_2c'

   PUBLIC :: build_2c_coulomb_matrix_kp, build_2c_coulomb_matrix_kp_small_cell, &
             add_ewald_reciprocal_fm, gaussian_basis_fourier

   ! Reciprocal-space cutoff of the Ewald sum: exp(-G^2/(4 omega^2)) is kept down to
   ! exp(-ewald_tail_exponent), with ewald_lmax_margin more e-folds per unit of angular momentum
   ! for the polynomial prefactors of the Fourier-transformed Gaussians (as in cp2k#6447).
   REAL(KIND=dp), PARAMETER, PRIVATE :: ewald_tail_exponent = 50.0_dp, ewald_lmax_margin = 4.0_dp

! **************************************************************************************************

   TYPE two_d_util_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)  :: block
   END TYPE two_d_util_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_kp ...
!> \param kpoints ...
!> \param basis_type ...
!> \param cell ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param size_lattice_sum ...
!> \param operator_type ...
!> \param ikp_start ...
!> \param ikp_end ...
!> \param omega range-separation parameter for the erf/erfc operators (operator_verf, operator_verfc)
! **************************************************************************************************
   SUBROUTINE build_2c_coulomb_matrix_kp(matrix_v_kp, kpoints, basis_type, cell, particle_set, qs_kind_set, &
                                         atomic_kind_set, size_lattice_sum, operator_type, ikp_start, ikp_end, &
                                         omega)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER                                            :: size_lattice_sum, operator_type, &
                                                            ikp_start, ikp_end
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_2c_coulomb_matrix_kp'

      INTEGER                                            :: handle
      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp

      CALL timeset(routineN, handle)

      CALL allocate_tmp(matrix_v_L_tmp, matrix_v_kp, ikp_start)

      CALL lattice_sum(matrix_v_kp, kpoints, basis_type, cell, particle_set, &
                       qs_kind_set, atomic_kind_set, size_lattice_sum, matrix_v_L_tmp, &
                       operator_type, ikp_start, ikp_end, omega)

      CALL deallocate_tmp(matrix_v_L_tmp)

      CALL timestop(handle)

   END SUBROUTINE build_2c_coulomb_matrix_kp

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_kp ...
!> \param kpoints ...
!> \param basis_type ...
!> \param cell ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param size_lattice_sum ...
!> \param matrix_v_L_tmp ...
!> \param operator_type ...
!> \param ikp_start ...
!> \param ikp_end ...
! **************************************************************************************************
   SUBROUTINE lattice_sum(matrix_v_kp, kpoints, basis_type, cell, particle_set, &
                          qs_kind_set, atomic_kind_set, size_lattice_sum, matrix_v_L_tmp, &
                          operator_type, ikp_start, ikp_end, omega)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER                                            :: size_lattice_sum
      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp
      INTEGER                                            :: operator_type, ikp_start, ikp_end
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: omega

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lattice_sum'

      INTEGER :: factor, handle, handle2, i_block, i_x, i_x_inner, i_x_outer, ik, j_y, j_y_inner, &
         j_y_outer, k_z, k_z_inner, k_z_outer, x_max, x_min, y_max, y_min, z_max, z_min
      INTEGER, DIMENSION(3)                              :: nkp_grid
      REAL(KIND=dp)                                      :: coskl, sinkl
      REAL(KIND=dp), DIMENSION(3)                        :: vec_L, vec_s
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_v_L, blocks_v_L_store
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_v_kp

      CALL timeset(routineN, handle)

      CALL get_factor_and_xyz_min_max(cell, kpoints, size_lattice_sum, factor, hmat, &
                                      x_min, x_max, y_min, y_max, z_min, z_max, nkp_grid)

      CALL allocate_blocks_v_kp(blocks_v_kp, matrix_v_kp, ikp_start, ikp_end)
      CALL allocate_blocks_v_L(blocks_v_L, matrix_v_L_tmp)
      CALL allocate_blocks_v_L(blocks_v_L_store, matrix_v_L_tmp)

      DO i_x_inner = 0, 2*nkp_grid(1) - 1
         DO j_y_inner = 0, 2*nkp_grid(2) - 1
            DO k_z_inner = 0, 2*nkp_grid(3) - 1

               DO i_x_outer = x_min, x_max + nkp_grid(1), 2*nkp_grid(1)
                  DO j_y_outer = y_min, y_max + nkp_grid(2), 2*nkp_grid(2)
                     DO k_z_outer = z_min, z_max + nkp_grid(3), 2*nkp_grid(3)

                        i_x = i_x_inner + i_x_outer
                        j_y = j_y_inner + j_y_outer
                        k_z = k_z_inner + k_z_outer

                        IF (i_x > x_max .OR. i_x < x_min .OR. &
                            j_y > y_max .OR. j_y < y_min .OR. &
                            k_z > z_max .OR. k_z < z_min) CYCLE

                        vec_s = [REAL(i_x, dp), REAL(j_y, dp), REAL(k_z, dp)]

                        vec_L = MATMUL(hmat, vec_s)

                        ! Compute (P 0 | Q vec_L) and store it in matrix_v_L_tmp
                        CALL compute_v_transl(matrix_v_L_tmp, blocks_v_L, vec_L, particle_set, &
                                              qs_kind_set, atomic_kind_set, basis_type, cell, &
                                              operator_type, omega)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             SHARED(blocks_v_L, blocks_v_L_store) &
!$OMP             PRIVATE(i_block)
                        DO i_block = 1, SIZE(blocks_v_L)
                           blocks_v_L_store(i_block)%block(:, :) = blocks_v_L_store(i_block)%block(:, :) &
                                                                   + blocks_v_L(i_block)%block(:, :)
                        END DO
!$OMP END PARALLEL DO

                     END DO
                  END DO
               END DO

               CALL timeset(routineN//"_R_to_k", handle2)

               ! add exp(iq*vec_L) * (P 0 | Q vec_L) to V_PQ(q)
               DO ik = ikp_start, ikp_end

                  ! coskl and sinkl are identical for all i_x_outer, j_y_outer, k_z_outer
                  coskl = COS(twopi*DOT_PRODUCT(vec_s(1:3), kpoints%xkp(1:3, ik)))
                  sinkl = SIN(twopi*DOT_PRODUCT(vec_s(1:3), kpoints%xkp(1:3, ik)))

                  DO i_block = 1, SIZE(blocks_v_L)

                     blocks_v_kp(ik, 1, i_block)%block(:, :) = blocks_v_kp(ik, 1, i_block)%block(:, :) &
                                                               + coskl*blocks_v_L_store(i_block)%block(:, :)
                     blocks_v_kp(ik, 2, i_block)%block(:, :) = blocks_v_kp(ik, 2, i_block)%block(:, :) &
                                                               + sinkl*blocks_v_L_store(i_block)%block(:, :)

                  END DO

               END DO

               DO i_block = 1, SIZE(blocks_v_L)

                  blocks_v_L_store(i_block)%block(:, :) = 0.0_dp

               END DO

               CALL timestop(handle2)

            END DO
         END DO
      END DO

      CALL set_blocks_to_matrix_v_kp(matrix_v_kp, blocks_v_kp, ikp_start, ikp_end)

      CALL deallocate_blocks_v_kp(blocks_v_kp)
      CALL deallocate_blocks_v_L(blocks_v_L)
      CALL deallocate_blocks_v_L(blocks_v_L_store)

      CALL timestop(handle)

   END SUBROUTINE lattice_sum

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_kp ...
!> \param blocks_v_kp ...
!> \param ikp_start ...
!> \param ikp_end ...
! **************************************************************************************************
   SUBROUTINE set_blocks_to_matrix_v_kp(matrix_v_kp, blocks_v_kp, ikp_start, ikp_end)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_v_kp
      INTEGER                                            :: ikp_start, ikp_end

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_blocks_to_matrix_v_kp'

      INTEGER                                            :: col, handle, i_block, i_real_im, ik, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      DO ik = ikp_start, ikp_end

         DO i_real_im = 1, 2

            i_block = 1

            CALL dbcsr_iterator_start(iter, matrix_v_kp(ik, i_real_im)%matrix)

            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row, col, data_block)

               data_block(:, :) = blocks_v_kp(ik, i_real_im, i_block)%block(:, :)

               i_block = i_block + 1

            END DO

            CALL dbcsr_iterator_stop(iter)

         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE set_blocks_to_matrix_v_kp

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_L_tmp ...
!> \param blocks_v_L ...
!> \param vec_L ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param basis_type ...
!> \param cell ...
!> \param operator_type ...
! **************************************************************************************************
   SUBROUTINE compute_v_transl(matrix_v_L_tmp, blocks_v_L, vec_L, particle_set, &
                               qs_kind_set, atomic_kind_set, basis_type, cell, operator_type, omega)

      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_v_L
      REAL(KIND=dp), DIMENSION(3)                        :: vec_L
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(cell_type), POINTER                           :: cell
      INTEGER                                            :: operator_type
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: omega

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_v_transl'

      INTEGER                                            :: col, handle, i_block, kind_a, kind_b, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(dp), DIMENSION(3)                             :: ra, rab_L, rb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: contr_a, contr_b
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b

      CALL timeset(routineN, handle)

      NULLIFY (basis_set_a, basis_set_b, data_block)

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      CALL dbcsr_set(matrix_v_L_tmp, 0.0_dp)

      i_block = 1

      CALL dbcsr_iterator_start(iter, matrix_v_L_tmp)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         kind_a = kind_of(row)
         kind_b = kind_of(col)

         CALL get_qs_kind(qs_kind=qs_kind_set(kind_a), basis_set=basis_set_a, basis_type=basis_type)
         CALL get_qs_kind(qs_kind=qs_kind_set(kind_b), basis_set=basis_set_b, basis_type=basis_type)

         ra(1:3) = pbc(particle_set(row)%r(1:3), cell)
         rb(1:3) = pbc(particle_set(col)%r(1:3), cell)

         rab_L(1:3) = rb(1:3) - ra(1:3) + vec_L(1:3)

         CALL contraction_matrix_shg(basis_set_a, contr_a)
         CALL contraction_matrix_shg(basis_set_b, contr_b)

         blocks_v_L(i_block)%block = 0.0_dp

         CALL int_operators_r12_ab_shg(operator_type, blocks_v_L(i_block)%block, rab=rab_L, &
                                       fba=basis_set_a, fbb=basis_set_b, scona_shg=contr_a, sconb_shg=contr_b, &
                                       omega=omega, calculate_forces=.FALSE.)

         i_block = i_block + 1

         DEALLOCATE (contr_a, contr_b)

      END DO

      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (kind_of)

      CALL timestop(handle)

   END SUBROUTINE compute_v_transl

! **************************************************************************************************
!> \brief ...
!> \param blocks_v_kp ...
! **************************************************************************************************
   SUBROUTINE deallocate_blocks_v_kp(blocks_v_kp)
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_v_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_blocks_v_kp'

      INTEGER                                            :: handle, i_block, i_real_img, ik

      CALL timeset(routineN, handle)

      DO ik = LBOUND(blocks_v_kp, 1), UBOUND(blocks_v_kp, 1)
         DO i_real_img = 1, SIZE(blocks_v_kp, 2)
            DO i_block = 1, SIZE(blocks_v_kp, 3)
               DEALLOCATE (blocks_v_kp(ik, i_real_img, i_block)%block)
            END DO
         END DO
      END DO

      DEALLOCATE (blocks_v_kp)

      CALL timestop(handle)

   END SUBROUTINE deallocate_blocks_v_kp

! **************************************************************************************************
!> \brief ...
!> \param blocks_v_L ...
! **************************************************************************************************
   SUBROUTINE deallocate_blocks_v_L(blocks_v_L)
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_v_L

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_blocks_v_L'

      INTEGER                                            :: handle, i_block

      CALL timeset(routineN, handle)

      DO i_block = 1, SIZE(blocks_v_L, 1)
         DEALLOCATE (blocks_v_L(i_block)%block)
      END DO

      DEALLOCATE (blocks_v_L)

      CALL timestop(handle)

   END SUBROUTINE deallocate_blocks_v_L

! **************************************************************************************************
!> \brief ...
!> \param blocks_v_L ...
!> \param matrix_v_L_tmp ...
! **************************************************************************************************
   SUBROUTINE allocate_blocks_v_L(blocks_v_L, matrix_v_L_tmp)
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_v_L
      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_blocks_v_L'

      INTEGER                                            :: col, handle, i_block, nblocks, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nblocks = 0

      CALL dbcsr_iterator_start(iter, matrix_v_L_tmp)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         nblocks = nblocks + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (blocks_v_L(nblocks))

      i_block = 1

      CALL dbcsr_iterator_start(iter, matrix_v_L_tmp)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         ALLOCATE (blocks_v_L(i_block)%block(SIZE(data_block, 1), SIZE(data_block, 2)))
         blocks_v_L(i_block)%block = 0.0_dp

         i_block = i_block + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE allocate_blocks_v_L

! **************************************************************************************************
!> \brief ...
!> \param blocks_v_kp ...
!> \param matrix_v_kp ...
!> \param ikp_start ...
!> \param ikp_end ...
! **************************************************************************************************
   SUBROUTINE allocate_blocks_v_kp(blocks_v_kp, matrix_v_kp, ikp_start, ikp_end)
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_v_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      INTEGER                                            :: ikp_start, ikp_end

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_blocks_v_kp'

      INTEGER                                            :: col, handle, i_block, i_real_img, ik, &
                                                            nblocks, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nblocks = 0

      CALL dbcsr_iterator_start(iter, matrix_v_kp(ikp_start, 1)%matrix)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         nblocks = nblocks + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (blocks_v_kp(ikp_start:ikp_end, SIZE(matrix_v_kp, 2), nblocks))

      DO ik = ikp_start, ikp_end

         DO i_real_img = 1, SIZE(matrix_v_kp, 2)

            i_block = 1

            CALL dbcsr_iterator_start(iter, matrix_v_kp(ik, i_real_img)%matrix)

            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row, col, data_block)

               ALLOCATE (blocks_v_kp(ik, i_real_img, i_block)%block(SIZE(data_block, 1), &
                                                                    SIZE(data_block, 2)))
               blocks_v_kp(ik, i_real_img, i_block)%block = 0.0_dp

               i_block = i_block + 1

            END DO

            CALL dbcsr_iterator_stop(iter)

         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_blocks_v_kp

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param kpoints ...
!> \param size_lattice_sum ...
!> \param factor ...
!> \param hmat ...
!> \param x_min ...
!> \param x_max ...
!> \param y_min ...
!> \param y_max ...
!> \param z_min ...
!> \param z_max ...
!> \param nkp_grid ...
! **************************************************************************************************
   SUBROUTINE get_factor_and_xyz_min_max(cell, kpoints, size_lattice_sum, factor, hmat, &
                                         x_min, x_max, y_min, y_max, z_min, z_max, nkp_grid)

      TYPE(cell_type), POINTER                           :: cell
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: size_lattice_sum, factor
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      INTEGER                                            :: x_min, x_max, y_min, y_max, z_min, z_max
      INTEGER, DIMENSION(3)                              :: nkp_grid

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_factor_and_xyz_min_max'

      INTEGER                                            :: handle, nkp
      INTEGER, DIMENSION(3)                              :: periodic

      CALL timeset(routineN, handle)

      CALL get_kpoint_info(kpoints, nkp_grid=nkp_grid, nkp=nkp)
      CALL get_cell(cell=cell, h=hmat, periodic=periodic)

      IF (periodic(1) == 0) THEN
         CPASSERT(nkp_grid(1) == 1)
      END IF
      IF (periodic(2) == 0) THEN
         CPASSERT(nkp_grid(2) == 1)
      END IF
      IF (periodic(3) == 0) THEN
         CPASSERT(nkp_grid(3) == 1)
      END IF

      IF (MODULO(nkp_grid(1), 2) == 1) THEN
         factor = 3**(size_lattice_sum - 1)
      ELSE IF (MODULO(nkp_grid(1), 2) == 0) THEN
         factor = 2**(size_lattice_sum - 1)
      END IF

      IF (MODULO(nkp_grid(1), 2) == 1) THEN
         x_min = -(factor*nkp_grid(1) - 1)/2
         x_max = (factor*nkp_grid(1) - 1)/2
      ELSE IF (MODULO(nkp_grid(1), 2) == 0) THEN
         x_min = -factor*nkp_grid(1)/2
         x_max = factor*nkp_grid(1)/2 - 1
      END IF
      IF (periodic(1) == 0) THEN
         x_min = 0
         x_max = 0
      END IF

      IF (MODULO(nkp_grid(2), 2) == 1) THEN
         y_min = -(factor*nkp_grid(2) - 1)/2
         y_max = (factor*nkp_grid(2) - 1)/2
      ELSE IF (MODULO(nkp_grid(2), 2) == 0) THEN
         y_min = -factor*nkp_grid(2)/2
         y_max = factor*nkp_grid(2)/2 - 1
      END IF
      IF (periodic(2) == 0) THEN
         y_min = 0
         y_max = 0
      END IF

      IF (MODULO(nkp_grid(3), 2) == 1) THEN
         z_min = -(factor*nkp_grid(3) - 1)/2
         z_max = (factor*nkp_grid(3) - 1)/2
      ELSE IF (MODULO(nkp_grid(3), 2) == 0) THEN
         z_min = -factor*nkp_grid(3)/2
         z_max = factor*nkp_grid(3)/2 - 1
      END IF
      IF (periodic(3) == 0) THEN
         z_min = 0
         z_max = 0
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_factor_and_xyz_min_max

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_L_tmp ...
!> \param matrix_v_kp ...
!> \param ikp_start ...
! **************************************************************************************************
   SUBROUTINE allocate_tmp(matrix_v_L_tmp, matrix_v_kp, ikp_start)

      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      INTEGER                                            :: ikp_start

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'allocate_tmp'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY (matrix_v_L_tmp)
      CALL dbcsr_init_p(matrix_v_L_tmp)
      CALL dbcsr_create(matrix=matrix_v_L_tmp, &
                        template=matrix_v_kp(ikp_start, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_reserve_all_blocks(matrix_v_L_tmp)
      CALL dbcsr_set(matrix_v_L_tmp, 0.0_dp)

      CALL timestop(handle)

   END SUBROUTINE allocate_tmp

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_L_tmp ...
! **************************************************************************************************
   SUBROUTINE deallocate_tmp(matrix_v_L_tmp)

      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'deallocate_tmp'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbcsr_release_p(matrix_v_L_tmp)

      CALL timestop(handle)

   END SUBROUTINE deallocate_tmp

! **************************************************************************************************
!> \brief ...
!> \param V_k ...
!> \param qs_env ...
!> \param kpoints ...
!> \param size_lattice_sum ...
!> \param basis_type ...
!> \param ikp_start ...
!> \param ikp_end ...
! **************************************************************************************************
   SUBROUTINE build_2c_coulomb_matrix_kp_small_cell(V_k, qs_env, kpoints, size_lattice_sum, &
                                                    basis_type, ikp_start, ikp_end)
      COMPLEX(KIND=dp), DIMENSION(:, :, :)               :: V_k
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: size_lattice_sum
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      INTEGER                                            :: ikp_start, ikp_end

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_2c_coulomb_matrix_kp_small_cell'

      INTEGER :: factor, handle, handle2, i_cell, i_x, i_x_inner, i_x_outer, ik, ikp_local, j_y, &
         j_y_inner, j_y_outer, k_z, k_z_inner, k_z_outer, n_atom, n_bf, x_max, x_min, y_max, &
         y_min, z_max, z_min
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bf_end_from_atom, bf_start_from_atom
      INTEGER, DIMENSION(3)                              :: nkp_grid
      REAL(KIND=dp)                                      :: coskl, sinkl
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: V_L
      REAL(KIND=dp), DIMENSION(3)                        :: vec_L, vec_s
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      CALL get_factor_and_xyz_min_max(cell, kpoints, size_lattice_sum, factor, hmat, &
                                      x_min, x_max, y_min, y_max, z_min, z_max, nkp_grid)

      CALL get_basis_sizes(qs_env, n_atom, basis_type, bf_start_from_atom, bf_end_from_atom, n_bf)

      ALLOCATE (V_L(n_bf, n_bf))

      DO i_x_inner = 0, 2*nkp_grid(1) - 1
         DO j_y_inner = 0, 2*nkp_grid(2) - 1
            DO k_z_inner = 0, 2*nkp_grid(3) - 1

               V_L(:, :) = 0.0_dp
               i_cell = 0

               DO i_x_outer = x_min, x_max + nkp_grid(1), 2*nkp_grid(1)
                  DO j_y_outer = y_min, y_max + nkp_grid(2), 2*nkp_grid(2)
                     DO k_z_outer = z_min, z_max + nkp_grid(3), 2*nkp_grid(3)

                        i_x = i_x_inner + i_x_outer
                        j_y = j_y_inner + j_y_outer
                        k_z = k_z_inner + k_z_outer

                        IF (i_x > x_max .OR. i_x < x_min .OR. &
                            j_y > y_max .OR. j_y < y_min .OR. &
                            k_z > z_max .OR. k_z < z_min) CYCLE

                        i_cell = i_cell + 1

                        vec_s = [REAL(i_x, dp), REAL(j_y, dp), REAL(k_z, dp)]

                        IF (MODULO(i_cell, para_env%num_pe) /= para_env%mepos) CYCLE

                        vec_L = MATMUL(hmat, vec_s)

                        ! Compute (P 0 | Q vec_L) and add it to V_R
                        CALL add_V_L(V_L, vec_L, n_atom, bf_start_from_atom, bf_end_from_atom, &
                                     particle_set, qs_kind_set, atomic_kind_set, basis_type, cell)

                     END DO
                  END DO
               END DO

               CALL para_env%sync()
               CALL para_env%sum(V_L)

               CALL timeset(routineN//"_R_to_k", handle2)

               ikp_local = 0

               ! add exp(iq*vec_L) * (P 0 | Q vec_L) to V_PQ(q)
               DO ik = 1, ikp_end

                  IF (MODULO(ik, para_env%num_pe) /= para_env%mepos) CYCLE

                  ikp_local = ikp_local + 1

                  IF (ik < ikp_start) CYCLE

                  ! coskl and sinkl are identical for all i_x_outer, j_y_outer, k_z_outer
                  coskl = COS(twopi*DOT_PRODUCT(vec_s(1:3), kpoints%xkp(1:3, ik)))
                  sinkl = SIN(twopi*DOT_PRODUCT(vec_s(1:3), kpoints%xkp(1:3, ik)))

                  V_k(:, :, ikp_local) = V_k(:, :, ikp_local) + z_one*coskl*V_L(:, :) + &
                                         gaussi*sinkl*V_L(:, :)

               END DO

               CALL timestop(handle2)

            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_2c_coulomb_matrix_kp_small_cell

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param n_atom ...
!> \param basis_type ...
!> \param bf_start_from_atom ...
!> \param bf_end_from_atom ...
!> \param n_bf ...
! **************************************************************************************************
   SUBROUTINE get_basis_sizes(qs_env, n_atom, basis_type, bf_start_from_atom, bf_end_from_atom, n_bf)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: n_atom
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bf_start_from_atom, bf_end_from_atom
      INTEGER                                            :: n_bf

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_basis_sizes'

      INTEGER                                            :: handle, iatom, ikind, n_kind, nsgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      n_atom = SIZE(particle_set)
      n_kind = SIZE(qs_kind_set)

      DO ikind = 1, n_kind
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), basis_set=basis_set_a, &
                          basis_type=basis_type)
         CPASSERT(ASSOCIATED(basis_set_a))
      END DO

      ALLOCATE (bf_start_from_atom(n_atom), bf_end_from_atom(n_atom))

      n_bf = 0
      DO iatom = 1, n_atom
         bf_start_from_atom(iatom) = n_bf + 1
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type=basis_type)
         n_bf = n_bf + nsgf
         bf_end_from_atom(iatom) = n_bf
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_basis_sizes

! **************************************************************************************************
!> \brief ...
!> \param V_L ...
!> \param vec_L ...
!> \param n_atom ...
!> \param bf_start_from_atom ...
!> \param bf_end_from_atom ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param basis_type ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE add_V_L(V_L, vec_L, n_atom, bf_start_from_atom, bf_end_from_atom, &
                      particle_set, qs_kind_set, atomic_kind_set, basis_type, cell)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: V_L
      REAL(KIND=dp), DIMENSION(3)                        :: vec_L
      INTEGER                                            :: n_atom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bf_start_from_atom, bf_end_from_atom
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_V_L'

      INTEGER                                            :: a_1, a_2, atom_a, atom_b, b_1, b_2, &
                                                            handle, kind_a, kind_b
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(dp), DIMENSION(3)                             :: ra, rab_L, rb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: V_L_ab
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: contr_a, contr_b
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b

      CALL timeset(routineN, handle)

      NULLIFY (basis_set_a, basis_set_b)

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      DO atom_a = 1, n_atom

         DO atom_b = 1, n_atom

            kind_a = kind_of(atom_a)
            kind_b = kind_of(atom_b)

            CALL get_qs_kind(qs_kind=qs_kind_set(kind_a), basis_set=basis_set_a, &
                             basis_type=basis_type)
            CALL get_qs_kind(qs_kind=qs_kind_set(kind_b), basis_set=basis_set_b, &
                             basis_type=basis_type)

            ra(1:3) = pbc(particle_set(atom_a)%r(1:3), cell)
            rb(1:3) = pbc(particle_set(atom_b)%r(1:3), cell)

            rab_L(1:3) = rb(1:3) - ra(1:3) + vec_L(1:3)

            CALL contraction_matrix_shg(basis_set_a, contr_a)
            CALL contraction_matrix_shg(basis_set_b, contr_b)

            a_1 = bf_start_from_atom(atom_a)
            a_2 = bf_end_from_atom(atom_a)
            b_1 = bf_start_from_atom(atom_b)
            b_2 = bf_end_from_atom(atom_b)

            ALLOCATE (V_L_ab(a_2 - a_1 + 1, b_2 - b_1 + 1))

            CALL int_operators_r12_ab_shg(operator_coulomb, V_L_ab, rab=rab_L, &
                                          fba=basis_set_a, fbb=basis_set_b, &
                                          scona_shg=contr_a, sconb_shg=contr_b, &
                                          calculate_forces=.FALSE.)

            V_L(a_1:a_2, b_1:b_2) = V_L(a_1:a_2, b_1:b_2) + V_L_ab(:, :)

            DEALLOCATE (contr_a, contr_b, V_L_ab)

         END DO

      END DO

      DEALLOCATE (kind_of)

      CALL timestop(handle)

   END SUBROUTINE add_V_L

! **************************************************************************************************
!> \brief Adds the reciprocal-space part of the Ewald sum for the 3D periodic Coulomb matrix at
!>        q = 0 to a distributed full matrix, with the G = 0 Fourier coefficient of 1/r omitted:
!>
!>          V_PQ(0) <- V_PQ(0) + (4π/Ω) Σ_{G≠0} exp(−G²/4ω²)/G² ρ_P(G) ρ_Q(G)* − π/(Ω ω²) ρ_P(0) ρ_Q(0)
!>
!>        with ρ_P(G) = ∫dr ϕ_P(r) exp(iG·r) (gaussian_basis_fourier). If the matrix holds the
!>        real-space lattice sum Σ_L (P0|erfc(ωr)/r|QL) on entry, it holds the periodic Coulomb
!>        matrix "with the singular G = 0 term omitted" on exit: 1/r = erfc(ωr)/r + erf(ωr)/r, the
!>        erf part has the Fourier transform (4π/G²) exp(−G²/4ω²) = 4π/G² − π/ω² + O(G²), and
!>        dropping its 4π/G² singularity leaves −π/ω² at G = 0 (same construction as the 3D
!>        Γ-point Coulomb matrix of cp2k#6447, there on a replicated matrix). The result is
!>        independent of ω, which is a check on both halves.
!>
!>        Every rank adds the sum over ALL G to its own local block, from ρ(G) evaluated for the
!>        atoms that own one of its local rows or columns only: no replicated n_bf x n_bf array and
!>        no communication. Σ_G over ±G makes the sum real, so Re(ρ_P ρ_Q*) is accumulated per G.
!> \param fm_V distributed n_bf x n_bf matrix of the basis `basis_type`, updated in place
!> \param qs_env QS environment, for cell, particles and the basis
!> \param basis_type basis of the matrix, e.g. "RI_AUX"
!> \param omega Ewald splitting parameter (1/bohr)
!> \param n_wave number of reciprocal lattice vectors inside the cutoff (diagnostic)
! **************************************************************************************************
   SUBROUTINE add_ewald_reciprocal_fm(fm_V, qs_env, basis_type, omega, n_wave)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_V
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      REAL(KIND=dp), INTENT(IN)                          :: omega
      INTEGER, INTENT(OUT)                               :: n_wave

      CHARACTER(LEN=*), PARAMETER :: routineN = 'add_ewald_reciprocal_fm'
      INTEGER, PARAMETER                                 :: n_block = 128

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: fourier
      INTEGER                                            :: handle, i_atom, i_col, i_g, i_row, &
                                                            idim, ix, iy, iz, n_atom, n_bf, &
                                                            n_col, ncol_local, nrow_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bf_end, bf_start
      INTEGER, DIMENSION(3)                              :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: atom_needed
      REAL(KIND=dp)                                      :: cutoff, p2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: prefactor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cols_ri, p_batch, rows_ri
      REAL(KIND=dp), DIMENSION(3, 3)                     :: reciprocal
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds

      CALL timeset(routineN, handle)

      NULLIFY (cell, particles, qs_kinds)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particles, qs_kind_set=qs_kinds)
      CPASSERT(COUNT(cell%perd /= 0) == 3)
      CALL get_basis_sizes(qs_env, n_atom, basis_type, bf_start, bf_end, n_bf)
      CALL cp_fm_get_info(fm_V, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      ! Only the atoms carrying a local row or column need their ρ(G).
      ALLOCATE (atom_needed(n_atom))
      atom_needed(:) = .FALSE.
      DO i_atom = 1, n_atom
         IF (ANY(row_indices(1:nrow_local) >= bf_start(i_atom) .AND. &
                 row_indices(1:nrow_local) <= bf_end(i_atom))) atom_needed(i_atom) = .TRUE.
         IF (ANY(col_indices(1:ncol_local) >= bf_start(i_atom) .AND. &
                 col_indices(1:ncol_local) <= bf_end(i_atom))) atom_needed(i_atom) = .TRUE.
      END DO

      ! b_i = 2π (h⁻¹)ᵀ columns; |G| ≤ cutoff with |n_i| ≤ cutoff |a_i|/2π from n_i = G·a_i/2π.
      reciprocal = twopi*TRANSPOSE(cell%h_inv)
      cutoff = 2.0_dp*omega*SQRT(ewald_tail_exponent + ewald_lmax_margin*basis_lmax(qs_kinds, basis_type))
      DO idim = 1, 3
         bounds(idim) = CEILING(cutoff*NORM2(cell%hmat(:, idim))/twopi)
      END DO

      ! Batches of n_block wave vectors: ρ(G) for all needed basis functions, then the local
      ! rank-2n_block update  V_loc += [Re ρ_rows, Im ρ_rows] diag(w) [Re ρ_cols, Im ρ_cols]ᵀ.
      ALLOCATE (fourier(n_bf, n_block), p_batch(3, n_block), prefactor(n_block))
      ALLOCATE (rows_ri(nrow_local, 2*n_block), cols_ri(ncol_local, 2*n_block))
      n_wave = 0
      n_col = 0
      DO ix = -bounds(1), bounds(1)
         DO iy = -bounds(2), bounds(2)
            DO iz = -bounds(3), bounds(3)
               p_batch(:, n_col + 1) = MATMUL(reciprocal, REAL([ix, iy, iz], dp))
               p2 = SUM(p_batch(:, n_col + 1)**2)
               IF (p2 > cutoff**2) CYCLE
               n_col = n_col + 1
               n_wave = n_wave + 1
               IF (p2 <= 0.0_dp) THEN
                  ! G = 0: the finite remainder −π/(Ω ω²) of the erf part, its 4π/G² omitted.
                  prefactor(n_col) = -pi/(cell%deth*omega**2)
               ELSE
                  ! G ≠ 0: (4π/Ω) exp(−G²/4ω²)/G².
                  prefactor(n_col) = fourpi*EXP(-p2/(4.0_dp*omega**2))/(cell%deth*p2)
               END IF
               IF (n_col == n_block) THEN
                  CALL ewald_flush_batch(n_col)
                  n_col = 0
               END IF
            END DO
         END DO
      END DO
      IF (n_col > 0) CALL ewald_flush_batch(n_col)

      DEALLOCATE (fourier, p_batch, prefactor, rows_ri, cols_ri, atom_needed, bf_start, bf_end)

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief ρ(G) of the batch for the needed atoms and the local update of fm_V with it.
!> \param n number of wave vectors in the batch
! **************************************************************************************************
      SUBROUTINE ewald_flush_batch(n)
         INTEGER, INTENT(IN)                                :: n

         INTEGER                                            :: i, j

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(n, p_batch, particles, qs_kinds, cell, basis_type, bf_start, &
!$OMP                                  bf_end, fourier, atom_needed) PRIVATE(i)
         DO i = 1, n
            CALL gaussian_basis_fourier(p_batch(:, i), particles, qs_kinds, cell, basis_type, &
                                        bf_start, bf_end, fourier(:, i), atom_mask=atom_needed)
         END DO
!$OMP END PARALLEL DO
         ! rows_ri(:, i) = w_i Re ρ_row(G_i), rows_ri(:, n + i) = w_i Im ρ_row(G_i); cols without w.
         DO i = 1, n
            DO j = 1, nrow_local
               rows_ri(j, i) = prefactor(i)*REAL(fourier(row_indices(j), i), dp)
               rows_ri(j, n + i) = prefactor(i)*AIMAG(fourier(row_indices(j), i))
            END DO
            DO j = 1, ncol_local
               cols_ri(j, i) = REAL(fourier(col_indices(j), i), dp)
               cols_ri(j, n + i) = AIMAG(fourier(col_indices(j), i))
            END DO
         END DO
         ! Re(ρ_P ρ_Q*) = Re ρ_P Re ρ_Q + Im ρ_P Im ρ_Q, accumulated into the local block.
         IF (nrow_local > 0 .AND. ncol_local > 0) THEN
            CALL dgemm('N', 'T', nrow_local, ncol_local, 2*n, 1.0_dp, rows_ri, nrow_local, &
                       cols_ri, ncol_local, 1.0_dp, fm_V%local_data, SIZE(fm_V%local_data, 1))
         END IF

      END SUBROUTINE ewald_flush_batch

   END SUBROUTINE add_ewald_reciprocal_fm

! **************************************************************************************************
!> \brief Largest angular momentum of a basis type over all kinds.
!> \param qs_kinds ...
!> \param basis_type ...
!> \return ...
! **************************************************************************************************
   FUNCTION basis_lmax(qs_kinds, basis_type) RESULT(lmax)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      INTEGER                                            :: lmax

      INTEGER                                            :: ikind
      TYPE(gto_basis_set_type), POINTER                  :: basis

      lmax = 0
      DO ikind = 1, SIZE(qs_kinds)
         CALL get_qs_kind(qs_kinds(ikind), basis_set=basis, basis_type=basis_type)
         lmax = MAX(lmax, MAXVAL(basis%lmax))
      END DO
   END FUNCTION basis_lmax

! **************************************************************************************************
!> \brief Fourier transform ρ_P(p) = ∫dr ϕ_P(r) exp(ip·r) of every contracted spherical basis
!>        function of one basis type (ported from cp2k#6447).
!>
!>        For a Cartesian primitive x^i y^j z^k exp(−α|r − R|²) the integral factorizes into the
!>        one-dimensional moments m_l(p_x) = ∫dx x^l exp(−αx²) exp(ip_x x), with
!>          m_0 = √(π/α) exp(−p²/4α),  m_1 = ip m_0/(2α),  m_l = (ip m_{l−1} + (l−1) m_{l−2})/(2α),
!>        times the centre phase exp(ip·R); the contraction matrix sphi carries the primitives to
!>        the spherical functions. At p = 0 this is the charge ∫ϕ_P.
!> \param p wave vector (1/bohr)
!> \param particles ...
!> \param qs_kinds ...
!> \param cell ...
!> \param basis_type ...
!> \param bf_start first basis function of every atom
!> \param bf_end last basis function of every atom
!> \param values ρ_P(p), one entry per basis function
!> \param atom_mask if present, atoms with .FALSE. are skipped (their values are left zero)
! **************************************************************************************************
   SUBROUTINE gaussian_basis_fourier(p, particles, qs_kinds, cell, basis_type, &
                                     bf_start, bf_end, values, atom_mask)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: p
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds
      TYPE(cell_type), POINTER                           :: cell
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      INTEGER, DIMENSION(:), INTENT(IN)                  :: bf_start, bf_end
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(OUT)        :: values
      LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL        :: atom_mask

      COMPLEX(KIND=dp)                                   :: phase
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: primitive
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: moments
      INTEGER                                            :: iatom, ico, idim, ikind, il, ipgf, iset, &
                                                            isgf, lmax, ncart, offset
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), DIMENSION(3)                        :: centre
      TYPE(gto_basis_set_type), POINTER                  :: basis

      values = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO iatom = 1, SIZE(particles)
         IF (PRESENT(atom_mask)) THEN
            IF (.NOT. atom_mask(iatom)) CYCLE
         END IF
         ikind = particles(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kinds(ikind), basis_set=basis, basis_type=basis_type)
         centre = pbc(particles(iatom)%r, cell)
         phase = EXP(gaussi*SUM(p*centre))
         CPASSERT(bf_end(iatom) - bf_start(iatom) + 1 == basis%nsgf)
         DO iset = 1, basis%nset
            lmax = basis%lmax(iset)
            ncart = ncoset(lmax)
            ALLOCATE (moments(0:lmax, 3), primitive(ncart*basis%npgf(iset)))
            primitive = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            DO ipgf = 1, basis%npgf(iset)
               alpha = basis%zet(ipgf, iset)
               DO idim = 1, 3
                  moments(0, idim) = SQRT(pi/alpha)*EXP(-p(idim)**2/(4.0_dp*alpha))
                  IF (lmax > 0) moments(1, idim) = gaussi*p(idim)*moments(0, idim)/(2.0_dp*alpha)
                  DO il = 2, lmax
                     moments(il, idim) = (gaussi*p(idim)*moments(il - 1, idim) + &
                                          REAL(il - 1, dp)*moments(il - 2, idim))/(2.0_dp*alpha)
                  END DO
               END DO
               DO ico = ncoset(basis%lmin(iset) - 1) + 1, ncart
                  primitive((ipgf - 1)*ncart + ico) = &
                     moments(indco(1, ico), 1)*moments(indco(2, ico), 2)*moments(indco(3, ico), 3)
               END DO
            END DO
            offset = basis%first_sgf(1, iset)
            DO isgf = offset, offset + basis%nsgf_set(iset) - 1
               values(bf_start(iatom) + isgf - 1) = &
                  phase*SUM(basis%sphi(1:SIZE(primitive), isgf)*primitive)
            END DO
            DEALLOCATE (moments, primitive)
         END DO
      END DO
   END SUBROUTINE gaussian_basis_fourier

END MODULE kpoint_coulomb_2c
