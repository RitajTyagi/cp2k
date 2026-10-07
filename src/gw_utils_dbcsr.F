!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common DBCSR matrix operations used by GW modules.
!> \par History
!>      09.2026 created Jan Wilhelm
! **************************************************************************************************
MODULE gw_utils_dbcsr
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_release, dbcsr_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: dbcsr_contract_ABA, hadamard_product, hadamard_product_inplace, &
             release_dbcsr_topology_and_matrices, setup_square_topology

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils_dbcsr'

CONTAINS

! **************************************************************************************************
!> \brief Computes the scaled element-wise product C = factor (A ◦ B) while preserving the
!>        block structure of A. Blocks absent from B are retained in C with zero values.
!> \param matrix_A First factor and source of the block structure
!> \param matrix_B Second factor
!> \param matrix_C Scaled element-wise product
!> \param factor Scaling factor
! **************************************************************************************************
   SUBROUTINE hadamard_product(matrix_A, matrix_B, matrix_C, factor)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_C
      REAL(KIND=dp), INTENT(IN)                          :: factor

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hadamard_product'

      INTEGER                                            :: handle, icol, irow
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_B, block_C
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL timeset(routineN, handle)

      CALL dbcsr_copy(matrix_C, matrix_A)
      CALL dbcsr_iterator_start(iterator, matrix_C)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, irow, icol, block_C)
         CALL dbcsr_get_block_p(matrix_B, irow, icol, block_B, found)
         IF (found) THEN
            block_C(:, :) = factor*block_C(:, :)*block_B(:, :)
         ELSE
            block_C(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iterator)

      CALL timestop(handle)

   END SUBROUTINE hadamard_product

! **************************************************************************************************
!> \brief Form A = factor * (A element-wise B) without changing A's block structure.
!> \param matrix_A First factor, overwritten by the product.
!> \param matrix_B Second factor; a missing block represents zero.
!> \param factor Product scale factor.
! **************************************************************************************************
   SUBROUTINE hadamard_product_inplace(matrix_A, matrix_B, factor)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B
      REAL(KIND=dp), INTENT(IN)                          :: factor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hadamard_product_inplace'

      INTEGER                                            :: handle, icol, irow
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_A, block_B
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iterator, matrix_A)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, irow, icol, block_A)
         CALL dbcsr_get_block_p(matrix_B, irow, icol, block_B, found)
         IF (found) THEN
            block_A(:, :) = factor*block_A(:, :)*block_B(:, :)
         ELSE
            block_A(:, :) = 0.0_dp
         END IF
      END DO
      CALL dbcsr_iterator_stop(iterator)

      CALL timestop(handle)
   END SUBROUTINE hadamard_product_inplace

! **************************************************************************************************
!> \brief Computes C=A B A^T or C=A^T B A for DBCSR matrices.
!> \param trans_A_left transposition applied to the left occurrence of A
!> \param trans_A_right transposition applied to the right occurrence of A
!> \param matrix_A left and right matrix A
!> \param matrix_B input matrix B
!> \param matrix_C output matrix C
!> \param eps_filter filtering threshold for both matrix multiplications
!> \param retain_sparsity if true, only existing blocks of C are filled
! **************************************************************************************************
   SUBROUTINE dbcsr_contract_ABA(trans_A_left, trans_A_right, matrix_A, matrix_B, matrix_C, &
                                 eps_filter, retain_sparsity)
      CHARACTER(LEN=1), INTENT(IN)                       :: trans_A_left, trans_A_right
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_A, matrix_B, matrix_C
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_contract_ABA'

      INTEGER                                            :: handle
      LOGICAL                                            :: my_retain_sparsity
      TYPE(dbcsr_type)                                   :: work

      CALL timeset(routineN, handle)

      my_retain_sparsity = .FALSE.
      IF (PRESENT(retain_sparsity)) my_retain_sparsity = retain_sparsity

      CALL dbcsr_create(work, template=matrix_A)

      IF (trans_A_left == "N" .AND. trans_A_right == "T") THEN
         ! C = A B A^T
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_A, matrix_B, &
                             0.0_dp, work, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "T", 1.0_dp, work, matrix_A, &
                             0.0_dp, matrix_C, filter_eps=eps_filter, &
                             retain_sparsity=my_retain_sparsity)
      ELSE IF (trans_A_left == "T" .AND. trans_A_right == "N") THEN
         ! C = A^T B A
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_B, matrix_A, &
                             0.0_dp, work, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_A, work, &
                             0.0_dp, matrix_C, filter_eps=eps_filter, &
                             retain_sparsity=my_retain_sparsity)
      ELSE
         CPABORT("Unsupported transposition pair in dbcsr_contract_ABA")
      END IF
      CALL dbcsr_release(work)

      CALL timestop(handle)

   END SUBROUTINE dbcsr_contract_ABA

! **************************************************************************************************

   SUBROUTINE setup_square_topology(matrix_template, dim_type, square_dist, blk_sizes, mapped_dist)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      CHARACTER(LEN=*), INTENT(IN)                       :: dim_type
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: square_dist
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: blk_sizes, mapped_dist

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_square_topology'

      INTEGER                                            :: handle, i, np, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist
      TYPE(dbcsr_distribution_type)                      :: dist_template

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_template, distribution=dist_template, &
                          row_blk_size=row_blk, col_blk_size=col_blk)
      CALL dbcsr_distribution_get(dist_template, row_dist=row_dist, col_dist=col_dist, &
                                  nprows=nprows, npcols=npcols)

      IF (TRIM(dim_type) == 'ROW') THEN
         ! Creates ROW x ROW (e.g., Grid x Grid from mat_phi_mu_l)
         blk_sizes => row_blk
         np = npcols
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=row_dist, col_dist=mapped_dist)

      ELSE IF (TRIM(dim_type) == 'COL') THEN
         ! Creates COL x COL (e.g., Aux x Aux from mat_Z_lP)
         blk_sizes => col_blk
         np = nprows
         ALLOCATE (mapped_dist(SIZE(blk_sizes)))
         DO i = 1, SIZE(blk_sizes)
            mapped_dist(i) = MOD(i - 1, np)
         END DO
         CALL dbcsr_distribution_new(square_dist, template=dist_template, &
                                     row_dist=mapped_dist, col_dist=col_dist)
      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_square_topology

! **************************************************************************************************
!> \brief DBCSR matrices deallocation
!> \param dist ...
!> \param mapped_dist    ...
!> \param m1 ...
!> \param m2 ...
!> \param m3 ...
!> \param m4 ...
! **************************************************************************************************

   SUBROUTINE release_dbcsr_topology_and_matrices(dist, mapped_dist, m1, m2, m3, m4)

      TYPE(dbcsr_distribution_type), INTENT(INOUT), &
         OPTIONAL                                        :: dist
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL, &
         POINTER                                         :: mapped_dist
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: m1, m2, m3, m4

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_dbcsr_topology_and_matrices'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (PRESENT(dist)) CALL dbcsr_distribution_release(dist)
      IF (PRESENT(mapped_dist)) THEN
         IF (ASSOCIATED(mapped_dist)) THEN
            DEALLOCATE (mapped_dist)
            NULLIFY (mapped_dist)
         END IF
      END IF
      IF (PRESENT(m1)) CALL dbcsr_release(m1)
      IF (PRESENT(m2)) CALL dbcsr_release(m2)
      IF (PRESENT(m3)) CALL dbcsr_release(m3)
      IF (PRESENT(m4)) CALL dbcsr_release(m4)

      CALL timestop(handle)

   END SUBROUTINE release_dbcsr_topology_and_matrices

END MODULE gw_utils_dbcsr
