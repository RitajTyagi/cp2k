!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Full-matrix operations not provided by the CP2K FM packages.
!> \par History
!>      09.2026 created Jan Wilhelm
! **************************************************************************************************
MODULE gw_utils_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_diag,                      ONLY: cp_fm_power,&
                                              cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_equivalent,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE mathlib,                         ONLY: diag_complex,&
                                              gemm_square
   USE parallel_gemm_api,               ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_utils_fm'

   PUBLIC :: cfm_contract_ABA, fm_contract_ABA, fm_invert, fm_pseudo_inverse_relative, &
             fm_sqrt, local_add_on_diag, local_complex_power

   INTERFACE fm_contract_ABA
      MODULE PROCEDURE fm_contract_ABA, fm_contract_ABA_set
   END INTERFACE fm_contract_ABA

CONTAINS

! **************************************************************************************************
!> \brief Computes A^H B A for complex full matrices.
!> \param matrix_A left and right matrix A; H denotes its Hermitian transpose
!> \param matrix_B input matrix B; overwritten by A^H B A if matrix_C is absent
!> \param matrix_C optional output matrix C=A^H B A
!> \param workspace optional reusable scratch matrix with the same layout as A;
!>        its contents are overwritten and it must not alias the input or output matrices
! **************************************************************************************************
   SUBROUTINE cfm_contract_ABA(matrix_A, matrix_B, matrix_C, workspace)
      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix_A
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrix_B
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL         :: matrix_C, workspace

      INTEGER                                            :: ncol, nrow
      LOGICAL                                            :: compatible_workspace
      TYPE(cp_cfm_type)                                  :: work

      CPASSERT(ASSOCIATED(matrix_A%matrix_struct))
      CPASSERT(ASSOCIATED(matrix_B%matrix_struct))
      CALL cp_cfm_get_info(matrix_A, nrow_global=nrow, ncol_global=ncol)
      CPASSERT(nrow == ncol)
      IF (PRESENT(workspace)) THEN
         compatible_workspace = cp_fm_struct_equivalent(matrix_A%matrix_struct, workspace%matrix_struct)
         compatible_workspace = compatible_workspace .AND. &
                                ALL(matrix_A%matrix_struct%descriptor == workspace%matrix_struct%descriptor)
         CPASSERT(compatible_workspace)
         CPASSERT(ASSOCIATED(workspace%local_data))
         CPASSERT(.NOT. ASSOCIATED(workspace%local_data, matrix_A%local_data))
         CPASSERT(.NOT. ASSOCIATED(workspace%local_data, matrix_B%local_data))
         IF (PRESENT(matrix_C)) THEN
            CPASSERT(.NOT. ASSOCIATED(workspace%local_data, matrix_C%local_data))
         END IF
         work = workspace
      ELSE
         CALL cp_cfm_create(work, matrix_A%matrix_struct)
      END IF

      CALL parallel_gemm("N", "N", nrow, nrow, nrow, z_one, matrix_B, matrix_A, &
                         z_zero, work)
      IF (PRESENT(matrix_C)) THEN
         CPASSERT(ASSOCIATED(matrix_C%matrix_struct))
         CALL parallel_gemm("C", "N", nrow, nrow, nrow, z_one, matrix_A, work, &
                            z_zero, matrix_C)
      ELSE
         CALL parallel_gemm("C", "N", nrow, nrow, nrow, z_one, matrix_A, work, &
                            z_zero, matrix_B)
      END IF

      IF (.NOT. PRESENT(workspace)) CALL cp_cfm_release(work)

   END SUBROUTINE cfm_contract_ABA

! **************************************************************************************************
!> \brief Add a scalar to the diagonal of a local complex matrix.
!> \param matrix input/output matrix
!> \param alpha scalar added to each diagonal element
! **************************************************************************************************
   SUBROUTINE local_add_on_diag(matrix, alpha)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: matrix
      COMPLEX(KIND=dp), INTENT(IN)                       :: alpha

      INTEGER                                            :: i, n

      n = SIZE(matrix, 1)
      CPASSERT(n == SIZE(matrix, 2))

      DO i = 1, n
         matrix(i, i) = matrix(i, i) + alpha
      END DO

   END SUBROUTINE local_add_on_diag

! **************************************************************************************************
!> \brief Compute a spectral power of a local complex Hermitian matrix.
!>        Eigenvalues not larger than eps (or relative_eps times the largest one) are discarded.
!> \param matrix input/output matrix
!> \param exponent matrix exponent
!> \param eps eigenvalue threshold
!> \param cond_nr optional condition number
!> \param min_ev optional minimum absolute eigenvalue
!> \param max_ev optional maximum absolute eigenvalue
!> \param relative_eps optional threshold relative to the largest eigenvalue
! **************************************************************************************************
   SUBROUTINE local_complex_power(matrix, exponent, eps, cond_nr, min_ev, max_ev, relative_eps)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: exponent, eps
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: cond_nr, min_ev, max_ev
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: relative_eps

      CHARACTER(len=*), PARAMETER :: routineN = 'local_complex_power'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: eigenvectors
      INTEGER                                            :: handle, i, n
      REAL(KIND=dp)                                      :: eigenvalue_cutoff, sqrt_power_eval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues

      CALL timeset(routineN, handle)

      n = SIZE(matrix, 1)
      CPASSERT(n == SIZE(matrix, 2))
      matrix(:, :) = 0.5_dp*(matrix(:, :) + CONJG(TRANSPOSE(matrix(:, :))))

      ALLOCATE (eigenvalues(n), eigenvectors(n, n))
      CALL diag_complex(matrix, eigenvectors, eigenvalues)

      IF (PRESENT(cond_nr)) cond_nr = MAXVAL(ABS(eigenvalues))/MINVAL(ABS(eigenvalues))
      IF (PRESENT(min_ev)) min_ev = MINVAL(ABS(eigenvalues))
      IF (PRESENT(max_ev)) max_ev = MAXVAL(ABS(eigenvalues))
      eigenvalue_cutoff = eps
      IF (PRESENT(relative_eps)) eigenvalue_cutoff = MAX(eps, relative_eps*MAXVAL(eigenvalues))

      ! Reconstruct matrix**exponent as X*X^H, where
      ! X = U*diag(eigenvalues**(exponent/2)).
      DO i = 1, n
         IF (eigenvalue_cutoff < eigenvalues(i)) THEN
            sqrt_power_eval = eigenvalues(i)**(0.5_dp*exponent)
         ELSE
            sqrt_power_eval = 0.0_dp
         END IF
         eigenvectors(:, i) = eigenvectors(:, i)*sqrt_power_eval
      END DO

      ! X*X^H = U*diag(eigenvalues**exponent)*U^H.
      CALL gemm_square(eigenvectors, 'N', eigenvectors, 'C', matrix)

      DEALLOCATE (eigenvalues, eigenvectors)

      CALL timestop(handle)

   END SUBROUTINE local_complex_power

! **************************************************************************************************
!> \brief Inverts a symmetric matrix. First, Cholesky decomposition is tried.
!>        If it fails, the matrix is diagonalized and inverted by taking
!>        inverse eigenvalues.
!> \param matrix_A input A; output inverse or filtered pseudoinverse
!> \param eigenvalue_threshold optional relative eigenvalue threshold
!> \param unit_nr optional output unit for a failed Cholesky decomposition
! **************************************************************************************************
   SUBROUTINE fm_invert(matrix_A, eigenvalue_threshold, unit_nr)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix_A
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eigenvalue_threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      INTEGER                                            :: info, n_dependent
      LOGICAL                                            :: cholesky_successful
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_fm_type)                                   :: work

      threshold = 0.0_dp
      IF (PRESENT(eigenvalue_threshold)) threshold = eigenvalue_threshold
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      ! Preserve A because a failed Cholesky decomposition overwrites it.
      CALL cp_fm_to_fm(matrix_A, work)
      CALL cp_fm_cholesky_decompose(matrix_A, info_out=info)
      cholesky_successful = info == 0
      n_dependent = 0
      IF (cholesky_successful) THEN
         CALL cp_fm_cholesky_invert(matrix_A)
         CALL cp_fm_uplo_to_full(matrix_A, work)
      ELSE
         CALL cp_fm_to_fm(work, matrix_A)
         CALL cp_fm_power(matrix_A, work, -1.0_dp, threshold, n_dependent)
         IF (PRESENT(unit_nr)) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') &
                  'Cholesky decomposition failed; matrix inverted by diagonalization.'
               WRITE (unit_nr, '(T2,A,T72,I9)') 'Discarded eigenmodes:', n_dependent
            END IF
         END IF
      END IF
      CALL cp_fm_release(work)

   END SUBROUTINE fm_invert

! **************************************************************************************************
!> \brief For input A, computes B such that B^T B=A.
!>        First, Cholesky decomposition is tried. If it fails, B is computed
!>        by diagonalizing A.
!> \param matrix_A symmetric input matrix A, retained unchanged
!> \param matrix_B output matrix B
!> \param eigenvalue_threshold optional relative eigenvalue threshold
!> \param unit_nr optional output unit for a failed Cholesky decomposition
! **************************************************************************************************
   SUBROUTINE fm_sqrt(matrix_A, matrix_B, eigenvalue_threshold, unit_nr)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_A
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix_B
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eigenvalue_threshold
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      INTEGER                                            :: i_row, info, j_col, n_dependent, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: cholesky_successful
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_fm_type)                                   :: work

      threshold = 0.0_dp
      IF (PRESENT(eigenvalue_threshold)) threshold = eigenvalue_threshold
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      CALL cp_fm_to_fm(matrix_A, matrix_B)
      CALL cp_fm_cholesky_decompose(matrix_B, info_out=info)
      cholesky_successful = info == 0
      IF (cholesky_successful) THEN
         n_dependent = 0
         CALL cp_fm_get_info(matrix=matrix_B, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices)
         DO j_col = 1, ncol_local
            DO i_row = 1, nrow_local
               IF (row_indices(i_row) > col_indices(j_col)) THEN
                  matrix_B%local_data(i_row, j_col) = 0.0_dp
               END IF
            END DO
         END DO
      ELSE
         CALL cp_fm_to_fm(matrix_A, matrix_B)
         CALL cp_fm_power(matrix_B, work, 0.5_dp, threshold, n_dependent)
         IF (PRESENT(unit_nr)) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') &
                  'Cholesky decomposition failed, using diagonalization.'
               WRITE (unit_nr, '(T2,A,T72,I9)') 'Discarded eigenmodes:', n_dependent
            END IF
         END IF
      END IF
      CALL cp_fm_release(work)

   END SUBROUTINE fm_sqrt

! **************************************************************************************************
!> \brief Computes A^T B A.
!> \param matrix_A left and right matrix A
!> \param matrix_B input matrix B; overwritten by A^T B A if matrix_C is absent
!> \param matrix_C optional output matrix C=A^T B A
! **************************************************************************************************
   SUBROUTINE fm_contract_ABA(matrix_A, matrix_B, matrix_C)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_A
      TYPE(cp_fm_type), INTENT(INOUT)                    :: matrix_B
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL          :: matrix_C

      INTEGER                                            :: ncol, nrow
      TYPE(cp_fm_type)                                   :: work

      CPASSERT(ASSOCIATED(matrix_A%matrix_struct))
      CPASSERT(ASSOCIATED(matrix_B%matrix_struct))
      CALL cp_fm_get_info(matrix_A, nrow_global=nrow, ncol_global=ncol)
      CPASSERT(nrow == ncol)
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      CALL parallel_gemm("N", "N", nrow, nrow, nrow, 1.0_dp, matrix_B, matrix_A, &
                         0.0_dp, work)
      IF (PRESENT(matrix_C)) THEN
         CPASSERT(ASSOCIATED(matrix_C%matrix_struct))
         CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                            0.0_dp, matrix_C)
      ELSE
         CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                            0.0_dp, matrix_B)
      END IF

      CALL cp_fm_release(work)

   END SUBROUTINE fm_contract_ABA

! **************************************************************************************************
!> \brief Computes A^T B_i A for a set of matrices.
!> \param matrix_A left and right matrix A
!> \param matrix_B input matrices B_i; overwritten if matrix_C is absent
!> \param matrix_C optional output matrices C_i=A^T B_i A
! **************************************************************************************************
   SUBROUTINE fm_contract_ABA_set(matrix_A, matrix_B, matrix_C)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix_A
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: matrix_B
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: matrix_C

      INTEGER                                            :: i, ncol, nrow
      TYPE(cp_fm_type)                                   :: work

      CPASSERT(ASSOCIATED(matrix_A%matrix_struct))
      CALL cp_fm_get_info(matrix_A, nrow_global=nrow, ncol_global=ncol)
      CPASSERT(nrow == ncol)
      IF (PRESENT(matrix_C)) THEN
         CPASSERT(SIZE(matrix_C) == SIZE(matrix_B))
      END IF
      CALL cp_fm_create(work, matrix_A%matrix_struct)

      DO i = 1, SIZE(matrix_B)
         CPASSERT(ASSOCIATED(matrix_B(i)%matrix_struct))
         CALL parallel_gemm("N", "N", nrow, nrow, nrow, 1.0_dp, matrix_B(i), matrix_A, &
                            0.0_dp, work)
         IF (PRESENT(matrix_C)) THEN
            CPASSERT(ASSOCIATED(matrix_C(i)%matrix_struct))
            CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                               0.0_dp, matrix_C(i))
         ELSE
            CALL parallel_gemm("T", "N", nrow, nrow, nrow, 1.0_dp, matrix_A, work, &
                               0.0_dp, matrix_B(i))
         END IF
      END DO

      CALL cp_fm_release(work)

   END SUBROUTINE fm_contract_ABA_set

! **************************************************************************************************
!> \brief Relative-cut pseudo-inverse of a real symmetric matrix M = U diag(m_i) Uᵀ, in place:
!>          M⁺ = U diag(f_i) Uᵀ,   f_i = 1/m_i  for m_i ≥ x m_max,   f_i = 0  otherwise.
!>        Optionally the retained modes L = {i : x m_max ≤ m_i < x_c m_max} are returned as U_L
!>        and m_L, e.g. for M⁺_c = M⁺ − U_L diag(1/m_L) U_Lᵀ.
!> \param matrix symmetric matrix M on input, its pseudo-inverse M⁺ on output
!> \param rel_cut relative cut x (> 0)
!> \param n_dropped number of eigenvalues m_i < x m_max
!> \param m_kept smallest retained eigenvalue
!> \param m_max largest eigenvalue
!> \param rel_cut_low optional second cut x_c > x defining the set L above
!> \param n_low number of modes in L (0 if rel_cut_low is absent or L is empty)
!> \param fm_U_low the eigenvectors U_L, n x n_low, distributed like matrix; created only if
!>        n_low > 0
!> \param m_low their eigenvalues m_L (allocated with n_low entries)
! **************************************************************************************************
   SUBROUTINE fm_pseudo_inverse_relative(matrix, rel_cut, n_dropped, m_kept, m_max, &
                                         rel_cut_low, n_low, fm_U_low, m_low)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: rel_cut
      INTEGER, INTENT(OUT)                               :: n_dropped
      REAL(KIND=dp), INTENT(OUT)                         :: m_kept, m_max
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: rel_cut_low
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_low
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL          :: fm_U_low
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: m_low

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_pseudo_inverse_relative'

      INTEGER                                            :: handle, i, i_first, n, n_l
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval, factor
      TYPE(cp_fm_struct_type), POINTER                   :: struct_low
      TYPE(cp_fm_type)                                   :: fm_U, fm_U_scaled

      CALL timeset(routineN, handle)
      CPASSERT(rel_cut > 0.0_dp)

      CALL cp_fm_get_info(matrix, nrow_global=n)
      ALLOCATE (eigval(n), factor(n))
      CALL cp_fm_create(fm_U, matrix%matrix_struct)
      CALL cp_fm_create(fm_U_scaled, matrix%matrix_struct)

      ! M = U diag(m_i) Uᵀ, eigenvalues ascending (cp_fm_syevd destroys M). ScaLAPACK, not ELPA:
      ! rounding errors in the small retained m_i are amplified by 1/m_i.
      CALL cp_fm_syevd(matrix, fm_U, eigval)
      m_max = eigval(n)
      ! f_i = 1/m_i on the retained modes m_i ≥ x m_max, 0 on the dropped ones
      n_dropped = 0
      m_kept = m_max
      DO i = 1, n
         IF (eigval(i) < rel_cut*m_max) THEN
            factor(i) = 0.0_dp
            n_dropped = n_dropped + 1
         ELSE
            factor(i) = 1.0_dp/eigval(i)
            m_kept = MIN(m_kept, eigval(i))
         END IF
      END DO
      ! M⁺ = (U diag(f)) Uᵀ
      CALL cp_fm_to_fm(fm_U, fm_U_scaled)
      CALL cp_fm_column_scale(fm_U_scaled, factor)
      CALL parallel_gemm("N", "T", n, n, n, 1.0_dp, fm_U_scaled, fm_U, 0.0_dp, matrix)

      ! Optional: the retained modes below the second cut, L = {i : x m_max ≤ m_i < x_c m_max},
      ! i.e. the columns i_first .. i_first+n_l-1 of U (eigenvalues ascending)
      IF (PRESENT(rel_cut_low)) THEN
         CPASSERT(PRESENT(n_low) .AND. PRESENT(fm_U_low) .AND. PRESENT(m_low))
         i_first = n_dropped + 1
         n_l = 0
         DO i = i_first, n
            IF (eigval(i) < rel_cut_low*m_max) n_l = n_l + 1
         END DO
         n_low = n_l
         IF (ALLOCATED(m_low)) DEALLOCATE (m_low)
         ALLOCATE (m_low(n_l))
         IF (n_l > 0) THEN
            m_low(:) = eigval(i_first:i_first + n_l - 1)
            NULLIFY (struct_low)
            CALL cp_fm_struct_create(struct_low, template_fmstruct=matrix%matrix_struct, &
                                     ncol_global=n_l)
            CALL cp_fm_create(fm_U_low, struct_low)
            CALL cp_fm_struct_release(struct_low)
            CALL cp_fm_to_fm_submat(fm_U, fm_U_low, n, n_l, 1, i_first, 1, 1)
         END IF
      END IF

      CALL cp_fm_release(fm_U)
      CALL cp_fm_release(fm_U_scaled)
      DEALLOCATE (eigval, factor)

      CALL timestop(handle)

   END SUBROUTINE fm_pseudo_inverse_relative

END MODULE gw_utils_fm
