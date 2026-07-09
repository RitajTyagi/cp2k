!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Static COHSEX (st-COHSEX) quasiparticle energies in the molecular O(N^4) RPA/GW path.
!>
!>        The static screened Coulomb interaction W^C(omega=0) is built from the RPA
!>        polarizability Q(0). The correlation self-energy in the KS-MO basis is
!>
!>           Sigma^c_mn = 1/2 sum_p W^C_mppn(0) sgn(eps_p - eps_F)
!>
!>        and the (off-diagonal) quasiparticle Hamiltonian
!>
!>           H_mn = eps_KS_n delta_mn + (Sigma_x - v_xc)_mn + Sigma^c_mn
!>
!>        is diagonalized separately in the occupied-occupied and virtual-virtual blocks
!>        (one-shot @KS). W^C(0), the st-COHSEX eigenvalues and eigenvectors are stored on
!>        mp2_env%ri_g0w0 for a subsequent BSE calculation. No analytic continuation is needed.
!> \author Ritaj Tyagi
!> \date   2026
! **************************************************************************************************
MODULE rpa_gw_stcohsex
   USE cp_dbcsr_api,                    ONLY: dbcsr_release
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_trace,&
                                              cp_fm_uplo_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_stcohsex'

   PUBLIC :: gw_stcohsex_energies

CONTAINS

! **************************************************************************************************
!> \brief Compute static COHSEX quasiparticle energies and eigenvectors and store them (together
!>        with W^C(0)) on mp2_env%ri_g0w0 for later BSE use.
!> \param fm_mat_Q on entry the spin-summed 1 + Q(omega=0) in the RI basis; on exit it holds the
!>        correlation kernel W^C(0) = (1 + Q(0))^-1 - 1 (inverted in place).
!> \param fm_mat_S_gw per-spin RI 3-center integrals B, size dimen_RI x dimen_nm_gw.
!> \param fm_mat_S_gw_work per-spin work matrix (same shape as fm_mat_S_gw), used to hold W^C.B.
!> \param dimen_RI number of (reduced) RI basis functions.
!> \param dimen_nm_gw compound dimension = gw_corr_lev_tot * nmo.
!> \param nmo number of molecular orbitals.
!> \param homo per-spin HOMO index.
!> \param gw_corr_lev_occ per-spin number of corrected occupied levels.
!> \param gw_corr_lev_virt per-spin number of corrected virtual levels.
!> \param Eigenval KS eigenvalues, size nmo x nspins.
!> \param mo_coeff per-spin MO coefficients (AO x nmo).
!> \param para_env global parallel environment.
!> \param qs_env ...
!> \param mp2_env ...
!> \param unit_nr output unit (only > 0 on the ionode).
! **************************************************************************************************
   SUBROUTINE gw_stcohsex_energies(fm_mat_Q, fm_mat_S_gw, fm_mat_S_gw_work, &
                                   dimen_RI, dimen_nm_gw, nmo, homo, &
                                   gw_corr_lev_occ, gw_corr_lev_virt, &
                                   Eigenval, mo_coeff, para_env, qs_env, mp2_env, unit_nr)

      ! fm_mat_Q follows the cp_fm INTENT(IN) idiom: the matrix data is overwritten
      ! (inverted in place -> W^C(0)) but the descriptor itself is not reassigned
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S_gw, fm_mat_S_gw_work
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_nm_gw, nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, gw_corr_lev_occ, gw_corr_lev_virt
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type)                                     :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_stcohsex_energies'

      INTEGER                                            :: handle, i_global, i_win, iiB, ispin, &
                                                            j_global, jjB, m, n, n_ao, ncol_local, &
                                                            no, nrow_local, nspins, nv, nwin, p, &
                                                            win_start
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eig_occ, eig_virt, scale_vec
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_occ, h_virt, sigma_c, sxmvxc
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_block, fm_struct_mo_mo, &
                                                            fm_struct_win
      TYPE(cp_fm_type)                                   :: fm_ao, fm_mo_mo, fm_tmp_ao_mo, fm_wbblk, &
                                                            fm_work_ri
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_bblk

      CALL timeset(routineN, handle)

      nspins = SIZE(fm_mat_S_gw)

      ! ================================================================================
      ! 1. Build the static correlation kernel W^C(0) = (1 + Q(0))^-1 - 1 in the RI basis
      ! ================================================================================
      CALL cp_fm_create(fm_work_ri, fm_mat_Q%matrix_struct)
      CALL cp_fm_cholesky_invert(fm_mat_Q)
      CALL cp_fm_uplo_to_full(fm_mat_Q, fm_work_ri)

      CALL cp_fm_get_info(matrix=fm_mat_Q, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! subtract 1 on the diagonal to remove the bare exchange part -> correlation kernel W^C
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (i_global == j_global .AND. i_global <= dimen_RI) THEN
               fm_mat_Q%local_data(iiB, jjB) = fm_mat_Q%local_data(iiB, jjB) - 1.0_dp
            END IF
         END DO
      END DO
      CALL cp_fm_release(fm_work_ri)

      ! store W^C(0) for later BSE use
      CALL cp_fm_create(mp2_env%ri_g0w0%fm_W_c_zero, fm_mat_Q%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_Q, mp2_env%ri_g0w0%fm_W_c_zero)

      ! storage for QP eigenvalues (initialized with KS values) and eigenvectors
      ALLOCATE (mp2_env%ri_g0w0%stcohsex_eigenval(nmo, nspins))
      mp2_env%ri_g0w0%stcohsex_eigenval(:, :) = Eigenval(:, :)
      ALLOCATE (mp2_env%ri_g0w0%fm_stcohsex_evec(nspins))

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(/,T3,A)") &
            "STATIC COHSEX| Static screened Coulomb interaction W^C(0) built"
      END IF

      DO ispin = 1, nspins

         no = gw_corr_lev_occ(ispin)
         nv = gw_corr_lev_virt(ispin)
         nwin = no + nv
         win_start = homo(ispin) - no + 1

         ! ------------------------------------------------------------------------------
         ! 2. WB = W^C(0) . B, then column-scale by 1/2 sgn(eps_p - eps_F)
         !    (eps_p - eps_F < 0 for occupied p <= homo, > 0 for virtual p > homo)
         ! ------------------------------------------------------------------------------
         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_nm_gw, k=dimen_RI, &
                            alpha=1.0_dp, matrix_a=fm_mat_Q, matrix_b=fm_mat_S_gw(ispin), &
                            beta=0.0_dp, matrix_c=fm_mat_S_gw_work(ispin))

         ALLOCATE (scale_vec(dimen_nm_gw))
         DO jjB = 1, dimen_nm_gw
            i_win = MAX(1, jjB - 1)/nmo + 1
            p = jjB - (i_win - 1)*nmo
            IF (p <= homo(ispin)) THEN
               scale_vec(jjB) = -0.5_dp
            ELSE
               scale_vec(jjB) = 0.5_dp
            END IF
         END DO
         CALL cp_fm_column_scale(fm_mat_S_gw_work(ispin), scale_vec)
         DEALLOCATE (scale_vec)

         ! ------------------------------------------------------------------------------
         ! 3. Sigma^c_mn = sum_p [1/2 sgn(eps_p-eps_F)] sum_P WB_(m,p),P B_(n,p),P
         !    formed as the Frobenius product of the per-window (RI x nmo) column blocks
         ! ------------------------------------------------------------------------------
         NULLIFY (fm_struct_block)
         CALL cp_fm_struct_create(fm_struct_block, template_fmstruct=fm_mat_S_gw(ispin)%matrix_struct, &
                                  nrow_global=dimen_RI, ncol_global=nmo)

         ALLOCATE (fm_bblk(nwin))
         DO m = 1, nwin
            CALL cp_fm_create(fm_bblk(m), fm_struct_block)
            CALL cp_fm_to_fm_submat(fm_mat_S_gw(ispin), fm_bblk(m), nrow=dimen_RI, ncol=nmo, &
                                    s_firstrow=1, s_firstcol=(m - 1)*nmo + 1, t_firstrow=1, t_firstcol=1)
         END DO
         CALL cp_fm_create(fm_wbblk, fm_struct_block)

         ALLOCATE (sigma_c(nwin, nwin))
         sigma_c(:, :) = 0.0_dp
         DO m = 1, nwin
            CALL cp_fm_to_fm_submat(fm_mat_S_gw_work(ispin), fm_wbblk, nrow=dimen_RI, ncol=nmo, &
                                    s_firstrow=1, s_firstcol=(m - 1)*nmo + 1, t_firstrow=1, t_firstcol=1)
            DO n = 1, nwin
               CALL cp_fm_trace(fm_wbblk, fm_bblk(n), sigma_c(m, n))
            END DO
         END DO

         CALL cp_fm_release(fm_wbblk)
         DO m = 1, nwin
            CALL cp_fm_release(fm_bblk(m))
         END DO
         DEALLOCATE (fm_bblk)
         CALL cp_fm_struct_release(fm_struct_block)

         ! ------------------------------------------------------------------------------
         ! 4. (Sigma_x - v_xc)_mn in the MO window: C^T (Sigma_x - v_xc)_AO C
         ! ------------------------------------------------------------------------------
         CALL cp_fm_get_info(matrix=mo_coeff(ispin), nrow_global=n_ao)

         NULLIFY (fm_struct_mo_mo)
         CALL cp_fm_struct_create(fm_struct_mo_mo, template_fmstruct=mo_coeff(ispin)%matrix_struct, &
                                  nrow_global=n_ao, ncol_global=n_ao)
         CALL cp_fm_create(fm_ao, fm_struct_mo_mo)
         CALL cp_fm_struct_release(fm_struct_mo_mo)
         CALL copy_dbcsr_to_fm(mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc(ispin)%matrix, fm_ao)

         ! tmp(n_ao x nmo) = (Sigma_x - v_xc)_AO . C
         CALL cp_fm_create(fm_tmp_ao_mo, mo_coeff(ispin)%matrix_struct)
         CALL parallel_gemm(transa="N", transb="N", m=n_ao, n=nmo, k=n_ao, alpha=1.0_dp, &
                            matrix_a=fm_ao, matrix_b=mo_coeff(ispin), beta=0.0_dp, matrix_c=fm_tmp_ao_mo)
         CALL cp_fm_release(fm_ao)

         ! mo_mo(nmo x nmo) = C^T . tmp
         NULLIFY (fm_struct_mo_mo)
         CALL cp_fm_struct_create(fm_struct_mo_mo, template_fmstruct=mo_coeff(ispin)%matrix_struct, &
                                  nrow_global=nmo, ncol_global=nmo)
         CALL cp_fm_create(fm_mo_mo, fm_struct_mo_mo)
         CALL cp_fm_struct_release(fm_struct_mo_mo)
         CALL parallel_gemm(transa="T", transb="N", m=nmo, n=nmo, k=n_ao, alpha=1.0_dp, &
                            matrix_a=mo_coeff(ispin), matrix_b=fm_tmp_ao_mo, beta=0.0_dp, matrix_c=fm_mo_mo)
         CALL cp_fm_release(fm_tmp_ao_mo)

         ALLOCATE (sxmvxc(nwin, nwin))
         sxmvxc(:, :) = 0.0_dp
         CALL cp_fm_get_submatrix(fm_mo_mo, sxmvxc, start_row=win_start, start_col=win_start, &
                                  n_rows=nwin, n_cols=nwin)
         CALL cp_fm_release(fm_mo_mo)

         ! ------------------------------------------------------------------------------
         ! 5. Build and diagonalize the occ-occ and virt-virt quasiparticle blocks
         ! ------------------------------------------------------------------------------
         ! occupied block
         ALLOCATE (h_occ(no, no), eig_occ(no))
         DO n = 1, no
            DO m = 1, no
               h_occ(m, n) = sxmvxc(m, n) + sigma_c(m, n)
            END DO
            h_occ(n, n) = h_occ(n, n) + Eigenval(win_start + n - 1, ispin)
         END DO
         h_occ(:, :) = 0.5_dp*(h_occ(:, :) + TRANSPOSE(h_occ(:, :)))
         CALL diamat_all(h_occ, eig_occ)

         ! virtual block
         ALLOCATE (h_virt(nv, nv), eig_virt(nv))
         DO n = 1, nv
            DO m = 1, nv
               h_virt(m, n) = sxmvxc(no + m, no + n) + sigma_c(no + m, no + n)
            END DO
            h_virt(n, n) = h_virt(n, n) + Eigenval(win_start + no + n - 1, ispin)
         END DO
         h_virt(:, :) = 0.5_dp*(h_virt(:, :) + TRANSPOSE(h_virt(:, :)))
         CALL diamat_all(h_virt, eig_virt)

         ! ------------------------------------------------------------------------------
         ! 6. Store eigenvalues (in the window slots) and the block-diagonal eigenvectors
         ! ------------------------------------------------------------------------------
         DO n = 1, no
            mp2_env%ri_g0w0%stcohsex_eigenval(win_start + n - 1, ispin) = eig_occ(n)
         END DO
         DO n = 1, nv
            mp2_env%ri_g0w0%stcohsex_eigenval(win_start + no + n - 1, ispin) = eig_virt(n)
         END DO

         NULLIFY (fm_struct_win)
         CALL cp_fm_struct_create(fm_struct_win, template_fmstruct=fm_mat_S_gw(ispin)%matrix_struct, &
                                  nrow_global=nwin, ncol_global=nwin)
         CALL cp_fm_create(mp2_env%ri_g0w0%fm_stcohsex_evec(ispin), fm_struct_win)
         CALL cp_fm_struct_release(fm_struct_win)
         CALL cp_fm_set_all(mp2_env%ri_g0w0%fm_stcohsex_evec(ispin), 0.0_dp)
         CALL cp_fm_set_submatrix(mp2_env%ri_g0w0%fm_stcohsex_evec(ispin), h_occ, &
                                  start_row=1, start_col=1)
         CALL cp_fm_set_submatrix(mp2_env%ri_g0w0%fm_stcohsex_evec(ispin), h_virt, &
                                  start_row=no + 1, start_col=no + 1)

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,I2,A,F12.5,A,F12.5,A)") &
               "STATIC COHSEX| spin ", ispin, ":  HOMO = ", &
               mp2_env%ri_g0w0%stcohsex_eigenval(homo(ispin), ispin)*evolt, " eV,  LUMO = ", &
               mp2_env%ri_g0w0%stcohsex_eigenval(homo(ispin) + 1, ispin)*evolt, " eV"
            WRITE (UNIT=unit_nr, FMT="(T3,A,I2,A,F12.5,A)") &
               "STATIC COHSEX| spin ", ispin, ":  gap  = ", &
               (mp2_env%ri_g0w0%stcohsex_eigenval(homo(ispin) + 1, ispin) &
                - mp2_env%ri_g0w0%stcohsex_eigenval(homo(ispin), ispin))*evolt, " eV"
         END IF

         DEALLOCATE (sigma_c, sxmvxc, h_occ, h_virt, eig_occ, eig_virt)

      END DO ! ispin

      ! the persistent AO Sigma_x - v_xc matrix has been consumed (field is allocatable, release manually)
      IF (ALLOCATED(mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc)) THEN
         DO ispin = 1, SIZE(mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc)
            CALL dbcsr_release(mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc(ispin)%matrix)
            DEALLOCATE (mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc(ispin)%matrix)
         END DO
         DEALLOCATE (mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc)
      END IF

      MARK_USED(qs_env)
      MARK_USED(para_env)

      CALL timestop(handle)

   END SUBROUTINE gw_stcohsex_energies

END MODULE rpa_gw_stcohsex
