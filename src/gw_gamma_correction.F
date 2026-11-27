!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Γ-point (q → 0) head-and-wing correction of W for the Γ-only large-cell RI-RS GW code
!>        (2D slabs and 3D bulk cells), following J. Wilhelm, *Periodic real-space RI for full-k GW
!>        band structures*. The charge channel is split off the auxiliary basis and the dielectric
!>        head and wings are integrated over the Wigner-Seitz cell of the supercell Brillouin zone:
!>
!>          c_P   = ∫ dr ϕ_P0(r)                                     auxiliary charge
!>          Π_PQ  = δ_PQ − c_P c_Q / Σ_S |c_S|²                      neutral projector
!>          V_b   = Π V(0) Π                                          neutral body
!>          M⁺    = Σ_{m_i ≥ x_g m_max} u_i u_iᵀ/m_i                 RI-metric pseudo-inverse
!>          M⁺_c  = M⁺ − U_L diag(1/m_L) U_Lᵀ,   U_L: eigenvectors with x_g m_max ≤ m_i < x_c m_max
!>
!>        V(0) is the Ewald sum without its G = 0 term, (u_i, m_i) the eigenpairs of the RI metric,
!>        x_g = REGULARIZATION_RI and x_c = REGULARIZATION_RI_CORRECTION; M⁺_c is used by the
!>        correction. 2D: charge mode ϕ_0 = 1/√A, v(q) = 2π/q; ΔW is added to the folded W, M⁺ W M⁺.
!>        3D: ϕ_0 = 1/√Ω, v(q) = 4π/q², Π = 1; the sphere-tail term is added to W, head and wings to
!>        Σ^c in the orbital basis. Atomic units; P, Q: RI_AUX functions; i (a): occupied (empty)
!>        states; α, β: x, y, z; ℓ, ℓ': RI-RS grid points r_ℓ.
!> \author Ritaj Tyagi
! **************************************************************************************************
MODULE gw_gamma_correction
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE constants_operator,              ONLY: operator_verfc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks,&
                                              dbcsr_scale_by_vector
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE gw_ri_rs_small_cell_full_kp_lattice,ONLY: gamma_auxiliary_charges
   USE input_constants,                 ONLY: ri_rs_large_cell_Gamma
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_coulomb_2c,               ONLY: add_ewald_reciprocal_2d_fm,&
                                              add_ewald_reciprocal_fm,&
                                              build_2c_coulomb_matrix_kp
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE mathlib,                         ONLY: det_3x3,&
                                              gauss_legendre,&
                                              inv_3x3,&
                                              vector_product
   USE message_passing,                 ONLY: mp_comm_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_moments,                      ONLY: build_berry_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE ri_rs_types,                     ONLY: gamma_corr_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_gamma_correction'

   ! Relative eigenvalue cut for V_b^(±1/2), V_b = Π V(0) Π: removes the charge null direction
   ! and the linear dependence of the RI basis.
   REAL(KIND=dp), PARAMETER, PRIVATE :: eps_null_V_b = 1.0E-10_dp

   ! Out-of-plane form-factor classes (gamma_corr_classify_atoms): per kind, the atom heights d are
   ! split into planes at gaps > class_plane_gap and binned within a plane at d̄ + k·class_bin_width
   ! (d̄ the plane mean), so |d − d_class| ≤ class_bin_width/2.
   REAL(KIND=dp), PARAMETER, PRIVATE :: class_bin_width = 0.1_dp, &
                                        class_plane_gap = 3.0_dp*class_bin_width

   ! Points per direction of the slab zone-remainder tables (gamma_corr_tabulate); also used by
   ! the memory estimate of the Γ-only driver.
   INTEGER, PARAMETER, PUBLIC :: gamma_corr_n_tab_rem = 32

   PUBLIC :: gamma_corr_init, &
             gamma_corr_build_coulomb_maps, &
             gamma_corr_response, &
             gamma_corr_tabulate, &
             gamma_corr_add_folded_to_W, &
             gamma_corr_sigma_shifts, &
             gamma_corr_tail_operator, &
             gamma_corr_classify_atoms

CONTAINS

! **************************************************************************************************
!> \brief Madelung constant ξ of a 3D lattice in the Ewald convention of V(q = 0) (no G = 0 term),
!>
!>          v_E(x) = Σ_L erfc(ω|x+L|)/|x+L| + (4π/Ω) Σ_{G≠0} e^{−G²/4ω²} cos(G·x)/G² − π/(Ωω²)
!>          ξ      = lim_{x→0} [1/|x| − v_E(x)]
!>                 = 2ω/√π − Σ_{L≠0} erfc(ω|L|)/|L| − (4π/Ω) Σ_{G≠0} e^{−G²/4ω²}/G² + π/(Ωω²),
!>
!>        independent of ω (here ω = √π/Ω^{1/3}).
!> \param hmat cell matrix, columns are the lattice vectors (bohr)
!> \return ξ (1/bohr = Hartree for unit charges)
! **************************************************************************************************
   FUNCTION gamma_corr_madelung_3d(hmat) RESULT(xi)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      REAL(KIND=dp)                                      :: xi

      REAL(KIND=dp), PARAMETER                           :: arg_max = 8.0_dp

      INTEGER                                            :: i1, i2, i3
      INTEGER, DIMENSION(3)                              :: n_g, n_l
      REAL(KIND=dp)                                      :: g2, g_cut, l_cut, omega, r, sum_g, &
                                                            sum_l, volume
      REAL(KIND=dp), DIMENSION(3)                        :: g_vec, l_vec
      REAL(KIND=dp), DIMENSION(3, 3)                     :: b_mat, h_inv

      volume = ABS(det_3x3(hmat))
      h_inv = inv_3x3(hmat)
      ! reciprocal vectors b_i = 2π (h⁻¹)_i,: (rows), b_i · a_j = 2π δ_ij
      b_mat = twopi*TRANSPOSE(h_inv)
      omega = SQRT(pi)/volume**(1.0_dp/3.0_dp)
      ! erfc(ω L) < e^{−64} and e^{−G²/4ω²} < e^{−64} beyond these radii
      l_cut = arg_max/omega
      g_cut = 2.0_dp*omega*arg_max
      DO i1 = 1, 3
         ! lattice: |m_i| ≤ L / h_i with the cell height h_i = 2π/|b_i|
         n_l(i1) = CEILING(l_cut*NORM2(b_mat(:, i1))/twopi) + 1
         ! reciprocal: |m_i| = |G·a_i|/2π ≤ G |a_i|/2π
         n_g(i1) = CEILING(g_cut*NORM2(hmat(:, i1))/twopi) + 1
      END DO

      sum_l = 0.0_dp
      DO i1 = -n_l(1), n_l(1)
         DO i2 = -n_l(2), n_l(2)
            DO i3 = -n_l(3), n_l(3)
               IF (i1 == 0 .AND. i2 == 0 .AND. i3 == 0) CYCLE
               l_vec = i1*hmat(:, 1) + i2*hmat(:, 2) + i3*hmat(:, 3)
               r = NORM2(l_vec)
               IF (r > l_cut) CYCLE
               sum_l = sum_l + ERFC(omega*r)/r
            END DO
         END DO
      END DO

      sum_g = 0.0_dp
      DO i1 = -n_g(1), n_g(1)
         DO i2 = -n_g(2), n_g(2)
            DO i3 = -n_g(3), n_g(3)
               IF (i1 == 0 .AND. i2 == 0 .AND. i3 == 0) CYCLE
               g_vec = i1*b_mat(:, 1) + i2*b_mat(:, 2) + i3*b_mat(:, 3)
               g2 = DOT_PRODUCT(g_vec, g_vec)
               IF (g2 > g_cut**2) CYCLE
               sum_g = sum_g + EXP(-g2/(4.0_dp*omega**2))/g2
            END DO
         END DO
      END DO

      xi = 2.0_dp*omega/SQRT(pi) - sum_l - fourpi/volume*sum_g + pi/(volume*omega**2)

   END FUNCTION gamma_corr_madelung_3d

! **************************************************************************************************
!> \brief Sphere-tail operator of the 3D Γ-point correction in the RI basis,
!>
!>          K_PQ = V^tr_PQ − V^E_PQ − ξ c_P c_Q ,
!>
!>        V^tr = (P|θ(R_c − r)/r|Q) the truncated Coulomb operator of Σ^x, V^E = gc%fm_V_zero, ξ the
!>        Madelung constant (stored in gc%xi_madelung) and c_P = ∫ϕ_P the auxiliary charges.
!> \param qs_env QS environment, for the cell and the RI 2c integrals
!> \param bs_env band-structure environment; needs gc%fm_V_zero and the charges
!> \param fm_K K_PQ, created here in the RI x RI structure
! **************************************************************************************************
   SUBROUTINE gamma_corr_tail_operator(qs_env, bs_env, fm_K)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_K

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_tail_operator'

      INTEGER                                            :: handle, i_col, i_row, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: t_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)
      t_start = m_walltime()

      gc => bs_env%ri_rs%gamma_corr
      CPASSERT(.NOT. gc%slab)
      CPASSERT(gc%V_zero_cached)
      CPASSERT(.NOT. bs_env%auto_ri%enabled)

      CALL get_qs_env(qs_env, cell=cell)
      gc%xi_madelung = gamma_corr_madelung_3d(cell%hmat)

      ! V^tr_PQ = (P|v_tr|Q) at Γ, the raw integrals (no M⁻¹) of the Σ^x operator
      CALL RI_2c_integral_mat(qs_env, fm_Vtr, bs_env%fm_struct_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb)

      ! K = V^tr − V^E − ξ c cᵀ with the unnormalized charges c = |c| ĉ
      CALL cp_fm_create(fm_K, bs_env%fm_struct_RI_RI)
      ALLOCATE (c(bs_env%n_RI))
      c(:) = gc%charge_norm*gc%charge
      CALL cp_fm_get_info(fm_K, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      DO i_col = 1, ncol_local
         DO i_row = 1, nrow_local
            fm_K%local_data(i_row, i_col) = fm_Vtr(1, 1)%local_data(i_row, i_col) - &
                                            gc%fm_V_zero%local_data(i_row, i_col) - &
                                            gc%xi_madelung*c(row_indices(i_row))*c(col_indices(i_col))
         END DO
      END DO
      CALL cp_fm_release(fm_Vtr)
      DEALLOCATE (c)
      gc%wall_time = gc%wall_time + m_walltime() - t_start

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_tail_operator

! **************************************************************************************************
!> \brief Measure and reciprocal vectors b_i (b_i · a_j = 2π δ_ij) of a periodic cell of any shape:
!>        slab, the layer area A = |a_1 × a_2| and the in-plane b_1, b_2; bulk, the cell volume
!>        Ω = |det h| and b_i = 2π (h⁻¹)ᵀ e_i.
!> \param bs_env band-structure environment; only the cell matrix hmat is read
!> \param slab .TRUE. for PERIODIC XY, .FALSE. for PERIODIC XYZ
!> \param measure A in bohr² (slab) or Ω in bohr³ (bulk)
!> \param b_recip reciprocal vectors b_recip(1:3, i), Cartesian, 1/bohr; column 3 is zero for a slab
! **************************************************************************************************
   SUBROUTINE gamma_corr_cell_geometry(bs_env, slab, measure, b_recip)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL, INTENT(IN)                                :: slab
      REAL(KIND=dp), INTENT(OUT)                         :: measure
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: b_recip

      REAL(KIND=dp)                                      :: area
      REAL(KIND=dp), DIMENSION(3)                        :: a_1, a_2, normal
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv

      b_recip(:, :) = 0.0_dp

      IF (.NOT. slab) THEN
         ! Ω = |det h|; b_i are the rows of 2π h⁻¹ (b_i · a_j = 2π δ_ij by construction).
         measure = ABS(det_3x3(bs_env%hmat))
         CPASSERT(measure > 0.0_dp)
         h_inv = inv_3x3(bs_env%hmat)
         b_recip(:, :) = twopi*TRANSPOSE(h_inv)
         RETURN
      END IF

      ! The columns of hmat are the Cartesian cell vectors (bohr); a_1, a_2 span the layer.
      a_1 = bs_env%hmat(:, 1)
      a_2 = bs_env%hmat(:, 2)

      ! normal = a_1 × a_2; its length is the layer area A, its direction the layer normal n.
      normal = [a_1(2)*a_2(3) - a_1(3)*a_2(2), &
                a_1(3)*a_2(1) - a_1(1)*a_2(3), &
                a_1(1)*a_2(2) - a_1(2)*a_2(1)]
      area = NORM2(normal)
      CPASSERT(area > 0.0_dp)
      ! n = (a_1 × a_2)/A, the unit normal used to build the in-plane reciprocal vectors.
      normal = normal/area

      ! b_1 = 2π (a_2 × n)/A and b_2 = 2π (n × a_1)/A satisfy b_i·a_j = 2π δ_ij in the plane.
      b_recip(:, 1) = twopi*[a_2(2)*normal(3) - a_2(3)*normal(2), &
                             a_2(3)*normal(1) - a_2(1)*normal(3), &
                             a_2(1)*normal(2) - a_2(2)*normal(1)]/area
      b_recip(:, 2) = twopi*[normal(2)*a_1(3) - normal(3)*a_1(2), &
                             normal(3)*a_1(1) - normal(1)*a_1(3), &
                             normal(1)*a_1(2) - normal(2)*a_1(1)]/area
      measure = area

   END SUBROUTINE gamma_corr_cell_geometry

! **************************************************************************************************
!> \brief Coulomb side of the correction: V(0) (gamma_corr_coulomb_at_zero), the charges
!>        ĉ_P = c_P/|c|, the neutral body V_b = Π V(0) Π and the two maps
!>
!>          F_density = M⁺_c Π V_b^(1/2) ,    G_extern = V(0) Π V_b^(−1/2) ,
!>
!>        with Π_PQ = δ_PQ − ĉ_P ĉ_Q (slab) or Π = 1 (bulk) and V_b^(±1/2) restricted to the
!>        eigenvalues λ > eps_null_V_b λ_max.
!> \param qs_env QS environment, for the 2c Coulomb integrals and the RI_AUX basis
!> \param bs_env band-structure environment; fills bs_env%ri_rs%gamma_corr (charge, maps)
! **************************************************************************************************
   SUBROUTINE gamma_corr_build_coulomb_maps(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_build_coulomb_maps'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: charge_cplx
      INTEGER                                            :: col, handle, i_col, i_row, n_RI, &
                                                            ncol_local, nrow_local, row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: cut, norm_c, t_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, scaling
      TYPE(cp_fm_type)                                   :: fm_inv_sqrt, fm_Pi, fm_sqrt, fm_Vb, &
                                                            fm_vec, fm_work
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)
      t_start = m_walltime()

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI

      ! Same values gamma_corr_init already set; recomputed so this routine is self-contained.
      CALL gamma_corr_cell_geometry(bs_env, gc%slab, gc%measure, gc%b_recip)

      ! 1. V_PQ(0), Ewald sum without G = 0, distributed in gc%fm_V_zero.
      CALL gamma_corr_coulomb_at_zero(qs_env, bs_env)

      ! 2. ĉ_P = c_P/|c|, c_P = ∫ dr ϕ_P(r), replicated.
      ALLOCATE (gc%charge(n_RI))
      ! gamma_auxiliary_charges returns c_P as complex numbers; they are real.
      ALLOCATE (charge_cplx(SIZE(gc%charge)))
      CALL gamma_auxiliary_charges(qs_env, charge_cplx)
      gc%charge(:) = REAL(charge_cplx, dp)
      DEALLOCATE (charge_cplx)
      norm_c = NORM2(gc%charge)
      IF (norm_c <= 0.0_dp) THEN
         CPABORT('Gamma-point correction: the RI_AUX basis carries no charge at all')
      END IF
      gc%charge(:) = gc%charge/norm_c
      gc%charge_norm = norm_c

      ! Π_PQ = δ_PQ − ĉ_P ĉ_Q (slab), Π_PQ = δ_PQ (bulk), on the local block.
      CALL cp_fm_create(fm_Pi, bs_env%fm_struct_RI_RI)
      CALL cp_fm_get_info(fm_Pi, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      DO i_col = 1, ncol_local
         col = col_indices(i_col)
         DO i_row = 1, nrow_local
            row = row_indices(i_row)
            fm_Pi%local_data(i_row, i_col) = 0.0_dp
            IF (gc%slab) fm_Pi%local_data(i_row, i_col) = -gc%charge(row)*gc%charge(col)
            IF (row == col) fm_Pi%local_data(i_row, i_col) = fm_Pi%local_data(i_row, i_col) + 1.0_dp
         END DO
      END DO

      ! 3. V_b = Π V(0) Π
      CALL cp_fm_create(fm_work, bs_env%fm_struct_RI_RI)
      CALL cp_fm_create(fm_Vb, bs_env%fm_struct_RI_RI)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, gc%fm_V_zero, fm_Pi, 0.0_dp, fm_work)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_Pi, fm_work, 0.0_dp, fm_Vb)

      ! V_b = U diag(λ) Uᵀ; eigenvalues λ ≤ cut = eps_null_V_b λ_max are dropped.
      CALL cp_fm_create(fm_vec, bs_env%fm_struct_RI_RI)
      ALLOCATE (eigenvalues(n_RI), scaling(n_RI))
      CALL cp_fm_syevd(fm_Vb, fm_vec, eigenvalues)
      cut = MAX(eps_null_V_b*MAXVAL(eigenvalues), 0.0_dp)
      ! V_b^(1/2) = (U diag(λ^(1/2))) Uᵀ
      CALL cp_fm_create(fm_sqrt, bs_env%fm_struct_RI_RI)
      scaling(:) = 0.0_dp
      WHERE (eigenvalues > cut) scaling = SQRT(eigenvalues)
      CALL cp_fm_to_fm(fm_vec, fm_work)
      CALL cp_fm_column_scale(fm_work, scaling)
      CALL parallel_gemm('N', 'T', n_RI, n_RI, n_RI, 1.0_dp, fm_work, fm_vec, 0.0_dp, fm_sqrt)
      ! V_b^(-1/2) = (U diag(λ^(-1/2))) Uᵀ on the same retained subspace
      CALL cp_fm_create(fm_inv_sqrt, bs_env%fm_struct_RI_RI)
      scaling(:) = 0.0_dp
      WHERE (eigenvalues > cut) scaling = 1.0_dp/SQRT(eigenvalues)
      CALL cp_fm_to_fm(fm_vec, fm_work)
      CALL cp_fm_column_scale(fm_work, scaling)
      CALL parallel_gemm('N', 'T', n_RI, n_RI, n_RI, 1.0_dp, fm_work, fm_vec, 0.0_dp, fm_inv_sqrt)
      CALL cp_fm_release(fm_vec)
      CALL cp_fm_release(fm_Vb)
      DEALLOCATE (eigenvalues, scaling)

      ! 4. F_density = M⁺_c Π V_b^(1/2), with the correction's inverse metric M⁺_c.
      CALL cp_fm_create(gc%fm_F_density, bs_env%fm_struct_RI_RI)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_Pi, fm_sqrt, 0.0_dp, fm_work)
      ! F_density = M⁺_c (Π V_b^(1/2))
      CALL gamma_corr_minv_c_fm(bs_env, fm_work, gc%fm_F_density)
      ! G_extern = V(0) Π V_b^(−1/2), the full V(0) on the left, no metric inverse.
      CALL cp_fm_create(gc%fm_G_extern, bs_env%fm_struct_RI_RI)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_Pi, fm_inv_sqrt, 0.0_dp, fm_work)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, gc%fm_V_zero, fm_work, &
                         0.0_dp, gc%fm_G_extern)

      CALL cp_fm_release(fm_Pi)
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_sqrt)
      CALL cp_fm_release(fm_inv_sqrt)

      gc%wall_time = gc%wall_time + m_walltime() - t_start
      CALL timestop(handle)

   END SUBROUTINE gamma_corr_build_coulomb_maps

! **************************************************************************************************
!> \brief V_PQ(q = 0) of the RI basis, the Ewald sum with the singular G = 0 term omitted.
!>        2D slab (reciprocal half: add_ewald_reciprocal_2d_fm):
!>
!>          V_PQ(0) = Σ_{L∈box} (P0|erfc(ω_s r)/r|QL)
!>                  + Re Σ_p [F_p ρ_P(p) ρ_Q(p)* + Σ_i w_i(p) R_PQ(p, k_i)] ,
!>
!>        p the in-plane reciprocal vectors, F_p = 2π erfc(p/2ω_s)/(A p), F_0 = −2√π/(A ω_s),
!>        Σ_i the Gauss-Legendre quadrature of the k_z integral, ρ_P(k) = ∫ dr ϕ_P(r) e^(−ik·r).
!>        3D bulk (reciprocal half: add_ewald_reciprocal_fm):
!>
!>          V_PQ(0) = Σ_{L∈box} (P0|erfc(ωr)/r|QL)
!>                  + (4π/Ω) Σ_{G≠0} exp(−G²/4ω²)/G² ρ_P(G) ρ_Q(G)* − π/(Ωω²) c_P c_Q .
!> \param qs_env QS environment, for cell, particles and the RI_AUX basis
!> \param bs_env band-structure environment (SIZE_LATTICE_SUM); V(0) is stored in gc%fm_V_zero
! **************************************************************************************************
   SUBROUTINE gamma_corr_coulomb_at_zero(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_coulomb_at_zero'
      REAL(KIND=dp), PARAMETER                           :: eps_imag_V0 = 1.0E-10_dp, &
                                                            omega_height_2d = 5.0_dp, &
                                                            omega_rc_3d = 8.0_dp

      INTEGER                                            :: handle, ikind, ipgf, iset, n_images, &
                                                            n_wave, re_im, size_box
      REAL(KIND=dp)                                      :: alpha_min, height_min, max_imag, &
                                                            max_real, mu_min, omega_eff, r_c
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(kpoint_type), POINTER                         :: kpoints_zero
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      ! Private one-point k-point set at q = 0, leaving kpoints_chi_eps_W untouched.
      NULLIFY (kpoints_zero)
      CALL kpoint_create(kpoints_zero)
      ALLOCATE (kpoints_zero%xkp(3, 1))
      kpoints_zero%xkp(:, 1) = 0.0_dp
      kpoints_zero%nkp = 1
      kpoints_zero%kp_scheme = "GENERAL"
      ! The box of the lattice sum comes from nkp_grid, so keep the production mesh here.
      kpoints_zero%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ! mat_V(1, 1) = Re V(0), mat_V(1, 2) = Im V(0), all blocks reserved and zeroed.
      NULLIFY (mat_V)
      ALLOCATE (mat_V(1:1, 2))
      DO re_im = 1, 2
         NULLIFY (mat_V(1, re_im)%matrix)
         ALLOCATE (mat_V(1, re_im)%matrix)
         CALL dbcsr_create(mat_V(1, re_im)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V(1, re_im)%matrix)
         CALL dbcsr_set(mat_V(1, re_im)%matrix, 0.0_dp)
      END DO

      gc => bs_env%ri_rs%gamma_corr
      IF (gc%slab) THEN
         ! Ewald parameter ω_s = 5/h_min, h_min = min_α 2π/|b_α| over the in-plane b_α.
         height_min = twopi/MAXVAL(NORM2(gc%b_recip(:, 1:2), DIM=1))
         gc%ewald_omega = omega_height_2d/height_min
         ! Most diffuse pair of RI_AUX primitives: μ_min = α_min/2 (a = b = α_min).
         alpha_min = HUGE(1.0_dp)
         DO ikind = 1, SIZE(qs_kind_set)
            NULLIFY (basis)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis, basis_type="RI_AUX")
            IF (.NOT. ASSOCIATED(basis)) CYCLE
            DO iset = 1, basis%nset
               DO ipgf = 1, basis%npgf(iset)
                  alpha_min = MIN(alpha_min, basis%zet(ipgf, iset))
               END DO
            END DO
         END DO
         mu_min = 0.5_dp*alpha_min
         ! Effective erfc decay rate of that pair, ω_eff² = μ ω_s²/(μ + ω_s²).
         omega_eff = SQRT(mu_min*gc%ewald_omega**2/(mu_min + gc%ewald_omega**2))
         ! Image box (nkp_grid 1 1 1): size_lattice_sum 2 → 3x3 if ω_eff h_min ≥ 4, else 3 → 9x9.
         IF (omega_eff*height_min >= 4.0_dp) THEN
            size_box = 2
         ELSE
            size_box = 3
         END IF
         kpoints_zero%nkp_grid = [1, 1, 1]
         ! Real-space half: Σ_{L ∈ box} (P0|erfc(ω_s r)/r|QL), omitted tail ≤ erfc(4).
         CALL build_2c_coulomb_matrix_kp(mat_V, kpoints_zero, basis_type="RI_AUX", cell=cell, &
                                         particle_set=particle_set, qs_kind_set=qs_kind_set, &
                                         atomic_kind_set=atomic_kind_set, &
                                         size_lattice_sum=size_box, &
                                         operator_type=operator_verfc, &
                                         ikp_start=1, ikp_end=1, omega=gc%ewald_omega)
      ELSE
         ! Ewald parameter ω = 8/r_c, r_c = ((N−1)/2) h_min the image-box radius, N = 3^(S−1) images
         ! per direction (S = size_lattice_sum_V), h_min = min_α 2π/|b_α|.
         n_images = 3**(bs_env%size_lattice_sum_V - 1)
         IF (n_images < 3) CPABORT('3D Gamma-point correction: SIZE_LATTICE_SUM must be at least 2')
         height_min = twopi/MAXVAL(NORM2(gc%b_recip, DIM=1))
         r_c = 0.5_dp*REAL(n_images - 1, dp)*height_min
         gc%ewald_omega = omega_rc_3d/r_c
         ! Real-space half: Σ_{R ∈ box} (P0|erfc(ωr)/r|QR), absolutely convergent.
         CALL build_2c_coulomb_matrix_kp(mat_V, kpoints_zero, basis_type="RI_AUX", cell=cell, &
                                         particle_set=particle_set, qs_kind_set=qs_kind_set, &
                                         atomic_kind_set=atomic_kind_set, &
                                         size_lattice_sum=bs_env%size_lattice_sum_V, &
                                         operator_type=operator_verfc, &
                                         ikp_start=1, ikp_end=1, omega=gc%ewald_omega)
      END IF

      ! Re V(0) is kept in gc%fm_V_zero, which the regular W path reuses.
      CALL cp_fm_create(gc%fm_V_zero, bs_env%fm_struct_RI_RI)
      CALL cp_fm_set_all(gc%fm_V_zero, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_V(1, 1)%matrix, gc%fm_V_zero)
      IF (gc%slab) THEN
         ! Reciprocal half at q = 0, singular G = 0 term omitted, added to the local block.
         CALL add_ewald_reciprocal_2d_fm(gc%fm_V_zero, qs_env, "RI_AUX", gc%ewald_omega, n_wave)
      ELSE
         ! Reciprocal half, G = 0 omitted, added to the local block of the distributed matrix.
         CALL add_ewald_reciprocal_fm(gc%fm_V_zero, qs_env, "RI_AUX", gc%ewald_omega, n_wave)
      END IF
      gc%V_zero_cached = .TRUE.
      CALL cp_fm_maxabsval(gc%fm_V_zero, max_real)

      ! Im V(0), checked on the distributed matrix.
      CALL cp_fm_create(fm_work, bs_env%fm_struct_RI_RI)
      CALL cp_fm_set_all(fm_work, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_V(1, 2)%matrix, fm_work)
      CALL cp_fm_maxabsval(fm_work, max_imag)
      CALL cp_fm_release(fm_work)
      CALL dbcsr_deallocate_matrix_set(mat_V)
      DEALLOCATE (kpoints_zero%xkp)
      CALL kpoint_release(kpoints_zero)

      ! Im V(0) = 0 for a centrosymmetric image box.
      IF (max_imag > eps_imag_V0*MAX(1.0_dp, max_real)) THEN
         CPABORT('Gamma-point correction: V(q=0) came back complex')
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_coulomb_at_zero

! **************************************************************************************************
!> \brief Decide whether the Γ-point correction runs and set slab or bulk, the 3D sphere-tail
!>        switch, the measure (layer area A or cell volume Ω) and the reciprocal vectors b_i.
!> \param bs_env band-structure environment; sets gamma_corr%enabled, %slab, %tail_3d, %measure
!>        and %b_recip of bs_env%ri_rs
! **************************************************************************************************
   SUBROUTINE gamma_corr_init(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      TYPE(gamma_corr_type), POINTER                     :: gc

      gc => bs_env%ri_rs%gamma_corr

      ! On for the periodic Γ-only RI-RS driver (PERIODIC XY or XYZ), off for molecules.
      gc%enabled = bs_env%gw_implementation == ri_rs_large_cell_Gamma
      IF (.NOT. gc%enabled) RETURN

      ! Slab (2π/q head over the layer area) or bulk (4π/q² head over the cell volume).
      gc%slab = ALL(bs_env%periodic == [1, 1, 0])
      ! Sphere-tail term of the 3D head (gamma_corr_tail_operator): off for slabs and AUTO_RI; an
      ! anisotropic head switches it off later (gamma_corr_direction_coefficients).
      gc%tail_3d = .NOT. (gc%slab .OR. bs_env%auto_ri%enabled)
      ! A or Ω and the reciprocal vectors.
      CALL gamma_corr_cell_geometry(bs_env, gc%slab, gc%measure, gc%b_recip)

   END SUBROUTINE gamma_corr_init

! **************************************************************************************************
!> \brief AO matrices of the commutator [H, r^α] = −∂^α + [V_nl, r^α], α = x, y, z, built as in
!>        the COMMUTATOR_HR print path of qs_scf_post_gpw.
!> \param qs_env QS environment, for the AO overlap template, neighbor lists and pseudopotential
!> \param moments the three antisymmetric AO matrices, allocated here (caller deallocates)
! **************************************************************************************************
   SUBROUTINE gamma_corr_commutator_matrices(qs_env, moments)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments

      INTEGER                                            :: i_dir
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (matrix_s, particle_set, qs_kind_set, sab_orb, sap_ppnl, cell, dft_control)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, sab_orb=sab_orb, sap_ppnl=sap_ppnl, &
                      cell=cell, dft_control=dft_control)

      ! Antisymmetric (real, anti-Hermitian operator), sparsity of sab_orb, zero-initialized.
      CALL dbcsr_allocate_matrix_set(moments, 3)
      DO i_dir = 1, 3
         ALLOCATE (moments(i_dir)%matrix)
         CALL dbcsr_create(moments(i_dir)%matrix, template=matrix_s(1)%matrix, &
                           name='COMMUTATOR_HR', matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(moments(i_dir)%matrix, sab_orb)
         CALL dbcsr_set(moments(i_dir)%matrix, 0.0_dp)
      END DO
      ! Accumulate [r^α, V_nl] (build_com_mom_nl, matrix_rv) ...
      CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                            particle_set, cell, matrix_rv=moments)
      ! ... plus ⟨µ|∂^α|ν⟩ (build_lin_mom_matrix), ...
      CALL build_lin_mom_matrix(qs_env, moments)
      ! ... and negate: −(∂^α + [r^α, V_nl]) = −∂^α + [V_nl, r^α] = [H, r^α].
      DO i_dir = 1, 3
         CALL dbcsr_scale(moments(i_dir)%matrix, alpha_scalar=-1.0_dp)
      END DO

   END SUBROUTINE gamma_corr_commutator_matrices

! **************************************************************************************************
!> \brief Optical transition dipoles between occupied i and empty a at Γ, in the velocity gauge,
!>
!>          µ^α_ia = ⟨a|[H, r^α]|i⟩/(ε_a − ε_i),   [H, r^α] = −∂^α + [V_nl, r^α] ,
!>
!>        built as distributed n_occ x n_vir matrices (gamma_corr_ov_block).
!> \param qs_env QS environment, for the AO overlap template, neighbor lists and pseudopotential
!> \param bs_env band-structure environment, for the Γ MO coefficients and eigenvalues
!> \param n_dir number of Cartesian directions α = 1 .. n_dir (2 for a slab in the xy plane)
!> \param fm_mu fm_mu(α, spin) = µ^α_ia, n_occ x n_vir (row i, column a − n_occ), created here
! **************************************************************************************************
   SUBROUTINE gamma_corr_transition_dipoles(qs_env, bs_env, n_dir, fm_mu)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: n_dir
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: fm_mu

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_transition_dipoles'

      INTEGER                                            :: handle, i_col, i_dir, i_row, i_spin, &
                                                            n_occ, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: delta
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments

      CALL timeset(routineN, handle)

      ! ⟨a|r^α|i⟩ = ⟨a|[H, r^α]|i⟩/(ε_a − ε_i): the position operator itself is not used under PBC.

      ! moments(α) = AO matrix of [H, r^α], α = x, y, z (gamma_corr_commutator_matrices).
      NULLIFY (moments)
      CALL gamma_corr_commutator_matrices(qs_env, moments)

      ALLOCATE (fm_mu(n_dir, bs_env%n_spin))
      DO i_spin = 1, bs_env%n_spin
         n_occ = bs_env%n_occ(i_spin)
         DO i_dir = 1, n_dir
            ! ⟨a|[H, r^α]|i⟩ for occupied i (rows) and empty a (columns)
            CALL gamma_corr_ov_block(bs_env, i_spin, moments(i_dir)%matrix, fm_mu(i_dir, i_spin))
            ! µ^α_ia = ⟨a|[H, r^α]|i⟩ / Δ_ia on the local elements, Δ_ia = ε_a − ε_i > 0
            CALL cp_fm_get_info(fm_mu(i_dir, i_spin), nrow_local=nrow_local, &
                                ncol_local=ncol_local, row_indices=row_indices, &
                                col_indices=col_indices)
            DO i_col = 1, ncol_local
               DO i_row = 1, nrow_local
                  delta = bs_env%eigenval_scf_Gamma(n_occ + col_indices(i_col), i_spin) - &
                          bs_env%eigenval_scf_Gamma(row_indices(i_row), i_spin)
                  CPASSERT(delta > 0.0_dp)
                  fm_mu(i_dir, i_spin)%local_data(i_row, i_col) = &
                     fm_mu(i_dir, i_spin)%local_data(i_row, i_col)/delta
               END DO
            END DO
         END DO
      END DO

      CALL dbcsr_deallocate_matrix_set(moments)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_transition_dipoles

! **************************************************************************************************
!> \brief Occupied-empty block of an AO matrix X (any DBCSR symmetry) at Γ,
!>
!>          O_ia = ⟨a|X|i⟩ = Σ_µν C_µa X_µν C_νi ,
!>          i ≤ n_occ (rows),  a = n_occ + 1 .. n_mo (columns),
!>
!>        computed as G = X C_occ and O = Gᵀ C_vir.
!> \param bs_env band-structure environment, for the Γ MO coefficients
!> \param i_spin spin channel
!> \param mat_ao the AO matrix X (DBCSR, any symmetry)
!> \param fm_ov O, n_occ x n_vir, distributed; created here
! **************************************************************************************************
   SUBROUTINE gamma_corr_ov_block(bs_env, i_spin, mat_ao, fm_ov)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: i_spin
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_ao
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_ov

      INTEGER                                            :: n_ao, n_mo, n_occ, n_vir
      TYPE(cp_fm_struct_type), POINTER                   :: struct_ao_occ, struct_ov
      TYPE(cp_fm_type)                                   :: fm_C_occ, fm_G
      TYPE(dbcsr_type)                                   :: mat_full

      n_ao = bs_env%n_ao
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = n_ao
      n_occ = bs_env%n_occ(i_spin)
      n_vir = n_mo - n_occ
      CPASSERT(n_occ > 0 .AND. n_vir > 0)

      NULLIFY (struct_ao_occ, struct_ov)
      CALL cp_fm_struct_create(struct_ao_occ, &
                               template_fmstruct=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct, &
                               ncol_global=n_occ)
      CALL cp_fm_struct_create(struct_ov, &
                               template_fmstruct=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct, &
                               nrow_global=n_occ, ncol_global=n_vir)
      CALL cp_fm_create(fm_G, struct_ao_occ)
      CALL cp_fm_create(fm_C_occ, struct_ao_occ)
      CALL cp_fm_create(fm_ov, struct_ov)
      CALL cp_fm_struct_release(struct_ao_occ)
      CALL cp_fm_struct_release(struct_ov)

      ! X with full (non-symmetric) storage
      CALL dbcsr_create(mat_full, template=mat_ao, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(mat_ao, mat_full)
      ! G = X C_occ, C_occ the first n_occ MO columns (copied to its own matrix)
      CALL cp_fm_to_fm_submat(bs_env%fm_mo_coeff_Gamma(i_spin), fm_C_occ, n_ao, n_occ, 1, 1, 1, 1)
      CALL cp_dbcsr_sm_fm_multiply(mat_full, fm_C_occ, fm_G, ncol=n_occ)
      CALL dbcsr_release(mat_full)
      CALL cp_fm_release(fm_C_occ)
      ! O = Gᵀ C_vir, C_vir the MO columns n_occ+1 .. n_mo: O_ia = Σ_µ (X C)_µi C_µa = ⟨a|X|i⟩
      CALL parallel_gemm('T', 'N', n_occ, n_vir, n_ao, 1.0_dp, fm_G, &
                         bs_env%fm_mo_coeff_Gamma(i_spin), 0.0_dp, fm_ov, b_first_col=n_occ + 1)
      CALL cp_fm_release(fm_G)

   END SUBROUTINE gamma_corr_ov_block

! **************************************************************************************************
!> \brief FINITE_Q head (3D bulk): replace the velocity-gauge head S(iω) by one extrapolated from
!>        the exact head without local fields at q = g_d, 2g_d along six directions g_d,
!>
!>          Y(g, iω)    = (2/n_s/Ω) Σ_σ Σ_ia |⟨a|e^(i g·r)|i⟩|²/|g|² · 2Δ_ia/(Δ_ia² + ω²) ,
!>          a(q̂, q, iω) = q̂ᵀ S(iω) q̂ / (1 + γ(iω) q²) ,
!>
!>        with Δ_ia = ε_a − ε_i, a fit 1/Y(q) = A + B q² per direction and the local-field ratio
!>        r_d = ĝ_dᵀ S_vel ĝ_d / ĝ_dᵀ Y_vel ĝ_d of the velocity-gauge heads.
!> \param qs_env QS environment, for the Berry-phase AO matrices
!> \param bs_env band-structure environment; gc%S_head is overwritten, gc%head_gamma is set
!> \param Y_vel velocity-gauge head WITHOUT local fields, Y(iω), 3 x 3 x n_freq
! **************************************************************************************************
   SUBROUTINE gamma_corr_finite_q_head(qs_env, bs_env, Y_vel)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Y_vel

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_finite_q_head'
      INTEGER, PARAMETER                                 :: n_fit_dir = 6

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: handle, i_a, i_col, i_d, i_i, i_row, &
                                                            i_spin, i_w, info, m, n_bad, n_good, &
                                                            n_mo, n_occ, n_t, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(n_fit_dir)                      :: ipiv
      LOGICAL, DIMENSION(n_fit_dir)                      :: good_d
      REAL(KIND=dp)                                      :: coef_A, coef_B, delta, g_norm2, measure, &
                                                            pair_w, r_lf, spin_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: lorentz
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: y_q
      REAL(KIND=dp), DIMENSION(3)                        :: g_hat, kvec
      REAL(KIND=dp), DIMENSION(3, 3)                     :: b_recip
      REAL(KIND=dp), DIMENSION(3, n_fit_dir)             :: g_dir
      REAL(KIND=dp), DIMENSION(n_fit_dir)                :: gamma_d, s_d, y0_d
      REAL(KIND=dp), DIMENSION(n_fit_dir, 1)             :: rhs
      REAL(KIND=dp), DIMENSION(n_fit_dir, n_fit_dir)     :: dyad, dyad_work
      TYPE(cp_fm_type)                                   :: fm_mo_cos, fm_mo_sin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_t = bs_env%num_time_freq_points
      CALL gamma_corr_cell_geometry(bs_env, .FALSE., measure, b_recip)
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = bs_env%n_ao
      spin_factor = 2.0_dp/REAL(bs_env%n_spin, dp)

      ! The six fit directions g_d (supercell reciprocal vectors, 1/bohr).
      g_dir(:, 1:3) = b_recip(:, 1:3)
      g_dir(:, 4) = b_recip(:, 1) + b_recip(:, 2)
      g_dir(:, 5) = b_recip(:, 2) + b_recip(:, 3)
      g_dir(:, 6) = b_recip(:, 3) + b_recip(:, 1)

      ! Berry-phase AO matrices ⟨µ|cos(k·r)|ν⟩, ⟨µ|sin(k·r)|ν⟩ on the overlap pattern.
      NULLIFY (matrix_s, sab_orb)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_create(cosmat, template=matrix_s(1)%matrix, name='FINITE_Q COS')
      CALL dbcsr_create(sinmat, template=matrix_s(1)%matrix, name='FINITE_Q SIN')
      CALL cp_dbcsr_alloc_block_from_nbl(cosmat, sab_orb)
      CALL cp_dbcsr_alloc_block_from_nbl(sinmat, sab_orb)
      ! ⟨a|Cos|i⟩ and ⟨a|Sin|i⟩ share one distributed layout; each rank sums Y(g, iω) over its
      ! local (i, a), then y_q is summed over ranks.
      ALLOCATE (lorentz(n_t))

      ! 1. y_q(i_d, m, i_w) = Y(m g_d, iω_j)
      !      = (2/n_s/Ω) Σ_ia [(Cᵀ Cos C)_ai² + (Cᵀ Sin C)_ai²]/|m g_d|² · 2Δ_ia/(Δ_ia² + ω_j²)
      ALLOCATE (y_q(n_fit_dir, 2, n_t))
      y_q(:, :, :) = 0.0_dp
      DO i_d = 1, n_fit_dir
         DO m = 1, 2
            kvec(:) = REAL(m, dp)*g_dir(:, i_d)
            g_norm2 = SUM(kvec**2)
            CALL dbcsr_set(cosmat, 0.0_dp)
            CALL dbcsr_set(sinmat, 0.0_dp)
            CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
            DO i_spin = 1, bs_env%n_spin
               n_occ = bs_env%n_occ(i_spin)
               CALL gamma_corr_ov_block(bs_env, i_spin, cosmat, fm_mo_cos)
               CALL gamma_corr_ov_block(bs_env, i_spin, sinmat, fm_mo_sin)
               CALL cp_fm_get_info(fm_mo_cos, nrow_local=nrow_local, ncol_local=ncol_local, &
                                   row_indices=row_indices, col_indices=col_indices)
               ! Local elements (row i occupied, column a empty)
               DO i_col = 1, ncol_local
                  i_a = n_occ + col_indices(i_col)
                  DO i_row = 1, nrow_local
                     i_i = row_indices(i_row)
                     ! Δ_ia = ε_a − ε_i (Hartree)
                     delta = bs_env%eigenval_scf_Gamma(i_a, i_spin) - &
                             bs_env%eigenval_scf_Gamma(i_i, i_spin)
                     ! (2/n_s/Ω) |M_ai|²/|q|², Y(q, iω) without its frequency factor
                     pair_w = spin_factor*(fm_mo_cos%local_data(i_row, i_col)**2 + &
                                           fm_mo_sin%local_data(i_row, i_col)**2)/ &
                              g_norm2/measure
                     ! 2Δ/(Δ² + ω_j²) for all frequencies
                     lorentz(:) = 2.0_dp*delta/ &
                                  (delta**2 + bs_env%time_frequency_grid%frequency(1:n_t)**2)
                     y_q(i_d, m, :) = y_q(i_d, m, :) + pair_w*lorentz(:)
                  END DO
               END DO
               CALL cp_fm_release(fm_mo_cos)
               CALL cp_fm_release(fm_mo_sin)
            END DO
         END DO
      END DO
      CALL dbcsr_release(cosmat)
      CALL dbcsr_release(sinmat)
      DEALLOCATE (cosmat, sinmat, lorentz)
      ! Complete Y(q, iω) over the (a, i) elements held by the other ranks.
      CALL bs_env%fm_mo_coeff_Gamma(1)%matrix_struct%para_env%sum(y_q)

      ! 2. 6 x 6 system ĝ_dᵀ S ĝ_d = s_d for the unknowns (S_xx, S_yy, S_zz, S_xy, S_yz, S_xz).
      DO i_d = 1, n_fit_dir
         g_hat(:) = g_dir(:, i_d)/NORM2(g_dir(:, i_d))
         dyad(i_d, :) = [g_hat(1)**2, g_hat(2)**2, g_hat(3)**2, 2.0_dp*g_hat(1)*g_hat(2), &
                         2.0_dp*g_hat(2)*g_hat(3), 2.0_dp*g_hat(1)*g_hat(3)]
      END DO

      IF (ALLOCATED(gc%head_gamma)) DEALLOCATE (gc%head_gamma)
      ALLOCATE (gc%head_gamma(n_t))
      n_bad = 0
      DO i_w = 1, n_t
         DO i_d = 1, n_fit_dir
            g_norm2 = SUM(g_dir(:, i_d)**2)
            g_hat(:) = g_dir(:, i_d)/SQRT(g_norm2)
            ! 3. 1/Y(q) = A + B q² through (g, Y_1) and (2g, Y_2):
            !      B = (1/Y_2 − 1/Y_1)/(3|g|²),  A = 1/Y_1 − B|g|²,  Y_0 = 1/A,  γ_d = B/A
            coef_B = (1.0_dp/y_q(i_d, 2, i_w) - 1.0_dp/y_q(i_d, 1, i_w))/(3.0_dp*g_norm2)
            coef_A = 1.0_dp/y_q(i_d, 1, i_w) - coef_B*g_norm2
            good_d(i_d) = coef_A > 0.0_dp .AND. coef_B >= 0.0_dp
            IF (good_d(i_d)) THEN
               y0_d(i_d) = 1.0_dp/coef_A
               gamma_d(i_d) = coef_B/coef_A
            ELSE
               ! Not of screening form (A ≤ 0 or B < 0): keep the q = g value, no dispersion.
               n_bad = n_bad + 1
               y0_d(i_d) = y_q(i_d, 1, i_w)
               gamma_d(i_d) = 0.0_dp
            END IF
            ! 4. s_d = r_d Y_0,d with the local-field ratio of the velocity-gauge head along ĝ_d,
            !      r_d = ĝ_dᵀ S_vel ĝ_d / ĝ_dᵀ Y_vel ĝ_d
            r_lf = DOT_PRODUCT(g_hat, MATMUL(gc%S_head(:, :, i_w), g_hat))/ &
                   DOT_PRODUCT(g_hat, MATMUL(Y_vel(:, :, i_w), g_hat))
            s_d(i_d) = r_lf*y0_d(i_d)
         END DO
         ! 5. S(iω) from the six directional values s_d, γ(iω) the average of the γ_d;
         !      a(q̂, q, iω) = q̂ᵀ S q̂/(1 + γ q²)
         rhs(:, 1) = s_d(:)
         ! dgesv overwrites its matrix with the LU factors; dyad is reused at every frequency.
         dyad_work(:, :) = dyad(:, :)
         CALL dgesv(n_fit_dir, 1, dyad_work, n_fit_dir, ipiv, rhs, n_fit_dir, info)
         IF (info /= 0) CPABORT('Gamma-point correction: finite-q head tensor fit failed')
         gc%S_head(:, :, i_w) = RESHAPE([rhs(1, 1), rhs(4, 1), rhs(6, 1), &
                                         rhs(4, 1), rhs(2, 1), rhs(5, 1), &
                                         rhs(6, 1), rhs(5, 1), rhs(3, 1)], [3, 3])
         ! γ: average over the directions of screening form only.
         n_good = COUNT(good_d)
         IF (n_good > 0) THEN
            gc%head_gamma(i_w) = SUM(gamma_d, MASK=good_d)/REAL(n_good, dp)
         ELSE
            gc%head_gamma(i_w) = 0.0_dp
         END IF
      END DO

      IF (n_bad > 0) THEN
         WRITE (message, '(A,I0,A,I0,A)') 'Gamma-point correction, finite-q head: ', n_bad, ' of ', &
            n_fit_dir*n_t, ' direction/frequency fits are not of screening form; the q = g value '// &
            'is kept there without dispersion.'
         CPWARN(TRIM(message))
      END IF

      DEALLOCATE (y_q)
      CALL timestop(handle)

   END SUBROUTINE gamma_corr_finite_q_head

! **************************************************************************************************
!> \brief Optical response of the Γ-only supercell: head tensor Y_αβ and Cartesian wing grid vectors
!>        X_ℓα in imaginary time (Δ_ia = ε_a − ε_i; A = layer area (slab) or cell volume Ω (bulk)),
!>
!>          Y_αβ(iτ) = (2/n_spin/A)  Σ_σ Σ_ia µ^α_ai µ^β_ai e^(−Δ_ia τ)
!>          X_ℓα(iτ) = (2/n_spin/√A) Σ_σ Σ_ia ψ_i(r_ℓ) ψ_a(r_ℓ) µ^α_ai e^(−Δ_ia τ)
!>
!>        cosine-transformed to iω and reduced to the head S (gc%S_head) and wing map T (gc%T_RI),
!>
!>          S_αβ(iω) = Y_αβ − [Eᵀ B⁻¹ E]_αβ ,   T_Pα(iω) = [G_extern B⁻¹ E]_Pα ,
!>          B = 1 + F_densityᵀ χ_PQ(Γ,iω) F_density ,   E = F_densityᵀ (Zᵀ X) ,
!>
!>        followed by the finite-q head: ρ(q, iω) in 2D, γ(iω) in 3D (where S_head is replaced).
!> \param bs_env band-structure environment; S_head and T_RI are stored in its gamma_corr
!> \param qs_env QS environment, for the transition dipoles
!> \param mat_chi_Gamma_tau χ_PQ(Γ,iτ), the regular polarizability the driver already built
! **************************************************************************************************
   SUBROUTINE gamma_corr_response(bs_env, qs_env, mat_chi_Gamma_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_response'

      INTEGER                                            :: handle, i_dir, i_spin, i_t, i_w, n_dir, &
                                                            n_grid, n_RI, n_t
      REAL(KIND=dp)                                      :: leak_cc, spin_factor, t_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: leak_v0
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: X_freq, X_time, Y_head_vel, Y_time
      TYPE(cp_fm_type)                                   :: fm_chi
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_mu
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)
      t_start = m_walltime()

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI
      n_grid = bs_env%ri_rs%n_grid_points
      n_t = bs_env%num_time_freq_points

      ! Spin factor 2/n_spin per channel (2 for a closed shell).
      spin_factor = 2.0_dp/REAL(bs_env%n_spin, dp)

      ! Cartesian directions: a slab in the xy plane needs only x and y (in-plane q̂); z stays zero.
      n_dir = 3
      IF (gc%slab .AND. ALL(gc%b_recip(3, 1:2) == 0.0_dp)) n_dir = 2

      ! µ^α_ia for all spins, distributed n_occ x n_vir per direction and spin.
      CALL gamma_corr_transition_dipoles(qs_env, bs_env, n_dir, fm_mu)

      ! Y_time(α, β, i_t): head, replicated; X_time(ℓ, α, i_t): wings, partial per rank.
      ALLOCATE (Y_time(3, 3, n_t), X_time(n_grid, 3, n_t))
      Y_time(:, :, :) = 0.0_dp
      X_time(:, :, :) = 0.0_dp

      ! Sum over the spin channels.
      DO i_spin = 1, bs_env%n_spin
         CALL gamma_corr_accumulate_response(bs_env, i_spin, n_dir, fm_mu, spin_factor, &
                                             Y_time, X_time)
      END DO
      DO i_spin = 1, bs_env%n_spin
         DO i_dir = 1, n_dir
            CALL cp_fm_release(fm_mu(i_dir, i_spin))
         END DO
      END DO
      DEALLOCATE (fm_mu)

      ! Complete X over the grid rows held by the other ranks (Y is complete already).
      CALL bs_env%para_env%sum(X_time)

      ! Cosine transform iτ → iω: F(iω_j) = Σ_k γ^cos_jk F(iτ_k), γ^cos = cosine_time_to_frequency.
      ! S_head holds Y(iω) here; −Eᵀ B⁻¹ E is subtracted per frequency below.
      ALLOCATE (gc%S_head(3, 3, n_t), X_freq(n_grid, 3, n_t))
      gc%S_head(:, :, :) = 0.0_dp
      X_freq(:, :, :) = 0.0_dp
      DO i_w = 1, n_t
         DO i_t = 1, n_t
            gc%S_head(:, :, i_w) = gc%S_head(:, :, i_w) + &
                                   bs_env%time_frequency_grid%cosine_time_to_frequency(i_w, i_t)* &
                                   Y_time(:, :, i_t)
         END DO
      END DO
      ! X_ℓα(iω_j) = Σ_k γ^cos_jk X_ℓα(iτ_k), threaded over the independent columns (α, j).
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(static) PRIVATE(i_w, i_dir, i_t) &
!$OMP          SHARED(n_t, X_freq, X_time, bs_env)
      DO i_w = 1, n_t
         DO i_dir = 1, 3
            DO i_t = 1, n_t
               X_freq(:, i_dir, i_w) = X_freq(:, i_dir, i_w) + &
                                       bs_env%time_frequency_grid%cosine_time_to_frequency(i_w, i_t)* &
                                       X_time(:, i_dir, i_t)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      DEALLOCATE (Y_time, X_time)

      ! Velocity-gauge Y(iω) without local fields, kept for the finite-q head.
      ALLOCATE (Y_head_vel(3, 3, n_t))
      Y_head_vel(:, :, :) = gc%S_head(:, :, :)

      ALLOCATE (gc%T_RI(n_RI, 3, n_t))
      gc%T_RI(:, :, :) = 0.0_dp
      CALL cp_fm_create(fm_chi, bs_env%fm_struct_RI_RI)

      ! Charge-leak diagnostic: v0 = M⁺_c(Γ) ĉ (χ_cc = v0ᵀ χ_PQ(iω_1) v0 in gamma_corr_response_fm).
      ALLOCATE (leak_v0(n_RI, 1))
      CALL gamma_corr_minv_c_rep(bs_env, RESHAPE(gc%charge, [n_RI, 1]), leak_v0)
      leak_cc = 0.0_dp

      ! Per frequency: body B, wings E, Schur complement S and wing vector T.
      CALL gamma_corr_response_fm(bs_env, fm_chi, mat_chi_Gamma_tau, X_freq, leak_v0, leak_cc)
      CALL cp_fm_release(fm_chi)
      DEALLOCATE (X_freq)
      ! The two maps are consumed; free them before W and Σ are built.
      CALL cp_fm_release(gc%fm_F_density)
      CALL cp_fm_release(gc%fm_G_extern)
      IF (gc%slab) THEN
         ! 2D: ρ(q, iω) fitted to the exact pair weights; S_head is kept.
         CALL gamma_corr_finite_q_head_2d(qs_env, bs_env, Y_head_vel)
      ELSE
         ! 3D: S(iω) and γ(iω) from the exact finite-q head (replaces S_head).
         CALL gamma_corr_finite_q_head(qs_env, bs_env, Y_head_vel)
      END IF
      DEALLOCATE (Y_head_vel)
      CALL gamma_corr_leak_diagnostic(bs_env, leak_cc)
      DEALLOCATE (leak_v0)

      gc%wall_time = gc%wall_time + m_walltime() - t_start
      CALL timestop(handle)

   END SUBROUTINE gamma_corr_response

! **************************************************************************************************
!> \brief 2D finite-q head (slabs): per frequency, a 1/q²-weighted least-squares fit of
!>        ρ(q) = c/(1 + γ q²) to ρ_k = y(g_k, iω)/ĝ_kᵀ Y_vel(iω) ĝ_k at the in-plane shells g_k = b₁,
!>        the longer of b₁ ± b₂ and 2b₁, with the exact no-local-field head
!>
!>          y(g, iω) = (2/n_s/A) Σ_σ Σ_ia |⟨a|e^(i g·r)|i⟩|²/|g|² · 2Δ_ia/(Δ_ia² + ω²) .
!>
!>        The head is then a(q, q̂, iω) = q̂ᵀ S(iω) q̂ ρ(q, iω); ρ = 1 where the fit is not of
!>        screening form.
!> \param qs_env QS environment, for the Berry-phase AO matrices
!> \param bs_env band-structure environment; sets gc%q_form and gc%q_par
!> \param Y_vel velocity-gauge head without local fields, Y(iω), 3 x 3 x n_freq
! **************************************************************************************************
   SUBROUTINE gamma_corr_finite_q_head_2d(qs_env, bs_env, Y_vel)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Y_vel

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_finite_q_head_2d'
      INTEGER, PARAMETER                                 :: n_shell = 3

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: handle, i_a, i_col, i_i, i_k, i_row, &
                                                            i_spin, i_w, n_bad, n_mo, n_occ, n_t, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: good
      REAL(KIND=dp)                                      :: delta, det, g_norm2, measure, p1, p2, &
                                                            pair_w, spin_factor, wk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: lorentz
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rho_q, y_q
      REAL(KIND=dp), DIMENSION(2)                        :: rhs
      REAL(KIND=dp), DIMENSION(2, 2)                     :: normal
      REAL(KIND=dp), DIMENSION(3)                        :: g_hat
      REAL(KIND=dp), DIMENSION(3, 3)                     :: b_recip
      REAL(KIND=dp), DIMENSION(3, n_shell)               :: g_shell
      REAL(KIND=dp), DIMENSION(n_shell)                  :: f1, f2, q_shell, t_k
      TYPE(cp_fm_type)                                   :: fm_mo_cos, fm_mo_sin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_t = bs_env%num_time_freq_points
      ! layer area A and the in-plane reciprocal vectors b₁, b₂ (third column zero)
      CALL gamma_corr_cell_geometry(bs_env, .TRUE., measure, b_recip)
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = bs_env%n_ao
      spin_factor = 2.0_dp/REAL(bs_env%n_spin, dp)

      ! the three shells: g₁ = b₁, the √3 shell (longer of b₁ ± b₂), 2g₁
      g_shell(:, 1) = b_recip(:, 1)
      IF (SUM((b_recip(:, 1) + b_recip(:, 2))**2) > SUM((b_recip(:, 1) - b_recip(:, 2))**2)) THEN
         g_shell(:, 2) = b_recip(:, 1) + b_recip(:, 2)
      ELSE
         g_shell(:, 2) = b_recip(:, 1) - b_recip(:, 2)
      END IF
      g_shell(:, 3) = 2.0_dp*b_recip(:, 1)
      DO i_k = 1, n_shell
         q_shell(i_k) = NORM2(g_shell(:, i_k))
      END DO

      ! 1. y_q(i_k, i_w) = y(g_k, iω_j), same pair-weight loop as gamma_corr_finite_q_head
      NULLIFY (matrix_s, sab_orb)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_create(cosmat, template=matrix_s(1)%matrix, name='FINITE_Q_2D COS')
      CALL dbcsr_create(sinmat, template=matrix_s(1)%matrix, name='FINITE_Q_2D SIN')
      CALL cp_dbcsr_alloc_block_from_nbl(cosmat, sab_orb)
      CALL cp_dbcsr_alloc_block_from_nbl(sinmat, sab_orb)
      ALLOCATE (lorentz(n_t), y_q(n_shell, n_t), rho_q(n_shell, n_t))
      y_q(:, :) = 0.0_dp
      DO i_k = 1, n_shell
         g_norm2 = q_shell(i_k)**2
         CALL dbcsr_set(cosmat, 0.0_dp)
         CALL dbcsr_set(sinmat, 0.0_dp)
         CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, g_shell(:, i_k))
         DO i_spin = 1, bs_env%n_spin
            n_occ = bs_env%n_occ(i_spin)
            ! occupied-empty blocks ⟨a|Cos|i⟩, ⟨a|Sin|i⟩ (gamma_corr_ov_block)
            CALL gamma_corr_ov_block(bs_env, i_spin, cosmat, fm_mo_cos)
            CALL gamma_corr_ov_block(bs_env, i_spin, sinmat, fm_mo_sin)
            CALL cp_fm_get_info(fm_mo_cos, nrow_local=nrow_local, ncol_local=ncol_local, &
                                row_indices=row_indices, col_indices=col_indices)
            ! local elements (row i occupied, column a empty)
            DO i_col = 1, ncol_local
               i_a = n_occ + col_indices(i_col)
               DO i_row = 1, nrow_local
                  i_i = row_indices(i_row)
                  ! Δ_ia = ε_a − ε_i
                  delta = bs_env%eigenval_scf_Gamma(i_a, i_spin) - bs_env%eigenval_scf_Gamma(i_i, i_spin)
                  ! (2/n_s/A) |M_ai(g)|²/|g|²
                  pair_w = spin_factor*(fm_mo_cos%local_data(i_row, i_col)**2 + &
                                        fm_mo_sin%local_data(i_row, i_col)**2)/g_norm2/measure
                  ! 2Δ/(Δ² + ω²) at all frequencies
                  lorentz(:) = 2.0_dp*delta/(delta**2 + bs_env%time_frequency_grid%frequency(1:n_t)**2)
                  y_q(i_k, :) = y_q(i_k, :) + pair_w*lorentz(:)
               END DO
            END DO
            CALL cp_fm_release(fm_mo_cos)
            CALL cp_fm_release(fm_mo_sin)
         END DO
      END DO
      CALL dbcsr_release(cosmat)
      CALL dbcsr_release(sinmat)
      DEALLOCATE (cosmat, sinmat, lorentz)
      CALL bs_env%fm_mo_coeff_Gamma(1)%matrix_struct%para_env%sum(y_q)

      ! 2. Per frequency: ρ_k = y_k / ĝ_kᵀ Y_vel ĝ_k and weighted least squares t_k = 1/ρ_k ≈
      !    p1 + p2 q_k² with weights w_k = (ρ_k/q_k)² (1/q²-weighted relative residuals of ρ).
      gc%q_form = 1
      IF (ALLOCATED(gc%q_par)) DEALLOCATE (gc%q_par)
      ALLOCATE (gc%q_par(2, n_t))
      n_bad = 0
      DO i_w = 1, n_t
         DO i_k = 1, n_shell
            g_hat(:) = g_shell(:, i_k)/q_shell(i_k)
            rho_q(i_k, i_w) = y_q(i_k, i_w)/DOT_PRODUCT(g_hat, MATMUL(Y_vel(:, :, i_w), g_hat))
            t_k(i_k) = 1.0_dp/rho_q(i_k, i_w)
            f1(i_k) = 1.0_dp
            f2(i_k) = q_shell(i_k)**2
         END DO
         normal(:, :) = 0.0_dp
         rhs(:) = 0.0_dp
         DO i_k = 1, n_shell
            wk = (rho_q(i_k, i_w)/q_shell(i_k))**2
            normal(1, 1) = normal(1, 1) + wk*f1(i_k)*f1(i_k)
            normal(1, 2) = normal(1, 2) + wk*f1(i_k)*f2(i_k)
            normal(2, 2) = normal(2, 2) + wk*f2(i_k)*f2(i_k)
            rhs(1) = rhs(1) + wk*f1(i_k)*t_k(i_k)
            rhs(2) = rhs(2) + wk*f2(i_k)*t_k(i_k)
         END DO
         normal(2, 1) = normal(1, 2)
         det = normal(1, 1)*normal(2, 2) - normal(1, 2)**2
         p1 = (normal(2, 2)*rhs(1) - normal(1, 2)*rhs(2))/det
         p2 = (normal(1, 1)*rhs(2) - normal(1, 2)*rhs(1))/det
         ! ρ = c/(1 + γq²) with c = 1/p1, γ = p2/p1; screening form needs p1 > 0, p2 ≥ 0
         good = p1 > 0.0_dp .AND. p2 >= 0.0_dp
         IF (good) THEN
            gc%q_par(:, i_w) = [1.0_dp/p1, p2/p1]
         ELSE
            gc%q_par(:, i_w) = [1.0_dp, 0.0_dp]
         END IF
         IF (.NOT. good) n_bad = n_bad + 1
      END DO

      IF (n_bad > 0) THEN
         WRITE (message, '(A,I0,A,I0,A)') 'Gamma-point correction, 2D finite-q head: the fit at ', &
            n_bad, ' of ', n_t, ' frequencies is not of screening form; rho = 1 is used there.'
         CPWARN(TRIM(message))
      END IF
      DEALLOCATE (y_q, rho_q)
      CALL timestop(handle)

   END SUBROUTINE gamma_corr_finite_q_head_2d

! **************************************************************************************************
!> \brief Ratio ρ(q, iω) = c(iω)/(1 + γ(iω) q²) of the in-plane head to q̂ᵀ S(iω) q̂
!>        (fitted in gamma_corr_finite_q_head_2d); 1 without the 2D finite-q head.
!> \param gc Γ-correction data
!> \param q in-plane |q| (1/bohr)
!> \param i_w frequency index
!> \return ρ(q, iω)
! **************************************************************************************************
   PURE FUNCTION gamma_corr_head_ratio(gc, q, i_w) RESULT(rho)
      TYPE(gamma_corr_type), INTENT(IN)                  :: gc
      REAL(KIND=dp), INTENT(IN)                          :: q
      INTEGER, INTENT(IN)                                :: i_w
      REAL(KIND=dp)                                      :: rho

      IF (gc%q_form == 1) THEN
         ! c/(1 + γ q²)
         rho = gc%q_par(1, i_w)/(1.0_dp + gc%q_par(2, i_w)*q*q)
      ELSE
         rho = 1.0_dp
      END IF
   END FUNCTION gamma_corr_head_ratio

! **************************************************************************************************
!> \brief Diagnostic only (no result changes): the RI-fitted occupied-virtual pair densities carry a
!>        spurious charge that adds 2πλ/q (2D) or 4πλ/q² (3D) to the head of ε, with
!>
!>          λ = |c|² χ_cc / A (2D)  or  |c|² χ_cc / Ω (3D),
!>          χ_cc = v0ᵀ χ_PQ(iω_1) v0,   v0 = M⁺_c(Γ) ĉ,   ĉ = c/|c| .
!>
!>        Warns if the estimated static gap reduction Δ from |q| > q_min (a lower bound) exceeds 5 meV.
!> \param bs_env band-structure environment; uses gc%charge_norm, measure, b_recip, slab, S_head
!> \param chi_cc χ_cc = v0ᵀ χ_PQ(iω_1) v0
! **************************************************************************************************
   SUBROUTINE gamma_corr_leak_diagnostic(bs_env, chi_cc)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: chi_cc

      REAL(KIND=dp), PARAMETER                           :: warn_ev = 0.005_dp

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: n1, n2, n3, n3_max
      REAL(KIND=dp)                                      :: alpha, delta, eps_inf, lambda, q_min, &
                                                            q_s, x
      REAL(KIND=dp), DIMENSION(3)                        :: g_vec
      TYPE(gamma_corr_type), POINTER                     :: gc

      gc => bs_env%ri_rs%gamma_corr

      ! λ from the leaked charge; χ_cc ≥ 0 up to round-off.
      lambda = gc%charge_norm**2*MAX(chi_cc, 0.0_dp)/gc%measure

      ! q_min: shortest nonzero supercell reciprocal vector, from the combinations n_i ∈ {−1,0,1}
      ! (in-plane only for a slab).
      n3_max = 1
      IF (gc%slab) n3_max = 0
      q_min = HUGE(1.0_dp)
      DO n1 = -1, 1
         DO n2 = -1, 1
            DO n3 = -n3_max, n3_max
               IF (n1 == 0 .AND. n2 == 0 .AND. n3 == 0) CYCLE
               g_vec = MATMUL(gc%b_recip, REAL([n1, n2, n3], dp))
               q_min = MIN(q_min, NORM2(g_vec))
            END DO
         END DO
      END DO

      delta = 0.0_dp
      IF (gc%slab) THEN
         ! 2D: Δ = 2π λ [ln((1+x)/x) − 1/(1+x)], x = 2π α q_min
         alpha = 0.5_dp*(gc%S_head(1, 1, 1) + gc%S_head(2, 2, 1))
         IF (alpha > 0.0_dp) THEN
            x = twopi*alpha*q_min
            delta = twopi*lambda*(LOG((1.0_dp + x)/x) - 1.0_dp/(1.0_dp + x))
         END IF
      ELSE
         ! 3D: Δ = (8 λ/ε_∞²)(1/q_s)[π/2 − arctan(q_min/q_s)], q_s² = 4π λ/ε_∞
         eps_inf = 1.0_dp + fourpi*(gc%S_head(1, 1, 1) + gc%S_head(2, 2, 1) + gc%S_head(3, 3, 1))/3.0_dp
         IF (lambda > 0.0_dp .AND. eps_inf > 0.0_dp) THEN
            q_s = SQRT(fourpi*lambda/eps_inf)
            delta = 8.0_dp*lambda/eps_inf**2/q_s*(0.5_dp*pi - ATAN(q_min/q_s))
         END IF
      END IF

      IF (delta*evolt > warn_ev) THEN
         WRITE (message, '(A,F7.3,A)') 'RI charge leak of the occupied-virtual pair densities: '// &
            'the estimated long-wavelength part of the static gap error at this cell size is ', &
            -delta*evolt, ' eV. The total error can be several times larger and grows with the '// &
            'cell size. Use a larger RI basis or a longer-range RI metric.'
         CALL cp_warn(__LOCATION__, TRIM(message))
      END IF

   END SUBROUTINE gamma_corr_leak_diagnostic

! **************************************************************************************************
!> \brief Per-frequency body, wings, Schur complement and wing vector of gamma_corr_response:
!>        E = F_dᵀ Zᵀ X (n_RI x 3), x = B⁻¹ E with the body B = 1 + F_dᵀ χ(iω) F_d (by conjugate
!>        gradients, B is not formed), S −= Eᵀ x and T = G_extern x.
!> \param bs_env band-structure environment; uses gc%fm_F_density, gc%fm_G_extern, adds to
!>        gc%S_head and fills gc%T_RI
!> \param fm_chi work matrix in the RI x RI layout
!> \param mat_chi_Gamma_tau χ_PQ(Γ, iτ)
!> \param X_freq X_ℓα(iω), replicated over grid rows
!> \param leak_v0 v0 = M⁺_c(Γ) ĉ, replicated (charge-leak diagnostic)
!> \param leak_cc χ_cc = v0ᵀ χ_PQ(iω_1) v0 (charge-leak diagnostic)
! **************************************************************************************************
   SUBROUTINE gamma_corr_response_fm(bs_env, fm_chi, mat_chi_Gamma_tau, X_freq, leak_v0, leak_cc)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_chi
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: X_freq
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: leak_v0
      REAL(KIND=dp), INTENT(OUT)                         :: leak_cc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_response_fm'

      INTEGER                                            :: handle, i_dir, i_w, j_dir, n_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: chi_v0, E_rep, x_rep, ZX_rep
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI

      ALLOCATE (E_rep(n_RI, 3), x_rep(n_RI, 3), ZX_rep(n_RI, 3))

      DO i_w = 1, bs_env%num_time_freq_points
         ! χ_PQ(Γ, iω) distributed in fm_chi.
         CALL gamma_corr_chi_PQ_at_freq(bs_env, fm_chi, i_w, mat_chi_Gamma_tau)
         ! Charge-leak diagnostic: χ_cc = v0ᵀ χ_PQ(iω_1) v0.
         IF (i_w == 1) THEN
            ALLOCATE (chi_v0(n_RI, 1))
            CALL gamma_corr_fm_times_rep('N', fm_chi, leak_v0, chi_v0)
            leak_cc = DOT_PRODUCT(leak_v0(:, 1), chi_v0(:, 1))
            DEALLOCATE (chi_v0)
         END IF

         ! E = F_dᵀ (Zᵀ X), with Zᵀ X (n_RI x 3) replicated.
         CALL gamma_corr_project_grid_to_RI(bs_env, X_freq(:, :, i_w), ZX_rep)
         CALL gamma_corr_fm_times_rep('T', gc%fm_F_density, ZX_rep, E_rep)

         ! x = B⁻¹ E, B = 1 + F_dᵀ χ F_d, by conjugate gradients
         CALL gamma_corr_body_cg(gc%fm_F_density, fm_chi, E_rep, x_rep)

         ! S_αβ = Y_αβ − Σ_P E_Pα x_Pβ
         DO j_dir = 1, 3
            DO i_dir = 1, 3
               gc%S_head(i_dir, j_dir, i_w) = gc%S_head(i_dir, j_dir, i_w) - &
                                              DOT_PRODUCT(E_rep(:, i_dir), x_rep(:, j_dir))
            END DO
         END DO

         ! T_Pα(iω) = [G_extern x]_Pα, the wing vector in the RI projection representation.
         CALL gamma_corr_fm_times_rep('N', gc%fm_G_extern, x_rep, gc%T_RI(:, :, i_w))
      END DO

      DEALLOCATE (E_rep, x_rep, ZX_rep)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_response_fm

! **************************************************************************************************
!> \brief x = B⁻¹ E for the dielectric body B = 1 + F_dᵀ χ F_d without forming B: conjugate
!>        gradients on the columns of E in lockstep with B p = p + F_dᵀ (χ (F_d p)), until
!>        |r_j| ≤ 1e-12 |E_j| for every column. Aborts after 300 iterations or if p·Bp ≤ 0.
!> \param fm_F F_d, n_RI x n_RI distributed
!> \param fm_chi χ(iω), n_RI x n_RI distributed
!> \param E right-hand sides, n_RI x n_k replicated
!> \param x solutions, n_RI x n_k replicated
! **************************************************************************************************
   SUBROUTINE gamma_corr_body_cg(fm_F, fm_chi, E, x)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_F, fm_chi
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: E
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: x

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_body_cg'
      INTEGER, PARAMETER                                 :: max_iter = 300
      REAL(KIND=dp), PARAMETER                           :: tol = 1.0E-12_dp

      CHARACTER(LEN=default_string_length)               :: msg
      INTEGER                                            :: handle, iter, j, n_k
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: active
      REAL(KIND=dp)                                      :: alpha, beta, pq, rr_new
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rr, rr0
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: p, q, r, y1, y2

      CALL timeset(routineN, handle)

      n_k = SIZE(E, 2)
      ALLOCATE (p(SIZE(E, 1), n_k), q(SIZE(E, 1), n_k), r(SIZE(E, 1), n_k), &
                y1(SIZE(E, 1), n_k), y2(SIZE(E, 1), n_k), rr(n_k), rr0(n_k), active(n_k))
      x(:, :) = 0.0_dp
      r(:, :) = E(:, :)
      p(:, :) = E(:, :)
      DO j = 1, n_k
         rr(j) = DOT_PRODUCT(r(:, j), r(:, j))
         rr0(j) = rr(j)
         active(j) = rr0(j) > 0.0_dp
      END DO

      iter = 0
      DO WHILE (ANY(active))
         iter = iter + 1
         IF (iter > max_iter) THEN
            WRITE (msg, '(A,I0,A,ES10.2)') 'Gamma-point correction: CG for the dielectric body did '// &
               'not converge in ', max_iter, ' iterations; largest relative residual ', &
               MAXVAL(SQRT(rr/MAX(rr0, TINY(1.0_dp))), MASK=active)
            CPABORT(TRIM(msg))
         END IF
         ! q = B p = p + F_dᵀ (χ (F_d p))
         CALL gamma_corr_fm_times_rep('N', fm_F, p, y1)
         CALL gamma_corr_fm_times_rep('N', fm_chi, y1, y2)
         CALL gamma_corr_fm_times_rep('T', fm_F, y2, q)
         q(:, :) = q(:, :) + p(:, :)
         DO j = 1, n_k
            IF (.NOT. active(j)) CYCLE
            pq = DOT_PRODUCT(p(:, j), q(:, j))
            IF (pq <= 0.0_dp) THEN
               CALL cp_abort(__LOCATION__, 'Gamma-point correction: dielectric body not positive definite')
            END IF
            alpha = rr(j)/pq
            x(:, j) = x(:, j) + alpha*p(:, j)
            r(:, j) = r(:, j) - alpha*q(:, j)
            rr_new = DOT_PRODUCT(r(:, j), r(:, j))
            IF (rr_new <= tol*tol*rr0(j)) THEN
               active(j) = .FALSE.
            ELSE
               beta = rr_new/rr(j)
               p(:, j) = r(:, j) + beta*p(:, j)
            END IF
            rr(j) = rr_new
         END DO
      END DO

      DEALLOCATE (p, q, r, y1, y2, rr, rr0, active)
      CALL timestop(handle)

   END SUBROUTINE gamma_corr_body_cg

! **************************************************************************************************
!> \brief Y = op(A) X for a distributed cp_fm A and a replicated skinny X (n x k), with Y replicated
!>        on every rank of A: one DGEMM on the local block of A and one sum over the ranks.
!>        Replaces parallel_gemm, which can abort in COSMA for such skinny operands.
!> \param trans 'N' for Y = A X, 'T' for Y = Aᵀ X
!> \param fm_A distributed matrix
!> \param X replicated right-hand side, SIZE(X, 1) = columns of op(A)
!> \param Y replicated result, SIZE(Y, 1) = rows of op(A)
! **************************************************************************************************
   SUBROUTINE gamma_corr_fm_times_rep(trans, fm_A, X, Y)
      CHARACTER(LEN=1), INTENT(IN)                       :: trans
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_A
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: X
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: Y

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_fm_times_rep'

      INTEGER                                            :: handle, i_col, i_row, k, lda, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: X_loc, Y_loc

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_A, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      k = SIZE(X, 2)
      CPASSERT(SIZE(Y, 2) == k)
      Y(:, :) = 0.0_dp

      IF (nrow_local > 0 .AND. ncol_local > 0) THEN
         lda = SIZE(fm_A%local_data, 1)
         IF (trans == 'N' .OR. trans == 'n') THEN
            ! Y_loc = A_loc X(cols of A_loc, :), then Y(rows of A_loc, :) += Y_loc
            ALLOCATE (X_loc(ncol_local, k), Y_loc(nrow_local, k))
            DO i_col = 1, ncol_local
               X_loc(i_col, :) = X(col_indices(i_col), :)
            END DO
            CALL dgemm('N', 'N', nrow_local, k, ncol_local, 1.0_dp, fm_A%local_data, lda, &
                       X_loc, ncol_local, 0.0_dp, Y_loc, nrow_local)
            DO i_row = 1, nrow_local
               Y(row_indices(i_row), :) = Y(row_indices(i_row), :) + Y_loc(i_row, :)
            END DO
         ELSE
            ! Y_loc = A_locᵀ X(rows of A_loc, :), then Y(cols of A_loc, :) += Y_loc
            ALLOCATE (X_loc(nrow_local, k), Y_loc(ncol_local, k))
            DO i_row = 1, nrow_local
               X_loc(i_row, :) = X(row_indices(i_row), :)
            END DO
            CALL dgemm('T', 'N', ncol_local, k, nrow_local, 1.0_dp, fm_A%local_data, lda, &
                       X_loc, nrow_local, 0.0_dp, Y_loc, ncol_local)
            DO i_col = 1, ncol_local
               Y(col_indices(i_col), :) = Y(col_indices(i_col), :) + Y_loc(i_col, :)
            END DO
         END IF
         DEALLOCATE (X_loc, Y_loc)
      END IF

      ! Each element of A contributed on exactly one rank: the sum over A's ranks completes Y.
      CALL fm_A%matrix_struct%para_env%sum(Y)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_fm_times_rep

! **************************************************************************************************
!> \brief Y = M⁺_c X for a replicated X (n_RI x k), with the inverse metric of the correction
!>        M⁺_c = M⁺ − U_L diag(1/m_L) U_Lᵀ; M⁺ = bs_env%fm_Minv_Gamma, and U_L, m_L are the modes
!>        between the global cut and REGULARIZATION_RI_CORRECTION (compute_Minv_Gamma).
!> \param bs_env band-structure environment
!> \param X replicated right-hand side, n_RI x k
!> \param Y replicated result, n_RI x k
! **************************************************************************************************
   SUBROUTINE gamma_corr_minv_c_rep(bs_env, X, Y)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: X
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: Y

      INTEGER                                            :: i_l
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: D, G
      TYPE(gamma_corr_type), POINTER                     :: gc

      gc => bs_env%ri_rs%gamma_corr
      ! Y = M⁺ X
      CALL gamma_corr_fm_times_rep('N', bs_env%fm_Minv_Gamma, X, Y)
      IF (gc%n_low == 0) RETURN
      ALLOCATE (G(gc%n_low, SIZE(X, 2)), D(SIZE(X, 1), SIZE(X, 2)))
      ! G = diag(1/m_L) U_Lᵀ X
      CALL gamma_corr_fm_times_rep('T', gc%fm_U_low, X, G)
      DO i_l = 1, gc%n_low
         G(i_l, :) = G(i_l, :)/gc%m_low(i_l)
      END DO
      ! Y −= U_L G
      CALL gamma_corr_fm_times_rep('N', gc%fm_U_low, G, D)
      Y(:, :) = Y - D
      DEALLOCATE (G, D)

   END SUBROUTINE gamma_corr_minv_c_rep

! **************************************************************************************************
!> \brief fm_Y = M⁺_c fm_X for a distributed n_RI x n_RI matrix X (see gamma_corr_minv_c_rep):
!>        Y = M⁺ X − [U_L diag(1/m_L)] (U_Lᵀ X).
!> \param bs_env band-structure environment
!> \param fm_X distributed right-hand side (RI x RI structure)
!> \param fm_Y distributed result, created by the caller
! **************************************************************************************************
   SUBROUTINE gamma_corr_minv_c_fm(bs_env, fm_X, fm_Y)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_X, fm_Y

      INTEGER                                            :: n_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: inv_m
      TYPE(cp_fm_struct_type), POINTER                   :: struct_lr
      TYPE(cp_fm_type)                                   :: fm_G, fm_Us
      TYPE(gamma_corr_type), POINTER                     :: gc

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI
      ! Y = M⁺ X
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, bs_env%fm_Minv_Gamma, fm_X, &
                         0.0_dp, fm_Y)
      IF (gc%n_low == 0) RETURN
      ! G = U_Lᵀ X (n_low x n_RI)
      NULLIFY (struct_lr)
      CALL cp_fm_struct_create(struct_lr, template_fmstruct=fm_X%matrix_struct, &
                               nrow_global=gc%n_low)
      CALL cp_fm_create(fm_G, struct_lr)
      CALL cp_fm_struct_release(struct_lr)
      CALL parallel_gemm('T', 'N', gc%n_low, n_RI, n_RI, 1.0_dp, gc%fm_U_low, fm_X, 0.0_dp, fm_G)
      ! U_s = U_L diag(1/m_L);  Y −= U_s G
      CALL cp_fm_create(fm_Us, gc%fm_U_low%matrix_struct)
      CALL cp_fm_to_fm(gc%fm_U_low, fm_Us)
      ALLOCATE (inv_m(gc%n_low))
      inv_m(:) = 1.0_dp/gc%m_low(:)
      CALL cp_fm_column_scale(fm_Us, inv_m)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, gc%n_low, -1.0_dp, fm_Us, fm_G, 1.0_dp, fm_Y)
      DEALLOCATE (inv_m)
      CALL cp_fm_release(fm_G)
      CALL cp_fm_release(fm_Us)

   END SUBROUTINE gamma_corr_minv_c_fm

! **************************************************************************************************
!> \brief Adds one spin channel to the head tensor Y_αβ(iτ) and the wing grid vectors X_ℓα(iτ)
!>        (A = gc%measure, e_i = e^(−|ε_i−ε_F|τ), e_a = e^(−|ε_a−ε_F|τ), so e^(−Δ_ia τ) = e_i e_a):
!>
!>          Y_αβ(iτ) = (2/n_spin/A) Σ_ia µ^α_ia µ^β_ia e_i e_a
!>          X_ℓα(iτ) = (2/n_spin/√A) Σ_ia ψ_i(r_ℓ) ψ_a(r_ℓ) µ^α_ia e_i e_a
!>                   = (2/n_spin/√A) Σ_i A_ℓi B^α_ℓi ,
!>          A = Φ C_occ ,   B^α = Φ (C_vir D^αᵀ) ,   D^α = diag(e_i) µ^α diag(e_a) .
!>
!>        X is partial over the grid rows of this rank; the caller sums it over the ranks.
!> \param bs_env band-structure environment, for eigenvalues, ε_F, MO coefficients, τ grid,
!>        Φ = mat_phi_mu_l and the measure A
!> \param i_spin spin channel
!> \param n_dir number of Cartesian directions in fm_mu
!> \param fm_mu µ^α_ia, fm_mu(α, spin) distributed n_occ x n_vir (gamma_corr_transition_dipoles)
!> \param spin_factor 2/n_spin
!> \param Y_time accumulated head tensor, Y(α,β,iτ)
!> \param X_time accumulated wing grid vectors, X(ℓ,α,iτ), partial over grid rows
! **************************************************************************************************
   SUBROUTINE gamma_corr_accumulate_response(bs_env, i_spin, n_dir, fm_mu, spin_factor, &
                                             Y_time, X_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: i_spin, n_dir
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: fm_mu
      REAL(KIND=dp), INTENT(IN)                          :: spin_factor
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Y_time, X_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_accumulate_response'

      INTEGER                                            :: handle, i_a, i_col, i_dir, i_i, i_row, &
                                                            i_t, j_dir, n_ao, n_grid, n_mo, n_occ, &
                                                            n_vir, ncol_ab, ncol_local, nrow_ab, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, col_indices_ab, &
                                                            row_indices, row_indices_ab
      REAL(KIND=dp)                                      :: e_fermi, prefactor, tau, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exp_occ, exp_vir, x_partial
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Y_loc
      TYPE(cp_fm_struct_type), POINTER                   :: struct_ao_occ, struct_grid_occ
      TYPE(cp_fm_type)                                   :: fm_A, fm_B, fm_D, fm_Q

      CALL timeset(routineN, handle)

      ! MO range: occupied 1..n_occ, empty n_occ+1..n_mo (retained MOs, all n_ao if no window).
      n_ao = bs_env%n_ao
      n_occ = bs_env%n_occ(i_spin)
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = n_ao
      n_vir = n_mo - n_occ
      n_grid = bs_env%ri_rs%n_grid_points
      ! Reference energy for the factorized exponential (Hartree).
      e_fermi = bs_env%e_fermi(i_spin)
      CPASSERT(n_occ > 0 .AND. n_vir > 0)

      ! Work matrices: A, B (n_grid x n_occ), Q (n_AO x n_occ) and D (n_occ x n_vir, layout of fm_mu).
      NULLIFY (struct_grid_occ, struct_ao_occ)
      CALL cp_fm_struct_create(struct_grid_occ, &
                               template_fmstruct=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct, &
                               nrow_global=n_grid, ncol_global=n_occ)
      CALL cp_fm_struct_create(struct_ao_occ, &
                               template_fmstruct=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct, &
                               ncol_global=n_occ)
      CALL cp_fm_create(fm_A, struct_grid_occ)
      CALL cp_fm_create(fm_B, struct_grid_occ)
      CALL cp_fm_create(fm_Q, struct_ao_occ)
      CALL cp_fm_create(fm_D, fm_mu(1, i_spin)%matrix_struct)
      CALL cp_fm_struct_release(struct_grid_occ)
      CALL cp_fm_struct_release(struct_ao_occ)
      CALL cp_fm_get_info(fm_D, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      CALL cp_fm_get_info(fm_A, nrow_local=nrow_ab, ncol_local=ncol_ab, &
                          row_indices=row_indices_ab, col_indices=col_indices_ab)

      ! A = Φ C_occ, the occupied orbitals on the grid (C_occ copied into fm_Q)
      CALL cp_fm_to_fm_submat(bs_env%fm_mo_coeff_Gamma(i_spin), fm_Q, n_ao, n_occ, 1, 1, 1, 1)
      CALL cp_dbcsr_sm_fm_multiply(bs_env%ri_rs%mat_phi_mu_l, fm_Q, fm_A, ncol=n_occ)

      ! prefactor = (2/n_spin)/√A of X
      ALLOCATE (exp_occ(n_occ), exp_vir(n_vir), x_partial(n_grid), Y_loc(3, 3, SIZE(Y_time, 3)))
      Y_loc(:, :, :) = 0.0_dp
      prefactor = spin_factor/SQRT(bs_env%ri_rs%gamma_corr%measure)

      DO i_t = 1, bs_env%num_time_freq_points
         ! τ_k of the minimax grid (Hartree⁻¹).
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         ! exp_occ(i) = e^(−|ε_i − ε_F| τ), exp_vir(a) = e^(−|ε_a − ε_F| τ); their product is
         ! e^(−Δ_ia τ) since ε_i < ε_F < ε_a.
         DO i_i = 1, n_occ
            exp_occ(i_i) = EXP(-ABS(bs_env%eigenval_scf_Gamma(i_i, i_spin) - e_fermi)*tau)
         END DO
         DO i_a = 1, n_vir
            exp_vir(i_a) = EXP(-ABS(bs_env%eigenval_scf_Gamma(n_occ + i_a, i_spin) - e_fermi)*tau)
         END DO

         ! Y_αβ(iτ) += (2/n_spin/A) Σ_ia µ^α_ia µ^β_ia e^(−Δ_ia τ) on the local elements
         DO i_col = 1, ncol_local
            i_a = col_indices(i_col)
            DO i_row = 1, nrow_local
               i_i = row_indices(i_row)
               weight = spin_factor*exp_occ(i_i)*exp_vir(i_a)/bs_env%ri_rs%gamma_corr%measure
               DO j_dir = 1, n_dir
                  DO i_dir = 1, n_dir
                     Y_loc(i_dir, j_dir, i_t) = Y_loc(i_dir, j_dir, i_t) + &
                                                weight*fm_mu(i_dir, i_spin)%local_data(i_row, i_col)* &
                                                fm_mu(j_dir, i_spin)%local_data(i_row, i_col)
                  END DO
               END DO
            END DO
         END DO

         DO i_dir = 1, n_dir
            ! D_ia = e_i µ^α_ia e_a on the local elements
            DO i_col = 1, ncol_local
               DO i_row = 1, nrow_local
                  fm_D%local_data(i_row, i_col) = exp_occ(row_indices(i_row))* &
                                                  fm_mu(i_dir, i_spin)%local_data(i_row, i_col)* &
                                                  exp_vir(col_indices(i_col))
               END DO
            END DO
            ! Q = C_vir Dᵀ (n_AO x n_occ), C_vir the MO columns n_occ+1 .. n_mo
            CALL parallel_gemm('N', 'T', n_ao, n_occ, n_vir, 1.0_dp, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), fm_D, 0.0_dp, fm_Q, &
                               a_first_col=n_occ + 1)
            ! B = Φ Q
            CALL cp_dbcsr_sm_fm_multiply(bs_env%ri_rs%mat_phi_mu_l, fm_Q, fm_B, ncol=n_occ)
            ! X_ℓα += (2/n_spin/√A) Σ_i A_ℓi B^α_ℓi on this rank's elements
            x_partial(:) = 0.0_dp
            DO i_col = 1, ncol_ab
               DO i_row = 1, nrow_ab
                  x_partial(row_indices_ab(i_row)) = x_partial(row_indices_ab(i_row)) + &
                                                     fm_A%local_data(i_row, i_col)* &
                                                     fm_B%local_data(i_row, i_col)
               END DO
            END DO
            X_time(:, i_dir, i_t) = X_time(:, i_dir, i_t) + prefactor*x_partial
         END DO
      END DO

      ! Complete Y over the (i, a) elements held by the other ranks.
      CALL fm_mu(1, i_spin)%matrix_struct%para_env%sum(Y_loc)
      Y_time(:, :, :) = Y_time + Y_loc

      DEALLOCATE (exp_occ, exp_vir, x_partial, Y_loc)
      CALL cp_fm_release(fm_A)
      CALL cp_fm_release(fm_B)
      CALL cp_fm_release(fm_Q)
      CALL cp_fm_release(fm_D)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_accumulate_response

! **************************************************************************************************
!> \brief χ_PQ(Γ,iω_j) = Σ_k γ^cos_jk χ_PQ(Γ,iτ_k), γ^cos = cosine_time_to_frequency, accumulated
!>        in DBCSR and copied to the distributed fm_chi.
!> \param bs_env band-structure environment; its mat_RI_RI is used as DBCSR scratch
!> \param fm_chi scratch full matrix with the RI x RI structure
!> \param j_w frequency index
!> \param mat_chi_Gamma_tau χ_PQ(Γ,iτ_k), one DBCSR matrix per time point
! **************************************************************************************************
   SUBROUTINE gamma_corr_chi_PQ_at_freq(bs_env, fm_chi, j_w, mat_chi_Gamma_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_chi
      INTEGER, INTENT(IN)                                :: j_w
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      INTEGER                                            :: i_t
      REAL(KIND=dp)                                      :: weight

      ! mat_RI_RI ← 0, then += γ^cos_jk χ(iτ_k) for every time point k.
      CALL dbcsr_set(bs_env%mat_RI_RI%matrix, 0.0_dp)
      DO i_t = 1, bs_env%num_time_freq_points
         weight = bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, i_t)
         CALL dbcsr_add(bs_env%mat_RI_RI%matrix, mat_chi_Gamma_tau(i_t)%matrix, 1.0_dp, weight)
      END DO
      ! DBCSR → block-cyclic full matrix.
      CALL copy_dbcsr_to_fm(bs_env%mat_RI_RI%matrix, fm_chi)

   END SUBROUTINE gamma_corr_chi_PQ_at_freq

! **************************************************************************************************
!> \brief (Zᵀ X)_Qα = Σ_ℓ Z_ℓQ X_ℓα for a three-column grid quantity X, from the local blocks of
!>        the sparse Z and one sum over the ranks.
!> \param bs_env band-structure environment, for mat_Z_lP and the communicator
!> \param X_grid X_ℓα, replicated over grid rows
!> \param X_RI result, n_RI x 3, replicated
! **************************************************************************************************
   SUBROUTINE gamma_corr_project_grid_to_RI(bs_env, X_grid, X_RI)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: X_grid
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: X_RI

      INTEGER                                            :: col, col_start, i_col, i_dir, row, &
                                                            row_start
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      X_RI(:, :) = 0.0_dp

      ! X_RI(Q, α) += Σ_{ℓ∈block} Z_ℓQ X_ℓα for each local block of Z, one dot product per column Q.
      CALL dbcsr_iterator_start(iter, bs_env%ri_rs%mat_Z_lP)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_offset=row_start, col_offset=col_start)
         DO i_dir = 1, 3
            DO i_col = 1, SIZE(block, 2)
               X_RI(col_start + i_col - 1, i_dir) = X_RI(col_start + i_col - 1, i_dir) + &
                                                    DOT_PRODUCT(block(:, i_col), &
                                                                X_grid(row_start:row_start + SIZE(block, 1) - 1, i_dir))
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! Each rank summed only its blocks of Z; complete the ℓ sum across ranks.
      CALL bs_env%para_env%sum(X_RI)

   END SUBROUTINE gamma_corr_project_grid_to_RI

! **************************************************************************************************
!> \brief Slab: tabulates the head and wing kernels, integrated over the Brillouin-zone cell around
!>        Γ, on a regular mesh of in-plane minimum-image displacements R, per frequency:
!>
!>          K_h(R,iω)   = (A/(2π)²) ∮dq̂ w(q̂) ∫_0^{q_max} dq q · [−4π² a/(1 + 2πq a)] cos(q q̂·R)
!>          K_w,α(R,iω) = (A/(2π)²) ∮dq̂ w(q̂) q̂_α ∫_0^{q_max} dq q · [2π/(1 + 2πq a)] sin(q q̂·R)
!>
!>        a = q̂ᵀ S q̂ · ρ(q, iω); K_h1, K_h2, K_w1,α have an extra factor q, q², q. Then calls
!>        gamma_corr_tabulate_remainder_2d. Bulk: only calls gamma_corr_direction_coefficients.
!> \param bs_env band-structure environment; the tables are stored in its gamma_corr
!> \param n_tab points per direction of the tabulation mesh
!> \param qs_env QS environment, for the atoms and the RI_AUX basis of the zone remainder
! **************************************************************************************************
   SUBROUTINE gamma_corr_tabulate(bs_env, n_tab, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: n_tab
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_tabulate'
      INTEGER, PARAMETER :: n_face = 24, n_radial = 32, n_tab_rem = gamma_corr_n_tab_rem

      INTEGER                                            :: handle, i_dir, i_node, i_tab, i_u, i_v, &
                                                            i_w, n_dir, n_t
      REAL(KIND=dp)                                      :: a_q, a_scalar, area_bz, denom, jac, &
                                                            proj, q, q_max, radial_w, start_time, &
                                                            term_h, term_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cos_q, dir_radius, dir_weight, kh, kh1, &
                                                            kh2, nodes, rule, sin_q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dir_hat, kw, kw1, q_node
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: coef_h, coef_w
      REAL(KIND=dp), DIMENSION(3)                        :: R_cart
      TYPE(gamma_corr_type), POINTER                     :: gc

      ! n_tab_rem: points per direction of the remainder tables (the remainder is smooth in R)

      CALL timeset(routineN, handle)

      start_time = m_walltime()
      gc => bs_env%ri_rs%gamma_corr
      n_t = bs_env%num_time_freq_points

      ! Bulk cell: no tables; head constant and A, B tensors by quadrature.
      IF (.NOT. gc%slab) THEN
         CALL gamma_corr_direction_coefficients(bs_env)
         gc%wall_time = gc%wall_time + m_walltime() - start_time
         CALL timestop(handle)
         RETURN
      END IF

      ! Area of the small Brillouin zone, (2π)²/A; dividing by it gives the A/(2π)² prefactor.
      area_bz = twopi**2/gc%measure

      ! Angular quadrature (q̂, w(q̂), q_max(q̂)) and the radial Gauss-Legendre rule on [−1, 1].
      CALL gamma_corr_ws_cell(bs_env, .TRUE., n_face, dir_hat, dir_weight, dir_radius)
      n_dir = SIZE(dir_weight)
      ALLOCATE (nodes(n_radial), rule(n_radial))
      CALL gauss_legendre(n_radial, nodes, rule)

      ! Square n_tab x n_tab table, flattened: K_head(i_R, i_w), K_wing(i_R, α, i_w).
      gc%n_tab = [n_tab, n_tab]
      ALLOCATE (gc%K_head(n_tab*n_tab, n_t), gc%K_wing(n_tab*n_tab, 3, n_t))
      gc%K_head(:, :) = 0.0_dp
      gc%K_wing(:, :, :) = 0.0_dp
      ! Out-of-plane expansion: the same integrands times q and q² (K_h1, K_h2, K_w1).
      ALLOCATE (gc%K_head1(n_tab*n_tab, n_t), gc%K_head2(n_tab*n_tab, n_t), &
                gc%K_wing1(n_tab*n_tab, 3, n_t))
      gc%K_head1(:, :) = 0.0_dp
      gc%K_head2(:, :) = 0.0_dp
      gc%K_wing1(:, :, :) = 0.0_dp

      ! R-independent radial factors per direction, frequency and node x ∈ [−1, 1]:
      !   q = q_max (x + 1)/2,   a = q̂ᵀ S q̂ · ρ(q, iω),   radial_w = (q_max/2) w_GL · q / (1 + 2πq a)
      !   coef_h = w(q̂) A/(2π)² · radial_w · 4π² a,   coef_w = w(q̂) A/(2π)² · radial_w · 2π
      ALLOCATE (q_node(n_radial, n_dir), coef_h(n_radial, n_t, n_dir), coef_w(n_radial, n_t, n_dir))
      DO i_dir = 1, n_dir
         q_max = dir_radius(i_dir)
         jac = dir_weight(i_dir)/area_bz
         DO i_node = 1, n_radial
            q_node(i_node, i_dir) = 0.5_dp*q_max*(nodes(i_node) + 1.0_dp)
         END DO
         DO i_w = 1, n_t
            a_scalar = DOT_PRODUCT(dir_hat(:, i_dir), &
                                   MATMUL(gc%S_head(:, :, i_w), dir_hat(:, i_dir)))
            DO i_node = 1, n_radial
               q = q_node(i_node, i_dir)
               a_q = a_scalar*gamma_corr_head_ratio(gc, q, i_w)
               denom = 1.0_dp + twopi*q*a_q
               radial_w = 0.5_dp*q_max*rule(i_node)*q/denom
               coef_h(i_node, i_w, i_dir) = jac*radial_w*4.0_dp*pi**2*a_q
               coef_w(i_node, i_w, i_dir) = jac*radial_w*twopi
            END DO
         END DO
      END DO

      ! Table points evenly spaced in fractional coordinates over [−1/2, 1/2]² (both ends included),
      ! round-robin over the ranks, then over the threads.
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          PRIVATE(i_tab, i_u, i_v, R_cart, i_dir, proj, i_w, i_node, q, term_h, term_w, &
!$OMP                  cos_q, sin_q, kh, kh1, kh2, kw, kw1) &
!$OMP          SHARED(gc, bs_env, n_tab, n_dir, n_t, dir_hat, q_node, coef_h, coef_w)
      ALLOCATE (cos_q(n_radial), sin_q(n_radial), kh(n_t), kh1(n_t), kh2(n_t), kw(3, n_t), &
                kw1(3, n_t))
!$OMP DO SCHEDULE(dynamic, 16)
      DO i_tab = bs_env%para_env%mepos + 1, n_tab*n_tab, bs_env%para_env%num_pe
         ! i_tab − 1 = i_u + n_tab·i_v, with i_u, i_v ∈ {0, …, n_tab − 1}.
         i_u = MODULO(i_tab - 1, n_tab)
         i_v = (i_tab - 1)/n_tab
         ! R = (i_u/(n_tab−1) − 1/2) a_1 + (i_v/(n_tab−1) − 1/2) a_2, Cartesian, bohr.
         R_cart = bs_env%hmat(:, 1)*(REAL(i_u, dp)/REAL(n_tab - 1, dp) - 0.5_dp) + &
                  bs_env%hmat(:, 2)*(REAL(i_v, dp)/REAL(n_tab - 1, dp) - 0.5_dp)
         kh(:) = 0.0_dp
         kh1(:) = 0.0_dp
         kh2(:) = 0.0_dp
         kw(:, :) = 0.0_dp
         kw1(:, :) = 0.0_dp

         DO i_dir = 1, n_dir
            ! projection t = q̂·R and the phases cos(q t), sin(q t) at the radial nodes
            proj = DOT_PRODUCT(dir_hat(:, i_dir), R_cart)
            DO i_node = 1, n_radial
               cos_q(i_node) = COS(q_node(i_node, i_dir)*proj)
               sin_q(i_node) = SIN(q_node(i_node, i_dir)*proj)
            END DO
            DO i_w = 1, n_t
               DO i_node = 1, n_radial
                  q = q_node(i_node, i_dir)
                  ! K_h += w(q̂) A/(2π)² · dq q · [−4π² a/(1 + 2πq a)] cos(q q̂·R)
                  term_h = -coef_h(i_node, i_w, i_dir)*cos_q(i_node)
                  kh(i_w) = kh(i_w) + term_h
                  ! K_w,α += w(q̂) A/(2π)² · dq q · [2π/(1 + 2πq a)] sin(q q̂·R) q̂_α
                  term_w = coef_w(i_node, i_w, i_dir)*sin_q(i_node)
                  kw(:, i_w) = kw(:, i_w) + term_w*dir_hat(:, i_dir)
                  ! K_h1, K_h2, K_w1,α: the integrands times q, q², q (expansion of e^{−q|z−z_s|})
                  kh1(i_w) = kh1(i_w) + term_h*q
                  kh2(i_w) = kh2(i_w) + term_h*q*q
                  kw1(:, i_w) = kw1(:, i_w) + term_w*q*dir_hat(:, i_dir)
               END DO
            END DO
         END DO
         gc%K_head(i_tab, :) = kh(:)
         gc%K_head1(i_tab, :) = kh1(:)
         gc%K_head2(i_tab, :) = kh2(:)
         gc%K_wing(i_tab, :, :) = kw(:, :)
         gc%K_wing1(i_tab, :, :) = kw1(:, :)
      END DO
!$OMP END DO
      DEALLOCATE (cos_q, sin_q, kh, kh1, kh2, kw, kw1)
!$OMP END PARALLEL
      DEALLOCATE (q_node, coef_h, coef_w)

      ! Every table point was computed by exactly one rank (zeros elsewhere): sum = gather.
      CALL bs_env%para_env%sum(gc%K_head)
      CALL bs_env%para_env%sum(gc%K_wing)
      CALL bs_env%para_env%sum(gc%K_head1)
      CALL bs_env%para_env%sum(gc%K_head2)
      CALL bs_env%para_env%sum(gc%K_wing1)

      ! remainder ∫_all − ∫_cell − (1/A) Σ_{g≠0} as separate tables, per class of atoms
      CALL gamma_corr_tabulate_remainder_2d(bs_env, qs_env, n_tab_rem)

      ! the tables are consumed by gamma_corr_add_folded_to_W; the head tensor is not needed any more
      DEALLOCATE (gc%S_head)

      DEALLOCATE (dir_hat, dir_weight, dir_radius, nodes, rule)

      gc%wall_time = gc%wall_time + m_walltime() - start_time

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_tabulate

! **************************************************************************************************
!> \brief Slab: tabulates the remainder of the head and wing model outside the Brillouin-zone cell
!>        for the form-factor classes A, B (form factors φ_A from gamma_corr_form_factors):
!>
!>          K^rem_AB(R, iω) = A [ ∫_all d²q/(2π)² f_h φ_A φ_B cos(q·R) − ∫_cell (…)
!>                                − (1/A) Σ_{g≠0} f_h(g) φ_A(g) φ_B(g) cos(g·R) ] ,
!>          K^rem_A,α(R, iω) = A [ same with f_w,α φ_A sin(q·R) ] ,
!>          f_h = −4π² a/(1 + 2πq a),   f_w,α = 2π q̂_α/(1 + 2πq a),   a = q̂ᵀ S q̂ · ρ(q, iω).
!>
!>        All terms are damped by e^{−σ²q²}, then extrapolated quadratically in σ² to σ = 0.
!> \param bs_env band-structure environment: gc%S_head, gc%measure, hmat; sets gc%K_rem_head/wing
!> \param qs_env QS environment, for the atoms and the RI_AUX basis of the form factors
!> \param n_tab points per direction of the tabulation mesh
! **************************************************************************************************
   SUBROUTINE gamma_corr_tabulate_remainder_2d(bs_env, qs_env, n_tab)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: n_tab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_tabulate_remainder_2d'
      INTEGER, PARAMETER :: n_angle = 512, n_chunk = 32, n_face = 24, n_grid_phi = 4096, &
                            n_radial = 192, n_radial_cell = 32, n_sigma = 3
      REAL(KIND=dp), DIMENSION(n_sigma), PARAMETER :: &
         sigma_frac = [1.0_dp/16.0_dp, 1.0_dp/20.0_dp, 1.0_dp/24.0_dp]
      REAL(KIND=dp), PARAMETER                           :: q_cut_sigma = 6.0_dp

      INTEGER :: i_item, last_pass, n_chunks
      INTEGER :: cls_a, cls_b, handle, i_a, i_c, i_dir, i_g, i_pair, i_pass, i_q, i_r, i_s, &
                 i_tab, i_u, i_v, i_w, j_s, m1, m2, m_max, n_c, n_cell, n_dir, n_g, n_half, n_pair, n_pass, &
                 n_q, n_t, n_w
      REAL(KIND=dp)                                      :: a_q, area, b_len, denom, dq_phi, &
                                                            lagrange, q, q_cut, q_cut_max, r_min, &
                                                            theta, w_ang, w_node
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dir_radius, dir_weight, extrap, nodes, &
                                                            phi_a, phi_b, q_abs, rule
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: base_h, base_wx, base_wy, cos_qr, &
                                                            dir_hat, phi, q_vec, R_chunk, res, &
                                                            sin_qr, wgt
      REAL(KIND=dp), DIMENSION(2)                        :: b1, b2, g2, q_hat
      REAL(KIND=dp), DIMENSION(3)                        :: R_cart
      REAL(KIND=dp), DIMENSION(n_sigma)                  :: s2
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)
      gc => bs_env%ri_rs%gamma_corr
      n_t = bs_env%num_time_freq_points
      area = gc%measure
      gc%n_tab_rem = n_tab
      ! one weight column per (frequency, σ)
      n_w = n_t*n_sigma

      ! shortest in-plane lattice vector and the in-plane reciprocal vectors b_i·a_j = 2π δ_ij
      r_min = HUGE(1.0_dp)
      DO m1 = -2, 2
         DO m2 = -2, 2
            IF (m1 == 0 .AND. m2 == 0) CYCLE
            R_cart = REAL(m1, dp)*bs_env%hmat(:, 1) + REAL(m2, dp)*bs_env%hmat(:, 2)
            r_min = MIN(r_min, NORM2(R_cart(1:2)))
         END DO
      END DO
      CALL gamma_corr_reciprocal_2d(bs_env%hmat, b1, b2)
      DO i_s = 1, n_sigma
         s2(i_s) = (sigma_frac(i_s)*r_min)**2
      END DO
      q_cut_max = q_cut_sigma/(sigma_frac(n_sigma)*r_min)

      ! form-factor classes and φ_A(q) on a uniform grid [0, q_cut_max]
      dq_phi = q_cut_max/REAL(n_grid_phi, dp)
      CALL gamma_corr_form_factors(qs_env, bs_env, n_grid_phi, dq_phi, phi)

      ! node set: (a) ∫_all, (b) −∫_cell, (c) −(1/A) Σ_{g≠0}; vectors in q_vec, the base weights
      ! (quadrature weight · damping · f, without form factors) per (frequency, σ) column
      b_len = MIN(NORM2(b1), NORM2(b2))
      m_max = CEILING(2.0_dp*q_cut_max/b_len) + 1
      n_g = 0
      DO m1 = -m_max, m_max
         DO m2 = -m_max, m_max
            IF (m1 == 0 .AND. m2 == 0) CYCLE
            IF (NORM2(REAL(m1, dp)*b1 + REAL(m2, dp)*b2) < q_cut_max) n_g = n_g + 1
         END DO
      END DO
      CALL gamma_corr_ws_cell(bs_env, .TRUE., n_face, dir_hat, dir_weight, dir_radius)
      n_dir = SIZE(dir_weight)
      n_cell = n_dir*n_radial_cell
      n_q = (n_angle/2)*n_radial*n_sigma + n_cell*n_sigma + n_g
      ALLOCATE (q_vec(2, n_q), q_abs(n_q), base_h(n_q, n_w), base_wx(n_q, n_w), base_wy(n_q, n_w))
      base_h(:, :) = 0.0_dp
      base_wx(:, :) = 0.0_dp
      base_wy(:, :) = 0.0_dp
      i_q = 0

      ! (a) ∫_all d²q/(2π)²: half circle × 2, weight (1/(2π)²)(2π/n_angle)·2·w_GL q, only in the
      !     columns of its own σ
      ALLOCATE (nodes(n_radial), rule(n_radial))
      CALL gauss_legendre(n_radial, nodes, rule)
      w_ang = 2.0_dp*(twopi/REAL(n_angle, dp))/twopi**2
      DO i_s = 1, n_sigma
         q_cut = q_cut_sigma/SQRT(s2(i_s))
         DO i_a = 1, n_angle/2
            theta = twopi*REAL(i_a - 1, dp)/REAL(n_angle, dp)
            q_hat = [COS(theta), SIN(theta)]
            DO i_r = 1, n_radial
               i_q = i_q + 1
               q = 0.5_dp*q_cut*(nodes(i_r) + 1.0_dp)
               q_vec(:, i_q) = q*q_hat
               w_node = w_ang*0.5_dp*q_cut*rule(i_r)*q*EXP(-s2(i_s)*q*q)
               CALL node_weights(i_q, q, q_hat, w_node, i_s)
            END DO
         END DO
      END DO
      DEALLOCATE (nodes, rule)

      ! (b) −∫_cell d²q/(2π)²: the cell quadrature of gamma_corr_tabulate,
      !     weight −w(q̂)/(2π)² (q_max/2) w_GL q
      ALLOCATE (nodes(n_radial_cell), rule(n_radial_cell))
      CALL gauss_legendre(n_radial_cell, nodes, rule)
      DO i_s = 1, n_sigma
         DO i_dir = 1, n_dir
            q_hat = dir_hat(1:2, i_dir)
            DO i_r = 1, n_radial_cell
               i_q = i_q + 1
               q = 0.5_dp*dir_radius(i_dir)*(nodes(i_r) + 1.0_dp)
               q_vec(:, i_q) = q*q_hat
               w_node = -dir_weight(i_dir)/twopi**2*0.5_dp*dir_radius(i_dir)*rule(i_r)*q* &
                        EXP(-s2(i_s)*q*q)
               CALL node_weights(i_q, q, q_hat, w_node, i_s)
            END DO
         END DO
      END DO
      DEALLOCATE (nodes, rule, dir_hat, dir_weight, dir_radius)

      ! (c) −(1/A) Σ_{0<|g|<6/σ_min} f(g) e^{−σ²g²}, in every σ column
      i_g = 0
      DO m1 = -m_max, m_max
         DO m2 = -m_max, m_max
            IF (m1 == 0 .AND. m2 == 0) CYCLE
            g2 = REAL(m1, dp)*b1 + REAL(m2, dp)*b2
            q = NORM2(g2)
            IF (q >= q_cut_max) CYCLE
            i_q = i_q + 1
            q_vec(:, i_q) = g2
            DO i_s = 1, n_sigma
               CALL node_weights(i_q, q, g2/q, -EXP(-s2(i_s)*q*q)/area, i_s)
            END DO
         END DO
      END DO
      CPASSERT(i_q == n_q)
      DO i_q = 1, n_q
         q_abs(i_q) = NORM2(q_vec(:, i_q))
      END DO

      ! output tables: head per class pair, wing per class
      n_pair = gc%n_class*(gc%n_class + 1)/2
      IF (ALLOCATED(gc%K_rem_head)) DEALLOCATE (gc%K_rem_head)
      IF (ALLOCATED(gc%K_rem_wing)) DEALLOCATE (gc%K_rem_wing)
      ALLOCATE (gc%K_rem_head(n_tab*n_tab, n_t, n_pair), gc%K_rem_wing(n_tab*n_tab, 2, n_t, gc%n_class))
      gc%K_rem_head(:, :, :) = 0.0_dp
      gc%K_rem_wing(:, :, :, :) = 0.0_dp

      ! σ² → 0 extrapolation factors, times A
      ALLOCATE (extrap(n_sigma))
      DO i_s = 1, n_sigma
         lagrange = 1.0_dp
         DO j_s = 1, n_sigma
            IF (j_s /= i_s) lagrange = lagrange*(0.0_dp - s2(j_s))/(s2(i_s) - s2(j_s))
         END DO
         extrap(i_s) = lagrange*area
      END DO

      ! first half of the table points (the rest follows by R → −R symmetry)
      n_half = (n_tab*n_tab + 1)/2
      ! Work items (pass, chunk of table points), round-robin over the ranks. Passes: every class
      ! pair (A ≤ B) for the head, then every class for the wing x and y.
      ALLOCATE (wgt(n_q, n_w), res(n_chunk, n_w), cos_qr(n_chunk, n_q), sin_qr(n_chunk, n_q), &
                R_chunk(2, n_chunk))
      ALLOCATE (phi_a(n_q), phi_b(n_q))
      n_pass = n_pair + 2*gc%n_class
      n_chunks = (n_half + n_chunk - 1)/n_chunk
      last_pass = 0
      DO i_item = bs_env%para_env%mepos + 1, n_pass*n_chunks, bs_env%para_env%num_pe
         i_pass = (i_item - 1)/n_chunks + 1
         i_tab = MODULO(i_item - 1, n_chunks)*n_chunk + 1
         IF (i_pass /= last_pass) THEN
            last_pass = i_pass
            ! weights of this pass: base × φ_A φ_B (head) or base × φ_A (wing)
            IF (i_pass <= n_pair) THEN
               ! pair index i_pair = A + B(B−1)/2 with A ≤ B
               i_pair = i_pass
               cls_b = CEILING((SQRT(8.0_dp*REAL(i_pair, dp) + 1.0_dp) - 1.0_dp)/2.0_dp)
               cls_a = i_pair - cls_b*(cls_b - 1)/2
               CALL interpolate_phi(cls_a, phi_a)
               CALL interpolate_phi(cls_b, phi_b)
               DO i_c = 1, n_w
                  wgt(:, i_c) = base_h(:, i_c)*phi_a(:)*phi_b(:)
               END DO
            ELSE
               cls_a = (i_pass - n_pair + 1)/2
               CALL interpolate_phi(cls_a, phi_a)
               IF (MODULO(i_pass - n_pair, 2) == 1) THEN
                  DO i_c = 1, n_w
                     wgt(:, i_c) = base_wx(:, i_c)*phi_a(:)
                  END DO
               ELSE
                  DO i_c = 1, n_w
                     wgt(:, i_c) = base_wy(:, i_c)*phi_a(:)
                  END DO
               END IF
            END IF
         END IF
         ! table points of the chunk: cos/sin(q·R) (chunk x n_q) times the weights (n_q x n_w)
         n_c = MIN(n_chunk, n_half - i_tab + 1)
         DO i_c = 1, n_c
            i_u = MODULO(i_tab + i_c - 2, n_tab)
            i_v = (i_tab + i_c - 2)/n_tab
            R_cart = bs_env%hmat(:, 1)*(REAL(i_u, dp)/REAL(n_tab - 1, dp) - 0.5_dp) + &
                     bs_env%hmat(:, 2)*(REAL(i_v, dp)/REAL(n_tab - 1, dp) - 0.5_dp)
            R_chunk(:, i_c) = R_cart(1:2)
         END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(static) PRIVATE(i_q, i_c, theta) &
!$OMP             SHARED(n_q, n_c, q_vec, R_chunk, cos_qr, sin_qr)
         DO i_q = 1, n_q
            DO i_c = 1, n_c
               theta = q_vec(1, i_q)*R_chunk(1, i_c) + q_vec(2, i_q)*R_chunk(2, i_c)
               cos_qr(i_c, i_q) = COS(theta)
               sin_qr(i_c, i_q) = SIN(theta)
            END DO
         END DO
!$OMP END PARALLEL DO
         IF (i_pass <= n_pair) THEN
            CALL dgemm('N', 'N', n_c, n_w, n_q, 1.0_dp, cos_qr, n_chunk, wgt, n_q, 0.0_dp, res, n_chunk)
         ELSE
            CALL dgemm('N', 'N', n_c, n_w, n_q, 1.0_dp, sin_qr, n_chunk, wgt, n_q, 0.0_dp, res, n_chunk)
         END IF
         ! D(σ) = D(0) + c₁σ² + c₂σ⁴ → D(0), times A
         DO i_s = 1, n_sigma
            DO i_w = 1, n_t
               i_c = i_w + (i_s - 1)*n_t
               IF (i_pass <= n_pair) THEN
                  gc%K_rem_head(i_tab:i_tab + n_c - 1, i_w, i_pair) = &
                     gc%K_rem_head(i_tab:i_tab + n_c - 1, i_w, i_pair) + extrap(i_s)*res(1:n_c, i_c)
               ELSE IF (MODULO(i_pass - n_pair, 2) == 1) THEN
                  gc%K_rem_wing(i_tab:i_tab + n_c - 1, 1, i_w, cls_a) = &
                     gc%K_rem_wing(i_tab:i_tab + n_c - 1, 1, i_w, cls_a) + extrap(i_s)*res(1:n_c, i_c)
               ELSE
                  gc%K_rem_wing(i_tab:i_tab + n_c - 1, 2, i_w, cls_a) = &
                     gc%K_rem_wing(i_tab:i_tab + n_c - 1, 2, i_w, cls_a) + extrap(i_s)*res(1:n_c, i_c)
               END IF
            END DO
         END DO
      END DO
      CALL bs_env%para_env%sum(gc%K_rem_head)
      CALL bs_env%para_env%sum(gc%K_rem_wing)
      ! mirror of the second half: K(−R) = K(R) for the head, K(−R) = −K(R) for the wing
      DO i_tab = n_half + 1, n_tab*n_tab
         gc%K_rem_head(i_tab, :, :) = gc%K_rem_head(n_tab*n_tab + 1 - i_tab, :, :)
         gc%K_rem_wing(i_tab, :, :, :) = -gc%K_rem_wing(n_tab*n_tab + 1 - i_tab, :, :, :)
      END DO
      DEALLOCATE (wgt, res, cos_qr, sin_qr, R_chunk, phi_a, phi_b, base_h, base_wx, base_wy, q_vec, q_abs, &
                  extrap, phi)

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief Base weights of node i_q in the σ columns i_s: w · f_h and w · f_w,α with
!>        f_h = −4π² a/(1 + 2πq a), f_w,α = 2π q̂_α/(1 + 2πq a), a = q̂ᵀSq̂ ρ(q, iω).
!> \param i_q node index
!> \param qq |q|
!> \param qh unit vector q̂ (in-plane)
!> \param w quadrature weight times damping (sign included)
!> \param is σ index
! **************************************************************************************************
      SUBROUTINE node_weights(i_q, qq, qh, w, is)
      INTEGER, INTENT(IN)                                :: i_q
      REAL(KIND=dp), INTENT(IN)                          :: qq
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: qh
      REAL(KIND=dp), INTENT(IN)                          :: w
      INTEGER, INTENT(IN)                                :: is

      INTEGER                                            :: ic, iw

         DO iw = 1, n_t
            a_q = DOT_PRODUCT(qh, MATMUL(gc%S_head(1:2, 1:2, iw), qh))*gamma_corr_head_ratio(gc, qq, iw)
            denom = 1.0_dp + twopi*qq*a_q
            ic = iw + (is - 1)*n_t
            base_h(i_q, ic) = w*(-4.0_dp*pi**2*a_q/denom)
            base_wx(i_q, ic) = w*twopi*qh(1)/denom
            base_wy(i_q, ic) = w*twopi*qh(2)/denom
         END DO
      END SUBROUTINE node_weights

! **************************************************************************************************
!> \brief φ_class(|q|) at every node by linear interpolation on the uniform grid.
!> \param cls form-factor class
!> \param out φ at the nodes
! **************************************************************************************************
      SUBROUTINE interpolate_phi(cls, out)
      INTEGER, INTENT(IN)                                :: cls
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: out

      INTEGER                                            :: ig, iq
      REAL(KIND=dp)                                      :: x

         DO iq = 1, n_q
            x = q_abs(iq)/dq_phi
            ig = MIN(INT(x), n_grid_phi - 1)
            x = x - REAL(ig, dp)
            out(iq) = (1.0_dp - x)*phi(ig, cls) + x*phi(ig + 1, cls)
         END DO
      END SUBROUTINE interpolate_phi

   END SUBROUTINE gamma_corr_tabulate_remainder_2d

! **************************************************************************************************
!> \brief Form-factor classes of the atoms of a slab: per kind, the heights d = Z_A − z_s above the
!>        sheet plane z_s (mean atom height) are split into planes at gaps > class_plane_gap and
!>        binned around the plane mean d̄: class = (kind, plane, k), k = NINT((d − d̄)/w),
!>        w = class_bin_width, class height d_class = d̄ + k w. The first atom of a class represents it.
!> \param qs_env QS environment, for the particles and the cell
!> \param n_class number of classes
!> \param atom_class class of every atom (n_atom)
!> \param d_class height d of every class (first n_class entries used)
!> \param rep_atom representative atom of every class (first n_class entries used)
! **************************************************************************************************
   SUBROUTINE gamma_corr_classify_atoms(qs_env, n_class, atom_class, d_class, rep_atom)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: n_class
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: atom_class
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: d_class
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: rep_atom

      INTEGER                                            :: i_cls, i_first, i_kind, i_s, iatom, j_s, &
                                                            n_atom, n_k, n_kind, n_plane
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of, bin_of, kind_class, order, &
                                                            plane_class, plane_of
      LOGICAL                                            :: found, new_plane
      REAL(KIND=dp)                                      :: d, l_z, z_sheet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: d_sorted, plane_mean, z_atom
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles

      CALL get_qs_env(qs_env, particle_set=particles, cell=cell)
      n_atom = SIZE(particles)
      IF (ALLOCATED(atom_class)) DEALLOCATE (atom_class)
      ALLOCATE (atom_class(n_atom), z_atom(n_atom), kind_class(n_atom), d_class(n_atom), &
                rep_atom(n_atom), plane_of(n_atom), plane_mean(n_atom), bin_of(n_atom), &
                order(n_atom), d_sorted(n_atom), plane_class(n_atom), atom_of(n_atom))

      ! atom heights, unwrapped along z relative to atom 1, and the sheet plane z_s
      l_z = cell%hmat(3, 3)
      DO iatom = 1, n_atom
         d = particles(iatom)%r(3) - particles(1)%r(3)
         IF (l_z > 0.0_dp) d = d - l_z*ANINT(d/l_z)
         z_atom(iatom) = particles(1)%r(3) + d
      END DO
      z_sheet = SUM(z_atom)/REAL(n_atom, dp)

      ! planes per kind: sorted heights split at gaps > class_plane_gap; d̄ = d_1 + mean(d − d_1)
      n_plane = 0
      n_kind = 0
      DO iatom = 1, n_atom
         n_kind = MAX(n_kind, particles(iatom)%atomic_kind%kind_number)
      END DO
      DO i_kind = 1, n_kind
         n_k = 0
         DO iatom = 1, n_atom
            IF (particles(iatom)%atomic_kind%kind_number /= i_kind) CYCLE
            n_k = n_k + 1
            d_sorted(n_k) = z_atom(iatom) - z_sheet
            atom_of(n_k) = iatom
         END DO
         IF (n_k == 0) CYCLE
         ! sort overwrites order with the original positions: sorted entry i is atom atom_of(order(i))
         CALL sort(d_sorted, n_k, order)
         i_first = 1
         DO i_s = 1, n_k
            IF (i_s == n_k) THEN
               new_plane = .TRUE.
            ELSE
               new_plane = d_sorted(i_s + 1) - d_sorted(i_s) > class_plane_gap
            END IF
            IF (.NOT. new_plane) CYCLE
            ! atoms i_first .. i_s of the sorted list form one plane
            n_plane = n_plane + 1
            plane_mean(n_plane) = d_sorted(i_first) + &
                                  SUM(d_sorted(i_first:i_s) - d_sorted(i_first))/REAL(i_s - i_first + 1, dp)
            DO j_s = i_first, i_s
               plane_of(atom_of(order(j_s))) = n_plane
               bin_of(atom_of(order(j_s))) = NINT((d_sorted(j_s) - plane_mean(n_plane))/class_bin_width)
            END DO
            i_first = i_s + 1
         END DO
      END DO

      ! classes (plane, bin), numbered in the order of their first atom
      n_class = 0
      DO iatom = 1, n_atom
         found = .FALSE.
         DO i_cls = 1, n_class
            IF (plane_class(i_cls) == plane_of(iatom) .AND. kind_class(i_cls) == bin_of(iatom)) THEN
               atom_class(iatom) = i_cls
               found = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. found) THEN
            n_class = n_class + 1
            ! kind_class holds the bin index here (the plane already fixes the kind)
            plane_class(n_class) = plane_of(iatom)
            kind_class(n_class) = bin_of(iatom)
            d_class(n_class) = plane_mean(plane_of(iatom)) + REAL(bin_of(iatom), dp)*class_bin_width
            rep_atom(n_class) = iatom
            atom_class(iatom) = n_class
         END IF
      END DO
      DEALLOCATE (plane_of, plane_mean, bin_of, order, d_sorted, plane_class, atom_of)
      DEALLOCATE (z_atom, kind_class)

   END SUBROUTINE gamma_corr_classify_atoms

! **************************************************************************************************
!> \brief Out-of-plane form factors of the auxiliary functions per class of atoms
!>        (gamma_corr_classify_atoms) on a uniform q grid, at the class height d:
!>
!>          φ(q) = Σ_P ∫dr ϕ_P(r) e^{−q|z−z_s|} / Σ_P ∫dr ϕ_P(r) ,
!>
!>        where a Cartesian primitive x^{n_x} y^{n_y} u^n e^{−αr²} (u = z − Z_A, n even) adds
!>        μ_{n_x} μ_{n_y} G_n(d, q) to the numerator, G_n(d, q) = ∫du u^n e^{−αu²} e^{−q|u+d|}.
!> \param qs_env QS environment, for the particles, the cell and the RI_AUX basis
!> \param bs_env band-structure environment; sets gc%n_class and gc%atom_class
!> \param n_grid number of grid intervals
!> \param dq grid spacing (1/bohr)
!> \param phi φ(0:n_grid, n_class), allocated here
! **************************************************************************************************
   SUBROUTINE gamma_corr_form_factors(qs_env, bs_env, n_grid, dq, phi)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: n_grid
      REAL(KIND=dp), INTENT(IN)                          :: dq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: phi

      INTEGER, PARAMETER                                 :: n_gl = 64
      REAL(KIND=dp), PARAMETER                           :: gauss_exp_cut = 40.0_dp

      INTEGER                                            :: i_cls, i_gl, i_piece, i_q, ico, il, &
                                                            ipgf, iset, isgf, jatom, lmax, ncart, &
                                                            nz, offset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rep_atom
      REAL(KIND=dp)                                      :: alpha, c_sum, d, hi, lo, u, w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c_prim, d_class, gl_w, gl_x, mom, num
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: g_prim
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds

      gc => bs_env%ri_rs%gamma_corr
      CALL get_qs_env(qs_env, particle_set=particles, qs_kind_set=qs_kinds)
      ! classes (kind, height d above the sheet plane) and their representative atoms
      CALL gamma_corr_classify_atoms(qs_env, gc%n_class, gc%atom_class, d_class, rep_atom)

      ALLOCATE (gl_x(n_gl), gl_w(n_gl))
      CALL gauss_legendre(n_gl, gl_x, gl_w)
      ALLOCATE (phi(0:n_grid, gc%n_class))
      ! Classes are independent and each fills its own column of phi: one class per thread.
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          PRIVATE(i_cls, jatom, d, basis, num, c_sum, iset, lmax, ncart, mom, c_prim, g_prim, &
!$OMP                  ipgf, alpha, il, ico, nz, i_piece, lo, hi, i_gl, u, w, i_q, offset, isgf) &
!$OMP          SHARED(gc, rep_atom, d_class, qs_kinds, particles, n_grid, dq, gl_x, gl_w, phi, &
!$OMP                 ncoset, indco)
      ALLOCATE (num(0:n_grid))
!$OMP DO SCHEDULE(dynamic)
      DO i_cls = 1, gc%n_class
         jatom = rep_atom(i_cls)
         d = d_class(i_cls)
         CALL get_qs_kind(qs_kinds(particles(jatom)%atomic_kind%kind_number), &
                          basis_set=basis, basis_type='RI_AUX')
         num(:) = 0.0_dp
         c_sum = 0.0_dp
         DO iset = 1, basis%nset
            lmax = basis%lmax(iset)
            ncart = ncoset(lmax)
            ALLOCATE (mom(0:lmax), c_prim(ncart*basis%npgf(iset)), &
                      g_prim(ncart*basis%npgf(iset), 0:n_grid))
            c_prim(:) = 0.0_dp
            g_prim(:, :) = 0.0_dp
            DO ipgf = 1, basis%npgf(iset)
               alpha = basis%zet(ipgf, iset)
               ! μ_k = ∫du u^k e^{−αu²}
               mom(:) = 0.0_dp
               mom(0) = SQRT(pi/alpha)
               DO il = 2, lmax, 2
                  mom(il) = REAL(il - 1, dp)*mom(il - 2)/(2.0_dp*alpha)
               END DO
               DO ico = ncoset(basis%lmin(iset) - 1) + 1, ncart
                  nz = indco(3, ico)
                  ! only the charge-carrying components (even n_z) enter the class average
                  IF (MODULO(nz, 2) == 1) CYCLE
                  ! charge μ_{n_x} μ_{n_y} μ_{n_z}
                  c_prim((ipgf - 1)*ncart + ico) = mom(indco(1, ico))*mom(indco(2, ico))*mom(nz)
                  ! z factor G_n(d, q) on the grid, two Gauss-Legendre pieces split at u = −d
                  DO i_piece = 1, 2
                     IF (i_piece == 1) THEN
                        lo = -SQRT(gauss_exp_cut/alpha)
                        hi = MIN(-d, SQRT(gauss_exp_cut/alpha))
                     ELSE
                        lo = MAX(-d, -SQRT(gauss_exp_cut/alpha))
                        hi = SQRT(gauss_exp_cut/alpha)
                     END IF
                     IF (hi <= lo) CYCLE
                     DO i_gl = 1, n_gl
                        u = 0.5_dp*(hi - lo)*gl_x(i_gl) + 0.5_dp*(hi + lo)
                        w = 0.5_dp*(hi - lo)*gl_w(i_gl)*u**nz*EXP(-alpha*u*u)* &
                            mom(indco(1, ico))*mom(indco(2, ico))
                        DO i_q = 0, n_grid
                           g_prim((ipgf - 1)*ncart + ico, i_q) = g_prim((ipgf - 1)*ncart + ico, i_q) + &
                                                                 w*EXP(-REAL(i_q, dp)*dq*ABS(u + d))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            ! contract to the spherical functions of the set and sum over them
            offset = basis%first_sgf(1, iset)
            DO isgf = offset, offset + basis%nsgf_set(iset) - 1
               c_sum = c_sum + SUM(basis%sphi(1:SIZE(c_prim), isgf)*c_prim)
               DO i_q = 0, n_grid
                  num(i_q) = num(i_q) + SUM(basis%sphi(1:SIZE(c_prim), isgf)*g_prim(:, i_q))
               END DO
            END DO
            DEALLOCATE (mom, c_prim, g_prim)
         END DO
         ! φ(q) = Σ_P ρ̂_P(q)/Σ_P c_P
         phi(:, i_cls) = num(:)/c_sum
      END DO
!$OMP END DO
      DEALLOCATE (num)
!$OMP END PARALLEL

      DEALLOCATE (d_class, rep_atom, gl_x, gl_w)

   END SUBROUTINE gamma_corr_form_factors

! **************************************************************************************************
!> \brief In-plane reciprocal vectors b_1, b_2 of the cell columns a_1, a_2: b_i·a_j = 2π δ_ij.
!> \param hmat cell matrix, columns a_i (bohr)
!> \param b1 first reciprocal vector (x, y)
!> \param b2 second reciprocal vector (x, y)
! **************************************************************************************************
   PURE SUBROUTINE gamma_corr_reciprocal_2d(hmat, b1, b2)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: b1, b2

      REAL(KIND=dp)                                      :: det

      ! [b1 b2]ᵀ [a1 a2] = 2π 1 for the 2x2 in-plane block
      det = hmat(1, 1)*hmat(2, 2) - hmat(1, 2)*hmat(2, 1)
      b1 = twopi/det*[hmat(2, 2), -hmat(1, 2)]
      b2 = twopi/det*[-hmat(2, 1), hmat(1, 1)]
   END SUBROUTINE gamma_corr_reciprocal_2d

! **************************************************************************************************
!> \brief Adds the Γ-point correction to the folded W, fm_W_folded_time(i_t) = M⁺ W(iτ_t) M⁺.
!>        3D: W_f(iτ_t) += [Σ_k γ̃^cos_tk (1/ε_0(iω_k) − 1)] M⁺ K M⁺ (K: gamma_corr_tail_operator).
!>        2D: with ∫dr ϕ_P e^{−q|z−z_s|} ≈ c_P − q m_P + (q²/2) s_P (gamma_corr_z_moments) and all
!>        kernels at the in-plane minimum-image displacement R_PQ of the atoms of P, Q (classes A, B):
!>
!>          ΔW_PQ(iω) = (1/A) [ c_P c_Q (K_h + K^rem_h,AB) − (m_P c_Q + c_P m_Q) K_h1
!>                              + (m_P m_Q + (s_P c_Q + c_P s_Q)/2) K_h2 ]
!>                    + (1/√A) Σ_α [ T_Pα (c_Q (K_w,α + K^rem_w,α,B) − m_Q K_w1,α)
!>                                   − (c_P (K_w,α + K^rem_w,α,A) − m_P K_w1,α) T_Qα ] ,
!>
!>        written as ΔW = E K Eᵀ with per-atom coefficient columns E = [E_a | E_T(iω)] and folded as
!>        M⁺_c ΔW M⁺_c = F K Fᵀ, F = M⁺_c E, M⁺_c = M⁺ on the modes above REGULARIZATION_RI_CORRECTION.
!> \param bs_env band-structure environment; uses and releases the correction data of gc
!> \param qs_env QS environment, for the atom positions and the z moments of the RI_AUX basis
!> \param fm_W_folded_time M⁺ W(iτ_t) M⁺ per time point; the folded correction is added
! **************************************************************************************************
   SUBROUTINE gamma_corr_add_folded_to_W(bs_env, qs_env, fm_W_folded_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: fm_W_folded_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_add_folded_to_W'

      INTEGER :: A_at, B_at, handle, i_atom, i_col, i_row, i_t, I_glob, J_glob, k, l, n_atom, n_ea, &
                 n_RI, n_t, ncol_local, nrow_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: inv_A, inv_sqrt_A, start_time, w_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c, m_z, s_z, v_ab, v_ba
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coef, frac_atom
      TYPE(cp_fm_struct_type), POINTER                   :: struct_aa, struct_F, struct_Ta
      TYPE(cp_fm_type)                                   :: fm_F_a, fm_F_T, fm_KT, fm_tmp
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_H, fm_K_aa, fm_K_Ta
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      start_time = m_walltime()

      gc => bs_env%ri_rs%gamma_corr
      CPASSERT(gc%ready)
      n_RI = bs_env%n_RI
      n_t = bs_env%num_time_freq_points

      ! ------------------------------------------------------------------------------------------
      ! 3D: W_f(iτ_t) += w_t M⁺ K M⁺,  w_t = Σ_k γ̃^cos_tk (1/ε_0(iω_k) − 1)
      ! ------------------------------------------------------------------------------------------
      IF (.NOT. gc%slab) THEN
         IF (gc%tail_3d) THEN
            CPASSERT(ALLOCATED(gc%tail_factor))
            CALL cp_fm_create(fm_tmp, bs_env%fm_struct_RI_RI)
            CALL cp_fm_create(fm_KT, bs_env%fm_struct_RI_RI)
            CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, gc%fm_K_tail, bs_env%fm_Minv_Gamma, &
                               0.0_dp, fm_tmp)
            CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, bs_env%fm_Minv_Gamma, fm_tmp, &
                               0.0_dp, fm_KT)
            CALL cp_fm_release(fm_tmp)
            CALL cp_fm_release(gc%fm_K_tail)
            DO i_t = 1, n_t
               w_t = SUM(bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, :)*gc%tail_factor)
               CALL cp_fm_scale_and_add(1.0_dp, fm_W_folded_time(i_t), w_t, fm_KT)
            END DO
            CALL cp_fm_release(fm_KT)
         END IF
         gc%wall_time = gc%wall_time + m_walltime() - start_time
         CALL timestop(handle)
         RETURN
      END IF

      ! ------------------------------------------------------------------------------------------
      ! 2D: 1. charges c_P = |c| ĉ_P, fractional atom positions, atom of every RI function and
      !        the out-of-plane moments m_P, s_P
      ! ------------------------------------------------------------------------------------------
      inv_A = 1.0_dp/gc%measure
      inv_sqrt_A = 1.0_dp/SQRT(gc%measure)
      ALLOCATE (c(n_RI))
      c(:) = gc%charge_norm*gc%charge
      NULLIFY (particle_set)
      CALL get_qs_env(qs_env, particle_set=particle_set)
      n_atom = SIZE(particle_set)
      ALLOCATE (frac_atom(2, n_atom), atom_of_RI(n_RI))
      DO i_atom = 1, n_atom
         frac_atom(1, i_atom) = DOT_PRODUCT(bs_env%ri_rs%cell%h_inv(1, :), particle_set(i_atom)%r)
         frac_atom(2, i_atom) = DOT_PRODUCT(bs_env%ri_rs%cell%h_inv(2, :), particle_set(i_atom)%r)
         atom_of_RI(bs_env%i_RI_start_from_atom(i_atom):bs_env%i_RI_end_from_atom(i_atom)) = i_atom
      END DO
      ALLOCATE (m_z(n_RI), s_z(n_RI))
      CALL gamma_corr_z_moments(qs_env, m_z, s_z)
      n_ea = 3*n_atom

      ! layouts: F (n_RI x n_ea, rows as M⁺), K_aa (n_ea x n_ea), K_Ta (n_ea x n_ea, rows (A,α))
      NULLIFY (struct_F, struct_aa, struct_Ta)
      CALL cp_fm_struct_create(struct_F, template_fmstruct=bs_env%fm_struct_RI_RI, ncol_global=n_ea)
      CALL cp_fm_struct_create(struct_aa, template_fmstruct=bs_env%fm_struct_RI_RI, &
                               nrow_global=n_ea, ncol_global=n_ea)
      struct_Ta => struct_aa

      ! 2. F_a = M⁺_c E_a, E_a(P, (A,k)) = {c_P, m_P, s_P}_k on the atom of P
      ALLOCATE (coef(n_RI, 3))
      coef(:, 1) = c(:)
      coef(:, 2) = m_z(:)
      coef(:, 3) = s_z(:)
      CALL cp_fm_create(fm_F_a, struct_F)
      CALL gamma_corr_minv_c_atom_columns(bs_env, coef, atom_of_RI, n_atom, fm_F_a)

      ! 3. K_aa(iτ_t) = Σ_j c_tj K_aa(iω_j) (c_tj = γ̃^cos_tj) and K_Ta(iω_j) on their local elements
      ALLOCATE (fm_K_aa(n_t), fm_K_Ta(n_t), v_ab(n_t), v_ba(n_t))
      DO i_t = 1, n_t
         CALL cp_fm_create(fm_K_aa(i_t), struct_aa)
         CALL cp_fm_create(fm_K_Ta(i_t), struct_Ta)
      END DO
      CALL cp_fm_get_info(fm_K_aa(1), nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic, 4) &
!$OMP          SHARED(ncol_local, nrow_local, row_indices, col_indices, fm_K_aa, fm_K_Ta, bs_env, gc, &
!$OMP                 frac_atom, n_t, inv_A, inv_sqrt_A) &
!$OMP          PRIVATE(i_col, i_row, I_glob, J_glob, A_at, B_at, k, l, v_ab, v_ba, i_t)
      DO i_col = 1, ncol_local
         J_glob = col_indices(i_col)
         B_at = (J_glob - 1)/3 + 1
         l = MOD(J_glob - 1, 3) + 1
         DO i_row = 1, nrow_local
            I_glob = row_indices(i_row)
            A_at = (I_glob - 1)/3 + 1
            k = MOD(I_glob - 1, 3) + 1
            ! head block, symmetrized over the two orderings of the atom pair
            CALL head_element(A_at, B_at, k, l, v_ab)
            CALL head_element(B_at, A_at, l, k, v_ba)
            DO i_t = 1, n_t
               fm_K_aa(i_t)%local_data(i_row, i_col) = 0.5_dp*inv_A* &
                                                      SUM(bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, :)*(v_ab + v_ba))
            END DO
            ! wing block: row (A,α = k), column (B,l); average of the Ta and the transposed aT form
            CALL wing_element(A_at, B_at, k, l, v_ab, v_ba)
            DO i_t = 1, n_t
               fm_K_Ta(i_t)%local_data(i_row, i_col) = 0.5_dp*inv_sqrt_A*(v_ab(i_t) + v_ba(i_t))
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      ! 4. H_t = Σ_j c_tj F_T(iω_j) K_Ta(iω_j)
      ALLOCATE (fm_H(n_t))
      DO i_t = 1, n_t
         CALL cp_fm_create(fm_H(i_t), struct_F, set_zero=.TRUE.)
      END DO
      CALL cp_fm_create(fm_F_T, struct_F)
      CALL cp_fm_create(fm_tmp, struct_F)
      DO k = 1, n_t
         ! F_T(iω_k) = M⁺_c E_T(iω_k), E_T(P, (A,α)) = T_Pα(iω_k) on the atom of P
         coef(:, 1:3) = gc%T_RI(:, 1:3, k)
         CALL gamma_corr_minv_c_atom_columns(bs_env, coef, atom_of_RI, n_atom, fm_F_T)
         CALL parallel_gemm('N', 'N', n_RI, n_ea, n_ea, 1.0_dp, fm_F_T, fm_K_Ta(k), 0.0_dp, fm_tmp)
         DO i_t = 1, n_t
            CALL cp_fm_scale_and_add(1.0_dp, fm_H(i_t), &
                                     bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, k), fm_tmp)
         END DO
      END DO
      CALL cp_fm_release(fm_F_T)
      CALL cp_fm_release(fm_tmp)
      DO i_t = 1, n_t
         CALL cp_fm_release(fm_K_Ta(i_t))
      END DO

      ! 5. Y_t = ½ F_a K_aa(iτ_t) + H_t,  W_f(iτ_t) += Y_t F_aᵀ + F_a Y_tᵀ
      DO i_t = 1, n_t
         CALL parallel_gemm('N', 'N', n_RI, n_ea, n_ea, 0.5_dp, fm_F_a, fm_K_aa(i_t), 1.0_dp, fm_H(i_t))
         CALL cp_fm_release(fm_K_aa(i_t))
         CALL parallel_gemm('N', 'T', n_RI, n_RI, n_ea, 1.0_dp, fm_H(i_t), fm_F_a, 1.0_dp, &
                            fm_W_folded_time(i_t))
         CALL parallel_gemm('N', 'T', n_RI, n_RI, n_ea, 1.0_dp, fm_F_a, fm_H(i_t), 1.0_dp, &
                            fm_W_folded_time(i_t))
         CALL cp_fm_release(fm_H(i_t))
      END DO
      DEALLOCATE (fm_H, fm_K_aa, fm_K_Ta)
      CALL cp_fm_release(fm_F_a)
      CALL cp_fm_struct_release(struct_F)
      CALL cp_fm_struct_release(struct_aa)

      ! The correction is now inside W; its ingredients are dead.
      DEALLOCATE (gc%K_head, gc%K_wing, gc%T_RI)
      DEALLOCATE (gc%K_head1, gc%K_head2, gc%K_wing1)
      DEALLOCATE (gc%K_rem_head, gc%K_rem_wing)
      DEALLOCATE (c, frac_atom, atom_of_RI, m_z, s_z, coef, v_ab, v_ba)

      gc%wall_time = gc%wall_time + m_walltime() - start_time
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A,T59,A,F7.1,A)') &
         'Computed Γ-point correction of W (2D slab),', 'Execution time', gc%wall_time, ' s'

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief Bilinear interpolation weights of a table over the minimum-image fractional displacement
!>        of atoms X − Y (layout of gamma_corr_tabulate: i_R = 1 + i_u + n i_v over [−1/2, 1/2]²).
!> \param X first atom
!> \param Y second atom
!> \param n table points per direction
!> \param k00 first table index
!> \param w weights (w00, w10, w01, w11)
! **************************************************************************************************
      SUBROUTINE pair_weights(X, Y, n, k00, w)
      INTEGER, INTENT(IN)                                :: X, Y, n
      INTEGER, INTENT(OUT)                               :: k00
      REAL(KIND=dp), DIMENSION(4), INTENT(OUT)           :: w

      INTEGER                                            :: i_u, i_v
      REAL(KIND=dp)                                      :: du, dv, f_u, f_v

         du = frac_atom(1, X) - frac_atom(1, Y)
         dv = frac_atom(2, X) - frac_atom(2, Y)
         du = du - ANINT(du)
         dv = dv - ANINT(dv)
         f_u = (du + 0.5_dp)*REAL(n - 1, dp)
         f_v = (dv + 0.5_dp)*REAL(n - 1, dp)
         i_u = MIN(MAX(INT(f_u), 0), n - 2)
         i_v = MIN(MAX(INT(f_v), 0), n - 2)
         f_u = f_u - REAL(i_u, dp)
         f_v = f_v - REAL(i_v, dp)
         w(1) = (1.0_dp - f_u)*(1.0_dp - f_v)
         w(2) = f_u*(1.0_dp - f_v)
         w(3) = (1.0_dp - f_u)*f_v
         w(4) = f_u*f_v
         k00 = 1 + i_u + n*i_v
      END SUBROUTINE pair_weights

! **************************************************************************************************
!> \brief Head kernel of coefficient types (k, l) for the atom pair (X, Y), at R_XY, per frequency
!>        (without the 1/A): c c: K_h + K^rem_h, c m / m c: −K_h1, m m: K_h2, c s / s c: K_h2/2.
!> \param X atom of the first function
!> \param Y atom of the second function
!> \param k coefficient type of the first (1 c, 2 m, 3 s)
!> \param l coefficient type of the second
!> \param v kernel per frequency
! **************************************************************************************************
      SUBROUTINE head_element(X, Y, k, l, v)
      INTEGER, INTENT(IN)                                :: X, Y, k, l
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: v

      INTEGER                                            :: cls_hi, cls_lo, i_pair, k00, kl, n
      REAL(KIND=dp), DIMENSION(4)                        :: w

         v(:) = 0.0_dp
         kl = 10*MIN(k, l) + MAX(k, l)
         IF (kl == 23 .OR. kl == 33) RETURN
         n = gc%n_tab(1)
         CALL pair_weights(X, Y, n, k00, w)
         SELECT CASE (kl)
         CASE (11)
            v(:) = w(1)*gc%K_head(k00, :) + w(2)*gc%K_head(k00 + 1, :) + &
                   w(3)*gc%K_head(k00 + n, :) + w(4)*gc%K_head(k00 + n + 1, :)
         CASE (12)
            v(:) = -(w(1)*gc%K_head1(k00, :) + w(2)*gc%K_head1(k00 + 1, :) + &
                     w(3)*gc%K_head1(k00 + n, :) + w(4)*gc%K_head1(k00 + n + 1, :))
         CASE (22, 13)
            v(:) = w(1)*gc%K_head2(k00, :) + w(2)*gc%K_head2(k00 + 1, :) + &
                   w(3)*gc%K_head2(k00 + n, :) + w(4)*gc%K_head2(k00 + n + 1, :)
         END SELECT
         ! (s_P c_Q + c_P s_Q)/2: half of K_h2 in each of the c s and s c entries
         IF (kl == 13) v(:) = 0.5_dp*v(:)
         ! remainder outside the cell, charge-charge only, per pair of form-factor classes
         IF (kl == 11) THEN
            n = gc%n_tab_rem
            CALL pair_weights(X, Y, n, k00, w)
            cls_lo = MIN(gc%atom_class(X), gc%atom_class(Y))
            cls_hi = MAX(gc%atom_class(X), gc%atom_class(Y))
            i_pair = cls_lo + cls_hi*(cls_hi - 1)/2
            v(:) = v(:) + w(1)*gc%K_rem_head(k00, :, i_pair) + w(2)*gc%K_rem_head(k00 + 1, :, i_pair) + &
                   w(3)*gc%K_rem_head(k00 + n, :, i_pair) + w(4)*gc%K_rem_head(k00 + n + 1, :, i_pair)
         END IF
      END SUBROUTINE head_element

! **************************************************************************************************
!> \brief Wing kernels (without 1/√A) of the element (T_X α, coefficient type l of atom Y), with the
!>        class of Y, in the Ta form and in the aT form of the transposed element (s: 0):
!>          Ta at R_XY: c: K_w,α + K^rem_w,α,     m: −K_w1,α
!>          aT at R_YX: c: −(K_w,α + K^rem_w,α),  m: +K_w1,α     (K^rem only for α = 1, 2)
!> \param X atom of T
!> \param Y atom of the coefficient
!> \param alpha Cartesian component of T
!> \param l coefficient type
!> \param v_ta Ta form per frequency
!> \param v_at aT form per frequency
! **************************************************************************************************
      SUBROUTINE wing_element(X, Y, alpha, l, v_ta, v_at)
      INTEGER, INTENT(IN)                                :: X, Y, alpha, l
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: v_ta, v_at

         v_ta(:) = 0.0_dp
         v_at(:) = 0.0_dp
         IF (l == 3) RETURN
         CALL wing_kernel(X, Y, Y, alpha, l, v_ta)
         CALL wing_kernel(Y, X, Y, alpha, l, v_at)
         v_at(:) = -v_at(:)
      END SUBROUTINE wing_element

! **************************************************************************************************
!> \brief Wing kernel at R_XY with the class of atom Z: c: K_w,α + K^rem_w,α, m: −K_w1,α.
!> \param X first atom of the displacement
!> \param Y second atom of the displacement
!> \param Z atom whose height / class selects the form factor
!> \param alpha Cartesian component
!> \param l coefficient type (1 c, 2 m)
!> \param v kernel per frequency
! **************************************************************************************************
      SUBROUTINE wing_kernel(X, Y, Z, alpha, l, v)
      INTEGER, INTENT(IN)                                :: X, Y, Z, alpha, l
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: v

      INTEGER                                            :: k00, n
      REAL(KIND=dp), DIMENSION(4)                        :: w

         n = gc%n_tab(1)
         CALL pair_weights(X, Y, n, k00, w)
         IF (l == 1) THEN
            v(:) = w(1)*gc%K_wing(k00, alpha, :) + w(2)*gc%K_wing(k00 + 1, alpha, :) + &
                   w(3)*gc%K_wing(k00 + n, alpha, :) + w(4)*gc%K_wing(k00 + n + 1, alpha, :)
         ELSE
            v(:) = -(w(1)*gc%K_wing1(k00, alpha, :) + w(2)*gc%K_wing1(k00 + 1, alpha, :) + &
                     w(3)*gc%K_wing1(k00 + n, alpha, :) + w(4)*gc%K_wing1(k00 + n + 1, alpha, :))
         END IF
         IF (l == 1 .AND. alpha <= 2) THEN
            n = gc%n_tab_rem
            CALL pair_weights(X, Y, n, k00, w)
            v(:) = v(:) + w(1)*gc%K_rem_wing(k00, alpha, :, gc%atom_class(Z)) + &
                   w(2)*gc%K_rem_wing(k00 + 1, alpha, :, gc%atom_class(Z)) + &
                   w(3)*gc%K_rem_wing(k00 + n, alpha, :, gc%atom_class(Z)) + &
                   w(4)*gc%K_rem_wing(k00 + n + 1, alpha, :, gc%atom_class(Z))
         END IF
      END SUBROUTINE wing_kernel

   END SUBROUTINE gamma_corr_add_folded_to_W

! **************************************************************************************************
!> \brief F = M⁺_c E for a block-sparse E whose column (A, k) (index n_k(A−1) + k) holds coef(P, k)
!>        on the rows P of atom A only; F must have the row layout of M⁺:
!>          F = M⁺ E − [U_L diag(1/m_L)] (U_Lᵀ E),   (M⁺ E)(:, (A,k)) = Σ_{P∈A} M⁺(:, P) coef(P, k)
!> \param bs_env band-structure environment, for M⁺ = fm_Minv_Gamma and U_L, m_L of gc
!> \param coef per-function coefficients, (n_RI, n_k), replicated
!> \param atom_of_RI atom of every RI function
!> \param n_atom number of atoms
!> \param fm_F result, n_RI x n_k n_atom
! **************************************************************************************************
   SUBROUTINE gamma_corr_minv_c_atom_columns(bs_env, coef, atom_of_RI, n_atom, fm_F)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: coef
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_RI
      INTEGER, INTENT(IN)                                :: n_atom
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_F

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_minv_c_atom_columns'
      INTEGER, PARAMETER                                 :: chunk = 2048, row_blk = 256

      INTEGER :: c0, c1, handle, i_blk, i_col, i_l, i_row, j_col, jc, k, n_blk, n_e, n_k, n_RI, &
                 ncol_F, ncol_local, nrow_F, nrow_local, P, r0, r1
      INTEGER, DIMENSION(:), POINTER                     :: col_F, col_indices, row_F, row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: inv_m
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buf, S
      TYPE(cp_fm_struct_type), POINTER                   :: struct_S
      TYPE(cp_fm_type)                                   :: fm_S, fm_Us
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(mp_comm_type)                                 :: comm_row

      CALL timeset(routineN, handle)
      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI
      n_k = SIZE(coef, 2)
      n_e = n_k*n_atom

      CALL cp_fm_get_info(bs_env%fm_Minv_Gamma, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      CALL cp_fm_get_info(fm_F, nrow_local=nrow_F, ncol_local=ncol_F, row_indices=row_F, &
                          col_indices=col_F)
      CPASSERT(nrow_F == nrow_local)
      CPASSERT(ALL(row_F(1:nrow_local) == row_indices(1:nrow_local)))

      ! ranks sharing a process row hold the same rows of M⁺ and of F
      CALL comm_row%from_split(bs_env%fm_Minv_Gamma%matrix_struct%para_env, &
                               bs_env%fm_Minv_Gamma%matrix_struct%context%mepos(1), &
                               bs_env%fm_Minv_Gamma%matrix_struct%context%mepos(2))

      ! (M⁺ E) by column chunks of F
      DO c0 = 1, n_e, chunk
         c1 = MIN(c0 + chunk - 1, n_e)
         ALLOCATE (buf(nrow_local, c0:c1))
         buf(:, :) = 0.0_dp
         ! buf(:, (A,k)) = Σ_{P∈A} coef(P, k) M⁺(:, P) on the local rows; the threads split the rows
         ! (several P add into the same column)
         n_blk = (nrow_local + row_blk - 1)/row_blk
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(static) PRIVATE(i_blk, r0, r1, i_col, P, k, jc) &
!$OMP          SHARED(n_blk, nrow_local, ncol_local, col_indices, n_k, atom_of_RI, c0, c1, coef, &
!$OMP                 buf, bs_env)
         DO i_blk = 1, n_blk
            r0 = (i_blk - 1)*row_blk + 1
            r1 = MIN(i_blk*row_blk, nrow_local)
            DO i_col = 1, ncol_local
               P = col_indices(i_col)
               DO k = 1, n_k
                  jc = n_k*(atom_of_RI(P) - 1) + k
                  IF (jc < c0 .OR. jc > c1) CYCLE
                  IF (coef(P, k) == 0.0_dp) CYCLE
                  buf(r0:r1, jc) = buf(r0:r1, jc) + &
                                   coef(P, k)*bs_env%fm_Minv_Gamma%local_data(r0:r1, i_col)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL comm_row%sum(buf)
         DO j_col = 1, ncol_F
            jc = col_F(j_col)
            IF (jc < c0 .OR. jc > c1) CYCLE
            fm_F%local_data(1:nrow_local, j_col) = buf(:, jc)
         END DO
         DEALLOCATE (buf)
      END DO
      CALL comm_row%free()

      IF (gc%n_low > 0) THEN
         ! S = U_Lᵀ E (n_low x n_e), replicated: local elements of U_L, summed over all ranks
         ALLOCATE (S(gc%n_low, n_e))
         S(:, :) = 0.0_dp
         CALL cp_fm_get_info(gc%fm_U_low, nrow_local=nrow_F, ncol_local=ncol_F, row_indices=row_F, &
                             col_indices=col_F)
         ! every local column i_col is a different mode i_l, so the threads write disjoint rows of S
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(static) PRIVATE(i_col, i_l, i_row, P, k, jc) &
!$OMP          SHARED(ncol_F, nrow_F, col_F, row_F, n_k, atom_of_RI, S, gc, coef)
         DO i_col = 1, ncol_F
            i_l = col_F(i_col)
            DO i_row = 1, nrow_F
               P = row_F(i_row)
               DO k = 1, n_k
                  jc = n_k*(atom_of_RI(P) - 1) + k
                  S(i_l, jc) = S(i_l, jc) + gc%fm_U_low%local_data(i_row, i_col)*coef(P, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL gc%fm_U_low%matrix_struct%para_env%sum(S)
         ! F −= [U_L diag(1/m_L)] S
         NULLIFY (struct_S)
         CALL cp_fm_struct_create(struct_S, template_fmstruct=fm_F%matrix_struct, nrow_global=gc%n_low)
         CALL cp_fm_create(fm_S, struct_S)
         CALL cp_fm_struct_release(struct_S)
         CALL cp_fm_get_info(fm_S, nrow_local=nrow_F, ncol_local=ncol_F, row_indices=row_F, &
                             col_indices=col_F)
         DO i_col = 1, ncol_F
            DO i_row = 1, nrow_F
               fm_S%local_data(i_row, i_col) = S(row_F(i_row), col_F(i_col))
            END DO
         END DO
         DEALLOCATE (S)
         CALL cp_fm_create(fm_Us, gc%fm_U_low%matrix_struct)
         CALL cp_fm_to_fm(gc%fm_U_low, fm_Us)
         ALLOCATE (inv_m(gc%n_low))
         inv_m(:) = 1.0_dp/gc%m_low(:)
         CALL cp_fm_column_scale(fm_Us, inv_m)
         CALL parallel_gemm('N', 'N', n_RI, n_e, gc%n_low, -1.0_dp, fm_Us, fm_S, 1.0_dp, fm_F)
         DEALLOCATE (inv_m)
         CALL cp_fm_release(fm_S)
         CALL cp_fm_release(fm_Us)
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_minv_c_atom_columns

! **************************************************************************************************
!> \brief Out-of-plane moments of the auxiliary functions about the sheet plane z_s (mean atom
!>        height, minimum image along z relative to atom 1),
!>
!>          m_P = ∫dr ϕ_P(r) |z − z_s|,     s_P = ∫dr ϕ_P(r) (z − z_s)².
!>
!>        For a Cartesian primitive x^{n_x} y^{n_y} u^n e^{−α r²} at height Z_P (u = z − Z_P,
!>        d = Z_P − z_s) they are μ_{n_x} μ_{n_y} F_n(d) and μ_{n_x} μ_{n_y} (μ_{n+2} + 2d μ_{n+1} + d² μ_n),
!>        with μ_k = ∫du u^k e^{−αu²} and F_n(d) of z_abs_moment, contracted with sphi to the RI
!>        functions (ordering of gamma_auxiliary_charges).
!> \param qs_env QS environment, for the particles and the RI_AUX basis
!> \param m_z m_P (n_RI)
!> \param s_z s_P (n_RI)
! **************************************************************************************************
   SUBROUTINE gamma_corr_z_moments(qs_env, m_z, s_z)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: m_z, s_z

      INTEGER                                            :: first, iatom, ico, il, ipgf, iset, isgf, &
                                                            lmax, n_atom, ncart, nx, ny, nz, offset
      REAL(KIND=dp)                                      :: alpha, d, f_abs, f_sq, l_z, z_sheet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mom, p_abs, p_sq, z_atom
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds

      CALL get_qs_env(qs_env, particle_set=particles, qs_kind_set=qs_kinds, cell=cell)
      n_atom = SIZE(particles)

      ! Atom heights, unwrapped along z relative to atom 1 (minimum image), and their mean z_s.
      ALLOCATE (z_atom(n_atom))
      l_z = cell%hmat(3, 3)
      DO iatom = 1, n_atom
         d = particles(iatom)%r(3) - particles(1)%r(3)
         IF (l_z > 0.0_dp) d = d - l_z*ANINT(d/l_z)
         z_atom(iatom) = particles(1)%r(3) + d
      END DO
      z_sheet = SUM(z_atom)/REAL(n_atom, dp)

      m_z(:) = 0.0_dp
      s_z(:) = 0.0_dp
      first = 0
      DO iatom = 1, n_atom
         ! d = Z_P − z_s for every function on this atom
         d = z_atom(iatom) - z_sheet
         CALL get_qs_kind(qs_kinds(particles(iatom)%atomic_kind%kind_number), &
                          basis_set=basis, basis_type='RI_AUX')
         DO iset = 1, basis%nset
            lmax = basis%lmax(iset)
            ncart = ncoset(lmax)
            ALLOCATE (mom(0:lmax + 2), p_abs(ncart*basis%npgf(iset)), p_sq(ncart*basis%npgf(iset)))
            p_abs(:) = 0.0_dp
            p_sq(:) = 0.0_dp
            DO ipgf = 1, basis%npgf(iset)
               alpha = basis%zet(ipgf, iset)
               ! μ_k = ∫du u^k e^{−αu²}: μ_0 = √(π/α), μ_k = (k−1)/(2α) μ_{k−2}, odd k zero
               mom(:) = 0.0_dp
               mom(0) = SQRT(pi/alpha)
               DO il = 2, lmax + 2, 2
                  mom(il) = REAL(il - 1, dp)*mom(il - 2)/(2.0_dp*alpha)
               END DO
               DO ico = ncoset(basis%lmin(iset) - 1) + 1, ncart
                  nx = indco(1, ico)
                  ny = indco(2, ico)
                  nz = indco(3, ico)
                  ! F_n(d) = ∫u^n e^{−αu²}|u + d|
                  f_abs = z_abs_moment(nz, alpha, d, mom)
                  ! ∫u^n e^{−αu²}(u + d)²
                  f_sq = mom(nz + 2) + 2.0_dp*d*mom(nz + 1) + d*d*mom(nz)
                  p_abs((ipgf - 1)*ncart + ico) = mom(nx)*mom(ny)*f_abs
                  p_sq((ipgf - 1)*ncart + ico) = mom(nx)*mom(ny)*f_sq
               END DO
            END DO
            offset = basis%first_sgf(1, iset)
            DO isgf = offset, offset + basis%nsgf_set(iset) - 1
               m_z(first + isgf) = SUM(basis%sphi(1:SIZE(p_abs), isgf)*p_abs)
               s_z(first + isgf) = SUM(basis%sphi(1:SIZE(p_sq), isgf)*p_sq)
            END DO
            DEALLOCATE (mom, p_abs, p_sq)
         END DO
         first = first + basis%nsgf
      END DO
      CPASSERT(first == SIZE(m_z))
      DEALLOCATE (z_atom)

   CONTAINS

! **************************************************************************************************
!> \brief F_n(d) = ∫du u^n e^{−αu²} |u + d| = J_{n+1}(d) + d J_n(d), with
!>
!>          J_k(d) = ∫du u^k e^{−αu²} sign(u + d) = 2 I_k(−d) − μ_k,   I_k(a) = ∫_a^∞ du u^k e^{−αu²},
!>          I_0 = (1/2)√(π/α) erfc(√α a),   I_1 = e^{−αa²}/(2α),
!>          I_k = a^{k−1} e^{−αa²}/(2α) + (k−1)/(2α) I_{k−2}.
!> \param n power of u
!> \param alpha Gaussian exponent
!> \param d offset Z_P − z_s
!> \param mu full moments μ_k, k = 0 .. ≥ n + 1
!> \return F_n(d)
! **************************************************************************************************
      PURE FUNCTION z_abs_moment(n, alpha, d, mu) RESULT(f)
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: alpha, d
      REAL(KIND=dp), DIMENSION(0:), INTENT(IN)           :: mu
      REAL(KIND=dp)                                      :: f

      INTEGER                                            :: k
      REAL(KIND=dp)                                      :: a, gauss
      REAL(KIND=dp), DIMENSION(0:n+1)                    :: i_tail, j_sign

         ! I_k(a) at a = −d by the upward recurrence
         a = -d
         gauss = EXP(-alpha*a*a)
         i_tail(0) = 0.5_dp*SQRT(pi/alpha)*ERFC(SQRT(alpha)*a)
         i_tail(1) = gauss/(2.0_dp*alpha)
         DO k = 2, n + 1
            i_tail(k) = a**(k - 1)*gauss/(2.0_dp*alpha) + &
                        REAL(k - 1, dp)/(2.0_dp*alpha)*i_tail(k - 2)
         END DO
         ! J_k = 2 I_k(−d) − μ_k
         DO k = 0, n + 1
            j_sign(k) = 2.0_dp*i_tail(k) - mu(k)
         END DO
         f = j_sign(n + 1) + d*j_sign(n)
      END FUNCTION z_abs_moment

   END SUBROUTINE gamma_corr_z_moments

! **************************************************************************************************
!> \brief Directions q̂, angular weights and radii L(q̂) of the Wigner-Seitz cell of the supercell
!>        reciprocal lattice around Γ, for a bulk cell (solid angle) or a slab (in-plane angle):
!>
!>          WS = { q : q·g ≤ |g|²/2 for all g ≠ 0 } ,   L(q̂) = min_{g: q̂·g > 0} (|g|²/2)/(q̂·g) .
!>
!>        Bulk: each face (clipped by the other g, Sutherland-Hodgman) is fanned into triangles
!>        (c, v_k, v_k+1) with Gauss-Legendre in the Duffy map q(s, t) = c + s [(v_k − c) + t (v_k+1 − v_k)]
!>        and dΩ = (q·ĝ) dA/|q|³. Slab: Gauss-Legendre along each edge with dθ = |q × dq|/|q|².
!>        The −g faces are added as exact negative images, so the set is symmetric under q̂ → −q̂.
!> \param bs_env band-structure environment, for the cell (via gamma_corr_cell_geometry)
!> \param slab .TRUE. for a slab (PERIODIC XY), .FALSE. for bulk
!> \param n_face Gauss-Legendre points per edge (slab) or per triangle side (bulk)
!> \param dir_hat unit directions, 3 x n_dir, Cartesian
!> \param dir_weight angular weights; they sum to 2π (slab) or 4π (bulk)
!> \param dir_radius distance to the cell boundary along each direction (1/bohr)
! **************************************************************************************************
   SUBROUTINE gamma_corr_ws_cell(bs_env, slab, n_face, dir_hat, dir_weight, dir_radius)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL, INTENT(IN)                                :: slab
      INTEGER, INTENT(IN)                                :: n_face
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: dir_hat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: dir_weight, dir_radius

      INTEGER, PARAMETER                                 :: max_vert = 32, n_shell = 3
      REAL(KIND=dp), PARAMETER :: eps_measure = 1.0E-8_dp, eps_parallel = 1.0E-12_dp, &
                                  eps_vertex = 1.0E-10_dp, min_perp_norm2 = 0.1_dp

      INTEGER                                            :: i_dir, i_g, i_node, i_sign, i_v, j_c, &
                                                            j_node, n1, n2, n3, n_cand, n_dir, &
                                                            n_half, n_per, n_v, n_v_new
      LOGICAL                                            :: in_a, in_b
      REAL(KIND=dp)                                      :: big, d_a, d_b, g_len, half_g2, measure, &
                                                            proj, rhs, s_hi, s_lo, scale_b, tol, &
                                                            triangle, w_dir
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: nodes, r_tmp, rule, w_tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cand, h_tmp
      REAL(KIND=dp), DIMENSION(3)                        :: centre, cross, dq, e_1, e_2, g_hat, &
                                                            normal, q, t_line
      REAL(KIND=dp), DIMENSION(3, 3)                     :: b
      REAL(KIND=dp), DIMENSION(3, max_vert)              :: poly, poly_new

      CALL gamma_corr_cell_geometry(bs_env, slab, measure, b)
      ALLOCATE (nodes(n_face), rule(n_face))
      CALL gauss_legendre(n_face, nodes, rule)
      n_per = 3
      IF (slab) n_per = 2
      scale_b = MAXVAL(NORM2(b(:, 1:n_per), DIM=1))
      tol = eps_vertex*scale_b
      big = 4.0_dp*REAL(n_shell, dp)*SUM(NORM2(b(:, 1:n_per), DIM=1))

      ! Candidate lattice vectors g = Σ n_i b_i, |n_i| ≤ n_shell, g ≠ 0. The first n_half are
      ! the "positive" half (first non-zero coefficient > 0); the rest are their negatives.
      n_cand = (2*n_shell + 1)**n_per - 1
      ALLOCATE (cand(3, n_cand))
      n_half = 0
      DO n1 = -n_shell, n_shell
         DO n2 = -n_shell, n_shell
            DO n3 = -n_shell, n_shell
               IF (slab .AND. n3 /= 0) CYCLE
               IF (n1 < 0 .OR. (n1 == 0 .AND. n2 < 0) .OR. (n1 == 0 .AND. n2 == 0 .AND. n3 <= 0)) CYCLE
               n_half = n_half + 1
               cand(:, n_half) = REAL(n1, dp)*b(:, 1) + REAL(n2, dp)*b(:, 2)
               IF (.NOT. slab) cand(:, n_half) = cand(:, n_half) + REAL(n3, dp)*b(:, 3)
            END DO
         END DO
      END DO
      cand(:, n_half + 1:2*n_half) = -cand(:, 1:n_half)
      n_cand = 2*n_half

      ! Upper bound on the number of directions, for the temporary arrays.
      IF (slab) THEN
         n_dir = 2*n_half*n_face
      ELSE
         n_dir = 2*n_half*max_vert*n_face**2
      END IF
      ALLOCATE (h_tmp(3, n_dir), w_tmp(n_dir), r_tmp(n_dir))
      i_dir = 0

      IF (slab) THEN
         ! Unit normal of the periodic plane.
         normal = [b(2, 1)*b(3, 2) - b(3, 1)*b(2, 2), b(3, 1)*b(1, 2) - b(1, 1)*b(3, 2), &
                   b(1, 1)*b(2, 2) - b(2, 1)*b(1, 2)]
         normal = normal/NORM2(normal)
      END IF

      DO i_g = 1, n_half
         half_g2 = 0.5_dp*SUM(cand(:, i_g)**2)
         g_len = SQRT(2.0_dp*half_g2)
         g_hat = cand(:, i_g)/g_len

         IF (slab) THEN
            ! Edge: q = g/2 + s t, t = n × ĝ; each other g' bounds s via q·g' ≤ |g'|²/2.
            t_line = vector_product(normal, g_hat)
            s_lo = -big
            s_hi = big
            DO j_c = 1, n_cand
               IF (j_c == i_g) CYCLE
               ! s (t·g') ≤ |g'|²/2 − (g/2)·g'
               proj = DOT_PRODUCT(t_line, cand(:, j_c))
               rhs = 0.5_dp*SUM(cand(:, j_c)**2) - 0.5_dp*DOT_PRODUCT(cand(:, i_g), cand(:, j_c))
               IF (ABS(proj) < eps_parallel*scale_b) THEN
                  IF (rhs < -tol) s_hi = s_lo - 1.0_dp
               ELSE IF (proj > 0.0_dp) THEN
                  s_hi = MIN(s_hi, rhs/proj)
               ELSE
                  s_lo = MAX(s_lo, rhs/proj)
               END IF
            END DO
            IF (s_hi - s_lo <= tol) CYCLE
            ! Gauss-Legendre on the edge, x ∈ [−1, 1]: q = g/2 + (s_mid + x s_half) t.
            dq = 0.5_dp*(s_hi - s_lo)*t_line
            DO i_sign = 1, -1, -2
               DO i_node = 1, n_face
                  q = REAL(i_sign, dp)*(0.5_dp*cand(:, i_g) + 0.5_dp*(s_hi + s_lo)*t_line + &
                                        nodes(i_node)*dq)
                  i_dir = i_dir + 1
                  r_tmp(i_dir) = NORM2(q)
                  h_tmp(:, i_dir) = q/r_tmp(i_dir)
                  ! dθ = |q × dq/dx| dx/|q|²
                  cross = vector_product(q, dq)
                  w_tmp(i_dir) = rule(i_node)*NORM2(cross)/SUM(q*q)
               END DO
            END DO
            CYCLE
         END IF

         ! Bulk face: start from a large square in the plane q·ĝ = |g|/2, clip by every other g'.
         e_1 = [g_hat(2), -g_hat(1), 0.0_dp]
         IF (SUM(e_1**2) < min_perp_norm2) e_1 = [0.0_dp, g_hat(3), -g_hat(2)]
         e_1 = e_1/NORM2(e_1)
         e_2 = vector_product(g_hat, e_1)
         centre = 0.5_dp*cand(:, i_g)
         n_v = 4
         poly(:, 1) = centre + big*(e_1 + e_2)
         poly(:, 2) = centre + big*(-e_1 + e_2)
         poly(:, 3) = centre + big*(-e_1 - e_2)
         poly(:, 4) = centre + big*(e_1 - e_2)
         DO j_c = 1, n_cand
            IF (j_c == i_g .OR. n_v == 0) CYCLE
            rhs = 0.5_dp*SUM(cand(:, j_c)**2)
            ! Sutherland-Hodgman against the half-space q·g' ≤ |g'|²/2.
            n_v_new = 0
            DO i_v = 1, n_v
               d_a = DOT_PRODUCT(poly(:, i_v), cand(:, j_c)) - rhs
               d_b = DOT_PRODUCT(poly(:, MODULO(i_v, n_v) + 1), cand(:, j_c)) - rhs
               in_a = d_a <= 0.0_dp
               in_b = d_b <= 0.0_dp
               IF (in_a) THEN
                  n_v_new = n_v_new + 1
                  poly_new(:, n_v_new) = poly(:, i_v)
               END IF
               IF (in_a .NEQV. in_b) THEN
                  ! Crossing point of the edge with the plane.
                  n_v_new = n_v_new + 1
                  poly_new(:, n_v_new) = poly(:, i_v) + d_a/(d_a - d_b)* &
                                         (poly(:, MODULO(i_v, n_v) + 1) - poly(:, i_v))
               END IF
               IF (n_v_new > max_vert - 1) THEN
                  CPABORT('Gamma-point correction: Wigner-Seitz face has too many vertices')
               END IF
            END DO
            ! Drop consecutive (near-)duplicate vertices left by cuts through a vertex.
            n_v = 0
            DO i_v = 1, n_v_new
               IF (n_v > 0) THEN
                  IF (NORM2(poly_new(:, i_v) - poly(:, n_v)) < tol) CYCLE
               END IF
               n_v = n_v + 1
               poly(:, n_v) = poly_new(:, i_v)
            END DO
            IF (n_v > 1) THEN
               IF (NORM2(poly(:, n_v) - poly(:, 1)) < tol) n_v = n_v - 1
            END IF
            IF (n_v < 3) n_v = 0
         END DO
         IF (n_v < 3) CYCLE
         ! Fan of triangles (c, v_k, v_k+1) from the vertex mean c; zero-area triangles are skipped
         centre = SUM(poly(:, 1:n_v), DIM=2)/REAL(n_v, dp)
         DO i_v = 1, n_v
            ! Triangle sides e_1 = v_k − c, e_2 = v_k+1 − v_k and twice its area |e_1 × e_2|.
            e_1 = poly(:, i_v) - centre
            e_2 = poly(:, MODULO(i_v, n_v) + 1) - poly(:, i_v)
            cross = vector_product(e_1, e_2)
            triangle = NORM2(cross)
            IF (triangle < tol*scale_b) CYCLE
            DO i_node = 1, n_face
               DO j_node = 1, n_face
                  ! (s, t) ∈ [0, 1]² from the Gauss-Legendre nodes on [−1, 1]; weight
                  ! (w_s/2)(w_t/2) s |e_1 × e_2| (Duffy Jacobian), times dΩ/dA = (|g|/2)/|q|³.
                  q = centre + 0.5_dp*(nodes(i_node) + 1.0_dp)* &
                      (e_1 + 0.5_dp*(nodes(j_node) + 1.0_dp)*e_2)
                  w_dir = 0.25_dp*rule(i_node)*rule(j_node)*0.5_dp*(nodes(i_node) + 1.0_dp)* &
                          triangle*0.5_dp*g_len/NORM2(q)**3
                  DO i_sign = 1, -1, -2
                     i_dir = i_dir + 1
                     r_tmp(i_dir) = NORM2(q)
                     h_tmp(:, i_dir) = REAL(i_sign, dp)*q/r_tmp(i_dir)
                     w_tmp(i_dir) = w_dir
                  END DO
               END DO
            END DO
         END DO
      END DO

      n_dir = i_dir
      ALLOCATE (dir_hat(3, n_dir), dir_weight(n_dir), dir_radius(n_dir))
      dir_hat(:, :) = h_tmp(:, 1:n_dir)
      dir_weight(:) = w_tmp(1:n_dir)
      dir_radius(:) = r_tmp(1:n_dir)
      DEALLOCATE (cand, h_tmp, w_tmp, r_tmp, nodes, rule)

      ! Closure: Σ w L²/2 = (2π)²/A (slab), Σ w L³/3 = (2π)³/Ω (bulk).
      IF (slab) THEN
         proj = 0.5_dp*SUM(dir_weight*dir_radius**2)*measure/twopi**2
      ELSE
         proj = SUM(dir_weight*dir_radius**3)/3.0_dp*measure/twopi**3
      END IF
      IF (ABS(proj - 1.0_dp) > eps_measure) THEN
         CPABORT('Gamma-point correction: Wigner-Seitz cell measure check failed')
      END IF

   END SUBROUTINE gamma_corr_ws_cell

! **************************************************************************************************
!> \brief Bulk cell: head constant w(iω) (gamma_corr_head_zone_remainder) and the wing (A) and
!>        head-dipole (B) tensors of the Γ-point correction, integrated over the Wigner-Seitz cell
!>        around Γ (gamma_corr_ws_cell: directions q̂, weights w(q̂), radii L(q̂)),
!>
!>          A_αβ(iω) = (1/(2π)³) Σ_q̂ w(q̂) ∫_0^{L(q̂)} dq q² 4π/(1 + 4πa) q̂_α q̂_β,
!>          B_αβ(iω) = (1/(2π)³) Σ_q̂ w(q̂) ∫_0^{L(q̂)} dq q² 4π (1/(1 + 4πa) − 1) q̂_α q̂_β,
!>
!>        with the finite-q head a = q̂ᵀ S q̂/(1 + γ(iω) q²) (gamma_corr_finite_q_head), and the
!>        TAIL_3D screening factor 1/ε_0(iω) − 1 = −4π ā/(1 + 4π ā), ā = tr S(iω)/3.
!> \param bs_env band-structure environment; consumes gc%S_head, gc%head_gamma and sets gc%head_pot_freq,
!>        gc%A_wing, gc%B_hd, gc%q_max, gc%tail_factor (TAIL_3D is switched off for an anisotropic head)
! **************************************************************************************************
   SUBROUTINE gamma_corr_direction_coefficients(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_direction_coefficients'
      INTEGER, PARAMETER                                 :: n_face = 8, n_radial = 32
      REAL(KIND=dp), PARAMETER                           :: max_anisotropy = 0.01_dp

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: handle, i_dir, i_node, i_w, n_dir, n_t
      REAL(KIND=dp)                                      :: a_q, a_scalar, anisotropy, denom, &
                                                            int_inv, int_vol, L, q, w_q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dir_weight, nodes, rule, w_zone
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_t = bs_env%num_time_freq_points

      CALL gamma_corr_ws_cell(bs_env, .FALSE., n_face, gc%dir_hat, dir_weight, gc%dir_radius)
      n_dir = SIZE(dir_weight)
      ALLOCATE (gc%A_wing(3, 3, n_t), gc%B_hd(3, 3, n_t))
      gc%A_wing(:, :, :) = 0.0_dp
      gc%B_hd(:, :, :) = 0.0_dp

      ! Gauss-Legendre rule for the radial integrals on [0, L(q̂)]
      CPASSERT(ALLOCATED(gc%head_gamma))
      ALLOCATE (nodes(n_radial), rule(n_radial))
      CALL gauss_legendre(n_radial, nodes, rule)

      ! TAIL_3D screening factor 1/ε_0(iω) − 1 = −4π ā/(1 + 4π ā), ā = tr S(iω)/3. It needs an
      ! isotropic head (diagonal of S within 1 % of ā); otherwise fall back to the head-only treatment.
      IF (gc%tail_3d) THEN
         DO i_w = 1, n_t
            a_scalar = (gc%S_head(1, 1, i_w) + gc%S_head(2, 2, i_w) + gc%S_head(3, 3, i_w))/3.0_dp
            anisotropy = MAXVAL(ABS([gc%S_head(1, 1, i_w), gc%S_head(2, 2, i_w), &
                                     gc%S_head(3, 3, i_w)] - a_scalar))/ABS(a_scalar)
            IF (anisotropy > max_anisotropy) THEN
               WRITE (message, '(A,F0.3,A,I0,A)') 'Gamma-point correction: the dielectric head '// &
                  'S(i omega) is anisotropic (max |S_aa - tr S/3|/(tr S/3) = ', anisotropy, &
                  ' at frequency ', i_w, '); the sphere-tail term needs isotropic screening '// &
                  '(cubic cells). Continuing with the head-only treatment.'
               CPWARN(TRIM(message))
               gc%tail_3d = .FALSE.
               IF (ALLOCATED(gc%tail_factor)) DEALLOCATE (gc%tail_factor)
               CALL cp_fm_release(gc%fm_K_tail)
               EXIT
            END IF
         END DO
      END IF
      IF (gc%tail_3d) THEN
         IF (ALLOCATED(gc%tail_factor)) DEALLOCATE (gc%tail_factor)
         ALLOCATE (gc%tail_factor(n_t))
         DO i_w = 1, n_t
            a_scalar = (gc%S_head(1, 1, i_w) + gc%S_head(2, 2, i_w) + gc%S_head(3, 3, i_w))/3.0_dp
            gc%tail_factor(i_w) = -fourpi*a_scalar/(1.0_dp + fourpi*a_scalar)
         END DO
      END IF

      DO i_w = 1, n_t
         DO i_dir = 1, n_dir
            ! a(q̂, iω) = q̂ᵀ S q̂ and the macroscopic denominator 1 + 4π a.
            a_scalar = DOT_PRODUCT(gc%dir_hat(:, i_dir), &
                                   MATMUL(gc%S_head(:, :, i_w), gc%dir_hat(:, i_dir)))
            denom = 1.0_dp + fourpi*a_scalar
            ! ε_M(q̂, iω) = 1 + 4π q̂ᵀ S q̂ must be positive in every direction
            IF (denom <= 0.0_dp) THEN
               CALL cp_abort(__LOCATION__, "Gamma-point correction: the fitted dielectric "// &
                             "head 1 + 4 pi q^T S q is not positive along some direction. The "// &
                             "cell is too small for the finite-q head (smallest |g| too large); "// &
                             "use a larger supercell.")
            END IF
            ! Radial integrals along q̂ with a(q) = a/(1 + γ q²), Gauss-Legendre on [0, L]:
            !   int_inv  = ∫_0^L dq q² /(1 + 4π a(q))            (wing tensor A)
            !   int_vol  = ∫_0^L dq q² [1/(1 + 4π a(q)) − 1]     (head-dipole tensor B)
            L = gc%dir_radius(i_dir)
            int_inv = 0.0_dp
            int_vol = 0.0_dp
            DO i_node = 1, n_radial
               ! node map x ∈ [−1, 1] → q = L (x + 1)/2, weight L/2 · w_GL
               q = 0.5_dp*L*(nodes(i_node) + 1.0_dp)
               w_q = 0.5_dp*L*rule(i_node)
               a_q = a_scalar/(1.0_dp + gc%head_gamma(i_w)*q*q)
               int_inv = int_inv + w_q*q*q/(1.0_dp + fourpi*a_q)
               int_vol = int_vol + w_q*q*q*(1.0_dp/(1.0_dp + fourpi*a_q) - 1.0_dp)
            END DO
            ! A_αβ += w(q̂) 4π int_inv q̂_α q̂_β/(2π)³,   B_αβ += w(q̂) 4π int_vol q̂_α q̂_β/(2π)³
            gc%A_wing(:, :, i_w) = gc%A_wing(:, :, i_w) + &
                                   dir_weight(i_dir)*int_inv/twopi**3*fourpi* &
                                   SPREAD(gc%dir_hat(:, i_dir), 2, 3)*SPREAD(gc%dir_hat(:, i_dir), 1, 3)
            gc%B_hd(:, :, i_w) = gc%B_hd(:, :, i_w) + &
                                 dir_weight(i_dir)*int_vol/twopi**3*fourpi* &
                                 SPREAD(gc%dir_hat(:, i_dir), 2, 3)*SPREAD(gc%dir_hat(:, i_dir), 1, 3)
         END DO
      END DO
      ! Head constant w(iω) = ∫ d³q/(2π)³ W_h(q) − (1/Ω) Σ_{g≠0} W_h(g) (gamma_corr_head_zone_remainder)
      ALLOCATE (w_zone(n_t))
      CALL gamma_corr_head_zone_remainder(bs_env, dir_weight, gc%head_gamma, w_zone)
      DEALLOCATE (nodes, rule, gc%head_gamma)

      ! Head potential w(iω) = w_zone.
      ALLOCATE (gc%head_pot_freq(n_t))
      gc%head_pot_freq(:) = w_zone(:)

      ! Largest radius of the BZ cell; gamma_corr_sigma_shifts clips |d_nm| at 1/q_max
      gc%q_max = MAXVAL(gc%dir_radius)
      ! The direction data are consumed; T_RI, A_wing, B_hd and head_pot_freq are kept for gamma_corr_sigma_shifts
      DEALLOCATE (gc%S_head, dir_weight, gc%dir_hat, gc%dir_radius)
      DEALLOCATE (w_zone)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_direction_coefficients

! **************************************************************************************************
!> \brief Bulk cell: head constant w(iω) as the remainder of the zone sum of the head model
!>        (Gygi and Baldereschi, PRB 34, 4405 (1986)),
!>
!>          w(iω) = ∫ d³q/(2π)³ W_h(q) − (1/Ω) Σ_{g≠0} W_h(g),
!>          W_h(q) = −(4π/q²) 4π a_q/(1 + 4π a_q),   a_q = q̂ᵀ S q̂/(1 + γ q²).
!>
!>        Both terms are damped by e^{−σ²q²}; the damped difference D(σ) = w + c₁σ² + c₂σ⁴ +
!>        O(e^{−R²/4σ²}) (R = shortest supercell lattice vector) is extrapolated to σ = 0 by a
!>        quadratic fit in σ² through σ = R/8, R/10, R/12.
!> \param bs_env band-structure environment: hmat and gc%S_head, gc%b_recip, gc%measure, gc%dir_hat
!> \param dir_weight solid-angle weights of the BZ-cell cubature directions gc%dir_hat
!> \param head_gamma γ(iω) of the finite-q head (zero where no direction fit has screening form)
!> \param w_zone w(iω) per frequency (Hartree)
! **************************************************************************************************
   SUBROUTINE gamma_corr_head_zone_remainder(bs_env, dir_weight, head_gamma, w_zone)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: dir_weight, head_gamma
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: w_zone

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_head_zone_remainder'
      INTEGER, PARAMETER                                 :: n_radial = 64, n_sigma = 3
      REAL(KIND=dp), DIMENSION(n_sigma), PARAMETER :: &
         sigma_frac = [0.125_dp, 0.1_dp, 1.0_dp/12.0_dp]
      REAL(KIND=dp), PARAMETER                           :: q_cut_sigma = 6.0_dp

      INTEGER                                            :: handle, i_dim, i_dir, i_node, i_s, i_w, &
                                                            j_s, n1, n2, n3, n_dir, n_t
      INTEGER, DIMENSION(3)                              :: n_max
      REAL(KIND=dp)                                      :: a_dir, a_q, damp, g2, lagrange, q, &
                                                            q_cut, r_min, sigma, w_q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: nodes, rule
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: D
      REAL(KIND=dp), DIMENSION(3)                        :: g_hat, g_vec, r_vec
      REAL(KIND=dp), DIMENSION(n_sigma)                  :: s2
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_t = SIZE(w_zone)
      n_dir = SIZE(dir_weight)

      ! R = shortest supercell lattice vector Σ_i n_i a_i, n ∈ {−1, 0, 1}³ \ {0}, a_i = hmat(:, i)
      r_min = HUGE(1.0_dp)
      DO n1 = -1, 1
         DO n2 = -1, 1
            DO n3 = -1, 1
               IF (n1 == 0 .AND. n2 == 0 .AND. n3 == 0) CYCLE
               r_vec = MATMUL(bs_env%hmat, REAL([n1, n2, n3], dp))
               r_min = MIN(r_min, NORM2(r_vec))
            END DO
         END DO
      END DO

      ALLOCATE (nodes(n_radial), rule(n_radial))
      CALL gauss_legendre(n_radial, nodes, rule)
      ALLOCATE (D(n_t, n_sigma))
      D(:, :) = 0.0_dp

      DO i_s = 1, n_sigma
         sigma = sigma_frac(i_s)*r_min
         s2(i_s) = sigma**2
         ! e^{−σ²q²} < 3e-16 beyond q_cut = 6/σ
         q_cut = q_cut_sigma/sigma

         ! (a) ∫ d³q/(2π)³ W_h(q) e^{−σ²q²} = (1/(2π)³) Σ_q̂ w(q̂) ∫_0^{q_cut} dq q² W_h(q q̂) e^{−σ²q²},
         !     q² W_h = −16π² a_q/(1 + 4π a_q)
         DO i_dir = 1, n_dir
            DO i_w = 1, n_t
               ! a(q̂, iω) = q̂ᵀ S q̂
               a_dir = DOT_PRODUCT(gc%dir_hat(:, i_dir), &
                                   MATMUL(gc%S_head(:, :, i_w), gc%dir_hat(:, i_dir)))
               DO i_node = 1, n_radial
                  ! node map x ∈ [−1, 1] → q = q_cut (x + 1)/2, weight q_cut/2 · w_GL
                  q = 0.5_dp*q_cut*(nodes(i_node) + 1.0_dp)
                  w_q = 0.5_dp*q_cut*rule(i_node)
                  a_q = a_dir/(1.0_dp + head_gamma(i_w)*q*q)
                  D(i_w, i_s) = D(i_w, i_s) - dir_weight(i_dir)*w_q*16.0_dp*pi**2* &
                                a_q/(1.0_dp + fourpi*a_q)*EXP(-(sigma*q)**2)/twopi**3
               END DO
            END DO
         END DO

         ! (b) −(1/Ω) Σ_{0<|g|<q_cut} W_h(g) e^{−σ²g²} over g = Σ_i n_i b_i. n_i = a_i·g/(2π), so
         !     |n_i| ≤ q_cut |a_i|/(2π) bounds the integer box.
         DO i_dim = 1, 3
            n_max(i_dim) = CEILING(q_cut*NORM2(bs_env%hmat(:, i_dim))/twopi)
         END DO
         DO n1 = -n_max(1), n_max(1)
            DO n2 = -n_max(2), n_max(2)
               DO n3 = -n_max(3), n_max(3)
                  IF (n1 == 0 .AND. n2 == 0 .AND. n3 == 0) CYCLE
                  g_vec = MATMUL(gc%b_recip, REAL([n1, n2, n3], dp))
                  g2 = DOT_PRODUCT(g_vec, g_vec)
                  IF (g2 >= q_cut**2) CYCLE
                  g_hat = g_vec/SQRT(g2)
                  damp = EXP(-sigma*sigma*g2)
                  DO i_w = 1, n_t
                     ! a_g = ĝᵀ S ĝ/(1 + γ g²);  −W_h(g)/Ω = +(4π/g²) 4π a_g/(1 + 4π a_g)/Ω
                     a_q = DOT_PRODUCT(g_hat, MATMUL(gc%S_head(:, :, i_w), g_hat))/ &
                           (1.0_dp + head_gamma(i_w)*g2)
                     D(i_w, i_s) = D(i_w, i_s) + (fourpi/g2)*fourpi*a_q/(1.0_dp + fourpi*a_q)* &
                                   damp/gc%measure
                  END DO
               END DO
            END DO
         END DO
      END DO

      ! w = D(σ → 0): value at σ² = 0 of the quadratic through (σ_k², D_k), Lagrange form
      !   w = Σ_k D_k Π_{j≠k} (0 − σ_j²)/(σ_k² − σ_j²)
      w_zone(:) = 0.0_dp
      DO i_s = 1, n_sigma
         lagrange = 1.0_dp
         DO j_s = 1, n_sigma
            IF (j_s /= i_s) lagrange = lagrange*(-s2(j_s))/(s2(i_s) - s2(j_s))
         END DO
         w_zone(:) = w_zone(:) + lagrange*D(:, i_s)
      END DO

      DEALLOCATE (D, nodes, rule)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_head_zone_remainder

! **************************************************************************************************
!> \brief Bulk cell: head, head-dipole and wing corrections of the Γ point as orbital-diagonal
!>        shifts ΔΣ^c_n(±i|τ|), from the orbital-pair kernel
!>
!>          Ŵ_nm(iω) = w(iω) δ_nm + Σ_αβ B_αβ(iω) d_nm,α d_nm,β − 2√Ω Σ_αβ A_αβ(iω) V_nm,α(iω) d_nm,β
!>                   + Ω Σ_αβ A_αβ(iω) V_nm,α(iω) V_nm,β(iω),
!>
!>        with w, A and B of gamma_corr_direction_coefficients (the B term only without TAIL_3D),
!>        the dipoles d_nm = ⟨n|[H,r]|m⟩/(ε_n − ε_m) (zero in degenerate manifolds, clipped to
!>        |d_nm| ≤ 1/q_max) and the wing potential V_nm,α = Σ_ℓ ψ_n(r_ℓ) ψ_m(r_ℓ) w̃_ℓα, w̃ = Z M⁺_c T.
!>        With Ŵ(iτ) from the cosine transform of the minimax grid,
!>
!>          ΔΣ^c_n(−i|τ|) = −Σ_{m occ} e^{−|τ|(ε_F − ε_m)} Ŵ_nm(iτ),
!>          ΔΣ^c_n(+i|τ|) = +Σ_{m vir} e^{−|τ|(ε_m − ε_F)} Ŵ_nm(iτ).
!> \param bs_env band-structure environment: head potential w, tensors A and B, T_RI, q_max,
!>        MOs, levels and the grid maps Φ = mat_phi_mu_l, Z = mat_Z_lP
!> \param qs_env QS environment, for the commutator matrices of the dipoles
!> \param dSigma ΔΣ^c_n(±i|τ_j|), (n_mo, n_t, sector 1 = + / 2 = −, n_spin); allocated here
! **************************************************************************************************
   SUBROUTINE gamma_corr_sigma_shifts(bs_env, qs_env, dSigma)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: dSigma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_sigma_shifts'
      REAL(KIND=dp), PARAMETER                           :: eps_degenerate = 1.0E-5_dp

      INTEGER                                            :: handle, i_col, i_dir, i_row, i_sector, &
                                                            i_spin, i_t, i_w, j_dir, m, n, n_ao, &
                                                            n_grid, n_mo, n_occ, n_RI, n_t, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: d_max, d_norm, delta, sqrt_V, &
                                                            start_time, tau, tau_E, w_hat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: g_branch, y_loc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: u_RI, w_grid
      REAL(KIND=dp), DIMENSION(3)                        :: d_nm, v_nm
      TYPE(cp_fm_struct_type), POINTER                   :: struct_mo_mo
      TYPE(cp_fm_type)                                   :: fm_ao_ao, fm_ao_mo, fm_What
      TYPE(cp_fm_type), DIMENSION(3)                     :: fm_d, fm_V
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments
      TYPE(dbcsr_type)                                   :: mat_Omega, mat_phi_scaled
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)
      start_time = m_walltime()

      gc => bs_env%ri_rs%gamma_corr
      CPASSERT(gc%ready .AND. .NOT. gc%slab .AND. ALLOCATED(gc%head_pot_freq))
      CPASSERT(ALLOCATED(gc%T_RI) .AND. ALLOCATED(gc%A_wing))
      n_ao = bs_env%n_ao
      n_RI = bs_env%n_RI
      n_grid = bs_env%ri_rs%n_grid_points
      n_t = bs_env%num_time_freq_points
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = n_ao
      ! √Ω of the prefactor 2√Ω of the wing-dipole term (Ω = gc%measure)
      sqrt_V = SQRT(gc%measure)

      ! Work matrices: n_MO x n_MO for d, V and Ŵ; n_AO x n_AO and n_AO x n_MO intermediates
      NULLIFY (struct_mo_mo)
      CALL cp_fm_struct_create(struct_mo_mo, &
                               context=bs_env%fm_mo_coeff_Gamma(1)%matrix_struct%context, &
                               para_env=bs_env%fm_mo_coeff_Gamma(1)%matrix_struct%para_env, &
                               nrow_global=n_mo, ncol_global=n_mo)
      DO i_dir = 1, 3
         CALL cp_fm_create(fm_d(i_dir), struct_mo_mo)
         CALL cp_fm_create(fm_V(i_dir), struct_mo_mo)
      END DO
      CALL cp_fm_create(fm_What, struct_mo_mo)
      CALL cp_fm_create(fm_ao_ao, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_ao_mo, bs_env%fm_mo_coeff_Gamma(1)%matrix_struct)
      CALL cp_fm_get_info(fm_What, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      ALLOCATE (u_RI(n_RI, 3), w_grid(n_grid, 3), g_branch(n_mo), y_loc(n_mo), &
                dSigma(n_mo, n_t, 2, bs_env%n_spin))

      ! Scratch and result matrices of Φᵀ diag(w̃_α) Φ
      CALL dbcsr_create(mat_phi_scaled, template=bs_env%ri_rs%mat_phi_mu_l)
      CALL dbcsr_create(mat_Omega, template=bs_env%mat_ao_ao%matrix)

      ! [H, r^α] in the AO basis, shared by the spins.
      NULLIFY (moments)
      CALL gamma_corr_commutator_matrices(qs_env, moments)

      ! d_nm needs the Kohn-Sham eigenvalues of the H in the commutator: copy them on the first call,
      ! before evGW0 replaces eigenval_scf_Gamma by QP energies (g_m(τ) uses the current ones)
      IF (.NOT. ALLOCATED(gc%eigenval_dipole)) THEN
         ALLOCATE (gc%eigenval_dipole, SOURCE=bs_env%eigenval_scf_Gamma)
      END IF

      DO i_spin = 1, bs_env%n_spin
         n_occ = bs_env%n_occ(i_spin)

         ! ---------------------------------------------------------------------------------
         ! 1. Dipoles d_nm = ⟨n|[H,r^α]|m⟩/(ε_n − ε_m) of the retained MOs, clipped to |d_nm| ≤ 1/q_max
         ! ---------------------------------------------------------------------------------
         DO i_dir = 1, 3
            CALL cp_fm_set_all(fm_ao_ao, 0.0_dp)
            CALL copy_dbcsr_to_fm(moments(i_dir)%matrix, fm_ao_ao)
            CALL parallel_gemm('N', 'N', n_ao, n_mo, n_ao, 1.0_dp, fm_ao_ao, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), 0.0_dp, fm_ao_mo)
            CALL parallel_gemm('T', 'N', n_mo, n_mo, n_ao, 1.0_dp, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), fm_ao_mo, 0.0_dp, fm_d(i_dir))
            DO i_col = 1, ncol_local
               m = col_indices(i_col)
               DO i_row = 1, nrow_local
                  n = row_indices(i_row)
                  ! Kohn-Sham ε_n − ε_m (frozen at the first call); d_nm = 0 for degenerate pairs
                  delta = gc%eigenval_dipole(n, i_spin) - gc%eigenval_dipole(m, i_spin)
                  IF (ABS(delta) > eps_degenerate) THEN
                     fm_d(i_dir)%local_data(i_row, i_col) = fm_d(i_dir)%local_data(i_row, i_col)/delta
                  ELSE
                     fm_d(i_dir)%local_data(i_row, i_col) = 0.0_dp
                  END IF
               END DO
            END DO
         END DO
         ! Clip |d_nm| at d_max = 1/q_max along the direction of d_nm.
         d_max = 1.0_dp/gc%q_max
         DO i_col = 1, ncol_local
            DO i_row = 1, nrow_local
               d_norm = SQRT(fm_d(1)%local_data(i_row, i_col)**2 + fm_d(2)%local_data(i_row, i_col)**2 + &
                             fm_d(3)%local_data(i_row, i_col)**2)
               IF (d_norm > d_max) THEN
                  DO i_dir = 1, 3
                     fm_d(i_dir)%local_data(i_row, i_col) = fm_d(i_dir)%local_data(i_row, i_col)*d_max/d_norm
                  END DO
               END IF
            END DO
         END DO

         dSigma(:, :, :, i_spin) = 0.0_dp
         DO i_w = 1, n_t
            ! ------------------------------------------------------------------------------
            ! 2. Wing potential on the grid, w̃_ℓα = Σ_Q Z_ℓQ [M⁺_c T(iω)]_Qα, and its orbital
            !    matrix V_nm,α = Σ_ℓ ψ_n(r_ℓ) ψ_m(r_ℓ) w̃_ℓα = [Cᵀ (Φᵀ diag(w̃_α) Φ) C]_nm.
            ! ------------------------------------------------------------------------------
            ! u = M⁺_c T(iω), replicated n_RI x 3
            CALL gamma_corr_minv_c_rep(bs_env, gc%T_RI(:, :, i_w), u_RI)
            CALL gamma_corr_project_RI_to_grid(bs_env, u_RI, w_grid)
            DO i_dir = 1, 3
               CALL dbcsr_copy(mat_phi_scaled, bs_env%ri_rs%mat_phi_mu_l)
               CALL dbcsr_scale_by_vector(mat_phi_scaled, w_grid(:, i_dir), side='left')
               CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_phi_scaled, bs_env%ri_rs%mat_phi_mu_l, 0.0_dp, &
                                   mat_Omega, filter_eps=bs_env%eps_filter)
               CALL cp_fm_set_all(fm_ao_ao, 0.0_dp)
               CALL copy_dbcsr_to_fm(mat_Omega, fm_ao_ao)
               CALL parallel_gemm('N', 'N', n_ao, n_mo, n_ao, 1.0_dp, fm_ao_ao, &
                                  bs_env%fm_mo_coeff_Gamma(i_spin), 0.0_dp, fm_ao_mo)
               CALL parallel_gemm('T', 'N', n_mo, n_mo, n_ao, 1.0_dp, &
                                  bs_env%fm_mo_coeff_Gamma(i_spin), fm_ao_mo, 0.0_dp, fm_V(i_dir))
            END DO

            ! ------------------------------------------------------------------------------
            ! 3. Ŵ_nm(iω) on the local elements.
            ! ------------------------------------------------------------------------------
            DO i_col = 1, ncol_local
               m = col_indices(i_col)
               DO i_row = 1, nrow_local
                  n = row_indices(i_row)
                  DO i_dir = 1, 3
                     d_nm(i_dir) = fm_d(i_dir)%local_data(i_row, i_col)
                     v_nm(i_dir) = fm_V(i_dir)%local_data(i_row, i_col)
                  END DO
                  w_hat = 0.0_dp
                  DO j_dir = 1, 3
                     DO i_dir = 1, 3
                        ! B:dd − 2√Ω A:Vd + Ω A:VV (the body screens the pair dipole, d → d − √Ω V);
                        ! with TAIL_3D, B:dd is part of the sphere-tail term in W
                        IF (.NOT. gc%tail_3d) THEN
                           w_hat = w_hat + gc%B_hd(i_dir, j_dir, i_w)*d_nm(i_dir)*d_nm(j_dir)
                        END IF
                        w_hat = w_hat - &
                                2.0_dp*sqrt_V*gc%A_wing(i_dir, j_dir, i_w)*v_nm(i_dir)*d_nm(j_dir) + &
                                gc%measure*gc%A_wing(i_dir, j_dir, i_w)*v_nm(i_dir)*v_nm(j_dir)
                     END DO
                  END DO
                  IF (n == m) w_hat = w_hat + gc%head_pot_freq(i_w)
                  fm_What%local_data(i_row, i_col) = w_hat
               END DO
            END DO

            ! ------------------------------------------------------------------------------
            ! 4. Branch sums per time point, with g_m(τ) = e^{−|τ(ε_m − ε_F)|}:
            !      ΔΣ_n(+i|τ_t|) += γ̃^cos_tw Σ_{m vir} g_m Ŵ_nm(iω_w)   (sector 1)
            !      ΔΣ_n(−i|τ_t|) −= γ̃^cos_tw Σ_{m occ} g_m Ŵ_nm(iω_w)   (sector 2)
            ! ------------------------------------------------------------------------------
            DO i_t = 1, n_t
               tau = bs_env%time_frequency_grid%imaginary_time(i_t)
               DO i_sector = 1, 2
                  DO m = 1, n_mo
                     g_branch(m) = 0.0_dp
                     IF ((i_sector == 1 .AND. m > n_occ) .OR. (i_sector == 2 .AND. m <= n_occ)) THEN
                        tau_E = ABS(tau*(bs_env%eigenval_scf_Gamma(m, i_spin) - bs_env%e_fermi(i_spin)))
                        IF (tau_E < bs_env%stabilize_exp) g_branch(m) = EXP(-tau_E)
                     END IF
                  END DO
                  y_loc(:) = 0.0_dp
                  DO i_col = 1, ncol_local
                     m = col_indices(i_col)
                     IF (g_branch(m) == 0.0_dp) CYCLE
                     DO i_row = 1, nrow_local
                        n = row_indices(i_row)
                        y_loc(n) = y_loc(n) + fm_What%local_data(i_row, i_col)*g_branch(m)
                     END DO
                  END DO
                  CALL bs_env%para_env%sum(y_loc)
                  dSigma(:, i_t, i_sector, i_spin) = dSigma(:, i_t, i_sector, i_spin) + &
                                                     MERGE(1.0_dp, -1.0_dp, i_sector == 1)* &
                                                     bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, i_w)* &
                                                     y_loc(:)
               END DO
            END DO
         END DO ! i_w
      END DO ! i_spin

      CALL dbcsr_deallocate_matrix_set(moments)
      CALL dbcsr_release(mat_phi_scaled)
      CALL dbcsr_release(mat_Omega)
      DO i_dir = 1, 3
         CALL cp_fm_release(fm_d(i_dir))
         CALL cp_fm_release(fm_V(i_dir))
      END DO
      CALL cp_fm_release(fm_What)
      CALL cp_fm_release(fm_ao_ao)
      CALL cp_fm_release(fm_ao_mo)
      CALL cp_fm_struct_release(struct_mo_mo)
      DEALLOCATE (u_RI, w_grid, g_branch, y_loc)
      ! head_pot_*, A_wing, B_hd, T_RI and eigenval_dipole stay allocated for the next evGW0 iteration

      ! One timing line for the whole 3D correction (maps, response, tensors, tail, ΔΣ^c)
      gc%wall_time = gc%wall_time + m_walltime() - start_time
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A,T59,A,F7.1,A)') &
         'Computed Γ-point correction of W (3D bulk),', 'Execution time', gc%wall_time, ' s'

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_sigma_shifts

! **************************************************************************************************
!> \brief w_ℓα = Σ_Q Z_ℓQ u_Qα for a three-column RI quantity (the transpose of
!>        gamma_corr_project_grid_to_RI), by walking the local blocks of Z.
!> \param bs_env band-structure environment, for mat_Z_lP and the communicator
!> \param u_RI u_Qα, replicated, n_RI x 3
!> \param w_grid result, replicated over grid rows, n_grid x 3
! **************************************************************************************************
   SUBROUTINE gamma_corr_project_RI_to_grid(bs_env, u_RI, w_grid)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: u_RI
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: w_grid

      INTEGER                                            :: col, col_start, n_col, n_row, row, &
                                                            row_start
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      w_grid(:, :) = 0.0_dp
      CALL dbcsr_iterator_start(iter, bs_env%ri_rs%mat_Z_lP)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_offset=row_start, col_offset=col_start)
         n_row = SIZE(block, 1)
         n_col = SIZE(block, 2)
         w_grid(row_start:row_start + n_row - 1, :) = w_grid(row_start:row_start + n_row - 1, :) + &
                                                      MATMUL(block, u_RI(col_start:col_start + n_col - 1, :))
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL bs_env%para_env%sum(w_grid)

   END SUBROUTINE gamma_corr_project_RI_to_grid

END MODULE gw_gamma_correction
