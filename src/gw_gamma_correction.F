!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Γ-point head-and-wing correction of W for the Γ-only large-cell RI-RS GW code (2D slabs).
!>
!>        **Why this exists.** The large-cell code reconstructs χ(k), ε(k) and W(k) from Γ on a
!>        k-mesh of the supercell's own small Brillouin zone, and takes V(k) from a finite
!>        real-space image sum. An EVEN mesh is half-shifted and therefore never samples q = 0; an
!>        ODD mesh does, because compute_xkp gives xkp = (2·ix − n − 1)/(2n). At q = 0 the image
!>        sum of the bare Coulomb operator between two CHARGED auxiliary Gaussians diverges with
!>        the summation range, so `KPOINTS_W 1 1 1` does not merely sample badly, it produces a
!>        large spurious additive W. Measured on a 392-atom h-BN slab: Σ^c of the valence band
!>        maximum +2.195 → −24.043 eV and of the conduction band minimum −3.283 → +20.237 eV
!>        (opposite signs because Σ^c = −G·W and the two branches of G carry opposite sign),
!>        while Σ^x and v^xc stayed bit-identical.
!>
!>        **What it does.** The charge channel is split off the auxiliary basis and carried
!>        analytically instead of being summed in real space:
!>
!>          c_P   = ∫ dr ϕ_P0(r)                                     auxiliary charge
!>          Π_PQ  = δ_PQ − c_P c_Q / Σ_S |c_S|²                      neutral projector
!>          V_b   = Π V(0) Π                                          neutral body
!>          ϕ_0   = 1/√A,  v(q) = 2π/q                                separated charge mode
!>
!>        Between two neutral localized functions the interaction decays as R⁻³, so the 2D lattice
!>        sum of V_b converges absolutely where the unprojected one does not. V_b is the body of
!>        the CORRECTION only (Eqs. 66-67). The regular W path keeps the full V(0) with just its
!>        divergent G = 0 constant removed, V(0) − D_box c cᵀ (gwp42 Sec. II, "the singular G = 0
!>        term omitted"); see gamma_corr_subtract_divergent_head for why the projection must not
!>        be used there. The dielectric head and wings in the extended basis {0, 1, …, N_RI} are
!>        then integrated over the WHOLE small Brillouin zone of the supercell. Nothing is
!>        subtracted, unlike the full-k code of J. Wilhelm (cp2k#6447, gwp42 Eq. 78), because the
!>        regular path is Γ only and its sampled head was removed exactly.
!>
!>        **Γ only means Γ only.** setup_kpoints_chi_eps_W drops the extrapolation mesh for this
!>        case. With it, "KPOINTS_W 1 1 1" would be a 2x2 calculation in disguise: Γ enters the
!>        s–s elements of W with the extrapolation weight −2.414 and every other element with
!>        weight 0 (MIC_contribution_from_ikp), so the Γ-point W would be amplified and sign
!>        flipped where it matters and ignored everywhere else.
!>
!>        Equations follow J. Wilhelm, *Periodic real-space RI for full-k GW band structures*,
!>        Sec. VI; equation numbers below refer to that manuscript. The derivation of the Γ-only
!>        specialization is in GAMMA_W_CORRECTION_PLAN.md.
!>
!>        **Scope.** 2D slabs (PERIODIC XY) with the Γ-only RI-RS implementation. χ, ε and W stay
!>        at one real k-point, so none of this needs complex arithmetic, and the minimum-image
!>        clipping of χ (cp_cfm_power) has nothing to repair: χ(Γ) is a Gram matrix of
!>        occupied-virtual products and exactly positive semidefinite.
!>
!>        **Call sequence.** Read the routines in the order in which the GW driver calls them:
!>
!>          ri_rs_setup              gamma_corr_init                 decide, set A, b_1, b_2
!>          gw_ri_rs_large_cell_     gamma_corr_build_coulomb_maps   c, Π, V_b, F_density, F_extern
!>          gamma, step 4b           gamma_corr_response             S_αβ(iω), T_Pα(iω)
!>                                   gamma_corr_tabulate             K_h(R), K_w(R) tables
!>          gw_compute_V_and_W       gamma_corr_subtract_divergent_head   V(0) − D_box c cᵀ
!>          (compute_W_MIC)          gamma_corr_add_to_W             W^MIC_PQ(iτ) += ΔW_PQ(iτ)
!>
!>        **Where the correction is applied.** In the RI basis, on the auxiliary charges c_P and
!>        the atom-pair separations, and BEFORE the M⁻¹ W M⁻¹ contraction, so that it reaches
!>        the grid pairs through the same fit Z M⁻¹ as the regular W. An earlier version applied
!>        the same formula per grid PAIR (Z_ℓ0 Z_ℓ'0 K_h(R_ℓℓ') + wings, inside the Σ^c panel
!>        loop); that is numerically unusable, because the RI-RS interpolation functions far from
!>        the atoms carry fit weights of 10^4 .. 10^8 of alternating sign that cancel only in
!>        exact linear functionals, never under a kernel that varies from point to point. See
!>        gamma_corr_add_to_W.
!>
!>        **Units.** Atomic units throughout: lengths in bohr, areas in bohr², energies and
!>        frequencies in Hartree, imaginary times in Hartree⁻¹. Grid index ℓ (ℓ' for the second
!>        point of a pair) runs over the n_grid RI-RS grid points r_ℓ; P, Q run over the n_RI
!>        auxiliary (RI_AUX) functions; i over occupied and a over empty Kohn-Sham states; α, β
!>        over the Cartesian directions x, y, z.
!> \author Ritaj Tyagi
! **************************************************************************************************
MODULE gw_gamma_correction
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE input_constants,                 ONLY: do_potential_short,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE ri_rs_types,                     ONLY: gamma_corr_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_gamma_correction'

   ! Relative eigenvalue cut used to invert the neutral body V_b = Π V(0) Π on its range.
   ! This is a NUMERICAL null-space threshold and must not be confused with REGULARIZATION_RI,
   ! which is a Tikhonov shift of the RI metric and whose default is 1E-02: used as a relative
   ! eigenvalue cut that discards every mode below 1 % of λ_max, which on a 4704-function
   ! auxiliary basis throws away about 85 % of the neutral subspace. The only modes that must go
   ! are the exact charge direction and whatever linear dependence the RI basis itself carries,
   ! so the cut sits just above the level to which the projection Π can leave the (much larger)
   ! charge eigenvalue contaminating V_b, i.e. a few orders of magnitude above ε_machine.
   REAL(KIND=dp), PARAMETER, PRIVATE :: eps_null_V_b = 1.0E-10_dp

   PUBLIC :: gamma_corr_requested, &
             gamma_corr_init, &
             gamma_corr_build_coulomb_maps, &
             gamma_corr_transition_dipoles, &
             gamma_corr_response, &
             gamma_corr_tabulate, &
             gamma_corr_add_to_W, &
             gamma_corr_subtract_divergent_head, &
             gamma_corr_auxiliary_charges, &
             gamma_corr_metric_head, &
             gamma_corr_slab_geometry, &
             gamma_corr_symmetric_power

CONTAINS

! **************************************************************************************************
!> \brief Whether the Γ-point correction applies to this calculation.
!>
!>        It is meaningful only when the mesh for χ, ε, W actually samples q = 0, which with the
!>        ±k pairing of setup_kpoints_chi_eps_W happens for the 1x1x1 mesh alone (see below), and
!>        only for a 2D slab, which is the only geometry whose head/wing construction is
!>        implemented here.
!> \param bs_env band-structure environment; reads periodic and nkp_grid_chi_eps_W_orig only
!> \return TRUE if the correction should be built
! **************************************************************************************************
   FUNCTION gamma_corr_requested(bs_env) RESULT(do_corr)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL                                            :: do_corr

      INTEGER                                            :: i_dim

      ! Default: no correction. Each test below returns early with this value.
      do_corr = .FALSE.

      ! Only PERIODIC XY: the 3D construction needs the Ewald ΔV and the 1/q² head instead.
      IF (.NOT. ALL(bs_env%periodic == [1, 1, 0])) RETURN

      ! q = 0 is actually SAMPLED only for the 1x1x1 mesh. setup_kpoints_chi_eps_W keeps
      ! nkp_orig = n1*n2*n3/2 points of an n1 x n2 x n3 mesh, one per +-k pair, and for an odd
      ! n > 1 that integer division drops the self-paired Gamma point: a 3x3 mesh is sampled
      ! at its 8 non-Gamma points. Only MAX(1/2, 1) = 1 keeps Gamma, as the single point.
      DO i_dim = 1, 3
         IF (bs_env%nkp_grid_chi_eps_W_orig(i_dim) /= 1) RETURN
      END DO

      ! PERIODIC XY and the 1x1x1 mesh: the regular path samples q = 0 only, so correct it.
      do_corr = .TRUE.

   END FUNCTION gamma_corr_requested

! **************************************************************************************************
!> \brief Layer area A = |a_1 × a_2| and the two in-plane reciprocal vectors b_1, b_2 with
!>        b_i · a_j = 2π δ_ij, for an arbitrary (possibly non-orthogonal) XY cell.
!>
!>        The vacuum height never enters: a slab's Fourier transform averages over a layer cell,
!>        not over a 3D cell that mostly contains vacuum. This is the single substitution that
!>        turns the 3D formulas into the 2D ones, 1/√Ω → 1/√A (gwp42 Eqs. 61, 62).
!> \param bs_env band-structure environment; only the cell matrix hmat is read
!> \param area |a_1 × a_2| in bohr²
!> \param b_recip the two in-plane reciprocal vectors, b_recip(1:3, 1:2), Cartesian, in 1/bohr
! **************************************************************************************************
   SUBROUTINE gamma_corr_slab_geometry(bs_env, area, b_recip)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: area
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(OUT)        :: b_recip

      REAL(KIND=dp), DIMENSION(3)                        :: a_1, a_2, normal

      ! The columns of hmat are the Cartesian cell vectors (bohr); a_1, a_2 span the layer.
      a_1 = bs_env%hmat(:, 1)
      a_2 = bs_env%hmat(:, 2)

      ! normal = a_1 × a_2; its length is the layer area A, its direction the layer normal n.
      normal = [a_1(2)*a_2(3) - a_1(3)*a_2(2), &
                a_1(3)*a_2(1) - a_1(1)*a_2(3), &
                a_1(1)*a_2(2) - a_1(2)*a_2(1)]
      area = SQRT(SUM(normal**2))
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

   END SUBROUTINE gamma_corr_slab_geometry

! **************************************************************************************************
!> \brief Zero-wavevector Fourier value of the RI metric kernel,
!>
!>          ∫ dr v_ω(r) = π/ω²                for the erfc-attenuated metric
!>          ∫ dr v_trunc(r) = 2π r_c²         for the truncated Coulomb metric
!>
!>        The charge-channel right-hand side d_ℓ0 is multiplied by this before the fit and the
!>        fitted coefficients divided by it afterwards. That is a no-op for an exact solve but
!>        not for the regularized one: the fixed Tikhonov λ is not scale invariant, so the
!>        scaling puts the charge column on the same footing as the Gaussian columns d_ℓP.
!> \param bs_env band-structure environment; only the RI metric (ri_metric) is read
!> \return the metric head, in bohr² (∫ d³r of a kernel with units 1/bohr)
! **************************************************************************************************
   FUNCTION gamma_corr_metric_head(bs_env) RESULT(head)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: head

      SELECT CASE (bs_env%ri_metric%potential_type)
      CASE (do_potential_short)
         ! v_ω(r) = erfc(ω r)/r:  ∫ d³r erfc(ω r)/r = 4π ∫ dr r erfc(ω r) = π/ω².
         head = pi/bs_env%ri_metric%omega**2
      CASE (do_potential_truncated)
         ! v(r) = θ(r_c − r)/r:  ∫ d³r v(r) = 4π ∫_0^{r_c} dr r = 2π r_c².
         head = twopi*bs_env%ri_metric%cutoff_radius**2
      CASE DEFAULT
         ! Other metrics are not supported; the bare Coulomb one, for instance, has no finite
         ! q = 0 value.
         CPABORT('Gamma-point correction needs an erfc or truncated-Coulomb RI metric')
      END SELECT

      CPASSERT(head > 0.0_dp)

   END FUNCTION gamma_corr_metric_head

! **************************************************************************************************
!> \brief Charges of the RI_AUX basis functions, c_P = ∫ dr ϕ_P0(r).
!>
!>        For a Cartesian primitive x^i y^j z^k exp(−α r²) the integral factorizes into
!>        one-dimensional moments
!>
!>          m_0 = √(π/α),   m_l = (l−1) m_{l−2} / (2α)   for even l,   m_l = 0 for odd l,
!>
!>        so c for the primitive is m_i m_j m_k, and the spherical contraction coefficients
!>        (basis%sphi) carry it to the actual basis function. "Charge" here means the integral,
!>        not the norm: a function with positive and negative parts can integrate to zero, and
!>        ϕ_1 − (c_1/c_2) ϕ_2 is exactly neutral whenever c_2 ≠ 0.
!> \param qs_env QS environment, for the particle set and the RI_AUX basis of every kind
!> \param charge c_P, one entry per RI_AUX basis function, in the global RI ordering
! **************************************************************************************************
   SUBROUTINE gamma_corr_auxiliary_charges(qs_env, charge)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: charge

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_auxiliary_charges'

      INTEGER                                            :: first, handle, i_atom, i_co, i_l, i_pgf, &
                                                            i_set, i_sgf, l_max, n_cart, offset
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: moments, primitive
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set)

      charge(:) = 0.0_dp
      ! first = number of RI functions on all atoms before i_atom, i.e. the global RI offset of
      ! the current atom. Atoms are visited in particle order, which is the global RI ordering.
      first = 0

      DO i_atom = 1, SIZE(particle_set)
         CALL get_qs_kind(qs_kind_set(particle_set(i_atom)%atomic_kind%kind_number), &
                          basis_set=basis, basis_type='RI_AUX')

         DO i_set = 1, basis%nset
            ! Primitive layout of this set, as sphi expects it: for every primitive exponent
            ! i_pgf a block of n_cart = ncoset(l_max) Cartesian components (all l <= l_max),
            ! primitive((i_pgf − 1)·n_cart + i_co). Components with l < lmin stay zero.
            l_max = basis%lmax(i_set)
            n_cart = ncoset(l_max)
            ALLOCATE (moments(0:l_max), primitive(n_cart*basis%npgf(i_set)))
            primitive(:) = 0.0_dp

            DO i_pgf = 1, basis%npgf(i_set)
               ! Gaussian exponent α of this primitive (bohr⁻²).
               alpha = basis%zet(i_pgf, i_set)

               ! One-dimensional Gaussian moments; odd orders vanish and stay zero.
               moments(:) = 0.0_dp
               moments(0) = SQRT(pi/alpha)
               DO i_l = 2, l_max, 2
                  moments(i_l) = REAL(i_l - 1, dp)*moments(i_l - 2)/(2.0_dp*alpha)
               END DO

               ! c = m_i m_j m_k for the Cartesian component (i,j,k) = indco(:, i_co).
               DO i_co = ncoset(basis%lmin(i_set) - 1) + 1, n_cart
                  primitive((i_pgf - 1)*n_cart + i_co) = PRODUCT(moments(indco(:, i_co)))
               END DO
            END DO

            ! Contract to the spherical functions of the set: c_P = Σ_k sphi(k, P) c_k, where
            ! k runs over the primitive Cartesian layout above and P over the set's spherical
            ! functions, numbered within the atom from first_sgf(1, i_set).
            offset = basis%first_sgf(1, i_set)
            DO i_sgf = offset, offset + basis%nsgf_set(i_set) - 1
               charge(first + i_sgf) = SUM(basis%sphi(1:SIZE(primitive), i_sgf)*primitive)
            END DO

            DEALLOCATE (moments, primitive)
         END DO

         ! Advance the global RI offset by the number of RI functions on this atom.
         first = first + basis%nsgf
      END DO

      ! Every RI function must have been visited exactly once.
      CPASSERT(first == SIZE(charge))

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_auxiliary_charges

! **************************************************************************************************
!> \brief Symmetric matrix power A^p by eigendecomposition, discarding modes below a relative
!>        threshold.
!>
!>        V_b = Π V(0) Π is singular by construction: the charge direction c is an exact null
!>        vector, and the neutral functions ϕ̄_Q = Σ_P ϕ_P Π_PQ are linearly dependent. Matrix
!>        powers must therefore act only on the nonzero neutral subspace, which is what the
!>        relative eigenvalue cut does. The threshold is eps_null_V_b, NOT REGULARIZATION_RI; see
!>        the comment on the parameter for why the two must not be confused.
!> \param matrix on entry the symmetric matrix, on exit its p-th power; replicated
!> \param exponent p
!> \param relative_eps modes with λ < relative_eps · max λ are set to zero
!> \param n_dropped optional count of discarded modes
!> \param lambda_max optional largest eigenvalue, for reporting the conditioning
!> \param lambda_kept optional smallest retained eigenvalue, for reporting the conditioning
!> \param exponent_2 optional second exponent, evaluated from the same eigendecomposition
!> \param matrix_2 optional result for exponent_2
! **************************************************************************************************
   SUBROUTINE gamma_corr_symmetric_power(matrix, exponent, relative_eps, n_dropped, &
                                         lambda_max, lambda_kept, exponent_2, matrix_2)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: exponent, relative_eps
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_dropped
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: lambda_max, lambda_kept
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: exponent_2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: matrix_2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_symmetric_power'

      INTEGER                                            :: handle, i, info, n, n_drop
      REAL(KIND=dp)                                      :: cut
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scaled, vectors

      CALL timeset(routineN, handle)

      n = SIZE(matrix, 1)
      CPASSERT(SIZE(matrix, 2) == n)

      ! LAPACK workspace: 34·n is (block size 32 + 2)·n, the blocked optimum for dsyev.
      ALLOCATE (vectors(n, n), eigenvalues(n), work(MAX(1, 34*n)))
      ! Symmetrize first: a product like Π V Π is symmetric only up to rounding, and dsyev reads
      ! one triangle, so the explicit average makes the result independent of which one.
      vectors(:, :) = 0.5_dp*(matrix + TRANSPOSE(matrix))

      ! Full eigendecomposition A = U diag(λ) Uᵀ; U overwrites vectors, λ ascending.
      CALL dsyev('V', 'U', n, vectors, n, eigenvalues, work, SIZE(work), info)
      IF (info /= 0) CPABORT('Gamma-point correction: dsyev failed')

      ! Absolute cut λ_cut = relative_eps · λ_max. The MAX(…, 0) keeps it non-negative, so
      ! negative eigenvalues (rounding noise on a PSD matrix) are always dropped and
      ! λ^p is never evaluated for λ <= 0.
      cut = relative_eps*MAXVAL(eigenvalues)
      cut = MAX(cut, 0.0_dp)

      ! A^p = U diag(λ^p) U^T on the retained modes; discarded modes stay quenched for every p,
      ! including negative ones, which is what keeps V_b^(-1/2) finite on the null space.
      ALLOCATE (scaled(n, n))
      n_drop = 0
      DO i = 1, n
         IF (eigenvalues(i) <= cut) THEN
            ! Discarded mode: column of U·diag(λ^p) set to zero.
            scaled(:, i) = 0.0_dp
            n_drop = n_drop + 1
         ELSE
            ! Retained mode: column i of U·diag(λ^p), i.e. u_i λ_i^p.
            scaled(:, i) = vectors(:, i)*eigenvalues(i)**exponent
         END IF
      END DO

      ! matrix = [U diag(λ^p)] Uᵀ = Σ_{retained i} λ_i^p u_i u_iᵀ.
      CALL dgemm('N', 'T', n, n, n, 1.0_dp, scaled, n, vectors, n, 0.0_dp, matrix, n)

      ! A second power from the SAME eigendecomposition: the diagonalization is the expensive
      ! part, and the two powers always share the retained modes.
      IF (PRESENT(matrix_2)) THEN
         CPASSERT(PRESENT(exponent_2))
         DO i = 1, n
            IF (eigenvalues(i) <= cut) THEN
               scaled(:, i) = 0.0_dp
            ELSE
               scaled(:, i) = vectors(:, i)*eigenvalues(i)**exponent_2
            END IF
         END DO
         CALL dgemm('N', 'T', n, n, n, 1.0_dp, scaled, n, vectors, n, 0.0_dp, matrix_2, n)
      END IF

      IF (PRESENT(n_dropped)) n_dropped = n_drop
      ! The pair (λ_max, smallest retained λ) says whether the cut fell inside a genuine gap or
      ! sliced through a continuum of physical modes; eigenvalues come out of dsyev ascending.
      IF (PRESENT(lambda_max)) lambda_max = eigenvalues(n)
      ! Smallest retained λ: the first eigenvalue above the cut in the ascending list.
      IF (PRESENT(lambda_kept)) THEN
         lambda_kept = 0.0_dp
         DO i = 1, n
            IF (eigenvalues(i) > cut) THEN
               lambda_kept = eigenvalues(i)
               EXIT
            END IF
         END DO
      END IF

      DEALLOCATE (vectors, eigenvalues, work, scaled)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_symmetric_power

! **************************************************************************************************
!> \brief Build everything the correction needs from the Coulomb side: the charge vector, the
!>        neutral body V_b = Π V(0) Π, and the two grid maps.
!>
!>        Steps, with gwp42 equation numbers:
!>
!>        1. V(0)   from the real-space 2c lattice sum at q = 0 (Eq. 63). Unprojected this
!>                  diverges with SIZE_LATTICE_SUM for charged functions; only V_b is ever used.
!>        2. c_P    auxiliary charges, normalized (Eq. 64)
!>           Π      = Id − c cᵀ
!>        3. V_b    = Π V(0) Π, and V_b^(±1/2) on the neutral subspace (Eq. 65)
!>        4. M(0)   attenuated RI metric at Γ, and M^-1(0)
!>        5. F_density = M^-1(0) Π V_b^(1/2)          (Eq. 66, the response-side map T = Z F_density)
!>           F_extern  = M^-1(0) V(0) Π V_b^(-1/2)    (Eq. 67, the Coulomb-side map A = Z F_extern)
!>
!>        The two maps are deliberately distinct. The right projector in F_extern makes the
!>        SOURCE function neutral; its potential may still integrate to something nonzero against
!>        a charged test function on the left, and those finite charge/body couplings are
!>        physical. Multiplying by M^-1(0) turns potential integrals into fit coefficients and
!>        does not make a neutral source charged.
!>
!>        Everything here is real: at q = 0 the image sum of a real kernel over a centrosymmetric
!>        box has no imaginary part, and the metric is real too.
!>
!>        Memory: V(0), Π, V_b^(±1/2), M^-1 and the two maps are dense n_RI x n_RI arrays
!>        REPLICATED on every rank (about 8 of them at the peak), and every dgemm below is a
!>        serial per-rank product. This is affordable because it is done once per calculation.
!> \param qs_env QS environment, for the 2c Coulomb integrals and the RI_AUX basis
!> \param bs_env band-structure environment; fills bs_env%ri_rs%gamma_corr (charge, maps)
! **************************************************************************************************
   SUBROUTINE gamma_corr_build_coulomb_maps(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_build_coulomb_maps'

      INTEGER                                            :: handle, i, n_drop_sqrt, n_RI
      REAL(KIND=dp)                                      :: lambda_kept, lambda_max, norm_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: M_inv, projector, V_b_inv_sqrt, &
                                                            V_b_sqrt, V_zero, work
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI

      ! Same values gamma_corr_init already set; recomputed so this routine is self-contained.
      CALL gamma_corr_slab_geometry(bs_env, gc%area, gc%b_recip)
      gc%metric_head = gamma_corr_metric_head(bs_env)

      ! ---------------------------------------------------------------------------------------
      ! 1. V_PQ(0): the finite real-space 2D lattice sum of the bare Coulomb operator (Eq. 63).
      ! ---------------------------------------------------------------------------------------
      ! V_zero(P, Q) = V_PQ(0) in Hartree, dense n_RI x n_RI, replicated on every rank.
      ALLOCATE (V_zero(n_RI, n_RI))
      CALL gamma_corr_coulomb_at_zero(qs_env, bs_env, V_zero)

      ! ---------------------------------------------------------------------------------------
      ! 2. c_P and the neutral projector Π_PQ = δ_PQ − c_P c_Q/Σ_S|c_S|² (Eq. 64).
      ! ---------------------------------------------------------------------------------------
      ! Unnormalized charges c_P = ∫ dr ϕ_P0(r) (units of ∫ϕ), replicated.
      ALLOCATE (gc%charge(n_RI))
      CALL gamma_corr_auxiliary_charges(qs_env, gc%charge)

      ! From here on gc%charge holds the UNIT vector ĉ = c/|c| and gc%charge_norm holds |c|,
      ! so that Π = Id − ĉ ĉᵀ and the unnormalized c = charge_norm · charge when needed
      ! (gamma_corr_subtract_divergent_head).
      norm_c = SQRT(DOT_PRODUCT(gc%charge, gc%charge))
      IF (norm_c <= 0.0_dp) &
         CPABORT('Gamma-point correction: the RI_AUX basis carries no charge at all')
      gc%charge(:) = gc%charge/norm_c
      gc%charge_norm = norm_c

      ! Dense projector, column by column: Π(:, i) = e_i − ĉ ĉ_i. Symmetric and idempotent,
      ! with ĉ as its only null vector.
      ALLOCATE (projector(n_RI, n_RI))
      DO i = 1, n_RI
         projector(:, i) = -gc%charge(:)*gc%charge(i)
         projector(i, i) = projector(i, i) + 1.0_dp
      END DO

      ! ---------------------------------------------------------------------------------------
      ! 3. V_b = Π V(0) Π and its ±1/2 powers on the neutral subspace (Eq. 65).
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (work(n_RI, n_RI), V_b_sqrt(n_RI, n_RI), V_b_inv_sqrt(n_RI, n_RI))
      ! work = V(0) Π
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, V_zero, n_RI, projector, n_RI, &
                 0.0_dp, work, n_RI)
      ! V_b = Π work = Π V(0) Π, held for now in the V_b_sqrt buffer.
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, projector, n_RI, work, n_RI, &
                 0.0_dp, V_b_sqrt, n_RI)

      ! One diagonalization for both V_b^(+1/2) and V_b^(−1/2), with the same relative cut
      ! eps_null_V_b on the same spectrum: both powers act on one and the same neutral subspace,
      ! which the maps below rely on.
      CALL gamma_corr_symmetric_power(V_b_sqrt, 0.5_dp, eps_null_V_b, n_drop_sqrt, &
                                      lambda_max=lambda_max, lambda_kept=lambda_kept, &
                                      exponent_2=-0.5_dp, matrix_2=V_b_inv_sqrt)

      ! ---------------------------------------------------------------------------------------
      ! 4. M^-1(0): NOT rebuilt here. bs_env%fm_Minv_Gamma already holds the inverse metric at Γ
      !    that compute_W_MIC folds into W, built with regularization_RI deliberately switched
      !    off so that the RI-RS and tensor paths agree. Building a second, differently
      !    regularized M^-1 for the correction is exactly the mistake that cost 3.8 eV in Σ^c
      !    once before, so the correction reuses this one verbatim.
      ! ---------------------------------------------------------------------------------------
      ! Gather the distributed M^-1(0) into a dense replicated copy, M_inv(P, Q).
      ALLOCATE (M_inv(n_RI, n_RI))
      CALL cp_fm_get_submatrix(bs_env%fm_Minv_Gamma, M_inv)

      ! ---------------------------------------------------------------------------------------
      ! 5. F_density = M^-1 Π V_b^(1/2)  and  F_extern = M^-1 V(0) Π V_b^(-1/2)  (Eqs. 66, 67).
      !    Both are stored replicated as F(Q, P): row Q is the RI fit index that Z_ℓQ contracts,
      !    column P is the index of the symmetrized neutral body (P in Eqs. 66-67). Both have
      !    rank n_RI − (number of discarded modes of V_b).
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (gc%F_density(n_RI, n_RI), gc%F_extern(n_RI, n_RI))

      ! work = Π V_b^(1/2). The retained eigenvectors of V_b are orthogonal to ĉ, so this
      ! equals V_b^(1/2) up to rounding; Π is kept to follow Eq. 66 literally.
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, projector, n_RI, V_b_sqrt, n_RI, &
                 0.0_dp, work, n_RI)
      ! F_density = M^-1 work = M^-1 Π V_b^(1/2)
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, M_inv, n_RI, work, n_RI, &
                 0.0_dp, gc%F_density, n_RI)

      ! work = Π V_b^(-1/2)
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, projector, n_RI, V_b_inv_sqrt, n_RI, &
                 0.0_dp, work, n_RI)
      ! V_b_sqrt is no longer needed and is reused as a buffer: V(0) Π V_b^(-1/2). Note the
      ! FULL V(0) on the left (not V_b): the charged test-function couplings are kept here.
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, V_zero, n_RI, work, n_RI, &
                 0.0_dp, V_b_sqrt, n_RI)
      ! G_extern = V(0) Π V_b^(-1/2) is kept as well: it is F_extern in the RI projection
      ! representation, the one W^MIC_PQ is in before the M^-1 W M^-1 contraction, which is where
      ! gamma_corr_add_to_W adds the wing.
      ALLOCATE (gc%G_extern(n_RI, n_RI))
      gc%G_extern(:, :) = V_b_sqrt
      ! F_extern = M^-1 V(0) Π V_b^(-1/2)
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, M_inv, n_RI, V_b_sqrt, n_RI, &
                 0.0_dp, gc%F_extern, n_RI)

      ! Diagnostics: the null-mode count and the (λ_max, smallest kept λ) pair show whether
      ! eps_null_V_b cut inside a spectral gap (see gamma_corr_symmetric_power).
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'Gamma-point correction of W (2D): charge channel separated from the RI basis'
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Layer area A (bohr^2)', gc%area
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'RI metric head int dr v(r)', gc%metric_head
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Norm of the charge vector |c|', norm_c
         WRITE (bs_env%unit_nr, '(T4,A,T72,I10)') 'Null modes of the neutral body V_b', n_drop_sqrt
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Largest eigenvalue of V_b', lambda_max
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Smallest retained eigenvalue of V_b', lambda_kept
      END IF

      DEALLOCATE (V_zero, projector, work, V_b_sqrt, V_b_inv_sqrt, M_inv)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_build_coulomb_maps

! **************************************************************************************************
!> \brief V_PQ(q=0) from the real-space 2c Coulomb lattice sum (gwp42 Eq. 63).
!>
!>        Built with a one-point k-mesh at the origin but with nkp_grid kept at the mesh for
!>        χ, ε, W, because get_factor_and_xyz_min_max takes the summation box from nkp_grid, not
!>        from the number of k-points: the box is ±(factor·n−1)/2 cells with factor = 3^(S−1) for
!>        odd n. The result is real; the imaginary part returned by the sum is zero at q = 0 for a
!>        centrosymmetric box and is asserted to be so.
!> \param qs_env QS environment, for cell, particles and the RI_AUX basis
!> \param bs_env band-structure environment, for the k-mesh, SIZE_LATTICE_SUM and RI templates
!> \param V_zero V_PQ(0), replicated on every rank
! **************************************************************************************************
   SUBROUTINE gamma_corr_coulomb_at_zero(qs_env, bs_env, V_zero)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: V_zero

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_coulomb_at_zero'

      INTEGER                                            :: handle, re_im
      REAL(KIND=dp)                                      :: max_imag
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(kpoint_type), POINTER                         :: kpoints_zero
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      ! A private one-point k-point set holding q = 0 only (fractional coordinates xkp = 0),
      ! so that the production kpoints_chi_eps_W object is left untouched.
      NULLIFY (kpoints_zero)
      CALL kpoint_create(kpoints_zero)
      ALLOCATE (kpoints_zero%xkp(3, 1))
      kpoints_zero%xkp(:, 1) = 0.0_dp
      kpoints_zero%nkp = 1
      kpoints_zero%kp_scheme = "GENERAL"
      ! The box of the lattice sum comes from nkp_grid, so keep the production mesh here.
      kpoints_zero%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ! mat_V(1, 1) = Re V(0), mat_V(1, 2) = Im V(0): RI x RI DBCSR matrices with every block
      ! reserved and zeroed, which is the layout build_2c_coulomb_matrix_kp accumulates into.
      NULLIFY (mat_V)
      ALLOCATE (mat_V(1:1, 2))
      DO re_im = 1, 2
         NULLIFY (mat_V(1, re_im)%matrix)
         ALLOCATE (mat_V(1, re_im)%matrix)
         CALL dbcsr_create(mat_V(1, re_im)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V(1, re_im)%matrix)
         CALL dbcsr_set(mat_V(1, re_im)%matrix, 0.0_dp)
      END DO

      ! V_PQ(0) = Σ_{R ∈ box} (P0|QR), the bare Coulomb lattice sum with unit phases (Eq. 63).
      CALL build_2c_coulomb_matrix_kp(mat_V, kpoints_zero, basis_type="RI_AUX", cell=cell, &
                                      particle_set=particle_set, qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=bs_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=1, ikp_end=1)

      ! The real part is kept distributed in gc%fm_V_zero: compute_V_k_by_lattice_sum needs the
      ! identical matrix for the regular W path and reuses it instead of repeating the lattice sum.
      gc => bs_env%ri_rs%gamma_corr
      CALL cp_fm_create(gc%fm_V_zero, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(gc%fm_V_zero, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_V(1, 1)%matrix, gc%fm_V_zero)
      gc%V_zero_cached = .TRUE.
      CALL cp_fm_get_submatrix(gc%fm_V_zero, V_zero)

      ! The imaginary part, checked on the distributed matrix (no second replication).
      CALL cp_fm_create(fm_work, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(fm_work, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_V(1, 2)%matrix, fm_work)
      CALL cp_fm_maxabsval(fm_work, max_imag)
      CALL cp_fm_release(fm_work)
      CALL dbcsr_deallocate_matrix_set(mat_V)
      DEALLOCATE (kpoints_zero%xkp)
      CALL kpoint_release(kpoints_zero)

      ! Im V(0) = 0 for a centrosymmetric image box; a nonzero value means the box is lopsided
      ! and the whole real-arithmetic argument of this module would be void.
      IF (max_imag > 1.0E-10_dp*MAX(1.0_dp, MAXVAL(ABS(V_zero)))) THEN
         CPABORT('Gamma-point correction: V(q=0) came back complex')
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_coulomb_at_zero

! **************************************************************************************************
!> \brief Decide whether the Γ-point correction runs and set the geometry it needs, early enough
!>        that the Z_lP fit can already carry the extra charge-channel right-hand side.
!>
!>        Only the cheap, purely geometric quantities are set here (layer area, in-plane
!>        reciprocal vectors, metric head). The Coulomb side -- V(0), the projector and the two
!>        maps -- is built later by gamma_corr_build_coulomb_maps, after the RI metric exists.
!> \param bs_env band-structure environment; sets bs_env%ri_rs%gamma_corr%enabled and geometry
! **************************************************************************************************
   SUBROUTINE gamma_corr_init(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      TYPE(gamma_corr_type), POINTER                     :: gc

      gc => bs_env%ri_rs%gamma_corr

      ! The single switch every later Γ-correction call site tests (gamma_corr%enabled).
      gc%enabled = gamma_corr_requested(bs_env)
      IF (.NOT. gc%enabled) RETURN

      ! A, b_1, b_2 and ∫dr v_ω(r): needed by gamma_corr_charge_rhs inside compute_Z_lP.
      CALL gamma_corr_slab_geometry(bs_env, gc%area, gc%b_recip)
      gc%metric_head = gamma_corr_metric_head(bs_env)

      ! Tell the user why the correction is switched on.
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'The mesh for chi, eps, W contains q = 0, where the real-space Coulomb lattice sum'
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'over charged auxiliary functions diverges. Separating the charge channel and'
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'adding the Gamma-point head-and-wing correction of W.'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

   END SUBROUTINE gamma_corr_init

! **************************************************************************************************
!> \brief Diagonal quadratic form of the AO values on the grid, x_ℓ = Σ_µν Φ_ℓµ G_µν Φ_ℓν, as the
!>        row dot products of Y = Φ G with Φ.
!>
!>        Y is formed on Φ's OWN block pattern (retain_sparsity): only those blocks meet a nonzero
!>        Φ_ℓµ in the row quadratic form, and the identical pattern is what lets the two matrices
!>        be walked in lockstep. (dbcsr_get_block_p on Φ from inside a Y iteration returned wrong
!>        blocks for about half the grid points; never fetch by index here.) Linear in the blocks
!>        of Φ, hence correct under any distribution: each rank adds the blocks it holds and the
!>        caller reduces over ranks.
!> \param mat_phi Φ, redistributed onto its owner ranks
!> \param mat_G n_AO x n_AO operand
!> \param mat_Y work matrix carrying Φ's block pattern (a copy of Φ made by the caller)
!> \param x result, replicated array, zero outside the rows of this rank's blocks
! **************************************************************************************************
   SUBROUTINE gamma_corr_phi_quadratic_form(mat_phi, mat_G, mat_Y, x)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi, mat_G, mat_Y
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: x

      INTEGER                                            :: col, col_phi, i, i_row, row, row_phi
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block, block_phi
      TYPE(dbcsr_iterator_type)                          :: iter, iter_phi

      ! Y = Φ G restricted to the pattern of Y, i.e. of Φ.
      CALL dbcsr_set(mat_Y, 0.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, mat_phi, mat_G, 0.0_dp, mat_Y, &
                          retain_sparsity=.TRUE.)

      ! x_ℓ = Σ_µ Y_ℓµ Φ_ℓµ. Walk Y and Φ in lockstep: identical patterns and distributions mean
      ! the two iterators return the same (row, col) blocks in the same order; the assertions
      ! check that.
      x(:) = 0.0_dp
      CALL dbcsr_iterator_start(iter, mat_Y)
      CALL dbcsr_iterator_start(iter_phi, mat_phi)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_offset=i_row)
         CPASSERT(dbcsr_iterator_blocks_left(iter_phi))
         CALL dbcsr_iterator_next_block(iter_phi, row_phi, col_phi, block_phi)
         CPASSERT(row == row_phi .AND. col == col_phi)
         ! Row-wise dot product over the AO columns of this block, accumulated per grid point.
         DO i = 1, SIZE(block, 1)
            x(i_row + i - 1) = x(i_row + i - 1) + SUM(block(i, :)*block_phi(i, :))
         END DO
      END DO
      CPASSERT(.NOT. dbcsr_iterator_blocks_left(iter_phi))
      CALL dbcsr_iterator_stop(iter_phi)
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE gamma_corr_phi_quadratic_form

! **************************************************************************************************
!> \brief Optical transition dipoles µ^α_ai between occupied i and empty a at Γ.
!>
!>        gwp42 Eq. 41 writes the interband position matrix element as the AO position matrix
!>        plus momentum derivatives of H and S,
!>
!>          µ^α_ai = Σ_µν C*_µa C_νi [ r^α_µν − i (∂_kα H_µν − ε_a ∂_kα S_µν)/(ε_a − ε_i) ] ,
!>
!>        because position is not an observable in a crystal. In a Γ-ONLY supercell run those
!>        derivative terms are identically zero and the first term is the whole answer:
!>        qs_neighbor_lists sets dokp = (kpoints%nkp > 0), so without a &KPOINTS section sab_kp is
!>        never built and matrix_ks_kp/matrix_s_kp carry a single folded image at R = 0, whence
!>        ∂_kα H = Σ_R i R^α H(R) ≡ 0 and likewise for S. So
!>
!>          µ^α_ai = Σ_µν C_µa r^α_µν C_νi ,   r^α_µν = Σ_R ∫ dr ϕ_µ0(r) r^α ϕ_νR(r) .
!>
!>        This is legitimate exactly because the cell is large: the position operator is well
!>        defined here only because the orbitals are localized well inside the supercell, so the
!>        result must be checked against cell size rather than assumed.
!>
!>        The code does NOT evaluate r^α_µν, however. A minimum-image position matrix is not
!>        translation invariant on a supercell (see the comment in the body), so the same matrix
!>        element is taken in the velocity gauge, from the commutator, which for eigenstates of
!>        H is identical: ⟨a|r^α|i⟩ = ⟨a|[H, r^α]|i⟩/(ε_a − ε_i), [H, r^α] = −∂^α + [V_nl, r^α].
!> \param qs_env QS environment, for the AO overlap template, neighbor lists and pseudopotential
!> \param bs_env band-structure environment, for the Γ MO coefficients and eigenvalues
!> \param mu µ^α_ai, indices (α, i, a, spin), replicated on every rank; a = 1 is MO n_occ + 1
! **************************************************************************************************
   SUBROUTINE gamma_corr_transition_dipoles(qs_env, bs_env, mu)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), INTENT(OUT)              :: mu

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_transition_dipoles'

      INTEGER                                            :: handle, i_dir, i_occ, i_spin, i_vir, &
                                                            n_ao, n_mo, n_occ, n_vir
      REAL(KIND=dp)                                      :: delta, n_elec
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mu_mo
      REAL(KIND=dp), DIMENSION(3)                        :: f_sum
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: struct_mo_mo
      TYPE(cp_fm_type)                                   :: fm_mo_mo, fm_r_ao, fm_r_mo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, moments
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! MO window: the retained MOs if a window was set (n_mo_retained > 0), otherwise all n_ao.
      n_ao = bs_env%n_ao
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = n_ao

      ! The commutator [H, r^α], NOT the position operator itself.
      !
      !   [H, r^α] = -∂^α + [V_nl, r^α]     and      ⟨a|r^α|i⟩ = ⟨a|[H, r^α]|i⟩ / (ε_a - ε_i)
      !
      ! Under periodic boundary conditions r^α is not a legal operator: the Γ-cell Bloch states
      ! are delocalized and ⟨a|r^α|i⟩ built from a minimum-image-folded position matrix is both
      ! origin-dependent and spuriously large, because folding couples pairs that carry no
      ! crystal momentum transfer. Its two signatures on a 14x14 h-BN cell were an in-plane head
      ! three orders of magnitude above the known 2D polarizability and S_xx /= S_yy on a
      ! hexagonal lattice, while S_zz -- the one direction that is NOT periodic and therefore
      ! never folded -- came out the right size. The commutator is translation invariant, so it
      ! has neither defect, and it vanishes between states of different crystal momentum as it
      ! must.
      !
      ! Construction follows the COMMUTATOR_HR print path in qs_scf_post_gpw: both builders
      ! ACCUMULATE into the same antisymmetric matrix, then one overall sign makes it [H, r].
      ! build_com_mom_nl returns immediately when sap_ppnl is null, so all-electron runs get the
      ! bare -∂^α, which is the correct commutator for a purely local potential.
      NULLIFY (matrix_s, particle_set, qs_kind_set, sab_orb, sap_ppnl, cell, dft_control)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, sab_orb=sab_orb, sap_ppnl=sap_ppnl, &
                      cell=cell, dft_control=dft_control)

      ! moments(α) = AO matrix of [H, r^α], α = x, y, z: antisymmetric (real, anti-Hermitian
      ! operator), with the sparsity of the orbital neighbor list sab_orb, zero-initialized.
      NULLIFY (moments)
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

      ! mu is sized for the worst case over spins: the largest occupied count and the largest
      ! empty count n_mo − min(n_occ). Entries beyond a spin's own (n_occ, n_vir) stay zero.
      n_occ = MAXVAL(bs_env%n_occ(1:bs_env%n_spin))
      n_vir = n_mo - MINVAL(bs_env%n_occ(1:bs_env%n_spin))
      CPASSERT(n_occ > 0 .AND. n_vir > 0)
      ALLOCATE (mu(3, n_occ, n_vir, bs_env%n_spin))
      mu(:, :, :, :) = 0.0_dp

      ! Work matrices: n_mo x n_mo for Cᵀ R C, AO x AO for R, AO x MO for R C.
      NULLIFY (struct_mo_mo)
      CALL cp_fm_struct_create(struct_mo_mo, &
                               context=bs_env%fm_mo_coeff_Gamma(1)%matrix_struct%context, &
                               para_env=bs_env%fm_mo_coeff_Gamma(1)%matrix_struct%para_env, &
                               nrow_global=n_mo, ncol_global=n_mo)
      CALL cp_fm_create(fm_mo_mo, struct_mo_mo)
      CALL cp_fm_create(fm_r_ao, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_r_mo, bs_env%fm_mo_coeff_Gamma(1)%matrix_struct)
      ALLOCATE (mu_mo(n_mo, n_mo))

      DO i_spin = 1, bs_env%n_spin
         DO i_dir = 1, 3
            ! µ^α = Cᵀ r^α C, built as two gemms through the AO x MO intermediate.
            ! (Here R = moments(α) holds [H, r^α], not r^α; see above.)
            ! R as a full AO x AO matrix; copy_dbcsr_to_fm expands the antisymmetric storage.
            CALL cp_fm_set_all(fm_r_ao, 0.0_dp)
            CALL copy_dbcsr_to_fm(moments(i_dir)%matrix, fm_r_ao)
            ! fm_r_mo = R C(:, 1:n_mo)
            CALL parallel_gemm('N', 'N', n_ao, n_mo, n_ao, 1.0_dp, fm_r_ao, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), 0.0_dp, fm_r_mo)
            ! fm_mo_mo = C(:, 1:n_mo)ᵀ R C(:, 1:n_mo), i.e. element (m, n) = ⟨m|[H, r^α]|n⟩
            CALL parallel_gemm('T', 'N', n_mo, n_mo, n_ao, 1.0_dp, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), &
                               fm_r_mo, 0.0_dp, fm_mo_mo)
            ! Replicate the n_mo x n_mo MO matrix on every rank.
            CALL cp_fm_get_submatrix(fm_mo_mo, mu_mo)

            ! Keep only the occupied-empty block and divide by the energy denominator:
            ! µ^α_ai = ⟨a|[H,r^α]|i⟩/(ε_a - ε_i). i runs over occupied, a over empty, so the
            ! denominator is the positive excitation energy and cannot vanish in a gapped system.
            DO i_vir = 1, n_mo - bs_env%n_occ(i_spin)
               DO i_occ = 1, bs_env%n_occ(i_spin)
                  ! Δ_ia = ε_a − ε_i (Hartree) from the Γ SCF eigenvalues; a = n_occ + i_vir.
                  delta = bs_env%eigenval_scf_Gamma(bs_env%n_occ(i_spin) + i_vir, i_spin) - &
                          bs_env%eigenval_scf_Gamma(i_occ, i_spin)
                  CPASSERT(delta > 0.0_dp)
                  ! µ^α_ai = ⟨a|[H,r^α]|i⟩/(ε_a − ε_i) = ⟨a|r^α|i⟩: row = empty a, column =
                  ! occupied i. The transposed element ⟨i|[H,r^α]|a⟩ is its NEGATIVE (the
                  ! commutator matrix is antisymmetric), which would flip the sign of every
                  ! wing while leaving the head and the sum rule, both quadratic in µ, unchanged.
                  mu(i_dir, i_occ, i_vir, i_spin) = &
                     mu_mo(bs_env%n_occ(i_spin) + i_vir, i_occ)/delta
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (mu_mo)
      CALL cp_fm_release(fm_mo_mo)
      CALL cp_fm_release(fm_r_ao)
      CALL cp_fm_release(fm_r_mo)
      CALL cp_fm_struct_release(struct_mo_mo)
      CALL dbcsr_deallocate_matrix_set(moments)

      ! Thomas-Reiche-Kuhn sum rule, the one basis- and gauge-independent check on µ:
      !
      !   (2/n_spin) Σ_spin Σ_ai 2 (ε_a - ε_i) |µ^α_ai|²  =  N_electrons
      !
      ! A truncated virtual space can only make the left side too SMALL, so a ratio above one is
      ! unambiguous evidence of spurious dipoles -- which is exactly what the folded position
      ! operator produced. Printed per direction because the in-plane and out-of-plane
      ! directions fail differently under periodic boundary conditions.
      IF (bs_env%unit_nr > 0) THEN
         ! f_sum(α) = Σ_spin Σ_ia 2 Δ_ia (µ^α_ai)², the summed oscillator strengths.
         f_sum(:) = 0.0_dp
         DO i_spin = 1, bs_env%n_spin
            DO i_vir = 1, n_mo - bs_env%n_occ(i_spin)
               DO i_occ = 1, bs_env%n_occ(i_spin)
                  delta = bs_env%eigenval_scf_Gamma(bs_env%n_occ(i_spin) + i_vir, i_spin) - &
                          bs_env%eigenval_scf_Gamma(i_occ, i_spin)
                  f_sum(:) = f_sum(:) + 2.0_dp*delta*mu(:, i_occ, i_vir, i_spin)**2
               END DO
            END DO
         END DO
         ! Occupation factor 2/n_spin per state (2 for closed shell, 1 per spin channel), on
         ! the oscillator-strength sum and on the electron count alike.
         f_sum(:) = f_sum(:)*2.0_dp/REAL(bs_env%n_spin, dp)
         n_elec = 2.0_dp/REAL(bs_env%n_spin, dp)*REAL(SUM(bs_env%n_occ(1:bs_env%n_spin)), dp)
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'Gamma-point correction: Thomas-Reiche-Kuhn sum rule of the transition dipoles'
         WRITE (bs_env%unit_nr, '(T4,A,T60,3F12.4)') 'f-sum / N_electrons   (x, y, z)', &
            f_sum(:)/MAX(n_elec, 1.0_dp)
         WRITE (bs_env%unit_nr, '(T4,A,T72,F10.1)') 'N_electrons in the retained MO window', n_elec
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_transition_dipoles

! **************************************************************************************************
!> \brief The optical response of the Γ-only supercell: the Cartesian head tensor and the three
!>        Cartesian wing grid vectors, in imaginary time and then in imaginary frequency.
!>
!>        A Γ-only supercell has no k sum at all -- every transition is vertical at the single
!>        point -- so gwp42's (1/N_k) Σ_k collapses to a bare sum over the supercell's own
!>        occupied i and empty a. With Δ_ia = ε_a − ε_i and the 2D normalizations of Eqs. 61-62,
!>
!>          Y_αβ(iτ) = (2/n_spin/A)  Σ_σ Σ_ia µ^α_ai µ^β_ai e^(−Δ_ia τ)            [Eq. 49]
!>          X_ℓα(iτ) = (2/n_spin/√A) Σ_σ Σ_ia ψ_i(r_ℓ) ψ_a(r_ℓ) µ^α_ai e^(−Δ_ia τ) [Eq. 48]
!>
!>        so that χ^ΓC_αβ = Y_αβ and χ^ΓC_ℓα = −i X_ℓα. Splitting the explicit factor i out keeps
!>        every stored array real: the i's cancel against each other in the Schur complement
!>        below, and in W they reappear only as a quantity plus its own conjugate.
!>
!>        The exponential factorizes, e^(−Δ_ia τ) = e^(−|ε_a−ε_F|τ) e^(−|ε_i−ε_F|τ), so X is the
!>        diagonal of Φ P^α(iτ) Φᵀ with the dense AO matrix P^α = C_occ e µ^α e C_virᵀ, evaluated
!>        by gamma_corr_accumulate_response with the same sparse kernel as the fitting operator.
!>
!>        Both are cosine-transformed with the same minimax weights as the regular response, and
!>        the frequency-domain quantities the correction actually needs are then assembled:
!>
!>          B_PQ(iω)  = δ_PQ + [F_densityᵀ χ_PQ(Γ,iω) F_density]_PQ                 [Eq. 72]
!>          E_Pα(iω)  = [F_densityᵀ (Zᵀ X)]_Pα                                      [Eq. 71]
!>          S_αβ(iω)  = Y_αβ − [Eᵀ B⁻¹ E]_αβ
!>          T^w_ℓα(iω) = Σ_Q Z_ℓQ [F_extern B⁻¹ E]_Qα
!>
!>        B is the regular Gaussian dielectric body at Γ. It needs NO new grid contraction:
!>        T = Z F_density, so Tᵀ χ_grid T = F_densityᵀ (Zᵀ χ_grid Z) F_density and the middle
!>        factor is exactly the χ_PQ(Γ) the code already computes.
!>
!>        S and T^w are everything the cell integrals need, because the 2D block inversion
!>        (Eqs. 73-74) collapses to closed form. With a(q̂,iω) = q̂ᵀ S q̂ and v(q) = 2π/q,
!>
!>          h(q,iω) − 1 = 2πq q̂ᵀ Y q̂ ,   w_P = −i √(2πq) (E q̂)_P
!>          h − w†B⁻¹w = 1 + 2πq a
!>          W^c head  = v h̃ = −4π² a /(1 + 2πq a)
!>          W^c wing  = i 2π t_ℓ/(1 + 2πq a) ,  t_ℓ(q̂) = Σ_α q̂_α T^w_ℓα
!>
!>        Both are finite as q → 0: h̃ ∝ q cancels 1/q and w̃_P ∝ √q cancels 1/√q. Note t is
!>        LINEAR in q̂, which is what lets the angular integral be done once into a table.
!>
!>        Outputs (in bs_env%ri_rs%gamma_corr), all replicated on every rank:
!>          S_head(α, β, i_w)   S_αβ(iω_j), 3 x 3 x n_freq
!>          T_wing(ℓ, α, i_w)   T^w_ℓα(iω_j), n_grid x 3 x n_freq
!>        The time and frequency grids have the same number of points, num_time_freq_points,
!>        so n_t below counts both.
!> \param bs_env band-structure environment; S_head and T_wing are stored in its gamma_corr
!> \param qs_env QS environment, for the transition dipoles
!> \param mat_chi_Gamma_tau χ_PQ(Γ,iτ), the regular polarizability the driver already built
! **************************************************************************************************
   SUBROUTINE gamma_corr_response(bs_env, qs_env, mat_chi_Gamma_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_response'

      INTEGER                                            :: handle, i_dir, i_spin, i_t, i_w, j_dir, &
                                                            n_grid, n_RI, n_t
      REAL(KIND=dp)                                      :: spin_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: chi_PQ, E_freq, work_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: X_freq, X_time, Y_time
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: mu
      TYPE(cp_fm_type)                                   :: fm_chi
      TYPE(dbcsr_type)                                   :: mat_phi_own
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI
      n_grid = bs_env%ri_rs%n_grid_points
      n_t = bs_env%num_time_freq_points

      ! Eqs. 48/49 carry an explicit factor two for the two spin states of a closed shell; for a
      ! spin-polarized run that factor is the spin sum instead, hence 2/n_spin per channel.
      spin_factor = 2.0_dp/REAL(bs_env%n_spin, dp)

      ! µ^α_ai for all spins, replicated; dimensions (3, max n_occ, max n_vir, n_spin).
      CALL gamma_corr_transition_dipoles(qs_env, bs_env, mu)

      ! Y_time(α, β, i_t): head, replicated. X_time(ℓ, α, i_t): wings, one column per Cartesian
      ! direction, accumulated per rank over its own grid rows.
      ALLOCATE (Y_time(3, 3, n_t), X_time(n_grid, 3, n_t))
      Y_time(:, :, :) = 0.0_dp
      X_time(:, :, :) = 0.0_dp

      ! Φ redistributed onto its owner ranks: the retain_sparsity product in the wing response
      ! is undefined on the as-built storage pattern (see gamma_corr_fit_charge_channel).
      CALL dbcsr_create(mat_phi_own, template=bs_env%ri_rs%mat_phi_mu_l)
      CALL dbcsr_complete_redistribute(bs_env%ri_rs%mat_phi_mu_l, mat_phi_own)
      ! Sum the spin channels; each needs its own dipoles and MO coefficients.
      DO i_spin = 1, bs_env%n_spin
         CALL gamma_corr_accumulate_response(bs_env, i_spin, mu, mat_phi_own, spin_factor, &
                                             Y_time, X_time)
      END DO
      CALL dbcsr_release(mat_phi_own)
      DEALLOCATE (mu)

      ! Complete X over the grid rows held by the other ranks. Y is NOT summed: every rank
      ! already computed the full Y redundantly (gamma_corr_accumulate_response).
      CALL bs_env%para_env%sum(X_time)

      ! ---------------------------------------------------------------------------------------
      ! Cosine transform iτ -> iω with the same minimax weights as the regular response.
      !   F(iω_j) = Σ_k cosine_time_to_frequency(j, k) F(iτ_k)
      ! S_head receives Y(iω) here; the body term −Eᵀ B⁻¹ E is subtracted per frequency below.
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (gc%S_head(3, 3, n_t), X_freq(n_grid, 3, n_t))
      gc%S_head(:, :, :) = 0.0_dp
      X_freq(:, :, :) = 0.0_dp
      DO i_w = 1, n_t
         DO i_t = 1, n_t
            gc%S_head(:, :, i_w) = gc%S_head(:, :, i_w) + &
                                   bs_env%time_frequency_grid%cosine_time_to_frequency(i_w, i_t)* &
                                   Y_time(:, :, i_t)
            X_freq(:, :, i_w) = X_freq(:, :, i_w) + &
                                bs_env%time_frequency_grid%cosine_time_to_frequency(i_w, i_t)* &
                                X_time(:, :, i_t)
         END DO
      END DO
      DEALLOCATE (Y_time, X_time)

      ! ---------------------------------------------------------------------------------------
      ! Per frequency: body, wings, Schur complement and the grid wing map.
      ! ---------------------------------------------------------------------------------------
      ! Dense replicated work arrays: chi_PQ (first χ, then overwritten by B), work_RI (an
      ! n_RI x n_RI scratch whose first three columns later hold Zᵀ X and B⁻¹E), E_freq
      ! (first E, later reused for F_extern B⁻¹ E).
      ALLOCATE (gc%T_RI(n_RI, 3, n_t))
      gc%T_RI(:, :, :) = 0.0_dp
      ALLOCATE (chi_PQ(n_RI, n_RI), work_RI(n_RI, n_RI), E_freq(n_RI, 3))
      CALL cp_fm_create(fm_chi, bs_env%fm_RI_RI%matrix_struct)

      DO i_w = 1, n_t
         ! χ_PQ(Γ,iω) from the stored χ(Γ,iτ): the same transform the regular W path uses.
         CALL gamma_corr_chi_PQ_at_freq(bs_env, fm_chi, i_w, mat_chi_Gamma_tau, chi_PQ)

         ! B_PQ = δ_PQ + F_densityᵀ χ_PQ F_density (Eq. 72).
         ! work_RI = χ F_density
         CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, chi_PQ, n_RI, gc%F_density, n_RI, &
                    0.0_dp, work_RI, n_RI)
         ! chi_PQ <- F_densityᵀ work_RI = F_densityᵀ χ F_density
         CALL dgemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, gc%F_density, n_RI, work_RI, n_RI, &
                    0.0_dp, chi_PQ, n_RI)
         ! chi_PQ <- chi_PQ + Id, so chi_PQ now holds B. (i_dir is only a loop counter here.)
         DO i_dir = 1, n_RI
            chi_PQ(i_dir, i_dir) = chi_PQ(i_dir, i_dir) + 1.0_dp
         END DO

         ! E_Pα = F_densityᵀ (Zᵀ X)_Qα (Eq. 71 without its √(2πq), which the cell integral adds).
         ! work_RI(:, 1:3) = Zᵀ X(iω_j), i.e. Σ_ℓ Z_ℓQ X_ℓα
         CALL gamma_corr_project_grid_to_RI(bs_env, X_freq(:, :, i_w), work_RI(:, 1:3))
         ! E_freq(P, α) = Σ_Q F_density(Q, P) work_RI(Q, α)
         CALL dgemm('T', 'N', n_RI, 3, n_RI, 1.0_dp, gc%F_density, n_RI, work_RI, n_RI, &
                    0.0_dp, E_freq, n_RI)

         ! S_αβ = Y_αβ − Eᵀ B⁻¹ E, by solving B x = E rather than inverting B.
         ! work_RI(:, 1:3) <- x = B⁻¹ E (B in chi_PQ is overwritten by its LU factors).
         CALL gamma_corr_body_solve(chi_PQ, E_freq, work_RI(:, 1:3))
         ! S_αβ −= Σ_P E_Pα x_Pβ
         DO j_dir = 1, 3
            DO i_dir = 1, 3
               gc%S_head(i_dir, j_dir, i_w) = gc%S_head(i_dir, j_dir, i_w) - &
                                              DOT_PRODUCT(E_freq(:, i_dir), work_RI(:, j_dir))
            END DO
         END DO

         ! T_Pα(iω) = [G_extern B⁻¹ E]_Pα = [V(0) Π V_b^(-1/2) B⁻¹ E]_Pα: the wing vector in the
         ! RI projection representation. (Its grid image would be Z M⁻¹ T = Z F_extern B⁻¹ E,
         ! but the correction is applied in RI space, see gamma_corr_add_to_W.)
         CALL dgemm('N', 'N', n_RI, 3, n_RI, 1.0_dp, gc%G_extern, n_RI, work_RI, n_RI, &
                    0.0_dp, gc%T_RI(:, :, i_w), n_RI)
      END DO

      CALL cp_fm_release(fm_chi)
      DEALLOCATE (chi_PQ, work_RI, E_freq, X_freq)
      ! The three n_RI x n_RI maps are consumed; free them before W and Σ are built.
      DEALLOCATE (gc%F_density, gc%F_extern, gc%G_extern)

      ! Print the diagonal of S(iω) per frequency: S_xx = S_yy is expected for a hexagonal
      ! layer, and S_zz is the out-of-plane response, which does not enter a(q̂) = q̂ᵀ S q̂
      ! for the in-plane directions q̂ of the cell integrals.
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'Gamma-point correction: Cartesian dielectric head S_ab(i*omega) [a.u.]'
         WRITE (bs_env%unit_nr, '(T4,A)') 'i_w           S_xx           S_yy           S_zz'
         DO i_w = 1, n_t
            WRITE (bs_env%unit_nr, '(T4,I3,3ES15.6)') i_w, gc%S_head(1, 1, i_w), &
               gc%S_head(2, 2, i_w), gc%S_head(3, 3, i_w)
         END DO
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_response

! **************************************************************************************************
!> \brief Accumulate one spin channel's contribution to the head tensor Y_αβ(iτ) and the wing grid
!>        vectors X_ℓα(iτ).
!>
!>        Y, Eq. 49, is a pure band quantity and is computed redundantly on every rank from the
!>        replicated dipoles: cheaper than communicating nine numbers.
!>
!>        X, Eq. 48, is evaluated as a diagonal quadratic form of the AO values on the grid,
!>
!>          X_ℓα(iτ) = (2/n_spin/√A) Σ_µν Φ_µ(r_ℓ) P^α_µν(iτ) Φ_ν(r_ℓ),
!>          P^α(iτ)  = C_occ diag(e^(−|ε_i−ε_F|τ)) µ^α diag(e^(−|ε_a−ε_F|τ)) C_virᵀ,
!>
!>        which is exactly the construction gamma_corr_apply_D uses for the fitting operator, with
!>        the dense n_AO x n_AO matrix P^α(iτ) in place of G: P^α comes from two ScaLAPACK gemms
!>        on the MO coefficients, the product Φ P^α is restricted to Φ's own block pattern, and
!>        the row dot products with Φ give X on the grid rows this rank holds; the caller reduces
!>        over ranks.
!>
!>        Why not the band-resolved form Σ_ia ψ_i(r_ℓ) ψ_a(r_ℓ) µ^α_ai with ψ = Φ C on the grid: it
!>        needs every MO column of ψ(r_ℓ) on the rank that holds row ℓ, but Φ C inherits the 2D
!>        process grid of Φ, so a rank sees only the MO blocks of its own process column, and
!>        contracting i with a across blocks it does not own silently drops every cross-column
!>        pair. Measured on the 3x3 h-BN test: 1 rank 5.945 eV, 4 ranks on a 2x2 grid 6.560 eV.
!>        The quadratic form is linear in the blocks of Φ and therefore correct under any
!>        distribution, and no orbital is ever stored densely on the grid.
!> \param bs_env band-structure environment, for eigenvalues, E_F, MO coefficients, τ grid and A
!> \param i_spin spin channel
!> \param mu µ^α_ai for all spins, mu(α, i, a, spin), replicated; a = 1 is MO n_occ + 1
!> \param mat_phi Φ_µ(r_ℓ), redistributed onto its owner ranks (see gamma_corr_fit_charge_channel)
!> \param spin_factor 2/n_spin
!> \param Y_time accumulated head tensor, Y(α,β,iτ)
!> \param X_time accumulated wing grid vectors, X(ℓ,α,iτ), partial over grid rows
! **************************************************************************************************
   SUBROUTINE gamma_corr_accumulate_response(bs_env, i_spin, mu, mat_phi, spin_factor, &
                                             Y_time, X_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: i_spin
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: mu
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi
      REAL(KIND=dp), INTENT(IN)                          :: spin_factor
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Y_time, X_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_accumulate_response'

      INTEGER                                            :: handle, i_a, i_dir, i_i, i_t, j_dir, &
                                                            n_ao, n_grid, n_mo, n_occ, n_vir
      REAL(KIND=dp)                                      :: e_fermi, prefactor, tau, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exp_occ, exp_vir, x_partial
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_ov
      TYPE(cp_fm_struct_type), POINTER                   :: struct_ao_vir, struct_occ_vir
      TYPE(cp_fm_type)                                   :: fm_D, fm_P, fm_W
      TYPE(dbcsr_type)                                   :: mat_P, mat_Y

      CALL timeset(routineN, handle)

      ! MO range of this spin: occupied 1..n_occ, empty n_occ+1..n_mo, with n_mo the number of
      ! retained MOs (all n_ao of them if no window is set).
      n_ao = bs_env%n_ao
      n_occ = bs_env%n_occ(i_spin)
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = n_ao
      n_vir = n_mo - n_occ
      n_grid = bs_env%ri_rs%n_grid_points
      ! Reference energy for the factorized exponential (Hartree).
      e_fermi = bs_env%e_fermi(i_spin)
      CPASSERT(n_occ > 0 .AND. n_vir > 0)
      CPASSERT(SIZE(mu, 2) >= n_occ .AND. SIZE(mu, 3) >= n_vir)

      ! Distributed work matrices: D = diag(e_occ) µ^α diag(e_vir) (n_occ x n_vir),
      ! W = C_occ D (n_AO x n_vir), P = W C_virᵀ (n_AO x n_AO), and the DBCSR copy of P in the
      ! AO blocking of Φ's columns.
      NULLIFY (struct_occ_vir, struct_ao_vir)
      CALL cp_fm_struct_create(struct_occ_vir, &
                               context=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct%context, &
                               para_env=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct%para_env, &
                               nrow_global=n_occ, ncol_global=n_vir)
      CALL cp_fm_struct_create(struct_ao_vir, &
                               context=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct%context, &
                               para_env=bs_env%fm_mo_coeff_Gamma(i_spin)%matrix_struct%para_env, &
                               nrow_global=n_ao, ncol_global=n_vir)
      CALL cp_fm_create(fm_D, struct_occ_vir)
      CALL cp_fm_create(fm_W, struct_ao_vir)
      CALL cp_fm_create(fm_P, bs_env%fm_s_Gamma%matrix_struct)
      CALL dbcsr_create(mat_P, template=bs_env%mat_ao_ao%matrix)
      ! Work matrix with Φ's block pattern for the retain_sparsity product Φ P.
      CALL dbcsr_create(mat_Y, template=mat_phi)
      CALL dbcsr_copy(mat_Y, mat_phi)

      ! Per-τ work arrays: exponential factors per band, D for one direction, and the partial X
      ! of this rank; prefactor = (2/n_spin)/√A of Eq. 48.
      ALLOCATE (exp_occ(n_occ), exp_vir(n_vir), d_ov(n_occ, n_vir), x_partial(n_grid))
      prefactor = spin_factor/SQRT(bs_env%ri_rs%gamma_corr%area)

      DO i_t = 1, bs_env%num_time_freq_points
         ! τ_k of the minimax grid (Hartree⁻¹).
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         ! e^(-Δ_ia τ) = e^(-|ε_a - ε_F| τ) e^(-|ε_i - ε_F| τ)
         ! exact because ε_i < ε_F < ε_a, so |ε_a − ε_F| + |ε_i − ε_F| = ε_a − ε_i = Δ_ia.
         ! exp_occ(i) = e^(−|ε_i − ε_F| τ), exp_vir(a) = e^(−|ε_a − ε_F| τ)
         DO i_i = 1, n_occ
            exp_occ(i_i) = EXP(-ABS(bs_env%eigenval_scf_Gamma(i_i, i_spin) - e_fermi)*tau)
         END DO
         DO i_a = 1, n_vir
            exp_vir(i_a) = EXP(-ABS(bs_env%eigenval_scf_Gamma(n_occ + i_a, i_spin) - e_fermi)*tau)
         END DO

         ! Y_αβ(iτ) = (2/n_spin/A) Σ_ia µ^α_ai µ^β_ai e^(-Δ_ia τ)                      [Eq. 49]
         DO i_a = 1, n_vir
            DO i_i = 1, n_occ
               ! weight = (2/n_spin) e^(−Δ_ia τ) / A
               weight = spin_factor*exp_occ(i_i)*exp_vir(i_a)/bs_env%ri_rs%gamma_corr%area
               ! Y_αβ += weight · µ^α_ai µ^β_ai (all nine components, symmetric in α, β)
               DO j_dir = 1, 3
                  DO i_dir = 1, 3
                     Y_time(i_dir, j_dir, i_t) = Y_time(i_dir, j_dir, i_t) + &
                                                 weight*mu(i_dir, i_i, i_a, i_spin)* &
                                                 mu(j_dir, i_i, i_a, i_spin)
                  END DO
               END DO
            END DO
         END DO

         DO i_dir = 1, 3
            ! D_ia = e_i µ^α_ia e_a, assembled replicated and scattered into the distributed fm.
            DO i_a = 1, n_vir
               DO i_i = 1, n_occ
                  d_ov(i_i, i_a) = exp_occ(i_i)*mu(i_dir, i_i, i_a, i_spin)*exp_vir(i_a)
               END DO
            END DO
            CALL cp_fm_set_submatrix(fm_D, d_ov)
            ! W = C_occ D, with C_occ the first n_occ columns of the MO coefficients.
            CALL parallel_gemm('N', 'N', n_ao, n_vir, n_occ, 1.0_dp, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), fm_D, 0.0_dp, fm_W)
            ! P = W C_virᵀ, with C_vir the columns n_occ+1 .. n_occ+n_vir.
            CALL parallel_gemm('N', 'T', n_ao, n_ao, n_vir, 1.0_dp, fm_W, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), 0.0_dp, fm_P, &
                               b_first_col=n_occ + 1)
            CALL copy_fm_to_dbcsr(fm_P, mat_P)
            ! X_ℓ = Σ_µ (Φ P)_ℓµ Φ_ℓµ on this rank's grid rows.                        [Eq. 48]
            CALL gamma_corr_phi_quadratic_form(mat_phi, mat_P, mat_Y, x_partial)
            X_time(:, i_dir, i_t) = X_time(:, i_dir, i_t) + prefactor*x_partial
         END DO
      END DO

      DEALLOCATE (exp_occ, exp_vir, d_ov, x_partial)
      CALL dbcsr_release(mat_P)
      CALL dbcsr_release(mat_Y)
      CALL cp_fm_release(fm_D)
      CALL cp_fm_release(fm_W)
      CALL cp_fm_release(fm_P)
      CALL cp_fm_struct_release(struct_occ_vir)
      CALL cp_fm_struct_release(struct_ao_vir)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_accumulate_response

! **************************************************************************************************
!> \brief χ_PQ(Γ,iω_j) from the stored χ_PQ(Γ,iτ), replicated as a dense block.
!>
!>          χ_PQ(iω_j) = Σ_k cosine_time_to_frequency(j, k) χ_PQ(iτ_k)
!>
!>        accumulated in DBCSR, then gathered into a dense replicated n_RI x n_RI array.
!> \param bs_env band-structure environment; its mat_RI_RI is used as DBCSR scratch
!> \param fm_chi scratch full matrix with the RI x RI structure
!> \param j_w frequency index
!> \param mat_chi_Gamma_tau χ_PQ(Γ,iτ_k), one DBCSR matrix per time point
!> \param chi_PQ χ_PQ(Γ,iω_j), replicated
! **************************************************************************************************
   SUBROUTINE gamma_corr_chi_PQ_at_freq(bs_env, fm_chi, j_w, mat_chi_Gamma_tau, chi_PQ)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_chi
      INTEGER, INTENT(IN)                                :: j_w
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: chi_PQ

      INTEGER                                            :: i_t
      REAL(KIND=dp)                                      :: weight

      ! mat_RI_RI <- 0, then += w_jk χ(iτ_k) for every time point k.
      CALL dbcsr_set(bs_env%mat_RI_RI%matrix, 0.0_dp)
      DO i_t = 1, bs_env%num_time_freq_points
         weight = bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, i_t)
         CALL dbcsr_add(bs_env%mat_RI_RI%matrix, mat_chi_Gamma_tau(i_t)%matrix, 1.0_dp, weight)
      END DO
      ! DBCSR -> block-cyclic full matrix -> dense copy on every rank.
      CALL copy_dbcsr_to_fm(bs_env%mat_RI_RI%matrix, fm_chi)
      CALL cp_fm_get_submatrix(fm_chi, chi_PQ)

   END SUBROUTINE gamma_corr_chi_PQ_at_freq

! **************************************************************************************************
!> \brief (Zᵀ X)_Qα = Σ_ℓ Z_ℓQ X_ℓα for a three-column grid quantity.
!>
!>        Done by walking the local blocks of the sparse Z rather than by a DBCSR product, which
!>        would need X wrapped as a three-column matrix with matching block structure.
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

      ! For each local block of Z (grid rows row_start.., RI columns col_start..) add
      ! X_RI(Q, α) += Σ_{ℓ∈block} Z_ℓQ X_ℓα, one dot product per RI column Q of the block.
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
!> \brief Solve B x = E for the three Cartesian right-hand sides.
!>
!>        gwp42 is explicit that B's inverse is never formed: the body is the regular Gaussian
!>        dielectric matrix at Γ, it is well conditioned because ε ⪰ Id, and three solves are
!>        cheaper and better behaved than an inversion.
!> \param B body matrix, destroyed (overwritten by its LU factors)
!> \param E right-hand sides, n_RI x 3
!> \param x solution, n_RI x 3
! **************************************************************************************************
   SUBROUTINE gamma_corr_body_solve(B, E, x)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: B
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: E
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: x

      INTEGER                                            :: info, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv

      n = SIZE(B, 1)
      ALLOCATE (ipiv(n))
      ! dgesv overwrites the right-hand sides with the solution, so start x from E.
      x(:, 1:3) = E(:, 1:3)
      ! LU with partial pivoting (general solver: B = Id + F_densityᵀ χ F_density is symmetric
      ! in exact arithmetic, but no symmetry is assumed here).
      CALL dgesv(n, 3, B, n, ipiv, x, n, info)
      IF (info /= 0) CPABORT('Gamma-point correction: the dielectric body solve failed')
      DEALLOCATE (ipiv)

   END SUBROUTINE gamma_corr_body_solve

! **************************************************************************************************
!> \brief Tabulate the two Brillouin-zone cell integrals of the singular model, as functions of the
!>        in-plane minimum-image displacement R.
!>
!>        **What differs from the full-k code.** There, gwp42 Eq. 78 replaces the four discrete
!>        near-Γ mesh samples by their cell integral and SUBTRACTS the samples, and the whole
!>        R dependence is over lattice vectors of the small cell, so one scalar per cell suffices.
!>        Here the supercell has a single cell, nothing was ever sampled at q = 0 because the
!>        charge mode was removed from the basis, so there is nothing to subtract, and all of the
!>        R dependence is the intra-cell grid displacement R = MIC(r_ℓ − r_ℓ'), which spans the
!>        whole supercell. The per-cell scalar therefore becomes a function of a continuous
!>        in-plane 2-vector:
!>
!>          ΔW^c_ℓℓ'(iω) = Z_ℓ0 Z_ℓ'0 K_h(R,iω)
!>                       + Σ_α [T^w_ℓα Z_ℓ'0 − Z_ℓ0 T^w_ℓ'α] K_w,α(R,iω)
!>
!>          K_h(R,iω)   = (A/(2π)²) ∮dq̂ w(q̂) ∫_0^{q_max} dq q · [−4π² a/(1 + 2πq a)] cos(q q̂·R)
!>          K_w,α(R,iω) = (A/(2π)²) ∮dq̂ w(q̂) q̂_α ∫_0^{q_max} dq q · [2π/(1 + 2πq a)] sin(q q̂·R)
!>
!>        with a(q̂,iω) = q̂ᵀ S q̂. The polar Jacobian q is included, so the 1/q head integrand is
!>        flat and the 1/√q wing integrand integrable -- that is the 2D analogue of the Coulomb
!>        singularity being integrable in the plane while its discrete sampling is not.
!>
!>        **Why cos and sin.** a is even in q̂ (a quadratic form) and t_ℓ = Σ_α q̂_α T^w_ℓα is odd
!>        (linear). Integrating over a region symmetric under q → −q, the head keeps only the even
!>        part of e^(−i q·R), i.e. cos, and the wing -- which carries an explicit factor i -- keeps
!>        only i·(−i sin) = sin. Both kernels come out real, K_h even in R and K_w,α odd, so the
!>        antisymmetric bracket above multiplied by the odd kernel gives a symmetric ΔW, as it
!>        must be.
!>
!>        **Why a table.** The kernels depend on R only through its in-plane part and are smooth,
!>        so they are built once on a regular mesh in fractional coordinates over [−1/2,1/2]² and
!>        bilinearly interpolated afterwards. Evaluating them per grid pair instead would cost
!>        n_grid² · n_dir · n_radial transcendentals, about 10^11 on a 392-atom cell, which is not
!>        affordable even once. With the table the per-pair work is a lookup and a multiply, the
!>        same shape as the Hadamard product Σ^c = −G ∘ W that the code already performs.
!>
!>        **Quadrature.** Polar coordinates around Γ: n_dir = 4·n_face directions q̂ from
!>        gamma_corr_bz_boundary with angular weights w(q̂) and boundary radii q_max(q̂), and
!>        n_radial Gauss-Legendre nodes on [0, q_max(q̂)] along each. The tables are filled in
!>        imaginary FREQUENCY, K_head(i_R, i_w) and K_wing(i_R, α, i_w), with the flattened table
!>        index i_R = 1 + i_u + n_tab·i_v; gamma_corr_pack_hot_loop transforms them to time.
!>        The table points are distributed round-robin over the ranks and summed at the end.
!> \param bs_env band-structure environment; the tables are stored in its gamma_corr
!> \param n_tab points per direction of the tabulation mesh
! **************************************************************************************************
   SUBROUTINE gamma_corr_tabulate(bs_env, n_tab)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: n_tab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_tabulate'
      INTEGER, PARAMETER                                 :: n_face = 24, n_radial = 32

      INTEGER                                            :: handle, i_dir, i_node, i_tab, i_u, i_v, &
                                                            i_w, n_dir, n_t
      REAL(KIND=dp)                                      :: a_scalar, area_bz, denom, jac, proj, q, &
                                                            q_max, radial_w, start_time
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dir_radius, dir_weight, nodes, rule
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dir_hat
      REAL(KIND=dp), DIMENSION(3)                        :: R_cart
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      start_time = m_walltime()
      gc => bs_env%ri_rs%gamma_corr
      n_t = bs_env%num_time_freq_points
      ! Area of the small Brillouin zone, (2π)²/A; dividing by it gives the A/(2π)² prefactor.
      area_bz = twopi**2/gc%area

      ! Angular quadrature (q̂, w(q̂), q_max(q̂)) and the radial Gauss-Legendre rule on [−1, 1].
      CALL gamma_corr_bz_boundary(bs_env, n_face, dir_hat, dir_weight, dir_radius)
      n_dir = SIZE(dir_weight)
      CALL gamma_corr_gauss_legendre(n_radial, nodes, rule)

      ! Square n_tab x n_tab table, flattened: K_head(i_R, i_w), K_wing(i_R, α, i_w).
      gc%n_tab = [n_tab, n_tab]
      ALLOCATE (gc%K_head(n_tab*n_tab, n_t), gc%K_wing(n_tab*n_tab, 3, n_t))
      gc%K_head(:, :) = 0.0_dp
      gc%K_wing(:, :, :) = 0.0_dp

      ! Table points are evenly spaced in fractional coordinates over the full minimum-image
      ! range [-1/2, 1/2] along each in-plane lattice vector, inclusive of both ends so that the
      ! interpolation never has to extrapolate.
      ! Rank mepos takes table points mepos+1, mepos+1+num_pe, ... (round-robin).
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(dynamic) &
!$OMP          PRIVATE(i_tab, i_u, i_v, R_cart, i_dir, q_max, proj, jac, i_w, a_scalar, i_node, &
!$OMP                  q, denom, radial_w) &
!$OMP          SHARED(gc, bs_env, n_tab, n_dir, n_t, dir_hat, dir_weight, dir_radius, nodes, rule, &
!$OMP                 area_bz)
      DO i_tab = bs_env%para_env%mepos + 1, n_tab*n_tab, bs_env%para_env%num_pe
         ! i_tab − 1 = i_u + n_tab·i_v, with i_u, i_v ∈ {0, …, n_tab − 1}.
         i_u = MODULO(i_tab - 1, n_tab)
         i_v = (i_tab - 1)/n_tab
         ! R = (i_u/(n_tab−1) − 1/2) a_1 + (i_v/(n_tab−1) − 1/2) a_2, Cartesian, bohr.
         R_cart = bs_env%hmat(:, 1)*(REAL(i_u, dp)/REAL(n_tab - 1, dp) - 0.5_dp) + &
                  bs_env%hmat(:, 2)*(REAL(i_v, dp)/REAL(n_tab - 1, dp) - 0.5_dp)

         DO i_dir = 1, n_dir
            ! Along direction q̂: boundary radius q_max(q̂), projection t = q̂·R, and the
            ! angular weight times the prefactor, jac = w(q̂) A/(2π)².
            q_max = dir_radius(i_dir)
            proj = DOT_PRODUCT(dir_hat(:, i_dir), R_cart)
            jac = dir_weight(i_dir)/area_bz

            DO i_w = 1, n_t
               ! a(q_hat, i_w) = q_hat^T S q_hat
               ! q̂ lies in the layer plane (b_1 and b_2 are orthogonal to n), so S's
               ! out-of-plane row/column do not contribute.
               a_scalar = DOT_PRODUCT(dir_hat(:, i_dir), &
                                      MATMUL(gc%S_head(:, :, i_w), dir_hat(:, i_dir)))

               DO i_node = 1, n_radial
                  ! Gauss-Legendre on [0, q_max]; the extra q is the polar Jacobian.
                  ! Node map x ∈ [−1, 1] -> q = q_max (x + 1)/2, weight scaled by q_max/2.
                  q = 0.5_dp*q_max*(nodes(i_node) + 1.0_dp)
                  ! Screening denominator 1 + 2πq a(q̂, iω) of both head and wing.
                  denom = 1.0_dp + twopi*q*a_scalar
                  ! radial_w = (q_max/2) w_GL · q / (1 + 2πq a): weight, Jacobian, denominator.
                  radial_w = 0.5_dp*q_max*rule(i_node)*q/denom

                  ! K_h += w(q̂) A/(2π)² · dq q · [−4π² a/(1 + 2πq a)] cos(q q̂·R)
                  gc%K_head(i_tab, i_w) = gc%K_head(i_tab, i_w) - &
                                          jac*radial_w*4.0_dp*pi**2*a_scalar*COS(q*proj)
                  ! K_w,α += w(q̂) A/(2π)² · dq q · [2π/(1 + 2πq a)] sin(q q̂·R) q̂_α
                  gc%K_wing(i_tab, :, i_w) = gc%K_wing(i_tab, :, i_w) + &
                                             jac*radial_w*twopi*SIN(q*proj)*dir_hat(:, i_dir)
               END DO
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      ! Every table point was computed by exactly one rank (zeros elsewhere): sum = gather.
      CALL bs_env%para_env%sum(gc%K_head)
      CALL bs_env%para_env%sum(gc%K_wing)

      ! The tables K_head(R, iω), K_wing(R, α, iω) and T_RI are consumed by gamma_corr_add_to_W;
      ! the head tensor is not needed any more.
      DEALLOCATE (gc%S_head)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,I4,A,I4,A,T62,F10.1,A)') &
            'Gamma-point correction: cell-integral tables (', n_tab, ' x', n_tab, &
            '), Exec. time', m_walltime() - start_time, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      DEALLOCATE (dir_hat, dir_weight, dir_radius, nodes, rule)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_tabulate

! **************************************************************************************************
!> \brief Add the Γ-point head-and-wing correction to W^MIC_PQ(iτ) in the RI basis.
!>
!>        Eq. DeltaW evaluated on the auxiliary functions instead of the grid functions:
!>
!>          (ϕ_P|ΔW|ϕ_Q)(iω) = (c_P c_Q / A) K_h(R_PQ, iω)
!>                           + (1/√A) Σ_α [T_Pα c_Q − c_P T_Qα] K_w,α(R_PQ, iω) ,
!>
!>        with c_P = ∫dr ϕ_P the auxiliary charges, T_Pα = [V(0) Π V_b^(-1/2) B⁻¹ E]_Pα the wing
!>        vector of gamma_corr_response and R_PQ the in-plane minimum-image separation of the
!>        atoms carrying P and Q (the Gaussians are treated as point charges at their centres,
!>        valid for q_max σ << 1). The result is transformed to iτ with the cosine weights of the
!>        minimax grid and added to W^MIC_PQ(iτ) BEFORE the M⁻¹ W M⁻¹ contraction, so that it
!>        reaches the grid pairs through the same fit Z M⁻¹ as the regular W.
!>
!>        Why not on the grid: Z_ℓ0 Z_ℓ'0 K_h(R_ℓℓ') with the grid charges Z_ℓ0 = ∫ξ_ℓ/√A is the
!>        same expression, but the RI-RS interpolation functions ξ_ℓ far from the atoms carry
!>        fit weights of 10^4 .. 10^8 of alternating sign. They cancel exactly in linear
!>        functionals Σ_ℓ ρ(ℓ) Z_ℓ0 (checked: the charge of a transition density comes out as
!>        3e-6 against partial sums of 2e-2) but not once each point is weighted by a different
!>        kernel value K(R_ℓℓ'), which is what any pair kernel on the grid does. The wing, being
!>        odd and cancellation-dominated, was wrong by orders of magnitude that way; the head
!>        survived only on large cells where the kernel is nearly flat over neighbouring points.
!> \param bs_env band-structure environment: the correction data, the cell, the RI-atom map
!> \param qs_env QS environment, for the atomic positions
!> \param fm_W_MIC_time W^MIC_PQ(iτ) in the projection representation, one fm per time point
! **************************************************************************************************
   SUBROUTINE gamma_corr_add_to_W(bs_env, qs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_add_to_W'

      INTEGER                                            :: atom_P, atom_Q, handle, i_atom, i_col, &
                                                            i_row, i_t, i_u, i_v, i_w, k00, &
                                                            n_atom, n_RI, n_t, n_tab, ncol_local, &
                                                            nrow_local, P, Q
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: c_P, c_Q, du, dv, f_u, f_v, inv_A, &
                                                            inv_sqrt_A, scale, start_time, w00, &
                                                            w01, w10, w11
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c, dW_freq, K_h
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: frac_atom, K_w
      TYPE(gamma_corr_type), POINTER                     :: gc
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      start_time = m_walltime()

      gc => bs_env%ri_rs%gamma_corr
      CPASSERT(gc%ready)
      n_RI = bs_env%n_RI
      n_t = bs_env%num_time_freq_points
      n_tab = gc%n_tab(1)
      ! n_tab − 1 table intervals per unit fractional length (bilinear interpolation below).
      scale = REAL(n_tab - 1, dp)
      inv_A = 1.0_dp/gc%area
      inv_sqrt_A = 1.0_dp/SQRT(gc%area)

      ! Unnormalized auxiliary charges c_P = ∫ϕ_P (gc%charge is the unit vector ĉ).
      ALLOCATE (c(n_RI))
      c(:) = gc%charge_norm*gc%charge

      ! Atom of every auxiliary function and the in-plane fractional coordinates of the atoms
      ! (unfolded; the minimum image is taken on the difference, as pbc() does).
      NULLIFY (particle_set)
      CALL get_qs_env(qs_env, particle_set=particle_set)
      n_atom = SIZE(particle_set)
      ALLOCATE (frac_atom(2, n_atom), atom_of_RI(n_RI))
      DO i_atom = 1, n_atom
         frac_atom(1, i_atom) = DOT_PRODUCT(bs_env%ri_rs%cell%h_inv(1, :), particle_set(i_atom)%r)
         frac_atom(2, i_atom) = DOT_PRODUCT(bs_env%ri_rs%cell%h_inv(2, :), particle_set(i_atom)%r)
         atom_of_RI(bs_env%i_RI_start_from_atom(i_atom):bs_env%i_RI_end_from_atom(i_atom)) = i_atom
      END DO

      ! Local block of the distributed W^MIC(iτ): every time point shares one layout.
      CALL cp_fm_get_info(fm_W_MIC_time(1), nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      ALLOCATE (dW_freq(n_t), K_h(n_t), K_w(3, n_t))

      DO i_col = 1, ncol_local
         Q = col_indices(i_col)
         atom_Q = atom_of_RI(Q)
         c_Q = c(Q)
         DO i_row = 1, nrow_local
            P = row_indices(i_row)
            atom_P = atom_of_RI(P)
            c_P = c(P)

            ! Minimum-image fractional displacement of R_P − R_Q, mapped onto [−1/2, 1/2].
            du = frac_atom(1, atom_P) - frac_atom(1, atom_Q)
            dv = frac_atom(2, atom_P) - frac_atom(2, atom_Q)
            du = du - ANINT(du)
            dv = dv - ANINT(dv)

            ! Bilinear interpolation of the tables at this displacement (same scheme as the
            ! table layout of gamma_corr_tabulate: i_R = 1 + i_u + n_tab i_v over [−1/2, 1/2]²).
            f_u = (du + 0.5_dp)*scale
            f_v = (dv + 0.5_dp)*scale
            i_u = MIN(MAX(INT(f_u), 0), n_tab - 2)
            i_v = MIN(MAX(INT(f_v), 0), n_tab - 2)
            f_u = f_u - REAL(i_u, dp)
            f_v = f_v - REAL(i_v, dp)
            w00 = (1.0_dp - f_u)*(1.0_dp - f_v)
            w10 = f_u*(1.0_dp - f_v)
            w01 = (1.0_dp - f_u)*f_v
            w11 = f_u*f_v
            k00 = 1 + i_u + n_tab*i_v
            K_h(:) = w00*gc%K_head(k00, :) + w10*gc%K_head(k00 + 1, :) + &
                     w01*gc%K_head(k00 + n_tab, :) + w11*gc%K_head(k00 + n_tab + 1, :)
            K_w(:, :) = w00*gc%K_wing(k00, :, :) + w10*gc%K_wing(k00 + 1, :, :) + &
                        w01*gc%K_wing(k00 + n_tab, :, :) + w11*gc%K_wing(k00 + n_tab + 1, :, :)

            ! (ϕ_P|ΔW|ϕ_Q)(iω) per frequency: head c_P c_Q K_h / A and wing
            ! (1/√A) Σ_α [T_Pα c_Q − c_P T_Qα] K_w,α. Symmetric in P <-> Q: K_h is even and K_w
            ! odd in R, the wing bracket is antisymmetric.
            DO i_w = 1, n_t
               dW_freq(i_w) = inv_A*c_P*c_Q*K_h(i_w) + inv_sqrt_A* &
                              SUM(K_w(:, i_w)*(gc%T_RI(P, :, i_w)*c_Q - c_P*gc%T_RI(Q, :, i_w)))
            END DO

            ! iω -> iτ with the cosine weights of the minimax grid (the transform W^MIC uses),
            ! accumulated straight into the local element of every W^MIC(iτ).
            DO i_t = 1, n_t
               fm_W_MIC_time(i_t)%local_data(i_row, i_col) = &
                  fm_W_MIC_time(i_t)%local_data(i_row, i_col) + &
                  SUM(bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, :)*dW_freq)
            END DO
         END DO
      END DO

      DEALLOCATE (c, frac_atom, atom_of_RI, dW_freq, K_h, K_w)
      ! The correction is now inside W^MIC; its ingredients are dead.
      DEALLOCATE (gc%K_head, gc%K_wing, gc%T_RI)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T62,F10.1,A)') &
            'Gamma-point correction added to W^MIC_PQ(itau) in RI space, Exec. time', &
            m_walltime() - start_time, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_add_to_W

! **************************************************************************************************
!> \brief Polar quadrature of the Brillouin-zone cell around Γ: directions q̂, the radius q_max to
!>        the cell boundary along each, and the angular weights.
!>
!>        The cell is the parallelogram spanned by ±b_1/2 and ±b_2/2, i.e. the single BZ cell
!>        around Γ, not the four near-Γ mesh cells of the full-k construction. Each of its four
!>        edges is sampled by Gauss-Legendre; a point on edge i at parameter t gives
!>        q = ±b_i/2 + t·b_j/2, and the weight |q × dq|/|q|² is the dθ subtended. Both signs are
!>        included, which is what makes the direction set symmetric under q̂ → −q̂ and lets the
!>        cos/sin split in gamma_corr_tabulate be exact.
!>
!>        The area check Σ w·q_max²/2 = (2π)²/A verifies the quadrature closes on the cell.
!> \param bs_env band-structure environment, for the cell (via gamma_corr_slab_geometry)
!> \param n_face Gauss-Legendre points per edge
!> \param dir_hat unit directions, 3 x n_dir, Cartesian, n_dir = 4·n_face
!> \param dir_weight angular weights dθ (radians); they sum to 2π
!> \param dir_radius distance to the cell boundary along each direction (1/bohr)
! **************************************************************************************************
   SUBROUTINE gamma_corr_bz_boundary(bs_env, n_face, dir_hat, dir_weight, dir_radius)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: n_face
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: dir_hat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: dir_weight, dir_radius

      INTEGER                                            :: i_axis, i_dir, i_node, i_sign, j_axis
      REAL(KIND=dp)                                      :: area, area_check
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: nodes, rule
      REAL(KIND=dp), DIMENSION(3)                        :: cross, dq, q
      REAL(KIND=dp), DIMENSION(3, 2)                     :: b

      CALL gamma_corr_slab_geometry(bs_env, area, b)
      CALL gamma_corr_gauss_legendre(n_face, nodes, rule)

      ALLOCATE (dir_hat(3, 4*n_face), dir_weight(4*n_face), dir_radius(4*n_face))
      i_dir = 0
      ! Edge pair i_axis: the two edges at ±b_i/2, each running along b_j (j ≠ i).
      DO i_axis = 1, 2
         j_axis = 3 - i_axis
         ! dq = dq/dt: a point on the edge is q(t) = ±b_i/2 + t b_j/2, t ∈ [−1, 1].
         dq = 0.5_dp*b(:, j_axis)
         DO i_sign = -1, 1, 2
            DO i_node = 1, n_face
               i_dir = i_dir + 1
               ! Boundary point at the Gauss-Legendre node t = nodes(i_node).
               q = REAL(i_sign, dp)*0.5_dp*b(:, i_axis) + nodes(i_node)*dq
               ! q_max(q̂) = |q| and q̂ = q/|q|.
               dir_radius(i_dir) = SQRT(SUM(q*q))
               dir_hat(:, i_dir) = q/dir_radius(i_dir)
               ! dθ = |q × dq/dt| dt / |q|²: the angle the edge element subtends at Γ.
               ! The absolute value makes every weight positive regardless of orientation.
               cross = [q(2)*dq(3) - q(3)*dq(2), q(3)*dq(1) - q(1)*dq(3), &
                        q(1)*dq(2) - q(2)*dq(1)]
               dir_weight(i_dir) = rule(i_node)*SQRT(SUM(cross**2))/SUM(q*q)
            END DO
         END DO
      END DO

      ! Cell area from the polar quadrature, ∮ dθ q_max(θ)²/2, must equal (2π)²/A.
      area_check = 0.5_dp*SUM(dir_weight*dir_radius**2)
      IF (ABS(area_check*area/twopi**2 - 1.0_dp) > 1.0E-8_dp) &
         CPABORT('Gamma-point correction: the BZ cell quadrature failed its area check')

      DEALLOCATE (nodes, rule)

   END SUBROUTINE gamma_corr_bz_boundary

! **************************************************************************************************
!> \brief Gauss-Legendre nodes and weights on [-1, 1], by Newton iteration on the Legendre
!>        polynomial with the standard asymptotic starting guess.
!> \param n number of points
!> \param nodes abscissae
!> \param weights weights
! **************************************************************************************************
   SUBROUTINE gamma_corr_gauss_legendre(n, nodes, weights)
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: nodes, weights

      INTEGER                                            :: i, iter, j
      REAL(KIND=dp)                                      :: dp_n, p0, p1, p2, x, x_old

      ALLOCATE (nodes(n), weights(n))

      DO i = 1, n
         ! Chebyshev-like initial guess for the i-th root of P_n.
         x = COS(pi*(REAL(i, dp) - 0.25_dp)/(REAL(n, dp) + 0.5_dp))
         DO iter = 1, 100
            ! Three-term recurrence for P_n(x) and its derivative.
            ! j P_j = (2j − 1) x P_{j−1} − (j − 1) P_{j−2}, from P_0 = 1, P_{−1} = 0.
            ! On exit p1 = P_n(x) and p2 = P_{n−1}(x).
            p1 = 1.0_dp
            p2 = 0.0_dp
            DO j = 1, n
               p0 = p2
               p2 = p1
               p1 = (REAL(2*j - 1, dp)*x*p2 - REAL(j - 1, dp)*p0)/REAL(j, dp)
            END DO
            ! P_n'(x) = n (x P_n − P_{n−1})/(x² − 1)
            dp_n = REAL(n, dp)*(x*p1 - p2)/(x*x - 1.0_dp)
            ! Newton step x <- x − P_n(x)/P_n'(x), stopped at 1E-14 change.
            x_old = x
            x = x - p1/dp_n
            IF (ABS(x - x_old) <= 1.0E-14_dp) EXIT
         END DO
         ! w_i = 2/((1 − x_i²) P_n'(x_i)²), with P_n' from the last Newton iteration.
         nodes(i) = x
         weights(i) = 2.0_dp/((1.0_dp - x*x)*dp_n*dp_n)
      END DO

   END SUBROUTINE gamma_corr_gauss_legendre

! **************************************************************************************************
!> \brief Remove the divergent G = 0 constant from the real-space lattice sum V(q=0) used by the
!>        REGULAR W path, and nothing else: V(0) <- V(0) - D_box c cᵀ.
!>
!>        gwp42 Sec. II: the regular V_PQ(0) is Eq. 23 "with the singular G = 0 term omitted".
!>        That term is (2π/(A q)) c_P c_Q with c_P = ∫dr ϕ_P, so omitting it is a RANK-ONE
!>        subtraction along the unnormalized charge vector. It is NOT the projection Π V Π of
!>        Eq. 65: Π V Π also deletes the charged direction's G ≠ 0 interactions, which are finite,
!>        short-ranged, physical, and not restored by the head/wing model. Measured on the 6 Å
!>        water test with the projection in this path: Σ^c of every state 5-9x too large.
!>
!>        For the finite real-space sum over the image box of lattice_sum, Poisson summation gives
!>
!>          Σ_{R∈box} (P0|QR) = (1/A) ∫_box d²r (P0|Q r) + [G ≠ 0 terms → V^reg]
!>
!>        and for a box much larger than the Gaussians the first term is c_P c_Q D_box with
!>
!>          D_box = (1/A) ∫_box d²r / |r| = (1/A) ∮ (r × dr) / |r|
!>
!>        (the in-plane integral of 1/|r − d| over a large region does not depend on the
!>        displacement d, so there is no finite multipole remainder; the −2π|z−z'| part of the 2D
!>        kernel, which the ϕ_0 = 1/√A model does not carry, is deliberately left inside V^reg).
!>        The box is the parallelogram ±F_1 a_1 ± F_2 a_2 covering the R set of
!>        get_factor_and_xyz_min_max, F_α = 3^(s−1) n_α / 2 for the odd mesh, and the boundary
!>        integral is done by Gauss-Legendre on its four edges with an area self-check.
!>
!>        Consistency: ΔW adds the cell integral of the singular model over the whole small BZ with
!>        nothing subtracted, which is exactly the head that this routine takes out of V(0).
!>
!>        Requires gc%charge (ĉ) and gc%charge_norm (|c|) from gamma_corr_build_coulomb_maps;
!>        c cᵀ = |c|² ĉ ĉᵀ. The cfm is complex only because the shared k-point path stores V(k)
!>        that way; at q = 0 the subtracted term is real.
!> \param bs_env band-structure environment, for the cell, the image box and gamma_corr
!> \param cfm_V V(q=0), overwritten by V(0) − D_box c cᵀ
! **************************************************************************************************
   SUBROUTINE gamma_corr_subtract_divergent_head(bs_env, cfm_V)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: cfm_V

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_subtract_divergent_head'
      INTEGER, PARAMETER                                 :: n_gl = 32

      COMPLEX(KIND=dp)                                   :: c_dot_u
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: u
      INTEGER                                            :: col, handle, i_col, i_edge, i_node, &
                                                            i_row, n_RI, row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: area_check, d_box, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: nodes, rule
      REAL(KIND=dp), DIMENSION(2)                        :: f_half
      REAL(KIND=dp), DIMENSION(3)                        :: dr, r
      REAL(KIND=dp), DIMENSION(3, 4)                     :: corner
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      CPASSERT(ALLOCATED(gc%charge))
      CPASSERT(gc%charge_norm > 0.0_dp .AND. gc%area > 0.0_dp)
      n_RI = bs_env%n_RI

      ! ---------------------------------------------------------------------------------------
      ! 1. D_box = (1/A) ∮ (r × dr)/|r| around the image box, counterclockwise.
      ! ---------------------------------------------------------------------------------------
      ! Half-widths F_α = 3^(s−1) n_α / 2 in units of a_α (odd mesh, factor 3^(s−1) as in
      ! get_factor_and_xyz_min_max): the box holds the 3^(s−1) n_α cells centred on
      ! −(3^(s−1) n_α − 1)/2 … +(3^(s−1) n_α − 1)/2, whose union spans ±F_α a_α.
      f_half(1:2) = 0.5_dp*REAL(3**(bs_env%size_lattice_sum_V - 1), dp)* &
                    REAL(bs_env%nkp_grid_chi_eps_W_orig(1:2), dp)
      ! Corners of the box parallelogram, counterclockwise for a right-handed (a_1, a_2).
      corner(:, 1) = -f_half(1)*bs_env%hmat(:, 1) - f_half(2)*bs_env%hmat(:, 2)
      corner(:, 2) = +f_half(1)*bs_env%hmat(:, 1) - f_half(2)*bs_env%hmat(:, 2)
      corner(:, 3) = +f_half(1)*bs_env%hmat(:, 1) + f_half(2)*bs_env%hmat(:, 2)
      corner(:, 4) = -f_half(1)*bs_env%hmat(:, 1) + f_half(2)*bs_env%hmat(:, 2)

      ! Polar form of the area integral: ∫_box d²r f(|r|) with f = 1/|r| becomes the boundary
      ! integral ∮ (r × dr)/|r|, because ∫_0^{s(θ)} ds s (1/s) = s(θ) and s² dθ = r × dr,
      ! with s(θ) the distance from the origin to the box boundary in direction θ.
      ! Gauss-Legendre with n_gl = 32 nodes per edge; the integrand is smooth on every edge
      ! because the origin lies strictly inside the box.
      CALL gamma_corr_gauss_legendre(n_gl, nodes, rule)
      d_box = 0.0_dp
      area_check = 0.0_dp
      DO i_edge = 1, 4
         ! Edge from corner i_edge to the next one, r(t) = midpoint + t·dr, t ∈ [−1, 1].
         dr = 0.5_dp*(corner(:, MODULO(i_edge, 4) + 1) - corner(:, i_edge))
         DO i_node = 1, n_gl
            r = 0.5_dp*(corner(:, MODULO(i_edge, 4) + 1) + corner(:, i_edge)) + nodes(i_node)*dr
            ! weight = w_GL (r × dr/dt)_z: the signed in-plane cross product (x-y components,
            ! i.e. a layer in the xy plane is assumed here).
            weight = rule(i_node)*(r(1)*dr(2) - r(2)*dr(1))
            ! ∮ (r × dr)/|r| accumulates D_box·A; (1/2) ∮ r × dr accumulates the box area.
            d_box = d_box + weight/SQRT(SUM(r(1:2)**2))
            area_check = area_check + 0.5_dp*weight
         END DO
      END DO
      DEALLOCATE (nodes, rule)
      ! D_box = (1/A) ∫_box d²r/|r|, in bohr⁻¹ (same units as V_PQ/(c_P c_Q)).
      d_box = d_box/gc%area
      ! The box area must come out as (2F_1)(2F_2)·A; a sign or orientation error fails here.
      IF (ABS(area_check/(4.0_dp*f_half(1)*f_half(2)*gc%area) - 1.0_dp) > 1.0E-8_dp) &
         CPABORT('Gamma-point correction: the image-box boundary quadrature failed its area check')

      ! ---------------------------------------------------------------------------------------
      ! 2. Report the charge weight before and after, then V <- V − D_box ||c||² ĉ ĉᵀ.
      ! ---------------------------------------------------------------------------------------
      ! Global (row, col) indices of this rank's local block-cyclic elements.
      CALL cp_cfm_get_info(cfm_V, row_indices=row_indices, col_indices=col_indices)

      ! u = V ĉ, from the local elements and one reduction (diagnostic only).
      ALLOCATE (u(n_RI))
      u(:) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO i_col = 1, SIZE(col_indices)
         col = col_indices(i_col)
         DO i_row = 1, SIZE(row_indices)
            row = row_indices(i_row)
            u(row) = u(row) + cfm_V%local_data(i_row, i_col)*gc%charge(col)
         END DO
      END DO
      CALL bs_env%para_env%sum(u)
      ! c_dot_u = ĉᵀ V ĉ, the charge-channel weight of V(0) before the subtraction.
      c_dot_u = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO row = 1, n_RI
         c_dot_u = c_dot_u + gc%charge(row)*u(row)
      END DO
      DEALLOCATE (u)

      ! V_PQ <- V_PQ − D_box |c|² ĉ_P ĉ_Q = V_PQ − D_box c_P c_Q on the local elements.
      ! (weight is reused here: it now holds D_box |c|², no longer a quadrature weight.)
      weight = d_box*gc%charge_norm**2
      DO i_col = 1, SIZE(col_indices)
         col = col_indices(i_col)
         DO i_row = 1, SIZE(row_indices)
            row = row_indices(i_row)
            cfm_V%local_data(i_row, i_col) = cfm_V%local_data(i_row, i_col) &
                                             - weight*gc%charge(row)*gc%charge(col)
         END DO
      END DO

      ! Print once per calculation (gc%reported_projection), although the routine may be called
      ! again for later V builds. ĉᵀ V ĉ after the subtraction is ĉᵀ V ĉ − D_box |c|², since
      ! ĉᵀ ĉ = 1. Only the x half-width is printed; both are equal for a square mesh.
      IF (bs_env%unit_nr > 0 .AND. .NOT. gc%reported_projection) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'Gamma-point correction: divergent G = 0 constant removed from V(q=0)'
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Image box half-widths (cells)', f_half(1)
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'D_box = (1/A) int_box d^2r/|r|', d_box
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Charge weight c^T V c before', REAL(c_dot_u, dp)
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Charge weight c^T V c after', &
            REAL(c_dot_u, dp) - weight
         WRITE (bs_env%unit_nr, '(A)') ' '
         gc%reported_projection = .TRUE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_subtract_divergent_head

END MODULE gw_gamma_correction
