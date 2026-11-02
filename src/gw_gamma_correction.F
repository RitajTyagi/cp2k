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
!>        sum of V_b converges absolutely where the unprojected one does not. The dielectric head
!>        and wings in the extended basis {0, 1, …, N_RI} are then integrated over the WHOLE small
!>        Brillouin zone of the supercell. Nothing is subtracted, unlike the full-k code of
!>        J. Wilhelm (cp2k#6447, gwp42 Eq. 78), because the charge mode was never sampled here in
!>        the first place.
!>
!>        Equations follow J. Wilhelm, *Periodic real-space RI for full-k GW band structures*,
!>        Sec. VI; equation numbers below refer to that manuscript. The derivation of the Γ-only
!>        specialization is in GAMMA_W_CORRECTION_PLAN.md.
!>
!>        **Scope.** 2D slabs (PERIODIC XY) with the Γ-only RI-RS implementation. χ, ε and W stay
!>        at one real k-point, so none of this needs complex arithmetic, and the minimum-image
!>        clipping of χ (cp_cfm_power) has nothing to repair: χ(Γ) is a Gram matrix of
!>        occupied-virtual products and exactly positive semidefinite.
!> \author Ritaj Tyagi
! **************************************************************************************************
MODULE gw_gamma_correction
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set
   USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                              cp_cfm_type
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
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
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE ri_rs_types,                     ONLY: gamma_corr_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_gamma_correction'

   PUBLIC :: gamma_corr_requested, &
             gamma_corr_init, &
             gamma_corr_grid_owner, &
             gamma_corr_build_coulomb_maps, &
             gamma_corr_charge_rhs, &
             gamma_corr_transition_dipoles, &
             gamma_corr_grid_orbitals, &
             gamma_corr_response, &
             gamma_corr_tabulate, &
             gamma_corr_delta_W, &
             gamma_corr_project_coulomb, &
             gamma_corr_auxiliary_charges, &
             gamma_corr_metric_head, &
             gamma_corr_slab_geometry, &
             gamma_corr_symmetric_power

CONTAINS

! **************************************************************************************************
!> \brief Whether the Γ-point correction applies to this calculation.
!>
!>        It is meaningful only when the mesh for χ, ε, W actually contains q = 0, i.e. when every
!>        periodic direction has an ODD number of points (compute_xkp puts a point at the origin
!>        for odd n and half-shifts for even n), and only for a 2D slab, which is the only
!>        geometry whose head/wing construction is implemented here.
!> \param bs_env ...
!> \return TRUE if the correction should be built
! **************************************************************************************************
   FUNCTION gamma_corr_requested(bs_env) RESULT(do_corr)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL                                            :: do_corr

      INTEGER                                            :: i_dim

      do_corr = .FALSE.

      ! Only PERIODIC XY: the 3D construction needs the Ewald ΔV and the 1/q² head instead.
      IF (.NOT. ALL(bs_env%periodic == [1, 1, 0])) RETURN

      ! q = 0 is in the mesh only if every periodic direction is odd.
      DO i_dim = 1, 3
         IF (bs_env%periodic(i_dim) == 0) CYCLE
         IF (MODULO(bs_env%nkp_grid_chi_eps_W_orig(i_dim), 2) == 0) RETURN
      END DO

      do_corr = .TRUE.

   END FUNCTION gamma_corr_requested

! **************************************************************************************************
!> \brief Layer area A = |a_1 × a_2| and the two in-plane reciprocal vectors b_1, b_2 with
!>        b_i · a_j = 2π δ_ij, for an arbitrary (possibly non-orthogonal) XY cell.
!>
!>        The vacuum height never enters: a slab's Fourier transform averages over a layer cell,
!>        not over a 3D cell that mostly contains vacuum. This is the single substitution that
!>        turns the 3D formulas into the 2D ones, 1/√Ω → 1/√A (gwp42 Eqs. 61, 62).
!> \param bs_env ...
!> \param area |a_1 × a_2|
!> \param b_recip the two in-plane reciprocal vectors, b_recip(1:3, 1:2)
! **************************************************************************************************
   SUBROUTINE gamma_corr_slab_geometry(bs_env, area, b_recip)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: area
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(OUT)        :: b_recip

      REAL(KIND=dp), DIMENSION(3)                        :: a_1, a_2, normal

      a_1 = bs_env%hmat(:, 1)
      a_2 = bs_env%hmat(:, 2)

      normal = [a_1(2)*a_2(3) - a_1(3)*a_2(2), &
                a_1(3)*a_2(1) - a_1(1)*a_2(3), &
                a_1(1)*a_2(2) - a_1(2)*a_2(1)]
      area = SQRT(SUM(normal**2))
      CPASSERT(area > 0.0_dp)
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
!> \param bs_env ...
!> \return the metric head
! **************************************************************************************************
   FUNCTION gamma_corr_metric_head(bs_env) RESULT(head)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: head

      SELECT CASE (bs_env%ri_metric%potential_type)
      CASE (do_potential_short)
         head = pi/bs_env%ri_metric%omega**2
      CASE (do_potential_truncated)
         head = twopi*bs_env%ri_metric%cutoff_radius**2
      CASE DEFAULT
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
!> \param qs_env ...
!> \param charge c_P, one entry per RI_AUX basis function, in the global RI ordering
! **************************************************************************************************
   SUBROUTINE gamma_corr_auxiliary_charges(qs_env, charge)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: charge

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_auxiliary_charges'

      INTEGER                                            :: first, handle, i_atom, i_co, i_l, &
                                                            i_pgf, i_set, i_sgf, l_max, n_cart, &
                                                            offset
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: moments, primitive
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set)

      charge(:) = 0.0_dp
      first = 0

      DO i_atom = 1, SIZE(particle_set)
         CALL get_qs_kind(qs_kind_set(particle_set(i_atom)%atomic_kind%kind_number), &
                          basis_set=basis, basis_type='RI_AUX')

         DO i_set = 1, basis%nset
            l_max = basis%lmax(i_set)
            n_cart = ncoset(l_max)
            ALLOCATE (moments(0:l_max), primitive(n_cart*basis%npgf(i_set)))
            primitive(:) = 0.0_dp

            DO i_pgf = 1, basis%npgf(i_set)
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

            offset = basis%first_sgf(1, i_set)
            DO i_sgf = offset, offset + basis%nsgf_set(i_set) - 1
               charge(first + i_sgf) = SUM(basis%sphi(1:SIZE(primitive), i_sgf)*primitive)
            END DO

            DEALLOCATE (moments, primitive)
         END DO

         first = first + basis%nsgf
      END DO

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
!>        relative eigenvalue cut does. The caller passes REGULARIZATION_RI as the threshold, the
!>        same knob the rest of the RI metric uses.
!> \param matrix on entry the symmetric matrix, on exit its p-th power; replicated
!> \param exponent p
!> \param relative_eps modes with λ < relative_eps · max λ are set to zero
!> \param n_dropped optional count of discarded modes
! **************************************************************************************************
   SUBROUTINE gamma_corr_symmetric_power(matrix, exponent, relative_eps, n_dropped)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: exponent, relative_eps
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_dropped

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_symmetric_power'

      INTEGER                                            :: handle, i, info, n, n_drop
      REAL(KIND=dp)                                      :: cut
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: scaled, vectors

      CALL timeset(routineN, handle)

      n = SIZE(matrix, 1)
      CPASSERT(SIZE(matrix, 2) == n)

      ALLOCATE (vectors(n, n), eigenvalues(n), work(MAX(1, 34*n)))
      vectors(:, :) = 0.5_dp*(matrix + TRANSPOSE(matrix))

      CALL dsyev('V', 'U', n, vectors, n, eigenvalues, work, SIZE(work), info)
      IF (info /= 0) CPABORT('Gamma-point correction: dsyev failed')

      cut = relative_eps*MAXVAL(eigenvalues)
      cut = MAX(cut, 0.0_dp)

      ! A^p = U diag(λ^p) U^T on the retained modes; discarded modes stay quenched for every p,
      ! including negative ones, which is what keeps V_b^(-1/2) finite on the null space.
      ALLOCATE (scaled(n, n))
      n_drop = 0
      DO i = 1, n
         IF (eigenvalues(i) <= cut) THEN
            scaled(:, i) = 0.0_dp
            n_drop = n_drop + 1
         ELSE
            scaled(:, i) = vectors(:, i)*eigenvalues(i)**exponent
         END IF
      END DO

      CALL dgemm('N', 'T', n, n, n, 1.0_dp, scaled, n, vectors, n, 0.0_dp, matrix, n)

      IF (PRESENT(n_dropped)) n_dropped = n_drop

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
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE gamma_corr_build_coulomb_maps(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_build_coulomb_maps'

      INTEGER                                            :: handle, i, n_drop_inv, n_drop_sqrt, n_RI
      REAL(KIND=dp)                                      :: norm_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: M_inv, projector, V_b_inv_sqrt, &
                                                            V_b_sqrt, V_zero, work
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI

      CALL gamma_corr_slab_geometry(bs_env, gc%area, gc%b_recip)
      gc%metric_head = gamma_corr_metric_head(bs_env)

      ! ---------------------------------------------------------------------------------------
      ! 1. V_PQ(0): the finite real-space 2D lattice sum of the bare Coulomb operator (Eq. 63).
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (V_zero(n_RI, n_RI))
      CALL gamma_corr_coulomb_at_zero(qs_env, bs_env, V_zero)

      ! ---------------------------------------------------------------------------------------
      ! 2. c_P and the neutral projector Π_PQ = δ_PQ − c_P c_Q/Σ_S|c_S|² (Eq. 64).
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (gc%charge(n_RI))
      CALL gamma_corr_auxiliary_charges(qs_env, gc%charge)

      norm_c = SQRT(DOT_PRODUCT(gc%charge, gc%charge))
      IF (norm_c <= 0.0_dp) &
         CPABORT('Gamma-point correction: the RI_AUX basis carries no charge at all')
      gc%charge(:) = gc%charge/norm_c

      ALLOCATE (projector(n_RI, n_RI))
      DO i = 1, n_RI
         projector(:, i) = -gc%charge(:)*gc%charge(i)
         projector(i, i) = projector(i, i) + 1.0_dp
      END DO

      ! ---------------------------------------------------------------------------------------
      ! 3. V_b = Π V(0) Π and its ±1/2 powers on the neutral subspace (Eq. 65).
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (work(n_RI, n_RI), V_b_sqrt(n_RI, n_RI), V_b_inv_sqrt(n_RI, n_RI))
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, V_zero, n_RI, projector, n_RI, &
                 0.0_dp, work, n_RI)
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, projector, n_RI, work, n_RI, &
                 0.0_dp, V_b_sqrt, n_RI)

      V_b_inv_sqrt(:, :) = V_b_sqrt
      CALL gamma_corr_symmetric_power(V_b_sqrt, 0.5_dp, bs_env%regularization_RI, n_drop_sqrt)
      CALL gamma_corr_symmetric_power(V_b_inv_sqrt, -0.5_dp, bs_env%regularization_RI, n_drop_inv)
      CPASSERT(n_drop_sqrt == n_drop_inv)

      ! ---------------------------------------------------------------------------------------
      ! 4. M^-1(0): NOT rebuilt here. bs_env%fm_Minv_Gamma already holds the inverse metric at Γ
      !    that compute_W_MIC folds into W, built with regularization_RI deliberately switched
      !    off so that the RI-RS and tensor paths agree. Building a second, differently
      !    regularized M^-1 for the correction is exactly the mistake that cost 3.8 eV in Σ^c
      !    once before, so the correction reuses this one verbatim.
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (M_inv(n_RI, n_RI))
      CALL cp_fm_get_submatrix(bs_env%fm_Minv_Gamma, M_inv)

      ! ---------------------------------------------------------------------------------------
      ! 5. F_density = M^-1 Π V_b^(1/2)  and  F_extern = M^-1 V(0) Π V_b^(-1/2)  (Eqs. 66, 67).
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (gc%F_density(n_RI, n_RI), gc%F_extern(n_RI, n_RI))

      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, projector, n_RI, V_b_sqrt, n_RI, &
                 0.0_dp, work, n_RI)
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, M_inv, n_RI, work, n_RI, &
                 0.0_dp, gc%F_density, n_RI)

      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, projector, n_RI, V_b_inv_sqrt, n_RI, &
                 0.0_dp, work, n_RI)
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, V_zero, n_RI, work, n_RI, &
                 0.0_dp, V_b_sqrt, n_RI)
      CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, M_inv, n_RI, V_b_sqrt, n_RI, &
                 0.0_dp, gc%F_extern, n_RI)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A)') &
            'Gamma-point correction of W (2D): charge channel separated from the RI basis'
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Layer area A (bohr^2)', gc%area
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'RI metric head int dr v(r)', gc%metric_head
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'Norm of the charge vector |c|', norm_c
         WRITE (bs_env%unit_nr, '(T4,A,T72,I10)') 'Null modes of the neutral body V_b', n_drop_sqrt
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
!> \param qs_env ...
!> \param bs_env ...
!> \param V_zero V_PQ(0), replicated on every rank
! **************************************************************************************************
   SUBROUTINE gamma_corr_coulomb_at_zero(qs_env, bs_env, V_zero)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: V_zero

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_coulomb_at_zero'

      INTEGER                                            :: handle, re_im
      REAL(KIND=dp)                                      :: max_imag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: V_imag
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V
      TYPE(kpoint_type), POINTER                         :: kpoints_zero
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                      qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)

      NULLIFY (kpoints_zero)
      CALL kpoint_create(kpoints_zero)
      ALLOCATE (kpoints_zero%xkp(3, 1))
      kpoints_zero%xkp(:, 1) = 0.0_dp
      kpoints_zero%nkp = 1
      kpoints_zero%kp_scheme = "GENERAL"
      ! The box of the lattice sum comes from nkp_grid, so keep the production mesh here.
      kpoints_zero%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      NULLIFY (mat_V)
      ALLOCATE (mat_V(1:1, 2))
      DO re_im = 1, 2
         NULLIFY (mat_V(1, re_im)%matrix)
         ALLOCATE (mat_V(1, re_im)%matrix)
         CALL dbcsr_create(mat_V(1, re_im)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V(1, re_im)%matrix)
         CALL dbcsr_set(mat_V(1, re_im)%matrix, 0.0_dp)
      END DO

      CALL build_2c_coulomb_matrix_kp(mat_V, kpoints_zero, basis_type="RI_AUX", cell=cell, &
                                      particle_set=particle_set, qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=bs_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=1, ikp_end=1)

      CALL cp_fm_create(fm_work, bs_env%fm_RI_RI%matrix_struct)

      CALL cp_fm_set_all(fm_work, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_V(1, 1)%matrix, fm_work)
      CALL cp_fm_get_submatrix(fm_work, V_zero)

      ALLOCATE (V_imag(SIZE(V_zero, 1), SIZE(V_zero, 2)))
      CALL cp_fm_set_all(fm_work, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_V(1, 2)%matrix, fm_work)
      CALL cp_fm_get_submatrix(fm_work, V_imag)
      max_imag = MAXVAL(ABS(V_imag))

      CALL cp_fm_release(fm_work)
      CALL dbcsr_deallocate_matrix_set(mat_V)
      DEALLOCATE (kpoints_zero%xkp)
      CALL kpoint_release(kpoints_zero)

      ! Im V(0) = 0 for a centrosymmetric image box; a nonzero value means the box is lopsided
      ! and the whole real-arithmetic argument of this module would be void.
      IF (max_imag > 1.0E-10_dp*MAX(1.0_dp, MAXVAL(ABS(V_zero)))) THEN
         CPABORT('Gamma-point correction: V(q=0) came back complex')
      END IF
      DEALLOCATE (V_imag)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_coulomb_at_zero

! **************************************************************************************************
!> \brief Right-hand side of the extra normal equation that fits the normalized constant
!>        function, gwp42 Eq. 61 for 2D:
!>
!>          Σ_ℓ' D_ℓℓ'(0) Z_ℓ'0 = d_ℓ0 ,      d_ℓ0 = (1/√A) ∫ dr f_ℓ(r)²
!>
!>        with the same f_ℓ(r) = Σ_µR ϕ_µR(r_ℓ) ϕ_µR(r) that the regular d_ℓP quadrature uses.
!>        Written out, the integral is an AO overlap contraction at the grid point,
!>
!>          ∫ dr f_ℓ(r)² = Σ_µν ϕ_µ(r_ℓ) S_µν ϕ_ν(r_ℓ) = ϕ_ℓᵀ S ϕ_ℓ ,
!>
!>        where S is the Γ-summed AO overlap, consistent with ϕ_µ(r_ℓ) already being the Γ Bloch
!>        sum over images. This routine never forms S densely. With the Cholesky factor
!>        S = Uᵀ U from cp_fm_cholesky_decompose,
!>
!>          ϕ_ℓᵀ S ϕ_ℓ = |U ϕ_ℓ|² = Σ_k (Σ_µ ϕ_µ(r_ℓ) U_kµ)² ,
!>
!>        so one sparse product G = Φ Uᵀ and the row norms of G give the whole vector. Φ is
!>        mat_phi_mu_l, already built and distributed over grid rows.
!>
!>        The result is scaled by the metric head ∫dr v_ω(r), which is what actually comes out of
!>        integrating the constant function against the RI metric kernel:
!>
!>          d_ℓ0 = Σ_µν ϕ_µ(r_ℓ) ϕ_ν(r_ℓ) (µν|ϕ_0)_ω
!>               = (1/√A) [∫ dr'' v_ω(r'')] Σ_µν ϕ_µ(r_ℓ) S_µν ϕ_ν(r_ℓ) .
!>
!>        The caller divides the fitted coefficients by the same factor again. That is a no-op
!>        for an exact solve but not for the Tikhonov-regularized one: λ is fixed and not scale
!>        invariant, so keeping the charge column at its physical magnitude during the solve is
!>        what makes it regularized on the same footing as the Gaussian columns.
!> \param bs_env ...
!> \param d_l0 d_ℓ0 for every grid point, replicated on every rank
! **************************************************************************************************
   SUBROUTINE gamma_corr_charge_rhs(bs_env, d_l0)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: d_l0

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_charge_rhs'

      INTEGER                                            :: col, handle, i_col, i_row, info, &
                                                            n_grid, row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: prefactor
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), POINTER :: local_data
      TYPE(cp_fm_type)                                   :: fm_U
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: mat_G, mat_U
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_grid = bs_env%ri_rs%n_grid_points
      CPASSERT(n_grid > 0)
      CPASSERT(gc%area > 0.0_dp .AND. gc%metric_head > 0.0_dp)

      ! ---------------------------------------------------------------------------------------
      ! 1. S = Uᵀ U. cp_fm_cholesky_decompose calls pdpotrf('U'), so U is the upper triangle and
      !    the strict lower triangle is left as workspace garbage; it has to be zeroed before the
      !    matrix is used as a factor.
      ! ---------------------------------------------------------------------------------------
      CALL cp_fm_create(fm_U, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_to_fm(bs_env%fm_s_Gamma, fm_U)
      CALL cp_fm_cholesky_decompose(fm_U, info_out=info)
      IF (info /= 0) &
         CPABORT('Gamma-point correction: the AO overlap is not positive definite')

      CALL cp_fm_get_info(fm_U, row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data)
      DO i_col = 1, SIZE(col_indices)
         col = col_indices(i_col)
         DO i_row = 1, SIZE(row_indices)
            row = row_indices(i_row)
            IF (row > col) local_data(i_row, i_col) = 0.0_dp
         END DO
      END DO

      CALL dbcsr_create(mat_U, template=bs_env%mat_ao_ao%matrix)
      CALL copy_fm_to_dbcsr(fm_U, mat_U)
      CALL cp_fm_release(fm_U)

      ! ---------------------------------------------------------------------------------------
      ! 2. G_ℓk = Σ_µ ϕ_µ(r_l) U_kµ, i.e. G = Φ Uᵀ.
      ! ---------------------------------------------------------------------------------------
      CALL dbcsr_create(mat_G, template=bs_env%ri_rs%mat_phi_mu_l)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, bs_env%ri_rs%mat_phi_mu_l, mat_U, &
                          0.0_dp, mat_G, filter_eps=bs_env%eps_filter)
      CALL dbcsr_release(mat_U)

      ! ---------------------------------------------------------------------------------------
      ! 3. d_ℓ0 = (metric_head/√A) Σ_k G_ℓk². Each rank sees only its own blocks of G, so the
      !    row sums are completed by one reduction.
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (d_l0(n_grid))
      d_l0(:) = 0.0_dp

      CALL dbcsr_iterator_start(iter, mat_G)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_offset=i_row)
         DO i_col = 1, SIZE(block, 2)
            d_l0(i_row:i_row + SIZE(block, 1) - 1) = &
               d_l0(i_row:i_row + SIZE(block, 1) - 1) + block(:, i_col)**2
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_release(mat_G)

      CALL bs_env%para_env%sum(d_l0)

      prefactor = gc%metric_head/SQRT(gc%area)
      d_l0(:) = prefactor*d_l0

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_charge_rhs

! **************************************************************************************************
!> \brief Decide whether the Γ-point correction runs and set the geometry it needs, early enough
!>        that the Z_lP fit can already carry the extra charge-channel right-hand side.
!>
!>        Only the cheap, purely geometric quantities are set here (layer area, in-plane
!>        reciprocal vectors, metric head). The Coulomb side -- V(0), the projector and the two
!>        maps -- is built later by gamma_corr_build_coulomb_maps, after the RI metric exists.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE gamma_corr_init(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      TYPE(gamma_corr_type), POINTER                     :: gc

      gc => bs_env%ri_rs%gamma_corr

      gc%enabled = gamma_corr_requested(bs_env)
      IF (.NOT. gc%enabled) RETURN

      CALL gamma_corr_slab_geometry(bs_env, gc%area, gc%b_recip)
      gc%metric_head = gamma_corr_metric_head(bs_env)

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
!> \brief Nearest atom to every grid point, under the minimum image convention.
!>
!>        The RI-RS grid is atom centered, but the integration spheres overlap, so a grid point
!>        r_l lies inside the sphere of several atoms P and the local fit of the single global
!>        charge channel produces a value for it once per sphere. Exactly one of them has to be
!>        kept, and the natural choice is the sphere in which r_l sits deepest, i.e. the one
!>        centered on the nearest atom, because that is where truncating the couplings D_ll' at
!>        the sphere boundary does the least damage. Ties go to the lowest atom index, which
!>        makes the assignment deterministic and rank independent.
!> \param bs_env ...
!> \param owner owner(l) = index of the atom whose fit supplies Z_l0 for grid point l
! **************************************************************************************************
   SUBROUTINE gamma_corr_grid_owner(bs_env, owner)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: owner

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_grid_owner'

      INTEGER                                            :: handle, i_atom, i_grid, n_atom, n_grid
      REAL(KIND=dp)                                      :: best, dist
      REAL(KIND=dp), DIMENSION(3)                        :: delta
      TYPE(cell_type), POINTER                           :: cell

      CALL timeset(routineN, handle)

      cell => bs_env%ri_rs%cell
      n_grid = bs_env%ri_rs%n_grid_points
      n_atom = bs_env%n_atom
      CPASSERT(n_grid > 0 .AND. n_atom > 0)

      ALLOCATE (owner(n_grid))
      owner(:) = 0

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_grid, n_atom, bs_env, cell, owner) &
      !$OMP PRIVATE(i_grid, i_atom, delta, dist, best) &
      !$OMP SCHEDULE(STATIC)
      DO i_grid = 1, n_grid
         best = HUGE(1.0_dp)
         DO i_atom = 1, n_atom
            delta = pbc(bs_env%ri_rs%grid_points(:, i_grid), &
                        bs_env%ri_rs%particle_set(i_atom)%r(:), cell)
            dist = SUM(delta*delta)
            IF (dist < best) THEN
               best = dist
               owner(i_grid) = i_atom
            END IF
         END DO
      END DO
      !$OMP END PARALLEL DO

      CPASSERT(ALL(owner > 0))

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_grid_owner

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
!>        result must be checked against cell size rather than assumed. The dipole integrals are
!>        taken with minimum_image=.TRUE., which references every pair to its nearest image and
!>        is the only consistent choice on a supercell.
!> \param qs_env ...
!> \param bs_env ...
!> \param mu µ^α_ai, indices (α, i, a, spin), replicated on every rank
! **************************************************************************************************
   SUBROUTINE gamma_corr_transition_dipoles(qs_env, bs_env, mu)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :), &
         INTENT(OUT)                                     :: mu

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_transition_dipoles'

      INTEGER                                            :: handle, i_dir, i_occ, i_spin, i_vir, &
                                                            n_ao, n_mo, n_occ, n_vir
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mu_mo
      TYPE(cp_fm_struct_type), POINTER                   :: struct_mo_mo
      TYPE(cp_fm_type)                                   :: fm_mo_mo, fm_r_ao, fm_r_mo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments

      CALL timeset(routineN, handle)

      n_ao = bs_env%n_ao
      n_mo = bs_env%n_mo_retained
      IF (n_mo <= 0) n_mo = n_ao

      ! r^α_µν for α = x, y, z, referenced to the cell origin with minimum-image pairing.
      NULLIFY (moments)
      ALLOCATE (moments(3))
      DO i_dir = 1, 3
         NULLIFY (moments(i_dir)%matrix)
         ALLOCATE (moments(i_dir)%matrix)
         CALL dbcsr_create(moments(i_dir)%matrix, template=bs_env%mat_ao_ao%matrix)
         CALL dbcsr_reserve_all_blocks(moments(i_dir)%matrix)
         CALL dbcsr_set(moments(i_dir)%matrix, 0.0_dp)
      END DO
      CALL build_local_moment_matrix(qs_env, moments, 1, &
                                     ref_point=[0.0_dp, 0.0_dp, 0.0_dp], minimum_image=.TRUE.)

      n_occ = MAXVAL(bs_env%n_occ(1:bs_env%n_spin))
      n_vir = n_mo - MINVAL(bs_env%n_occ(1:bs_env%n_spin))
      CPASSERT(n_occ > 0 .AND. n_vir > 0)
      ALLOCATE (mu(3, n_occ, n_vir, bs_env%n_spin))
      mu(:, :, :, :) = 0.0_dp

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
            CALL cp_fm_set_all(fm_r_ao, 0.0_dp)
            CALL copy_dbcsr_to_fm(moments(i_dir)%matrix, fm_r_ao)
            CALL parallel_gemm('N', 'N', n_ao, n_mo, n_ao, 1.0_dp, fm_r_ao, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), 0.0_dp, fm_r_mo)
            CALL parallel_gemm('T', 'N', n_mo, n_mo, n_ao, 1.0_dp, &
                               bs_env%fm_mo_coeff_Gamma(i_spin), &
                               fm_r_mo, 0.0_dp, fm_mo_mo)
            CALL cp_fm_get_submatrix(fm_mo_mo, mu_mo)

            ! Keep only the occupied-empty block; i runs over occupied, a over empty.
            DO i_vir = 1, n_mo - bs_env%n_occ(i_spin)
               DO i_occ = 1, bs_env%n_occ(i_spin)
                  mu(i_dir, i_occ, i_vir, i_spin) = &
                     mu_mo(i_occ, bs_env%n_occ(i_spin) + i_vir)
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

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_transition_dipoles

! **************************************************************************************************
!> \brief Kohn-Sham orbitals on the RI-RS grid, ψ_n(r_l) = Σ_µ ϕ_µ(r_l) C_µn (gwp42 Eq. 17 at Γ).
!>
!>        Needed only by the wing response, Eq. 48, which unlike the body cannot be written as a
!>        product of two Green's functions: the dipole couples i and a, so the occupied and empty
!>        sums do not factor and the band-resolved orbitals have to appear explicitly.
!> \param bs_env ...
!> \param i_spin spin channel
!> \param mat_psi ψ_n(r_l), grid rows x MO columns, inheriting mat_phi_mu_l's row distribution
! **************************************************************************************************
   SUBROUTINE gamma_corr_grid_orbitals(bs_env, i_spin, mat_psi)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: i_spin
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_psi

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_grid_orbitals'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: mat_C

      CALL timeset(routineN, handle)

      CALL dbcsr_create(mat_C, template=bs_env%mat_ao_ao%matrix)
      ! Spin 1 only here; the caller loops spins and calls this once per spin channel.
      CALL copy_fm_to_dbcsr(bs_env%fm_mo_coeff_Gamma(i_spin), mat_C)

      CALL dbcsr_create(mat_psi, template=bs_env%ri_rs%mat_phi_mu_l)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, bs_env%ri_rs%mat_phi_mu_l, mat_C, &
                          0.0_dp, mat_psi, filter_eps=bs_env%eps_filter)

      CALL dbcsr_release(mat_C)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_grid_orbitals

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
!>        The exponential factorizes, e^(−Δ_ia τ) = e^(−|ε_a−ε_F|τ) e^(−|ε_i−ε_F|τ), which is what
!>        makes X_ℓα a pair of gemms per direction rather than a triple loop: first
!>        Λ^α_ℓi = Σ_a µ^α_ai ψ_a(r_ℓ) e^(−|ε_a−ε_F|τ), then X_ℓα = Σ_i ψ_i(r_ℓ) e^(−|ε_i−ε_F|τ)
!>        Λ^α_ℓi row by row. Cost 3·n_grid·n_occ·n_vir per time point.
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
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau χ_PQ(Γ,iτ), the regular polarizability the driver already built
! **************************************************************************************************
   SUBROUTINE gamma_corr_response(bs_env, qs_env, mat_chi_Gamma_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_response'

      INTEGER                                            :: handle, i_dir, i_spin, i_t, i_w, j_dir, &
                                                            n_grid, n_occ, n_RI, n_t, n_vir
      REAL(KIND=dp)                                      :: spin_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: chi_PQ, E_freq, work_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: X_freq, X_time, Y_time, mu_spin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: mu
      TYPE(cp_fm_type)                                   :: fm_chi
      TYPE(dbcsr_type)                                   :: mat_psi
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      n_RI = bs_env%n_RI
      n_grid = bs_env%ri_rs%n_grid_points
      n_t = bs_env%num_time_freq_points

      ! Eqs. 48/49 carry an explicit factor two for the two spin states of a closed shell; for a
      ! spin-polarized run that factor is the spin sum instead, hence 2/n_spin per channel.
      spin_factor = 2.0_dp/REAL(bs_env%n_spin, dp)

      CALL gamma_corr_transition_dipoles(qs_env, bs_env, mu)
      n_occ = SIZE(mu, 2)
      n_vir = SIZE(mu, 3)

      ALLOCATE (Y_time(3, 3, n_t), X_time(n_grid, 3, n_t))
      Y_time(:, :, :) = 0.0_dp
      X_time(:, :, :) = 0.0_dp

      DO i_spin = 1, bs_env%n_spin
         ALLOCATE (mu_spin(3, n_occ, n_vir))
         mu_spin(:, :, :) = mu(:, :, :, i_spin)

         CALL gamma_corr_grid_orbitals(bs_env, i_spin, mat_psi)
         CALL gamma_corr_accumulate_response(bs_env, i_spin, mu_spin, mat_psi, spin_factor, &
                                             Y_time, X_time)
         CALL dbcsr_release(mat_psi)
         DEALLOCATE (mu_spin)
      END DO
      DEALLOCATE (mu)

      CALL bs_env%para_env%sum(X_time)

      ! ---------------------------------------------------------------------------------------
      ! Cosine transform iτ -> iω with the same minimax weights as the regular response.
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
      ALLOCATE (gc%T_wing(n_grid, 3, n_t))
      gc%T_wing(:, :, :) = 0.0_dp
      ALLOCATE (chi_PQ(n_RI, n_RI), work_RI(n_RI, n_RI), E_freq(n_RI, 3))
      CALL cp_fm_create(fm_chi, bs_env%fm_RI_RI%matrix_struct)

      DO i_w = 1, n_t
         ! χ_PQ(Γ,iω) from the stored χ(Γ,iτ): the same transform the regular W path uses.
         CALL gamma_corr_chi_PQ_at_freq(bs_env, fm_chi, i_w, mat_chi_Gamma_tau, chi_PQ)

         ! B_PQ = δ_PQ + F_densityᵀ χ_PQ F_density (Eq. 72).
         CALL dgemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, chi_PQ, n_RI, gc%F_density, n_RI, &
                    0.0_dp, work_RI, n_RI)
         CALL dgemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, gc%F_density, n_RI, work_RI, n_RI, &
                    0.0_dp, chi_PQ, n_RI)
         DO i_dir = 1, n_RI
            chi_PQ(i_dir, i_dir) = chi_PQ(i_dir, i_dir) + 1.0_dp
         END DO

         ! E_Pα = F_densityᵀ (Zᵀ X)_Qα (Eq. 71 without its √(2πq), which the cell integral adds).
         CALL gamma_corr_project_grid_to_RI(bs_env, X_freq(:, :, i_w), work_RI(:, 1:3))
         CALL dgemm('T', 'N', n_RI, 3, n_RI, 1.0_dp, gc%F_density, n_RI, work_RI, n_RI, &
                    0.0_dp, E_freq, n_RI)

         ! S_αβ = Y_αβ − Eᵀ B⁻¹ E, by solving B x = E rather than inverting B.
         CALL gamma_corr_body_solve(chi_PQ, E_freq, work_RI(:, 1:3))
         DO j_dir = 1, 3
            DO i_dir = 1, 3
               gc%S_head(i_dir, j_dir, i_w) = gc%S_head(i_dir, j_dir, i_w) - &
                                              DOT_PRODUCT(E_freq(:, i_dir), work_RI(:, j_dir))
            END DO
         END DO

         ! T^w_ℓα = Σ_Q Z_ℓQ [F_extern B⁻¹ E]_Qα: push B⁻¹E through the Coulomb-side map and
         ! back onto the grid.
         CALL dgemm('N', 'N', n_RI, 3, n_RI, 1.0_dp, gc%F_extern, n_RI, work_RI, n_RI, &
                    0.0_dp, E_freq, n_RI)
         CALL gamma_corr_project_RI_to_grid(bs_env, E_freq, gc%T_wing(:, :, i_w))
      END DO

      CALL cp_fm_release(fm_chi)
      DEALLOCATE (chi_PQ, work_RI, E_freq, X_freq)

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
!>        vectors X_ℓα(iτ), using the factorized exponential described in gamma_corr_response.
!>
!>        Each rank owns a set of grid-row blocks of ψ, so X comes out partial and the caller
!>        reduces it. Y is a pure band quantity and is computed redundantly on every rank, which
!>        is cheaper than communicating nine numbers.
!> \param bs_env ...
!> \param i_spin ...
!> \param mu µ^α_ai for this spin
!> \param mat_psi ψ_n(r_l) for this spin
!> \param spin_factor 2/n_spin
!> \param Y_time accumulated head tensor, Y(α,β,iτ)
!> \param X_time accumulated wing grid vectors, X(l,α,iτ), partial over grid rows
! **************************************************************************************************
   SUBROUTINE gamma_corr_accumulate_response(bs_env, i_spin, mu, mat_psi, spin_factor, &
                                             Y_time, X_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: i_spin
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: mu
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_psi
      REAL(KIND=dp), INTENT(IN)                          :: spin_factor
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Y_time, X_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_accumulate_response'

      INTEGER :: col, col_start, handle, i_a, i_dir, i_i, i_local, i_t, j_dir, n_local, n_mo, &
         n_occ, n_vir, row, row_start
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_of_global, row_global
      REAL(KIND=dp)                                      :: e_fermi, tau, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exp_occ, exp_vir
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: lambda, psi_dense, psi_occ, psi_vir
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      n_occ = bs_env%n_occ(i_spin)
      n_vir = SIZE(mu, 3)
      n_mo = n_occ + n_vir
      e_fermi = bs_env%e_fermi(i_spin)

      ! ---------------------------------------------------------------------------------------
      ! mat_psi's COLUMNS are blocked by AO atom, so one block covers only part of the MO range,
      ! while the occupied-empty contraction below couples indices that generally sit in
      ! different column blocks of the same row. The blocks this rank owns are therefore first
      ! gathered into a dense (local grid rows) x (MO) array; everything after that is plain
      ! BLAS. Two passes: first map the global grid rows this rank holds onto local positions,
      ! then scatter the blocks in.
      ! ---------------------------------------------------------------------------------------
      ALLOCATE (local_of_global(bs_env%ri_rs%n_grid_points))
      local_of_global(:) = 0

      CALL dbcsr_iterator_start(iter, mat_psi)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_offset=row_start)
         DO i_i = 1, SIZE(block, 1)
            local_of_global(row_start + i_i - 1) = 1
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      n_local = 0
      DO i_i = 1, bs_env%ri_rs%n_grid_points
         IF (local_of_global(i_i) == 0) CYCLE
         n_local = n_local + 1
         local_of_global(i_i) = n_local
      END DO

      ALLOCATE (row_global(MAX(n_local, 1)))
      DO i_i = 1, bs_env%ri_rs%n_grid_points
         IF (local_of_global(i_i) > 0) row_global(local_of_global(i_i)) = i_i
      END DO

      ALLOCATE (psi_dense(MAX(n_local, 1), n_mo))
      psi_dense(:, :) = 0.0_dp

      CALL dbcsr_iterator_start(iter, mat_psi)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_offset=row_start, col_offset=col_start)
         DO i_i = 1, SIZE(block, 2)
            IF (col_start + i_i - 1 > n_mo) CYCLE
            DO i_a = 1, SIZE(block, 1)
               psi_dense(local_of_global(row_start + i_a - 1), col_start + i_i - 1) = &
                  block(i_a, i_i)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
      DEALLOCATE (local_of_global)

      ALLOCATE (exp_occ(n_occ), exp_vir(n_vir))
      ALLOCATE (psi_occ(MAX(n_local, 1), n_occ), psi_vir(MAX(n_local, 1), n_vir), &
                lambda(MAX(n_local, 1), n_occ))

      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)

         ! e^(-Delta_ia tau) = e^(-|eps_a - eps_F| tau) e^(-|eps_i - eps_F| tau)
         DO i_i = 1, n_occ
            exp_occ(i_i) = EXP(-ABS(bs_env%eigenval_scf_Gamma(i_i, i_spin) - e_fermi)*tau)
         END DO
         DO i_a = 1, n_vir
            exp_vir(i_a) = EXP(-ABS(bs_env%eigenval_scf_Gamma(n_occ + i_a, i_spin) - e_fermi)*tau)
         END DO

         ! Y_ab(i tau) = (2/n_spin/A) sum_ia mu^a_ai mu^b_ai e^(-Delta_ia tau)      [Eq. 49]
         ! A band quantity, so computed redundantly on every rank: cheaper than communicating
         ! nine numbers.
         DO i_a = 1, n_vir
            DO i_i = 1, n_occ
               weight = spin_factor*exp_occ(i_i)*exp_vir(i_a)/bs_env%ri_rs%gamma_corr%area
               DO j_dir = 1, 3
                  DO i_dir = 1, 3
                     Y_time(i_dir, j_dir, i_t) = Y_time(i_dir, j_dir, i_t) + &
                                                 weight*mu(i_dir, i_i, i_a)*mu(j_dir, i_i, i_a)
                  END DO
               END DO
            END DO
         END DO

         IF (n_local == 0) CYCLE

         ! X_la(i tau) = (2/n_spin/sqrt(A)) sum_ia psi_i(r_l) psi_a(r_l) mu^a_ai e^(-Delta_ia tau)
         ! Lambda^a_li = sum_a mu^a_ai psi_a(r_l) e^(-...) is one gemm per direction; the
         ! remaining occupied index is contracted row by row.                       [Eq. 48]
         DO i_i = 1, n_occ
            psi_occ(1:n_local, i_i) = psi_dense(1:n_local, i_i)*exp_occ(i_i)
         END DO
         DO i_a = 1, n_vir
            psi_vir(1:n_local, i_a) = psi_dense(1:n_local, n_occ + i_a)*exp_vir(i_a)
         END DO

         DO i_dir = 1, 3
            CALL dgemm('N', 'T', n_local, n_occ, n_vir, 1.0_dp, psi_vir, SIZE(psi_vir, 1), &
                       mu(i_dir, :, :), n_occ, 0.0_dp, lambda, SIZE(lambda, 1))
            DO i_i = 1, n_occ
               DO i_local = 1, n_local
                  X_time(row_global(i_local), i_dir, i_t) = &
                     X_time(row_global(i_local), i_dir, i_t) + &
                     spin_factor/SQRT(bs_env%ri_rs%gamma_corr%area)* &
                     psi_occ(i_local, i_i)*lambda(i_local, i_i)
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (exp_occ, exp_vir, psi_occ, psi_vir, lambda, psi_dense, row_global)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_accumulate_response

! **************************************************************************************************
!> \brief χ_PQ(Γ,iω_j) from the stored χ_PQ(Γ,iτ), replicated as a dense block.
!> \param bs_env ...
!> \param fm_chi scratch full matrix with the RI x RI structure
!> \param j_w frequency index
!> \param mat_chi_Gamma_tau ...
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

      CALL dbcsr_set(bs_env%mat_RI_RI%matrix, 0.0_dp)
      DO i_t = 1, bs_env%num_time_freq_points
         weight = bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, i_t)
         CALL dbcsr_add(bs_env%mat_RI_RI%matrix, mat_chi_Gamma_tau(i_t)%matrix, 1.0_dp, weight)
      END DO
      CALL copy_dbcsr_to_fm(bs_env%mat_RI_RI%matrix, fm_chi)
      CALL cp_fm_get_submatrix(fm_chi, chi_PQ)

   END SUBROUTINE gamma_corr_chi_PQ_at_freq

! **************************************************************************************************
!> \brief (Zᵀ X)_Qα = Σ_ℓ Z_ℓQ X_ℓα for a three-column grid quantity.
!>
!>        Done by walking the local blocks of the sparse Z rather than by a DBCSR product, which
!>        would need X wrapped as a three-column matrix with matching block structure.
!> \param bs_env ...
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

      CALL bs_env%para_env%sum(X_RI)

   END SUBROUTINE gamma_corr_project_grid_to_RI

! **************************************************************************************************
!> \brief G_ℓα = Σ_Q Z_ℓQ G_Qα, the transpose of gamma_corr_project_grid_to_RI.
!>
!>        The result is partial over grid rows on each rank and completed by one reduction, so
!>        that every rank afterwards holds the whole grid vector; the cell-integral application
!>        needs arbitrary grid rows, not just its own.
!> \param bs_env ...
!> \param G_RI n_RI x 3, replicated
!> \param G_grid result, grid x 3, replicated
! **************************************************************************************************
   SUBROUTINE gamma_corr_project_RI_to_grid(bs_env, G_RI, G_grid)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: G_RI
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: G_grid

      INTEGER                                            :: col, col_start, i_col, i_dir, row, &
                                                            row_start
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      G_grid(:, :) = 0.0_dp

      CALL dbcsr_iterator_start(iter, bs_env%ri_rs%mat_Z_lP)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_offset=row_start, col_offset=col_start)
         DO i_dir = 1, 3
            DO i_col = 1, SIZE(block, 2)
               G_grid(row_start:row_start + SIZE(block, 1) - 1, i_dir) = &
                  G_grid(row_start:row_start + SIZE(block, 1) - 1, i_dir) + &
                  block(:, i_col)*G_RI(col_start + i_col - 1, i_dir)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL bs_env%para_env%sum(G_grid)

   END SUBROUTINE gamma_corr_project_RI_to_grid

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
      x(:, 1:3) = E(:, 1:3)
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
!> \param bs_env ...
!> \param n_tab points per direction of the tabulation mesh
! **************************************************************************************************
   SUBROUTINE gamma_corr_tabulate(bs_env, n_tab)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: n_tab

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_tabulate'
      INTEGER, PARAMETER                                 :: n_face = 24, n_radial = 32

      INTEGER                                            :: handle, i_dir, i_node, i_tab, i_u, &
                                                            i_v, i_w, n_dir, n_t
      REAL(KIND=dp) :: a_scalar, area_bz, denom, jac, proj, q, q_max, radial_w, start_time, t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dir_radius, dir_weight, nodes, rule
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dir_hat
      REAL(KIND=dp), DIMENSION(3)                        :: R_cart
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      start_time = m_walltime()
      gc => bs_env%ri_rs%gamma_corr
      n_t = bs_env%num_time_freq_points
      area_bz = twopi**2/gc%area

      CALL gamma_corr_bz_boundary(bs_env, n_face, dir_hat, dir_weight, dir_radius)
      n_dir = SIZE(dir_weight)
      CALL gamma_corr_gauss_legendre(n_radial, nodes, rule)

      gc%n_tab = [n_tab, n_tab]
      ALLOCATE (gc%K_head(n_tab*n_tab, n_t), gc%K_wing(n_tab*n_tab, 3, n_t))
      gc%K_head(:, :) = 0.0_dp
      gc%K_wing(:, :, :) = 0.0_dp

      ! Table points are evenly spaced in fractional coordinates over the full minimum-image
      ! range [-1/2, 1/2] along each in-plane lattice vector, inclusive of both ends so that the
      ! interpolation never has to extrapolate.
      DO i_tab = bs_env%para_env%mepos + 1, n_tab*n_tab, bs_env%para_env%num_pe
         i_u = MODULO(i_tab - 1, n_tab)
         i_v = (i_tab - 1)/n_tab
         R_cart = bs_env%hmat(:, 1)*(REAL(i_u, dp)/REAL(n_tab - 1, dp) - 0.5_dp) + &
                  bs_env%hmat(:, 2)*(REAL(i_v, dp)/REAL(n_tab - 1, dp) - 0.5_dp)

         DO i_dir = 1, n_dir
            q_max = dir_radius(i_dir)
            proj = DOT_PRODUCT(dir_hat(:, i_dir), R_cart)
            jac = dir_weight(i_dir)/area_bz

            DO i_w = 1, n_t
               ! a(q_hat, i_w) = q_hat^T S q_hat
               a_scalar = DOT_PRODUCT(dir_hat(:, i_dir), &
                                      MATMUL(gc%S_head(:, :, i_w), dir_hat(:, i_dir)))

               DO i_node = 1, n_radial
                  ! Gauss-Legendre on [0, q_max]; the extra q is the polar Jacobian.
                  q = 0.5_dp*q_max*(nodes(i_node) + 1.0_dp)
                  denom = 1.0_dp + twopi*q*a_scalar
                  radial_w = 0.5_dp*q_max*rule(i_node)*q/denom

                  gc%K_head(i_tab, i_w) = gc%K_head(i_tab, i_w) - &
                                          jac*radial_w*4.0_dp*pi**2*a_scalar*COS(q*proj)
                  gc%K_wing(i_tab, :, i_w) = gc%K_wing(i_tab, :, i_w) + &
                                             jac*radial_w*twopi*SIN(q*proj)*dir_hat(:, i_dir)
               END DO
            END DO
         END DO
      END DO

      CALL bs_env%para_env%sum(gc%K_head)
      CALL bs_env%para_env%sum(gc%K_wing)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T60,I6,A,I6)') &
            'Gamma-point correction: cell-integral table', n_tab, ' x', n_tab
         WRITE (bs_env%unit_nr, '(T4,A,T72,ES10.3)') 'K_head at R = 0, lowest frequency', &
            gc%K_head(1 + (n_tab/2) + n_tab*(n_tab/2), 1)
         WRITE (bs_env%unit_nr, '(T4,A,T70,F12.1)') 'Wall time [s]', m_walltime() - start_time
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      DEALLOCATE (dir_hat, dir_weight, dir_radius, nodes, rule)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_tabulate

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
!> \param bs_env ...
!> \param n_face Gauss-Legendre points per edge
!> \param dir_hat unit directions, 3 x n_dir
!> \param dir_weight angular weights
!> \param dir_radius distance to the cell boundary along each direction
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
      DO i_axis = 1, 2
         j_axis = 3 - i_axis
         dq = 0.5_dp*b(:, j_axis)
         DO i_sign = -1, 1, 2
            DO i_node = 1, n_face
               i_dir = i_dir + 1
               q = REAL(i_sign, dp)*0.5_dp*b(:, i_axis) + nodes(i_node)*dq
               dir_radius(i_dir) = SQRT(SUM(q*q))
               dir_hat(:, i_dir) = q/dir_radius(i_dir)
               cross = [q(2)*dq(3) - q(3)*dq(2), q(3)*dq(1) - q(1)*dq(3), &
                        q(1)*dq(2) - q(2)*dq(1)]
               dir_weight(i_dir) = rule(i_node)*SQRT(SUM(cross**2))/SUM(q*q)
            END DO
         END DO
      END DO

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
            p1 = 1.0_dp
            p2 = 0.0_dp
            DO j = 1, n
               p0 = p2
               p2 = p1
               p1 = (REAL(2*j - 1, dp)*x*p2 - REAL(j - 1, dp)*p0)/REAL(j, dp)
            END DO
            dp_n = REAL(n, dp)*(x*p1 - p2)/(x*x - 1.0_dp)
            x_old = x
            x = x - p1/dp_n
            IF (ABS(x - x_old) <= 1.0E-14_dp) EXIT
         END DO
         nodes(i) = x
         weights(i) = 2.0_dp/((1.0_dp - x*x)*dp_n*dp_n)
      END DO

   END SUBROUTINE gamma_corr_gauss_legendre

! **************************************************************************************************
!> \brief ΔW^c_ℓℓ'(iτ_t) for one pair of grid points, by bilinear interpolation of the cell-integral
!>        tables at the minimum-image in-plane displacement R = MIC(r_ℓ − r_ℓ').
!>
!>          ΔW^c_ℓℓ' = Z_ℓ0 Z_ℓ'0 K_h(R)
!>                   + Σ_α [T^w_ℓα Z_ℓ'0 − Z_ℓ0 T^w_ℓ'α] K_w,α(R)
!>
!>        K_h is even in R and the bracket is antisymmetric under ℓ ↔ ℓ', while K_w,α is odd, so
!>        ΔW is symmetric, as a screened interaction must be. The tables are indexed in fractional
!>        coordinates of the in-plane lattice vectors over [−1/2, 1/2], which is exactly the
!>        minimum-image range, so no extrapolation can occur.
!>
!>        The tables hold the frequency-domain kernels; i_t selects an imaginary-TIME point and
!>        the cosine transform back to time is applied here, which keeps the stored tables
!>        frequency-indexed (as the response naturally produces them) while the caller, the Σ^c
!>        panel loop, works in time.
!> \param bs_env ...
!> \param i_t imaginary-time index
!> \param l_grid first grid index
!> \param lp_grid second grid index
!> \return ΔW^c_ℓℓ'(iτ_t)
! **************************************************************************************************
   FUNCTION gamma_corr_delta_W(bs_env, i_t, l_grid, lp_grid) RESULT(delta_W)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: i_t, l_grid, lp_grid
      REAL(KIND=dp)                                      :: delta_W

      INTEGER                                            :: i_dir, i_u, i_v, i_w, n_tab
      REAL(KIND=dp)                                      :: f_u, f_v, head, weight, wing(3)
      REAL(KIND=dp), DIMENSION(2)                        :: frac
      REAL(KIND=dp), DIMENSION(3)                        :: R_cart
      TYPE(gamma_corr_type), POINTER                     :: gc

      gc => bs_env%ri_rs%gamma_corr
      n_tab = gc%n_tab(1)

      ! Minimum-image displacement R = r_l - r_l', then its fractional in-plane coordinates in
      ! [-1/2, 1/2]. NOTE the argument order: pbc3(ra, rb) returns rb - ra, so the grid points go
      ! in reversed. K_h is even in R and would not notice, but K_w is ODD, so getting this
      ! backwards would flip the sign of the whole wing contribution.
      R_cart = pbc(bs_env%ri_rs%grid_points(:, lp_grid), &
                   bs_env%ri_rs%grid_points(:, l_grid), bs_env%ri_rs%cell)
      frac(1) = DOT_PRODUCT(bs_env%ri_rs%cell%h_inv(1, :), R_cart)
      frac(2) = DOT_PRODUCT(bs_env%ri_rs%cell%h_inv(2, :), R_cart)
      frac(:) = frac - ANINT(frac)

      ! Bilinear weights on the table mesh.
      f_u = (frac(1) + 0.5_dp)*REAL(n_tab - 1, dp)
      f_v = (frac(2) + 0.5_dp)*REAL(n_tab - 1, dp)
      i_u = MIN(MAX(INT(f_u), 0), n_tab - 2)
      i_v = MIN(MAX(INT(f_v), 0), n_tab - 2)
      f_u = f_u - REAL(i_u, dp)
      f_v = f_v - REAL(i_v, dp)

      head = 0.0_dp
      wing(:) = 0.0_dp
      DO i_w = 1, bs_env%num_time_freq_points
         weight = bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, i_w)
         head = head + weight*gamma_corr_bilinear(gc%K_head(:, i_w), n_tab, i_u, i_v, f_u, f_v)
         DO i_dir = 1, 3
            wing(i_dir) = wing(i_dir) + &
                          weight*gamma_corr_bilinear(gc%K_wing(:, i_dir, i_w), n_tab, &
                                                     i_u, i_v, f_u, f_v)
         END DO
      END DO

      delta_W = gc%Z_l0(l_grid)*gc%Z_l0(lp_grid)*head
      DO i_dir = 1, 3
         delta_W = delta_W + wing(i_dir)* &
                   (gc%T_wing(l_grid, i_dir, i_t)*gc%Z_l0(lp_grid) - &
                    gc%Z_l0(l_grid)*gc%T_wing(lp_grid, i_dir, i_t))
      END DO

   END FUNCTION gamma_corr_delta_W

! **************************************************************************************************
!> \brief Bilinear interpolation on the flattened n_tab x n_tab table.
!> \param table flattened table, index 1 + i_u + n_tab*i_v
!> \param n_tab points per direction
!> \param i_u lower u index
!> \param i_v lower v index
!> \param f_u fractional offset in u
!> \param f_v fractional offset in v
!> \return interpolated value
! **************************************************************************************************
   PURE FUNCTION gamma_corr_bilinear(table, n_tab, i_u, i_v, f_u, f_v) RESULT(value)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: table
      INTEGER, INTENT(IN)                                :: n_tab, i_u, i_v
      REAL(KIND=dp), INTENT(IN)                          :: f_u, f_v
      REAL(KIND=dp)                                      :: value

      value = (1.0_dp - f_u)*(1.0_dp - f_v)*table(1 + i_u + n_tab*i_v) + &
              f_u*(1.0_dp - f_v)*table(2 + i_u + n_tab*i_v) + &
              (1.0_dp - f_u)*f_v*table(1 + i_u + n_tab*(i_v + 1)) + &
              f_u*f_v*table(2 + i_u + n_tab*(i_v + 1))

   END FUNCTION gamma_corr_bilinear

! **************************************************************************************************
!> \brief Project the charge channel out of V(k) in the REGULAR W path, V(k) <- Π V(k) Π.
!>
!>        This is the other half of separating the charge channel, and without it the correction
!>        would double count: ΔW adds the analytically integrated head and wings back, so they
!>        must not also be present in the V(k) that ε and W are built from. At q = 0 that term is
!>        not merely double counted, it is the divergent one -- the whole reason the correction
!>        exists -- so the projection is what makes a Γ-containing mesh well posed at all.
!>
!>        The projector is rank one, Π = Id − c cᵀ with c normalized, so the triple product is a
!>        rank-two update rather than two matrix multiplications. With u = V c,
!>
!>          (Π V Π)_PQ = V_PQ − c_P u_Q − u_P c_Q + (cᵀu) c_P c_Q .
!>
!>        V is Hermitian, so u is formed from the distributed matrix by one reduction and the
!>        update is then purely local.
!> \param bs_env ...
!> \param cfm_V V(k), overwritten by Π V(k) Π
! **************************************************************************************************
   SUBROUTINE gamma_corr_project_coulomb(bs_env, cfm_V)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: cfm_V

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_project_coulomb'

      COMPLEX(KIND=dp)                                   :: c_dot_u
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: u
      INTEGER                                            :: col, handle, i_col, i_row, n_RI, row
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(gamma_corr_type), POINTER                     :: gc

      CALL timeset(routineN, handle)

      gc => bs_env%ri_rs%gamma_corr
      CPASSERT(ALLOCATED(gc%charge))
      n_RI = bs_env%n_RI

      CALL cp_cfm_get_info(cfm_V, row_indices=row_indices, col_indices=col_indices)

      ! u_P = sum_Q V_PQ c_Q, accumulated from the local block and completed by one reduction.
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

      c_dot_u = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO row = 1, n_RI
         c_dot_u = c_dot_u + gc%charge(row)*u(row)
      END DO

      ! cᵀV(k)c is the charge-channel weight this projection removes. At q = 0 it is the piece
      ! that diverges with SIZE_LATTICE_SUM, so reporting it once is the direct check that the
      ! separation happened and how large the offending term was.
      IF (bs_env%unit_nr > 0 .AND. .NOT. gc%reported_projection) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T68,ES14.6)') &
            'Gamma-point correction: charge weight c^T V c removed from V(k)', REAL(c_dot_u, dp)
         gc%reported_projection = .TRUE.
      END IF

      ! Rank-two update; purely local once u and c are replicated.
      DO i_col = 1, SIZE(col_indices)
         col = col_indices(i_col)
         DO i_row = 1, SIZE(row_indices)
            row = row_indices(i_row)
            cfm_V%local_data(i_row, i_col) = cfm_V%local_data(i_row, i_col) &
                                             - gc%charge(row)*u(col) &
                                             - u(row)*gc%charge(col) &
                                             + c_dot_u*gc%charge(row)*gc%charge(col)
         END DO
      END DO

      DEALLOCATE (u)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_project_coulomb

END MODULE gw_gamma_correction
