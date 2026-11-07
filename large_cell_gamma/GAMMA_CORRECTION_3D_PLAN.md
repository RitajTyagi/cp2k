# Γ-point head-and-wing correction for 3D periodic cells: plan (2026-11-05)

Status: the 2D (PERIODIC XY) correction is validated (8x8 h-BN 6.710 eV vs full-k 6.690).
3D is NOT supported: gamma_corr_requested returns .FALSE. unless periodic = (1,1,0).
The colleague's full-k code (PR 6447) has both; the Γ-only version below mirrors its 3D
branch (gamma_point_correction, non-slab path; gamma_ri_coulomb_g0; gamma_cell_face_directions).

## What changes from 2D to 3D (physics)

1. Coulomb head: v(q) = 4π/q² instead of 2π/q. The singularity is still integrable over the
   small Brillouin-zone cell of the supercell, because the Jacobian is q² dq.

2. V_PQ(q=0) of the RI basis. In 2D the finite-box lattice sum minus the box constant
   D_box c cᵀ is a sound regularization (validated). In 3D the finite-box sum of 1/r between
   charged functions is conditionally convergent with a SHAPE-dependent surface term, so the
   box construction is not usable. Use the GPW periodic Coulomb integrals with the G = 0
   Fourier coefficient omitted (mp2_eri_2c_integrate_gpw with the periodic Poisson solver,
   exactly the PR's gamma_ri_coulomb_g0). Consequences:
   - no divergent-head subtraction (gamma_corr_subtract_divergent_head is 2D only);
   - the charge channel of V(0) is finite; ε(q=0) in that channel is ≈ 1 (the charge mode is
     unscreened at exactly q = 0), which is right, because the correction supplies the whole
     Brillouin-zone average of the screened head.
   The regular Γ W path already reuses the cached gc%fm_V_zero, so it inherits this V(0).

3. Dielectric head and wing, per VOLUME instead of per area:
     Y = (2/V) Σ µµ e^{-Δτ},  X = (2/√V) Σ ψψ µ e^{-Δτ},  S = Y − EᵀB⁻¹E,  a(q̂) = q̂ᵀSq̂,
     ε̃_00 = 1 + 4π a(q̂)          (q-independent: the macroscopic dielectric constant),
     W^c_00(q) = (4π/q²)(1/ε̃ − 1) = −(4π/q²)·4πa/(1+4πa),
     W^c_P0(q) = ±i (4π/q) q̂·T_P/(1+4πa),  T = V(0) Π V_b^(-1/2) B⁻¹E   (same block inverse).
   Cell integrals (BZ averages, prefactor V/(2π)³):
     K_h(R) = (V/(2π)³) ∫_cell d³q (4π/q²)(1/ε̃(q̂) − 1) cos(q·R),
     K_w(R) = (V/(2π)³) ∫_cell d³q (4π/q) q̂ sin(q·R)/(1+4πa(q̂)),
   and ΔW_PQ(iω) = c_P c_Q K_h(R_PQ)/V + (1/√V) Σ_α [T_Pα c_Q − c_P T_Qα] K_w,α(R_PQ),
   added to W^c_PQ(iτ) before M⁻¹ W M⁻¹ exactly as in 2D (gamma_corr_add_to_W).

4. Quadrature of the cell around Γ: the parallelepiped |q·a_i| ≤ π with six faces; a
   Gauss–Legendre n x n rule on each face gives the directions q̂ and the radii q_max(q̂)
   (PR: gamma_cell_face_directions), Gauss–Legendre along q to the face. The direction set is
   inversion symmetric, so the cos/sin split stays exact.

5. Kernel evaluation: directly at every local atom pair (P, Q) of the distributed W (no 3D
   tables), with the 3D minimum image of R_P − R_Q. Cost n_atom²/n_ranks x n_q x n_t,
   negligible next to W itself.

6. Unchanged: velocity-gauge dipoles, TRK check, Σ^x with the truncated Coulomb operator,
   the 2D path (kept verbatim, it is validated).

## Implementation steps (files)

- gw_gamma_correction.F: gate (XYZ allowed), gamma_corr_cell_geometry (area/volume, b_1..b_3),
  3D V(0) via GPW (new routine, modelled on the PR's gamma_ri_coulomb_g0), 3D face quadrature
  (generalize gamma_corr_bz_boundary), 3D kernels evaluated per atom pair (replace the
  2D tables by direct evaluation, or keep tables for 2D), 1/A → 1/measure everywhere.
- gw_compute_V_and_W.F: skip the divergent-head subtraction in 3D.
- ri_rs_types.F: measure (A or V) and dimensionality in gamma_corr_type.
- Smoke test first: does the Γ-only large-cell driver run at all with PERIODIC XYZ
  (MIC truncation, chunk centroids, Σ^x cutoff = half the smallest cell dimension)?

## Validation

Bulk system with the PR 3D full-k code as reference, same basis and settings on both sides:
diamond or cubic BN (2 atoms, cubic cell, SZV basis). Our n x n x n supercell Γ-only
(4x4x4 = 128 atoms on noctua, ~like the 8x8 slab) against the PR's N x N x N meshes up to
convergence; the 3x3x3 = 54 atoms is feasible locally if the aug-free basis permits it.
Quantities: gap at the k-points the supercell folds onto Γ, state-by-state Σ^c at Γ.

## To confirm before coding

(a) V(0) via the GPW periodic integrator with G = 0 omitted (needs the GPW machinery of
    qs_env still alive at GW time; the PR re-creates the k-point cell index for it);
(b) no bare-head term is needed in W^c (the bare 4π/q² average belongs to V, which enters Σ
    only through ε at q = 0 and Σ^x; the PR adds it to its V_grid_R for its own W = V + W^c
    bookkeeping);
(c) test system and reference meshes.

Not done (both cells): Gaussian form factors of the auxiliary charges. Numerically, the
head kernel changes by 1.1 % at 8x8 (q_max = 0.096 bohr⁻¹, α_min = 0.199) and 0.4 % at
14x14, on a head contribution of ~30 meV to the gap: < 1 meV. Only on a 3x3 cell (7 %) it
would reach ~20 meV. Not worth the per-exponent kernels.
