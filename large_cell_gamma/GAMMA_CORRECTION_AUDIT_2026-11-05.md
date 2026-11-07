# Γ-point head-and-wing correction: cleanup, audit and fixes (2026-11-05)

Branch `rirs_periodic_gamma`. Commits (newest last):
- f8090ceda7  remove debug-era diagnostics (misplaced-block count, operator symmetry self-test,
              K_head(R=0) print) from the charge-channel fit and the tables
- 23be26a03a  CP2K formatter + doxify pass on all touched files, unused imports dropped
- 35e6178c4f  BUG FIXES S1 + S2 and bit-identical performance changes (see below)
- d8556d7662  developer comments for every routine and step (comment-only, code verified identical)
- 7b644c7511  BUG FIX S3: dipole sign (3x3 test: 6.229 -> 7.040 eV; the 3x3 cell is extremely
              wing-sensitive because its small Brillouin zone is huge, q_max ~ |b|/6)

## Bugs found by the audit (all in the wing part of the correction)

| id | what | effect | fix |
|----|------|--------|-----|
| S1 | wing X_ℓα accumulated from the mat_psi = Φ C blocks each rank owns; on a 2D DBCSR process grid (npcol > 1) only occupied/virtual pairs on the same process column were formed | 3x3 SZV h-BN: 1 rank 5.945 eV, 4 ranks 6.560 eV. ALL 16-rank noctua Γ-only runs (4x4 grid) were affected | X = diag(Φ P^α(iτ) Φᵀ), P^α = C_occ e µ^α e C_virᵀ via two ScaLAPACK gemms + one DBCSR product on Φ's pattern (same kernel as the fitting operator). 1 rank identical to before; 4 ranks = 1 rank |
| S2 | wing kernel K_w(R,iω) was cosine-transformed to iτ and then multiplied by T^w at the frequency index i_t | wrong wing in every run; 3x3: gap 5.945 -> 6.229 eV when fixed | tables carry one channel per (ω, α): K_time(1+3n_t, n_R, n_t), U_grid(1+3n_t, n_grid); the ω sum is done per pair in the kernel |
| S3 | transition dipole stored as mu_mo(i, a)/Δ = ⟨i|[H,r]|a⟩/(ε_a−ε_i) = −⟨a|r|i⟩; Eq. 48 and the tex need +⟨a|r|i⟩; the PR 6447 code also uses the ⟨a|…|i⟩ orientation (gw_gamma_transition_moments) | wing sign flipped (head and TRK sum rule unaffected) | mu_mo(a, i)/Δ |

Deterministic test: `large_cell_gamma/tests/local_szv/gamma_3x3_det.inp`, 1 rank x 1 thread,
EPS_FILTER 1e-12, ~3.5 min; serial runs are bit-identical run to run. Band files of every stage
are kept next to it (`bandstructure_3x3_det_{run1,S1,S2,perf,S3}_1rank`).

FINAL VALIDATION STILL REQUIRED: 8x8 aug-SZV h-BN on noctua with the fixed code,
`large_cell_gamma/tests/881/new_rirs_fixed/hbn.inp` (16 ranks as before, ~4 min), to be compared
with the PR 6447 full-k reference 6.69 eV at P=(11/32,11/32) (`881/kpoint_ref_pr6447/RESULTS.md`).
The old 16-rank number 6.706 carried S1+S2+S3. The branch with all fixes is pushed to
origin (RitajTyagi/cp2k, rirs_periodic_gamma, HEAD 7b644c7511).

## Performance (392-atom h-BN, 16x24: the correction was 472 s of 1072 s = 44 %)

Done, bit-identical (verified on the 3x3 test):
- V(q=0) lattice sum reused by the regular W path instead of being repeated (70 s of 472)
- one dsyev for V_b^(+1/2) and V_b^(-1/2) (26 s)
- block kernel: OpenMP hoisted from the 32x32 block to the DBCSR iterator of the panel (the
  fork/join per block dominated the 47 s); ANINT replaced by INT(x+SIGN(0.5,x))
- gamma_corr_tabulate threaded; K_head/K_wing/T_wing/Z_l0/S_head freed after packing,
  F_density/F_extern freed after the response (2.3 GB/rank dead weight at N = 1000 removed)
- V_imag replication replaced by cp_fm_maxabsval on the distributed matrix

Not done (proposals, ordered by value; none bit-identical except the restart):
1. Replace the replicated n_RI x n_RI dense algebra (build_coulomb_maps, response, body_solve,
   chi_PQ_at_freq, cp_fm_get_submatrix x 13) by cp_fm/ScaLAPACK: 46 GB per rank and ~1700 s of
   redundant per-rank work at N = 2000; the first OOM at 4 GB/rank is V_zero at N ~ 1000.
2. Write Z_l0 (and S_head, T_wing) to the Z_lP restart file: the 131 s CG fit reruns on every
   restart today.
3. CG fit: Tikhonov 1E-3 (within the band the code comment validated) cuts ~750 -> ~250
   iterations; block-Jacobi preconditioner would be the real lever.
4. Dipoles: build only the occ-vir block (C_virᵀ r C_occ) instead of the full n_mo x n_mo product.
5. Wing kernel: SVD-compress K_w(R, ω) over ω (rank ~4) to bring the per-pair cost from
   1+3n_t channels back to ~1+3r.
6. compute_MinvVsqrt_Vsqrt: V(0) − D_box ccᵀ is indefinite along c, so the Cholesky always fails
   and the complex eigensolver runs; a square-root-free W = V(1−KV)⁻¹ − V would avoid it.

## Addendum (same day, after the noctua rerun): S4, the correction moves to RI space

The 8x8 rerun with S1+S2+S3 (`881/new_rirs_fixed/output.out`) gave 5.451 eV instead of ~6.69,
split the six symmetry-equivalent folded states at the supercell Γ (Σ^c −3.22 .. −4.39 eV where
the old run had −3.21 for all six) and moved deep virtuals by up to 8.6 eV. Diagnostics with
the PR 6447 code (its slab wing zeroed by a local hack): the wing changes the 8x8 primitive gap
by 1 meV (6.621 -> 6.620) and its head by 30 meV; ours changed the 3x3 gap by ±0.4 eV
(wing −1/0/+1: 6.229/6.648/7.040) and the 8x8 by > 1 eV. The formulas of the two codes agree
term by term (same F maps, same E, B, same cell integrals, same dipole orientation).

Root cause (`local_szv/dump`, Python analysis of Z_l0, T_wing, ψ on the grid): the RI-RS
interpolation functions far from the atoms (vacuum points, density 1e-9) carry fit weights
|Z_l0| up to 2e4 and |T_wing| up to 8e7 of alternating sign. They cancel exactly in linear
functionals of a density (Σ_l ρ(l) Z_l0 = 1/√A = 0.0757 for a diagonal pair, 3e-6 for
HOMO-LUMO; Σ ρ T = 3e-3 against partial sums of 6-8) but NOT when each point is weighted by a
different kernel value K(R_ll'), which any pair kernel on the grid does. The odd wing is
cancellation-dominated and was garbage; the even head survived only where the kernel is flat
over neighbouring points (large cells).

Fix (`gamma_corr_add_to_W`, called in compute_W_MIC before M⁻¹ W M⁻¹): the same formula on the
auxiliary functions, (ϕ_P|ΔW|ϕ_Q) = c_P c_Q K_h(R_PQ)/A + (1/√A) Σ_α [T_Pα c_Q − c_P T_Qα]
K_w,α(R_PQ), with c_P the auxiliary charges, T = V(0) Π V_b^(-1/2) B⁻¹E (RI representation),
R_PQ the atom-pair minimum image, tables interpolated at atom pairs. It reaches the grid
through the same fit Z M⁻¹ as the regular W. Removed: the global CG fit of Z_l0 (131 s), the
grid wing map, the packed tables, the per-pair kernel and its panel hook (the 44 % of the
392-atom run is gone; the correction now costs 0.1 s on the 3x3).

3x3 deterministic test, RI-level: wing −1/0/+1 = 6.618/6.625/6.632 eV (the wing is a 7 meV
effect, as in the PR code); folded manifolds degenerate to 0.02 eV; 4 ranks = 1 rank to 1 meV.
Grid-level head-only was 6.648 (23 meV off the RI-level head).

The wing SIGN cannot be validated numerically (7 meV on 3x3, ~1 meV on 8x8); the derivation
(charge-channel ket e^{iqr}/√A, χ_P0 = +i q·X, block inverse, pair phase e^{+iq(R_P − R_Q)})
gives the sign the code has. Gaussian form factors e^{−q²/4α} of the auxiliary charges are
neglected (point-charge limit); a few % of the correction on the 3x3, less on larger cells.

FINAL VALIDATION STILL REQUIRED (unchanged): the 8x8 noctua run with the pushed branch,
`881/new_rirs_fixed/hbn.inp`, against the PR reference 6.69 eV at P.

## VALIDATED (2026-11-05, 16:10): 8x8 aug-SZV h-BN on noctua with abacf6ad50

`881/new_rirs_fixed/output.out`, 16 ranks, 176 s: G0W0 gap 6.710 eV at P (PR 6447 converged
reference 6.690; old grid-level run 6.706, which contained S1+S2+S3 and the S4 defect that
cancelled on this cell size). At Q (supercell Γ): VB −5.036/−5.045 eV (PR 24x24: −5.034),
CB 2.070/2.090 (PR 24x24: 2.067); the six folded states have Σ^c −3.220/−3.202/−3.209 in
exactly degenerate pairs (the S1-S3 run had −3.22 .. −4.39). Changes vs the old run are
< 50 meV for 923 of 1024 states; the rest are virtuals 39-45 eV above the VBM where the
grid-level head was unreliable. The correction costs 0.0 s of the 176 s.

## Dense n_RI x n_RI algebra: replicated or distributed (2026-11-05, evening)

The last OOM hotspot of the correction was the replicated dense algebra of
gamma_corr_build_coulomb_maps / gamma_corr_response (up to ten n_RI² arrays per rank; with
12 RI functions per atom that is 4 GB per rank at ~600 atoms, 46 GB at 2000). Same philosophy
as the Z_lP solves: gamma_corr_init estimates the replicated peak (10 x 8 n_RI² bytes), measures
the memory per rank (mp_mem_avail_per_rank_GB) and takes the cp_fm/ScaLAPACK path when the
estimate exceeds a quarter of it (automatic only; the forcing keyword was removed on the
user's request 2026-11-06, the distributed path was validated with it before). Distributed path: Π filled locally, V_b = Π V Π by pdgemm,
cp_fm_syevd for V_b^(±1/2) with the same relative cut, F_density and G_extern by pdgemm, per
frequency B = 1 + F_dᵀ χ F_d, Cholesky inverse of B, x = B⁻¹E, S −= Eᵀx, T = G_extern x; only
n_RI x 3 arrays are replicated. F_extern (unused since the RI-level correction) dropped.
3x3 det test: AUTO (replicated) bit-identical; DISTRIBUTED 1 rank identical to all printed
digits; DISTRIBUTED 4 ranks 6.633 (= the replicated 4-rank number, panel-loop noise).

Gaussian form factors of the auxiliary charges: NOT implemented. Head-kernel change with the
most diffuse s-exponent of the aug-SZV RI basis (α = 0.199): 7 % at 3x3, 1.8 % at 6x6,
1.1 % at 8x8, 0.4 % at 14x14, on a head contribution of ~30 meV to the gap => < 1 meV on
production cells. (An aug-SZV 3x3/4x4/6x6 test cell could not even run: the minimum-image
k-point overlap is not positive definite with that diffuse basis below 8x8.)

Remaining replicated arrays: µ (3 x n_occ x n_vir, 0.8 GB at 2000 atoms), X(iτ) on the grid
(n_grid x 3 x n_t, 0.16 GB at 2000 atoms). Both grow linearly in system size.
