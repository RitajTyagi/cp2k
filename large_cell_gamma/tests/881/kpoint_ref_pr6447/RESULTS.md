# Cross-check of the Gamma-only RI-RS GW (branch rirs_periodic_gamma) against the
# full-k RI-RS GW of cp2k PR #6447 (J. Wilhelm, charge-neutralised chi(q), commit a2a9da4c91)

System: monolayer h-BN, a = 2.50223 A, 15 A vacuum, PERIODIC XY, MT solver, GPW, CUTOFF 600,
aug-SZV-MOLOPT-GTH-tier-1 / RI_aug-SZV-MOLOPT-GTH-tier-1_1 (my_gth_basis), GTH-PBE, 10 minimax
points, EPS_FILTER 1e-10, CUTOFF_RADIUS_RI 5, SIZE_LATTICE_SUM 4 -- identical in both codes.

PR build: /localhome/ritaj/project/rirs/cp2k-pr6447 (worktree of pr6447, LIBXS off).
Run: ./run.sh hbn_kpN 2 6 (2 ranks x 6 threads, ~0.7 GB, ~4.5 min each; the GW part is ~6 s,
the rest is RI-RS grid optimisation + d_lP). RS_AO_RATIO 7 (112 grid points); the g20 run
(320 points) changes the gap by 3 meV -> grid-converged.

Comparison points. Our 8x8 supercell with its 2x2 DOS mesh samples the primitive k-points
(i/8 +- 1/32, j/8 +- 1/32). The band edges (gap 6.706) sit at P = (11/32, 11/32) ~ K; the
supercell Gamma point folds Q = (3/8, 3/8), whose states are listed in the supercell band file.

## PR full-k code, G0W0 band edges (eV); Sigma^x cutoff = a_min N_k / 2

| mesh  | Sx cutoff | gap(K) | gap(P) | gap(Q) | GW VB(Q) | GW CB(Q) | HF gap |
|-------|-----------|--------|--------|--------|----------|----------|--------|
| 6x6   | 6.5 A     | 6.539  | 6.590  | 7.051  | -4.936   | 2.115    | 11.765 |
| 8x8   | 8.7 A     | 6.621  | 6.658  | 7.083  | -4.974   | 2.109    | 11.933 |
| 12x12 | 13.0 A    | 6.644  | 6.681  | 7.096  | -5.005   | 2.091    | 12.025 |
| 16x16 | 17.4 A    | 6.648  | 6.684  | 7.096  | -5.018   | 2.078    | 12.061 |
| 20x20 | 21.7 A    | 6.651  | 6.687  | 7.098  | -5.027   | 2.071    | 12.080 |
| 24x24 | 26.1 A    | 6.654  | 6.690  | 7.101  | -5.034  | 2.067  | 12.092 |

## Our Gamma-only supercell code (KPOINTS_W 1 1 1 + head/wing correction)

| cell  | Sx cutoff | gap(P-type DOS point) | gap(Q=Gamma_super) | GW VB(Q) | GW CB(Q) | HF gap |
|-------|-----------|-----------------------|--------------------|----------|----------|--------|
| 8x8   | 8.7 A     | 6.706                 | 7.120              | -5.039   | 2.081    | 12.031 |
| 10x10 | 10.8 A    | 6.705                 |                    |          |          |        |
| 14x14 | 15.2 A    | 6.699                 |                    |          |          | 12.072 |
(old k-mesh scheme, same cells: 8.152 (8x8), 6.582 (14x14, 4x4 mesh), 6.595 (14x14, 2x2 mesh))

## State-by-state at Q (eV): our 8x8 Gamma-only vs PR 8x8 (same Sigma^x cutoff 8.7 A)

|            | DFT    | Sigma^c | Sigma^x | v_xc    | G0W0   |
|------------|--------|---------|---------|---------|--------|
| VB ours    | -5.952 |  2.131  | -17.665 | -16.447 | -5.039 |
| VB PR 8x8  | -5.954 |  2.160  | -17.628 | -16.447 | -4.974 |
| VB PR 20x20| -5.953 |  2.221  | -17.742 | -16.447 | -5.027 |
| VB PR 24x24| -5.953 |  2.225  | -17.753 | -16.447 | -5.034 |
| CB ours    | -0.906 | -3.210  |  -4.744 | -10.941 |  2.081 |
| CB PR 8x8  | -0.907 | -3.164  |  -4.761 | -10.941 |  2.109 |
| CB PR 20x20| -0.907 | -3.218  |  -4.745 | -10.941 |  2.071 |
| CB PR 24x24| -0.907 | -3.222  |  -4.744 | -10.941 |  2.067 |

## Conclusions
1. The PR code converges with its mesh from below onto the Gamma-only numbers: gap(P)
   6.590 -> 6.690 (24x24) vs 6.699-6.706; absolute G0W0 energies at Q agree to 5 meV (VB) and
   14 meV (CB) at 24x24. The Gamma-only implementation is confirmed at the ~15 meV level.
2. The old k-mesh large-cell scheme (6.58 at 14x14, 8.15 at 8x8) is NOT a reference: it samples
   W(q) at q /= 0 with the spurious pair-density charge that PR #6447 identified, and its
   MIC-truncated V makes that a finite, cell-size-dependent uniform Sigma^c shift.
3. Sigma^x and Sigma^c individually still move with the PR mesh (Sigma^x cutoff a_min N_k/2),
   but their sum converges; only totals should be compared between the codes.

## Update 2026-11-05: Gamma-only 8x8 with the RI-level correction (commit abacf6ad50)

gap(P) 6.710, VB(Q) -5.036/-5.045, CB(Q) 2.070/2.090, HF gap 12.031. The 6.706 row above was
the grid-level correction with the wing bugs S1-S3 and the S4 defect (audit document).

## One-build rerun (2026-11-08): master after the merge of PR 6447, averaged Padé on both sides

Same inputs (stage copies + my_gth_basis in tests/881/one_build_v2/), binary periodic_3d_test/cp2k-v2
on Noctua largemem, dirs benchmark_2d/hbn_augszv_v2/{gamma_8x8,gamma_10x10,gamma_14x14,pr_kp6..24}.
HF gaps reproduce the tables above to the meV (12.031/12.047/12.072 ours; 11.766 ... 12.091 PR),
so all changes are the averaged Padé continuation (#6474).

PR, gap(K) / gap(P) / gap(Q), G0W0 VB(Q) / CB(Q), Sigma^c VB(Q) / CB(Q):

| mesh  | gap(K) | gap(P) | gap(Q) | VB(Q)  | CB(Q) | Sc VB / CB      |
|-------|--------|--------|--------|--------|-------|-----------------|
| 6x6   | 6.662  | 6.714  | 7.182  | -5.013 | 2.169 | 2.050 / -3.088  |
| 8x8   | 6.741  | 6.780  | 7.212  | -5.051 | 2.161 | 2.082 / -3.113  |
| 12x12 | 6.762  | 6.801  | 7.225  | -5.082 | 2.143 | 2.118 / -3.146  |
| 16x16 | 6.767  | 6.804  | 7.225  | -5.096 | 2.129 | 2.134 / -3.161  |
| 20x20 | 6.770  | 6.807  | 7.227  | -5.105 | 2.122 | 2.142 / -3.168  |
| 24x24 | 6.773  | 6.810  | 7.231  | -5.112 | 2.119 | 2.146 / -3.172  |

Gamma-only (printed gap = minimum over the DOS points nearest K; that point is P only for 8x8,
13/40 or 14/40 for 10x10 and 19/56 for 14x14, so the cell series is NOT at one fixed point):

| cell  | gap (DOS point) | gap(Q) | VB(Q)  | CB(Q) | Sc VB / CB     | Sx VB / CB        | wall 8 ranks |
|-------|-----------------|--------|--------|-------|----------------|-------------------|--------------|
| 8x8   | 6.866           | 7.271  | -5.128 | 2.143 | 2.043 / -3.148 | -17.665 / -4.744  | 864 s        |
| 10x10 | 6.855           |        |        |       |                |                   | 1767 s       |
| 14x14 | 6.832           |        |        |       |                |                   | 4821 s       |

- 8x8 vs PR 24x24 at Q: VB -0.016, CB +0.024, gap(Q) +0.040; at P: 6.866 vs 6.810 (+0.056).
  With the standard Padé the same comparison gave 5 / 14 / +0.019 and +0.020 (6.710 vs 6.690).
  Sigma^x and v_xc identical to the meV in both continuations; the averaged Padé is more
  sensitive to the small Sigma^c differences (ours 2.043 vs PR-8x8 2.082 at the same 8.7 A
  cutoff) because the continuation of the folded supercell spectrum and of the primitive-cell
  spectrum differ in their Padé orders.
- The full-k series still converges from below towards the Gamma-only values (gap(P) 6.714 ->
  6.810 vs 6.832-6.866), i.e. the agreement with one build is 20-60 meV instead of 10-20 meV.
