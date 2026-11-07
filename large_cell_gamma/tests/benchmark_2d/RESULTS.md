# 2D benchmark: Gamma-only RI-RS GW (KPOINTS_W 1 1 1 + head-and-wing correction, 9x9 supercell)
# vs the full-k RI-RS GW of PR 6447 (head 4e0d03c43b, standard 16-parameter Padé), primitive cell

Settings identical in both codes: DZVP-MOLOPT-SR-GTH, RI sets at ~2e-3 (B RI_20_4, N RI_26,
Si RI_26, C RI_26, Ga RI_61, Al RI_20_4), GTH-PBE, CUTOFF 500 / REL 50, MT Poisson, 15 A
(hBN, SiC, GaN, AlN) vacuum, 10 minimax points, EPS_FILTER 1e-10, CUTOFF_RADIUS_RI 5 A,
SIZE_LATTICE_SUM 4. Ours: 9x9 supercell, DOS mesh 1x1 (the folded Gamma manifold contains K);
tabulated RI-RS grids (data/ri_rs_grid). PR: SCF mesh = W mesh = DOS mesh n x n x 1, path G-M-K,
RI_RS grid optimisation RS_AO_RATIO 7 (his default, required by his code). Indirect G0W0 gap
over the sampled points (K for all four materials). Runs: Noctua largemem,
/scratch/hpc-prf-metdyn/metdyn07_Ritaj/benchmark_2d (inputs mirrored here).

| system | ours 9x9 (162 at.) | PR 12x12 | PR 18x18 | PR 24x24 | ours - PR24 |
|---|---:|---:|---:|---:|---:|
| h-BN  (a = 2.502 A) | 7.435 | 7.434 | 7.444 | 7.451 | -0.016 |
| SiC   (a = 3.090 A) | 3.982 (12x12: 3.980) | 4.021 | 4.029 | 4.031 | -0.049 (-0.051) |
| h-GaN (a = 3.210 A) | 4.345 (12x12: 4.363) | 4.402 | 4.416 | 4.420 | -0.075 (-0.057) |
| h-AlN (a = 3.110 A) | 6.082 (12x12: 6.090) | 6.319 | 6.331 | 6.339 | -0.257 (-0.249) |

Notes
- The PR with its averaged-Padé analytic continuation (commit a56e7a351e, since moved to a
  separate PR) gave hBN 7.828 (12x12) / 7.839 (18x18): the continuation, not W, was the
  0.4 eV discrepancy seen first. State by state at K (ours vs PR 12x12, standard Padé to be
  re-extracted): Sigma^x and v_xc agree to 0.07 eV.
- PR odd meshes (9x9) abort in his gamma_point_correction ("Inconsistent Gamma-neighbor
  transfer points"); even meshes only.
- MoS2 dropped: no tabulated RI-RS grid for Mo and our GRID_OPTIMIZATION is nonperiodic-only.
- Grid check: hbn PR 12x12 with RS_AO_RATIO 20 gives 7.432 (vs 7.434 at ratio 7): his grid is
  converged to 2 meV.
- AlN (indirect K -> Gamma gap, nearly-free-electron CBM at Gamma): state by state ours vs PR
  24x24: VBM(K) -6.506 vs -6.601, CBM(Gamma) -0.424 vs -0.262 (Sigma^c -2.493 vs -2.391,
  Sigma^x -4.894 vs -4.834 with cutoffs 12.1 vs 32.3 A). NOT a cell-size effect on our side:
  the 12x12 cell (288 atoms, Sigma^x cutoff 16.2 A) gives 6.090: VBM(K) Sigma^c 2.622 /
  Sigma^x -19.060 (PR: 2.589 / -19.131), CBM(Gamma) Sigma^c -2.482 / Sigma^x -4.888 (PR: -2.391 /
  -4.834). Sigma^x scatter between the codes is +-0.07 eV in every system (hBN: ours more
  negative, AlN: less), i.e. grid/RI-level differences, not a trend. The AlN-specific part is
  Sigma^c of the nearly-free-electron CBM at Gamma (0.09 eV): origin not resolved; candidate:
  the strictly in-plane charge model (phi_0 = 1/sqrt(A), no exp(-q|z|)) for a state that
  extends into the vacuum, which enters our whole-small-BZ integral but only the near-Gamma
  cells of the PR.

## AlN follow-up (2026-11-07)

- PR grid is converged for AlN too: 12x12 with RS_AO_RATIO 20 gives 6.318 eV (direct 7.200)
  vs 6.319 (7.198) at ratio 7; CBM(Gamma) Sigma^c -2.376 / Sigma^x -4.837 at ratio 20. So the
  0.25 eV is not on the PR grid side (dir aln/pr_kp12_ratio20, 33 s on one largemem node).
- Decomposition of ours - PR24 (gap -0.257): Sigma^x VBM(K) +0.071 and CBM(Gamma) -0.054 (both
  close the gap, 0.125 together), Sigma^c CBM -0.09, Sigma^c VBM +0.033. The e^{-q|z|}
  thickness form factor would scale with the small-BZ radius and should shrink the Sigma^c part
  by ~half from 9x9 to 12x12; observed shift 8 meV. It is therefore unlikely to be the cause.
- Arbiter on our side: tensor (exact-RI) driver with the same Gamma-only W + 2D correction
  (commit 9263cd4925), AlN 6x6 (aln/tensor_6x6) vs RI-RS 6x6 (aln/gamma_6x6): isolates the
  RI-RS grid fit from the Gamma model. Results below when available.
- RESULT (2026-11-07): AlN 6x6 (72 atoms, DZVP) RI-RS vs tensor exact-RI, both with the Gamma-only
  W + 2D head/wing correction (same binary, 9263cd4925): gap 6.013 vs 6.013 eV; VBM(K) state 144
  eps_DFT -5.231, Sigma^c 2.641, Sigma^x -19.063, v_xc -15.173, GW -6.480 in BOTH; CBM(Gamma)
  state 145 -1.868 / -2.529 / -4.901 / -8.832 / -0.466 in BOTH; all occupied and the first ~40
  virtual states agree to <= 1 meV (band files in results/aln_6x6/). Cost: RI-RS 326 s / 7.7 GB,
  tensor 1133 s / 44 GB per rank (chi 3.1 s vs 33 s per time point). So the RI-RS grid fit is
  NOT the origin of the AlN difference to the PR; it lies in the Gamma-point model (whole-small-BZ
  analytic head+wing with the charge projector vs the PR's near-Gamma cell model) or in the
  supercell-Gamma vs full-k formulation (Sigma^x truncation 8-16 A vs 32 A, folded spectra on the
  minimax grid), not in numerics on either side.

## Rebase validation (branch rirs_periodic_gamma_v2 on upstream master after PR 6447, 2026-11-07)

Same inputs, binary from the rebased branch (Noctua periodic_3d_test/cp2k-v2):
- regtests QS/regtest-gw-realspace: 45/49 correct; the 4 "wrong" were our periodic H2O
  references, shifted by upstream #6474 (averaged Padé orders, unconditional); HF edges identical,
  references updated (23.312 -> 23.519, evGW0 -6.003/17.865/23.868 -> -6.011/17.943/23.953).
- AlN 6x6 RI-RS: gap 6.241 (was 6.013); HF 11.184 (identical), Sigma^x -19.063/-4.901 (identical),
  v_xc identical; Sigma^c 2.641 -> 2.460 (VBM), -2.529 -> -2.482 (CBM): continuation only.
- diamond 3x3x3 RI-RS: 9.210 (was 9.039); HF 14.833 (was 14.838); correction diagnostics identical.
NOTE: every GW number on master now carries the averaged Padé continuation; the 2D benchmark table
above and the PR comparisons were made with the standard Padé on both sides and must be rerun
with one build on both sides before being compared to master numbers.
