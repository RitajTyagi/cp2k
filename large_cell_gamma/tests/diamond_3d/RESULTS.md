# Diamond (primitive 2-atom fcc cell, a = 3.567 A), PERIODIC XYZ, GPW, SZV-MOLOPT-SR-GTH,
# RI_20_2_2_1_1_0_0_0_1.6e-02, GTH-PBE-q4, CUTOFF 400, 10 minimax points, EPS_FILTER 1e-12,
# SIZE_LATTICE_SUM 2, default RI cutoff (half the cell height). 2026-11-06.
#
# The sections up to "Exact-RI tensor driver vs RI-RS, 3x3x3" were made with the standard Padé
# continuation (branch rirs_periodic_gamma) and PR 6447 builds that have since changed; the
# definitive one-build comparison (master after the merge of PR 6447, averaged Padé on both
# sides) is the LAST section of this file (2026-11-08).

Runs on Noctua: /scratch/hpc-prf-metdyn/metdyn07_Ritaj/periodic_3d_test (cn_* = compute nodes).

## Lessons
- 2x2x2 primitive supercell (16 atoms) is below the method: Sigma^x = -395 keV, DOS-mesh
  Cholesky fails. CUTOFF_RADIUS_RI must stay below half the cell height (the 2D input's 5 A gave
  Sigma^x garbage on the 3x3x3 cell); the default (3.0 A here) is fine.
- Colleague's PR 6447: 2x2x2 mesh unusable (Sigma^x cutoff 2.1 A, positive Sigma^x); 3x3x3
  ABORTS in his gamma_point_correction ("Inconsistent Gamma-neighbor transfer points", odd 3D
  mesh); even meshes run.

## Gamma-only supercell, 54 atoms (3x3x3), G0W0 gap over the folded Gamma states (eV)

| scheme | VBM | CBM | gap | note |
|---|---:|---:|---:|---|
| no correction (E0) | 10.689 | 22.939 | 12.249 | Gamma-only W with the Ewald V(0) |
| RI-space DeltaW, wing only (E2) | 10.692 | 21.464 | 10.772 | |
| RI-space DeltaW, head only (E1) | 86.198 | -9.919 | -96.117 | FAILS: fit amplification, triplet split 0.4 eV |
| RI-space DeltaW, full | 28.059 | -26.044 | -54.103 | FAILS |
| **head in the AO basis (ed33bd6ced)** | 12.271 | 21.141 | **8.871** | VBM +1.58, CBM -1.80 = the uniform-head shifts |

## Gamma-only supercell, head in the AO basis, cell series

| cell | atoms | HF gap | G0W0 gap | note |
|---|---:|---:|---:|---|
| 3x3x3 | 54 | 14.838 | 8.871 | Sigma^x cutoff 3.1 A |
| 4x4x4 | 128 | 15.123 | 9.722 | Sigma^x cutoff 4.1 A; PR 4x4x4: HF 14.113, G0W0 9.314 |
| 5x5x5 | 250 | | (running, 2 nodes) | |
| 6x6x6 | 432 | | (running, 4 nodes) | |

State by state at 4x4x4 (ours vs PR, same mesh density, same nominal cutoff): Sigma^x of the
occupied states is 0.4-0.65 eV more negative in ours (VBM -18.957 vs -18.302), virtuals agree
to 0.03 eV; Sigma^c VBM 2.061 vs 1.773, CBM -3.340 vs the PR's split triplet -3.296/-2.980/-3.025
(the PR breaks the Gamma CBM degeneracy by 0.3 eV at 4x4x4; ours keeps it). Net: VBM 0.34 eV
lower in ours, CBM equal -> gap +0.37. The exchange difference is a truncation-scheme difference
(spherical cutoff vs his), not W; both must converge with the cell/mesh.

Head potential K_h(0)/Omega = -3.498 eV at omega_1 (offline check -3.39, eps_inf = 3.05);
in imaginary time -1.429 eV at tau_1. HF gap 14.838 (Sigma^x cutoff 3.1 A).

## Colleague's full-k PR 6447 (primitive cell, MP mesh = W mesh, path Gamma/L/X)

| mesh | Sigma^x cutoff | HF gap | G0W0 gap |
|---|---|---:|---:|
| 2x2x2 | 2.1 A | 26.2 | 24.4 (unusable) |
| 4x4x4 | 4.1 A | 14.113 | 9.314 |
| 6x6x6 | | | 9.594 |
| 8x8x8 | | | 9.663 |

## Exchange validated against the exact-RI tensor code (128-atom cell, same cutoff 4.1 A)

| | Sigma^x VBM (Gamma triplet) | Sigma^x CBM | HF gap |
|---|---|---|---|
| tensor large-cell Gamma (exact RI, no RI-RS grid) | -18.957 / -18.957 / -18.957 | -9.767 | 15.125 |
| ours RI-RS Gamma-only | -18.957 / -18.956 / -18.956 | -9.767 / -9.768 | 15.123 |
| PR 6447 4x4x4 (standard Padé build) | -18.535 / -18.581 / -18.372 | -9.739 | 14.275 |

Our exchange is exact to the meV; the PR's exchange at 4x4x4 is 0.4-0.6 eV off and breaks the
Gamma triplet by 0.2 eV (its HF gap reaches 15.17 only at 8x8x8). The 0.4 eV gap difference
at convergence is therefore entirely in Sigma^c (ours -5.3, PR -> -5.75 eV of correlation
contribution to the gap); the standard CP2K tensor k-point GW on the primitive cell
(cn_tensorkp6/8) is the arbiter.

## Correction series (full kernel: head + head-dipole + wing + wing-wing, dipoles clipped at 1/q_max)

| cell | head only (v2) | + head-dipole (v3, wing ~0 by a factor bug) | full (v5/v6) |
|---|---:|---:|---:|
| 3x3x3 | 8.871 | 9.066 | 9.041 (v6 with wing-wing 9.039) |
| 4x4x4 | 9.722 | 7.750 (one state +2 eV: divergent dipole -> clip) | 9.782 (v6 9.782) |
| 5x5x5 | 9.910 | | 9.934 |

## Pending
- ours 4x4x4 (128 atoms, largemem), 5x5x5 (250 atoms, largemem).
- 3D wings / head-dipole terms (O(q.d), ~1/Omega) not yet included in the AO-basis scheme.

## Standard tensor k-point references (cn_tensorkp6, cn_tensorkp8): not completed

Primitive 2-atom cell, SCF mesh = W mesh = 6x6x6 / 8x8x8 (gw_tensor_small_cell_full_kp.F, no
supercell, no Gamma correction; W from the k-mesh extrapolation). Both jobs hit the 24 h limit on
one largemem node (8 ranks x 16 threads) after the FIRST of 10 chi time points: 37870 s (6^3) and
77988 s (8^3) per point, 87 GB peak per rank. The cost is the real-space lattice sums (1607 cells
for the 3c integrals, 12071 cells Delta R from the 9.0 / 9.6 A orbital / RI radii), not the basis
(SZV, 40 RI functions). The k-point code has no chi restart, so the runs cannot be chained;
a complete 6^3 run needs ~5 node-days. Dropped. The arbiter used instead is the tensor
large-cell Gamma driver with the Gamma-only W and the same head-and-wing correction
(exact-RI chi, W body and Sigma^c; only the optical response on the RI-RS grid), see below.

## Exact-RI tensor driver vs RI-RS, 3x3x3 (54 atoms), same binary 9263cd4925, same Gamma-only W + correction

| | RI-RS (cn_rgc_3x3x3) | tensor exact RI (cn_tgc_3x3x3_lm) |
|---|---:|---:|
| G0W0 gap (eV) | 9.039 | 9.386 |
| HF gap (eV) | 14.838 | 15.074 |
| VBM: Sigma^c / Sigma^x (eV) | 2.192 / -19.149 | 2.188 / -19.090 |
| CBM: Sigma^c / Sigma^x (eV) | -3.491 / -10.021 | -3.455 / -9.943 |
| correction diagnostics (K_h(0)/Omega, A_xx, B_xx) | identical | identical |
| wall time, 8 ranks | 1680 s (login node, shared) | 5025 s (largemem) |

State by state (results/exact_ri_3x3x3/): occupied states Sigma^x differs by -0.063 eV on average
(max 0.19), lowest 32 virtuals by -0.144 (max 0.27); Sigma^c at the band edges 0.004 / 0.036 eV.
The difference is the RI-RS grid representation of the truncated Coulomb operator of Sigma^x at
the 3.1 A cutoff of this cell (= the 3.0 A RI-metric cutoff); at 4x4x4 (4.1 A) RI-RS and exact
RI Sigma^x agree to 1 meV (cn_tensor_4x4x4 vs cn_v6_full_4x4x4: -18.956/-9.767 both). So the
54-atom RI-RS gap is exchange-limited by ~0.3 eV, the 128/250-atom gaps are not, and the
orbital-basis correction is consistent between grid and exact RI.

## One-build comparison (2026-11-08): master build cp2k-v2 (PR 6447 merged, averaged Padé on both sides)

Same inputs, both codes from the rebased branch rirs_periodic_gamma_v2. Dirs: periodic_3d_test/
{cn_v2port_rgc_3x3x3, cn_v2port_full_4x4x4, cn_v2port_full_5x5x5, cn_tgc_4x4x4, cn_pr_kp{4,6,8,10,12}_v2}.
Indirect G0W0 gap over the folded Gamma manifold / HF gap, eV. Full correction kernel.

| Gamma-only cell | atoms | Sigma^x cutoff | G0W0 gap | HF gap | wall (8 ranks/node) |
|---|---:|---:|---:|---:|---|
| 3x3x3 RI-RS | 54 | 3.1 A | 9.210 | 14.833 | |
| 4x4x4 RI-RS | 128 | 4.1 A | 9.970 | 15.123 | |
| 4x4x4 tensor, exact RI | 128 | 4.1 A | 9.970 | 15.125 | 3130 s, 1 node |
| 5x5x5 RI-RS | 250 | 5.1 A | 10.100 | 15.216 | 8277 s, 2 nodes |
| 6x6x6 RI-RS | 432 | 6.2 A | not run: the earlier attempt was still in Z_lP after 6.9 h on 4 nodes (109 s per atom against 33 s per atom for the whole 250-atom run) | | |

| PR full-k mesh | G0W0 gap | HF gap |
|---|---:|---:|
| 4x4x4 | 9.239 | 14.116 |
| 6x6x6 | 9.519 | 14.784 |
| 8x8x8 | 9.593 | 15.011 |
| 10x10x10 | 9.553 | 15.077 |
| 12x12x12 | 9.565 | 15.120 |

- Exact-RI arbiter at 128 atoms (results/exact_ri_4x4x4/): tensor and RI-RS agree to 1 meV in
  EVERY column of every state near the edges. VBM triplet (254-256): eps_DFT 11.917, Sigma^c
  1.950/1.932/1.932, Sigma^x -18.957, v_xc -17.093, G0W0 12.003/11.985/11.985 in both; CBM (257):
  18.377 / -3.206 / -9.767 / -16.568 / 21.973 (tensor) vs 21.972 (RI-RS). The RI-RS grid is
  therefore not the origin of ANY part of the 3D difference to the PR (the 54-atom cell's 0.24 eV
  exchange deficit is the only grid effect, and it is gone at 4.1 A).
- The Gamma-only gap still rises with the cell: 9.210 -> 9.970 -> 10.100 (+0.76, +0.13); the HF
  gap 14.833 -> 15.123 -> 15.216 (+0.29, +0.09). The PR HF gap rises 14.116 -> 15.120 (12^3,
  +0.04 per step at the end), i.e. both exchange series head for the same ~15.2-15.3 eV. The PR
  G0W0 gap is converged to ~0.01 eV at 9.56 eV from 8^3 on. At 250 atoms the Gamma-only gap is
  0.53 eV above it, and the whole difference is Sigma^c (exchange agrees to 0.1 eV and is
  still converging the same way on both sides).
- So the 3D situation after the one-build rerun: the exchange and the numerics are validated
  (exact RI), the correction terms scale as derived, but the Gamma-only head-and-wing model
  gives a gap ~0.5 eV above the full-k RI-RS at convergence. This is a difference between the
  two Gamma-point models (whole-small-BZ analytic head/wing in the orbital basis vs the PR's
  near-Gamma cell model with its Gamma-neighbour transfer), not resolved by a third code (the
  tensor k-point GW is unaffordable, see above).
- Averaged Padé vs standard Padé (same inputs): 3x3x3 9.039 -> 9.210, 4x4x4 9.782 -> 9.970,
  5x5x5 9.934 -> 10.100; PR 4^3 9.200 -> 9.239, 6^3 9.444 -> 9.519, 8^3 9.511 -> 9.593 (the PR
  builds also differ: 4e0d03c43b vs merged master; its HF gaps moved by up to 0.2 eV at 4^3).
