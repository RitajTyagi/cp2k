# Diamond (primitive 2-atom fcc cell, a = 3.567 A), PERIODIC XYZ, GPW, SZV-MOLOPT-SR-GTH,
# RI_20_2_2_1_1_0_0_0_1.6e-02, GTH-PBE-q4, CUTOFF 400, 10 minimax points, EPS_FILTER 1e-12,
# SIZE_LATTICE_SUM 2, default RI cutoff (half the cell height). 2026-11-06.

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
