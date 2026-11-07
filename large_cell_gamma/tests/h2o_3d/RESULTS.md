# Water in a 6 A cube, PERIODIC XYZ, GAPW, 10 minimax points (regtest 09 settings)

First 3D test of the Gamma-only path with the head-and-wing correction (2026-11-06).

| run | W treatment | G0W0 VBM | CBM | gap (eV) |
|-----|-------------|---------:|----:|---------:|
| h2o_k2 (= regtest 09) | KPOINTS_W 2 2 2, old k-mesh path | -5.665 | 17.647 | 23.312 |
| h2o_k3 | KPOINTS_W 3 3 3 (odd: Gamma dropped, 8 non-Gamma points) | | | 23.365 |
| h2o_k4 | KPOINTS_W 4 4 4 | -5.651 | 17.627 | 23.278 |
| h2o_k6 | KPOINTS_W 6 6 6 | | | (running) |
| h2o_g1 | KPOINTS_W 1 1 1 + 3D correction (Ewald V(0), SIZE_LATTICE_SUM 3, omega 0.176, 461 G) | -5.581 | 17.658 | 23.239 |
| h2o_g1_s2 | same, SIZE_LATTICE_SUM 2 (omega 0.706, 30647 G) | -5.581 | 17.658 | 23.239 |

- Ewald omega-independence: identical to the meV between s=2 and s=3 -> both halves of the
  Ewald V(0) (erfc lattice sum + reciprocal sum with the G=0 term omitted) are consistent.
- The head is tiny for a molecule in a box (S_xx ~ 2e-3 at the first frequency), so this test
  exercises the Ewald V(0), the unprojected maps and the 3D plumbing, not the kernel magnitudes.
- Gamma-only sits 39 meV below the 4x4x4 mesh value; the mesh series alternates (odd/even).
