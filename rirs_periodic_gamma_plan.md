# Efficient RI-RS GW for large-cell Gamma-point calculations

**Branch:** `rirs_periodic_gamma`, forked from `upstream/master` @ `badc0e2145`.

## Goal

Bring the periodic Gamma-point RI-RS path (`src/gw_ri_rs_large_cell_gamma.F`) up to the efficiency
of the molecular path, reusing the molecular optimizations and nothing more exotic. Two targets:

1. **RI-RS should clearly beat the tensor code** (`gw_tensor_large_cell_gamma.F`) for periodic
   systems, as it already does for molecules.
2. **Smaller converged supercells.** The tensor code needs roughly a 20 Angstrom supercell before
   HOMO/LUMO/gap converge; the hope is RI-RS converges sooner. This is a wish, not something to
   chase at the cost of efficiency.

## Scope decisions

- **One W implementation only: the auxiliary-basis `get_W_MIC`.** A grid-basis W was prototyped in
  an earlier session and measured ~12x the cost of the auxiliary path for the W build, because the
  minimum-image phase e^{i2 pi k.n(l,l')} is not separable into n(l) - n(l'), so chi(k) cannot be
  pre-contracted and the whole grid->RI contraction is repeated at every k. It is dropped.
- No new approximations, no new keywords beyond what the molecular path already exposes.
- Coding style follows the upstream RI-RS molecular code.

## Starting point

Upstream has already split the shared machinery into modules:

| module | role |
|---|---|
| `gw_ri_rs_compute_Z_lP.F` | optimized Z_lP solve (shared) |
| `gw_ri_rs_utils.F` | radii, AO evaluation on points (shared) |
| `gw_ri_rs_grid_setup_main.F` | grid assembly (shared) |
| `gw_ri_rs_non_periodic.F` | molecular driver **plus the panel machinery** |
| `gw_ri_rs_large_cell_gamma.F` | periodic driver, **stale unoptimized duplicate** |

The periodic file still carries its own `atomic_basis_at_grid_point`, `fill_phi_for_atom`,
`compute_Z_lP`, `compute_d_lp`, and chi / Sigma^x / Sigma^c that materialize full grid x grid
matrices. The molecular file carries the panel machinery that the periodic path needs.

## Work items

### 1. New shared module `src/gw_ri_rs_panels.F`

Move the 15 panel routines (lines 764-2077 of `gw_ri_rs_non_periodic.F`, 1314 lines) into their own
module: `mask_grid_blocks_near_panel`, `panel_template_elems`, `panel_mem_estimate_GB`,
`plan_grid_panels`, `resolve_grid_panels`, `print_ri_rs_memory_estimate`,
`build_geo_template_panel`, `extract_grid_panel`, `collect_used_col_blocks`,
`extract_masked_blocks`, `reserve_blocks_within_radius`, `create_product_matrix`, `build_G_ao`,
`contract_grid_panels`, `contract_grid_panels_sigma_c`.

The molecular file shrinks rather than grows, and `gw_ri_rs_*` matches upstream naming.
Needs a `src/CMakeLists.txt` entry.

### 2. PBC-correct the panel screening

The geometric tests use plain Euclidean distances. Since `d_euclid >= d_MIC`, a Euclidean cutoff
test **under-includes** under PBC: it silently drops boundary-wrapping pairs that are genuinely
inside the cutoff. Measure through `bs_env%ri_rs%cell`; for a molecule `perd = 0`, `pbc()` is the
identity and molecular results stay bit-identical.

`mask_grid_blocks_near_panel` needs more than a substitution: it measures the distance to an
axis-aligned bounding box, which is meaningless modulo the lattice. Minimize over the
`perd`-gated images. Its invariant is that `used` must be a *superset* of what the exact per-pair
test keeps -- over-inclusion only widens the neighbourhood slice, under-inclusion is the bug.

### 3. One front end for both drivers

Two shared routines are molecular-only today; the periodic bodies are strict generalizations that
collapse to the molecular ones when `perd = 0`:

- `evaluate_ao_basis_on_points` / `evaluate_ao_on_points` (`gw_ri_rs_utils.F`) take a single
  minimum image. Add an **opt-in** periodic image sum. Opt-in, not a semantics change, because
  `gw_ri_rs_grid_initialization.F` and `gw_ri_rs_grid_optimization.F` evaluate isolated atoms and
  must not get image sums.
- the d_lP build inside `gw_ri_rs_compute_Z_lP.F` needs the `(cell_R, cell_S)` image sum that the
  periodic `compute_d_lp` has.

That deletes `fill_phi_for_atom`, `compute_Z_lP` and `compute_d_lp` from the periodic file.

### 4. Rewire the periodic driver

Use the shared grid setup, the shared Z_lP, the shared `atomic_basis_at_grid_point`, and
panel-streamed chi / Sigma^x / Sigma^c. Keep `get_W_MIC` for W. Delete `build_G_grid` and the three
grid x grid routines. evGW0 follows for free through `compute_Sigma_c_and_QP_energies`.

Expected: `gw_ri_rs_large_cell_gamma.F` 1442 -> ~300 lines; `gw_ri_rs_non_periodic.F` -1314 lines.

### 5. Tests

There is **no periodic RI-RS regtest upstream** -- all 8 in `tests/QS/regtest-gw-realspace` are
molecular. Add one, sized for CI (< 20-30 s on 1 rank).

Validation at every step:
- the 8 molecular regtests must stay exact (items 2 and 3 touch code the molecular path uses)
- the periodic result must match the tensor reference on the same system
- timing: RI-RS vs `gw_tensor_large_cell_gamma` on a periodic system

Heavier runs (supercell convergence) go to noctua.
