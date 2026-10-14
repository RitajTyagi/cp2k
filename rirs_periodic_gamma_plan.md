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

---

# Progress

## Done (branch `rirs_periodic_gamma`, pushed to `origin`)

| commit | what |
|---|---|
| `e33df85948` | `gw_ri_rs_panels` extraction; topology helpers unified into `gw_utils_dbcsr` |
| `86d94220bd` | minimum-image panel screening (`min_image_dist2`) |
| `47b851677e` | periodic Gamma-point regtest (none existed) |
| `94a543d27b` | opt-in periodic image sums in the shared front end |
| `be274bcfc4` | exact minimum image + Bloch-summed phi in the Z_lP sphere |
| `f98648e643` | periodic driver runs on the shared optimized code |

Net against upstream: **+1856 / -2785 lines**, i.e. 929 lines smaller while adding capability.

| file | upstream | now |
|---|---|---|
| `gw_ri_rs_large_cell_gamma.F` | 1442 | **235** |
| `gw_ri_rs_non_periodic.F` | 2793 | 1410 |
| `gw_ri_rs_panels.F` | – | 1437 |
| `gw_utils_dbcsr.F` | 155 | 246 |

The periodic path now has panel streaming (no grid x grid matrix is ever materialized), atom-aligned
grid blocking, the optimized Z_lP solve, the cutoff keywords, the memory estimate and evGW0. Only
`get_W_MIC` is still periodic-specific.

## Verification

All 9 regtests pass at every step; the 8 molecular ones are bit-identical throughout.

| test | value |
|---|---|
| 01..08 molecular | 23.847 / 24.168 / 21.699 / 23.324 / 23.517 / 23.724 / 23.674 / 21.416 |
| 09 periodic RI-RS | 23.329 |
| 09 system, tensor code | 23.312 |

So the RI-RS approximation error on the periodic path is 0.017 eV.

Note: on clean upstream, test 01 gives 23.847 where `TEST_FILES.toml` says 23.846. That 1 meV is
pre-existing upstream, not from this work.

## Efficiency benchmark: silicon

`si_benchmark/` holds `si{n}_{rirs,tensor}.inp`: bulk Si, diamond structure, conventional cubic cell
(8 atoms, a = 5.431 A), supercell set by `MULTIPLE_UNIT_CELL` in **both** `&CELL` and `&TOPOLOGY`.

    n=1:   8 atoms,  5.43 A    n=3: 216 atoms, 16.29 A
    n=2:  64 atoms, 10.86 A    n=4: 512 atoms, 21.72 A

**The DOS k-mesh has to be Gamma only.** `check_positive_definite_overlap_mat`
(`post_scf_bandstructure_utils.F:1002`) loops over `kpoints_DOS` and reconstructs S(k) from
S(Gamma) by the minimum image; that reconstruction is not positive definite for these cell sizes,
and the run aborts with "the cell of the calculation is too small" -- at n=1 (5.43 A) and still at
n=2 (10.86 A), for the tensor code as much as for RI-RS. Setting `&DOS KPOINTS 1 1 1` removes the
reconstruction entirely, and it is also the physically consistent choice: in a Gamma-only supercell
scheme the HOMO/LUMO/gap to converge against cell size ARE the Gamma-point values.

Reaching the ~20 A where the tensor code is said to converge needs n=4, i.e. 512 atoms -- a
production job, which is precisely the regime the RI-RS linear scaling is meant for.

Run on the noctua login node (too heavy for the laptop), in
`/scratch/hpc-prf-metdyn/metdyn07_Ritaj/claude_periodic_test/`.

### noctua setup notes

The branch could not be checked out in `implementation/github/cp2k-dev`: upstream now tracks
`tests/QS/regtest-cohsex/`, which collides with untracked WIP files of the same name there. Rather
than move those, the benchmark uses a **separate git worktree** so that checkout and its build stay
untouched:

    git worktree add /scratch/.../claude_periodic_test/cp2k-rirs rirs_periodic_gamma

That needs its own configure, and the first attempt silently produced a binary **without libint**,
which aborts at the start of any GW run. The working configure adds
`-DCP2K_USE_LIBINT2=ON` plus `-DLibint2_DIR=<spack libint-2.11.2>/lib/cmake/libint2`.

## Per-state comparison on bulk Si: what works and what does not

`si1_ri33_{rirs,tensor}.inp`, 8 atoms, a = 5.431 A, DZVP-MOLOPT-GTH / RI_33, DOS k = Gamma.
Both runs print `bandstructure_SCF_and_G0W0`, so every stage can be compared state by state.

### Timings (8 ranks, noctua login node)

| run | wall | chi per tau | W per k batch |
|---|---|---|---|
| RI-RS  |  404 s | 100.0 s (10 points) | 40.4 s |
| tensor | 5067 s | 175.8 s             | 39.7 s |

At 64 atoms the chi gap is much wider: RI-RS 55 s vs tensor 1478 s per time point, i.e. 27x. The
RI-RS speedup is real, but at these sizes `get_W_MIC` -- which is **shared** with the tensor code
and scales as n_RI^3 -- already costs more than everything RI-RS accelerates (10750 s per k batch
at 64 atoms with RI_TZ). Goal 1 is met; the next bottleneck is not in the RI-RS code.

### Sigma^x is fine, Sigma^c is not

eps_nk^DFT and v_nk^xc are exactly degenerate in both codes, so the SCF and the grid are consistent.

* **Sigma^x**: RI-RS agrees with tensor to ~0.1 eV and breaks the symmetry degeneracy by only
  0.03-0.08 eV. That is the RI-RS grid fit error -- the grid is not symmetry adapted, so some
  splitting is expected at that level.
* **Sigma^c**: RI-RS is uniformly **2-8x too negative** (mean ~3.8x) and breaks degeneracy by
  **1-3 eV**.

      eps^DFT  mult |  spread Sigma^x        |  spread Sigma^c
                    | tensor   RI-RS         | tensor   RI-RS
       -2.266    6  |  0.000   0.045         |  0.039   2.151
        2.670    6  |  0.000   0.068         |  0.020   1.157
        6.342    6  |  0.000   0.080         |  0.060   3.133
       15.843    6  |  0.000   0.328         |  0.029   1.873
       17.703    6  |  0.000   0.031         |  0.058   1.644

### This is pre-existing, not from this port

The same input run with the *other* periodic RI-RS implementation (the 1442-line
`gw_large_cell_gamma_ri_rs.F` in `implementation/github/cp2k-dev`) gives **the same numbers to
1 meV**: the only differences in the whole file are five high-lying virtuals between 43 and 96 eV,
each 0.001 eV. So

1. the Sigma^c defect is a property of the periodic RI-RS Gamma formulation, not of this branch, and
2. the agreement is an independent validation of the port -- two separately written drivers, one of
   them 235 lines on top of the shared molecular code, land on the same answer.

### Why Sigma^c and not Sigma^x: the cell is far outside the MIC window

The run prints, for this cell,

    r_AO = 8.696 A      r_RI = 4.757 A      a_min/2 = 2.72 A

The AUX-path requirement is `a_min/2 >~ r_G + 2 r_RI`, violated here by about a factor 6. Every RI-RS
grid operator is built from Bloch-summed phi on a grid confined to one cell, so

    Sigma^c_lambda,sigma = sum_{l,l'} [sum_R1 phi] [sum_R2 G] [sum_R3 W] [sum_R4 phi]

while the Gamma-point expression needs a *single* image sum tying the four factors together. The
product of Bloch sums is not the Bloch sum of products, and with r_AO/a = 1.6 many images are in
range, so each grid pair collects image contributions multiplicatively. That is consistent with both
observations: a systematic overcount of a few, and state-dependent scatter.

Sigma^x escapes it because its second grid operator is the **truncated** Coulomb metric
(`CUTOFF_RADIUS_RI`, 3.0 A by default) rather than W: short-ranged enough that essentially one image
contributes, so `Z V Z^T` stays clean.

This also explains why the periodic H2O regtest agrees with the tensor code to 0.017 eV in a 6 A box
despite formally violating the same condition: a molecule in vacuum has negligible density in the
wrapped tails, so only R = 0 contributes whatever the formal radii say.

### Consequence for goal 2

The evidence does not support the hope that RI-RS converges at *smaller* supercells. The grid-space
Sigma^c contraction needs `a_min/2` above the range of both G and W, which is a **stronger** cell
requirement than the tensor code's, not a weaker one. Making RI-RS work at small cells means keeping
explicit R dependence in the grid operators and contracting `sum_R G^R o W^R` instead of multiplying
Bloch sums -- a change to the formulation, well beyond porting the molecular optimizations.

Two things would settle it, both too heavy for a login node:

* n=2 (64 atoms, a = 10.86 A, `a_min/2` = 5.43 A): if the mechanism above is right, the Sigma^c
  overcount factor should fall and the degeneracy spread shrink. Cost is dominated by
  `get_W_MIC` at n_RI = 2112 (n_RI^3), order hours.
* n=4 (512 atoms, 21.72 A), the cell where the tensor code is said to converge.

## n=2 settles it: the 8-atom result was a domain-of-validity artefact, not a bug

`si2_ri33_rirs.inp`, 64 atoms, a = 10.86 A, so `a_min/2` doubles from 2.72 to 5.43 A. Everything
else identical (RI_33, `KPOINTS_W 2 2 2`, `&DOS KPOINTS 1 1 1`, 10 time points). One node,
16 MPI x 8 OMP: **1062 s**.

| quantity | n=1 (a_min/2 = 2.72 A) | n=2 (a_min/2 = 5.43 A) |
|---|---|---|
| Sigma^x spread over degenerate sets | 0.03 - 0.33 eV | **0.000 eV** |
| Sigma^c spread over degenerate sets | 1.2 - 3.1 eV | **0.006 - 0.018 eV** |

Sigma^x is now *exactly* degenerate and the Sigma^c spread has fallen by a factor of order 100-200.
Doubling the minimum image distance removed the symmetry breaking almost entirely, which is what
the "product of Bloch sums is not the Bloch sum of products" explanation predicts and what a genuine
bug would not do. The periodic RI-RS Gamma path is behaving correctly; the 8-atom cell was simply
outside the method's window, with r_AO = 8.70 A against `a_min/2` = 2.72 A.

The earlier conclusion drawn from n=1 alone -- that RI-RS needs *larger* cells than the tensor code
-- does not survive this. What n=1 showed was a cell six times too small for the AUX condition, not
a property of the formulation.

### Still open at n=2

Sigma^c comes out **positive**, about +2.7 to +3.3 eV on the occupied states, where the tensor code
at n=1 gave -0.4 to -1.5 eV. The degeneracy is clean, so this is a magnitude question, not a
symmetry one, and it needs the tensor reference at the same cell to judge. Absolute gaps at these
cell sizes are meaningless anyway (Si G0W0: 3.68 eV at n=1 RI-RS, 5.75 eV at n=1 tensor, 3.60 eV at
n=2 RI-RS, against ~1.2 eV converged), because a Gamma-only supercell with `KPOINTS_W 2 2 2`
converges slowly -- which is exactly the convergence this sweep is meant to map.

`si2_ri33_tensor` needs more than one hour: chi alone is 296 s per time point. It is chained over
six 1 h passes using `SCF_GUESS RESTART`, which resumes from the chi matrices written by the
previous pass; each pass exits immediately once `DONE.txt` exists.

### W cost: the n_RI^3 estimate was far too pessimistic

n=2 RI-RS completed in 1062 s where a naive n_RI^3 extrapolation from n=1 predicted ~13 h on a node.
`get_W_MIC` is not cubic-dominated in this range. The hBN sweep below is sized on measurement
rather than on that extrapolation.

## 2D hBN supercell sweep

`claude_periodic_test/hbn/`, built from the existing `hbn-rirs/gw/hbn.inp`:
aug-SZV-MOLOPT-GTH-tier-2 (13 AO/atom), RI 26 (B) / 45 (N) = 35.5 RI/atom, `PERIODIC XY` with MT
Poisson, `CUTOFF_RADIUS_RI 5`, `REGULARIZATION_RI 0.001`, 10 time points, `&DOS KPOINTS 1 1 1`.
n = 5, 6, 7, 8, 10, 12, 14, each with RI-RS and tensor, as a job array capped at two concurrent.

2D is a much better probe than bulk: `a_min/2` grows linearly in n while the atom count grows as
n^2, so the minimum-image margin improves far faster per unit cost. n=8 already reaches the ~20 A
where the tensor code is said to converge.

| n | atoms | in-plane a | a_min/2 | n_RI |
|---|---|---|---|---|
|  5 |  50 | 12.5 A |  6.3 A |  1775 |
|  8 | 128 | 20.0 A | 10.0 A |  4544 |
| 10 | 200 | 25.0 A | 12.5 A |  7100 |
| 14 | 392 | 35.1 A | 17.5 A | 13916 |

`EXTENDED_FFT_LENGTHS` had to be dropped from the old input -- the keyword no longer exists in
current CP2K and aborts the parser.

## hBN results: the 2D Coulomb treatment is where RI-RS and tensor part

### RI-RS convergence curve (2D, CUTOFF_RADIUS_RI 5, one node, 8 MPI x 16 OMP)

| n | atoms | in-plane a | a_min/2 | G0W0 gap | wall |
|---|---|---|---|---|---|
|  5 |  50 | 12.51 A |  6.26 A | 11.848 eV |   327 s |
|  6 |  72 | 15.01 A |  7.51 A |  9.737 eV |   693 s |
|  7 |  98 | 17.52 A |  8.76 A |  9.659 eV |  1519 s |
|  8 | 128 | 20.02 A | 10.01 A |  9.428 eV |  3437 s |
| 10 | 200 | 25.02 A | 12.51 A |  8.632 eV |  6050 s |
| 12 | 288 | 30.03 A | 15.01 A |  7.544 eV | 20280 s |

The Gamma-only supercell gap is **still falling at 30 A**, so 2D hBN is not converged even there --
a harder case than the ~20 A quoted for bulk. RI-RS put 288 atoms within 5.6 h on a single node.

### Which of the three differences causes the Sigma^c discrepancy

hBN differs from Si in three ways at once, so each was isolated at n=5 with both codes:

| variant | mean abs dSigma^c | mean abs dSigma^x | gap RI-RS | gap tensor |
|---|---|---|---|---|
| baseline: 2D, `CUTOFF_RADIUS_RI 5` | 3.784 eV | 0.006 eV | 11.848 | 11.425 |
| 2D, `CUTOFF_RADIUS_RI 3`           | 3.517 eV | 0.007 eV |  9.686 |  9.581 |
| **3D stack (`PERIODIC XYZ`), RI 5**| **0.585 eV** | 0.005 eV | 11.890 | 12.466 |

* **z-wrapping under `PERIODIC XY` -- ruled out.** `pbc1` and `pbc3` in `subsys/cell_types.F` both
  multiply each component by `cell%perd`, so nothing wraps along the vacuum axis.
* **`CUTOFF_RADIUS_RI` -- ruled out.** Dropping 5 -> 3 A moves the absolute gap by over 2 eV in
  *both* codes but leaves the disagreement at 3.5 eV. It is an accuracy knob, not the cause.
* **The 2D Coulomb / k-extrapolation treatment -- confirmed.** Making the same cell `PERIODIC XYZ`,
  changing nothing else, collapses the Sigma^c discrepancy by a factor of 6.5. The aug- basis is
  identical across all three variants, so it is not the basis either.

Sigma^x agrees to 0.005-0.007 eV in every variant, which keeps the Z_lP fit and phi out of the
frame: the discrepancy travels through chi -> W. In 2D the Coulomb kernel V(k) ~ 1/|k| diverges at
Gamma and is handled by the extrapolation weights `wkp_s_p` / `wkp_no_extra` selected on
`bs_env%l_RI` in `MIC_contribution_from_ikp`. Those weights are applied to a chi that the two codes
build differently, and the 1/|k| divergence amplifies whatever difference exists at the Gamma term.
That is the next thing to look at.

### Restart granularity: a trap worth recording

The GW restart is staged SCF (wfn) -> chi (per-tau matrices) -> **W^MIC (all or nothing)** -> Sigma^c.
`W^MIC` is written only once every k-batch has completed, so a job killed in the middle of W discards
the whole stage and the next pass repeats it from batch 1. Chained one-hour passes therefore make no
progress at all on any case whose W stage exceeds the walltime -- two consecutive `p3d_rirs` passes
stopped at the identical batch 15/35. Si n=2 tensor only resumed cleanly because it died *after* W,
during Sigma^c. Size the walltime to fit the entire W stage.
