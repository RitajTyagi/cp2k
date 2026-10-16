# RI-RS GW for large unit cells at the Gamma point

Branch `rirs_periodic_gamma`, off `upstream/master`. This note states the equations as they are
actually implemented, explains why the polarizability is projected back onto the RI auxiliary basis
instead of staying in the grid basis, and sets out what an explicitly `R`-resolved formulation would
cost.

Notation: `mu, nu, lambda, sigma` are AO indices, `P, Q` RI auxiliary indices, `l, l'` RI-RS grid
points, `R` a lattice translation of the simulation cell, `i tau` imaginary time and `i omega`
imaginary frequency.

---

## 1. The large-cell Gamma-point framework

One supercell is treated at `k = 0`. Every two-index operator is stored as a single real matrix in
the reference cell, understood as the **minimum image** (MIC) contribution

$$O^{\mathrm{MIC}}_{\mu\nu} \;=\; O^{R_{\mathrm{MIC}}(\mu,\nu)}_{\mu\nu},
\qquad R_{\mathrm{MIC}}(\mu,\nu)=\arg\min_{R}\bigl|\mathbf r_\mu-(\mathbf r_\nu+\mathbf R)\bigr| .$$

A `k`-dependent operator is reconstructed from that single matrix by attaching one phase per atom
pair (`fm_gamma_to_kp`, `cfm_ikp_from_fm_Gamma`):

$$O_{\mu\nu}(\mathbf k)\;=\;e^{\,i\mathbf k\cdot\mathbf R_{\mathrm{MIC}}(\mu,\nu)}\;O^{\mathrm{MIC}}_{\mu\nu}.$$

This is exact only if a single image contributes per pair, i.e. if the cell satisfies

$$a_{\min}/2 \;\gtrsim\; r_G + 2 r_{\mathrm{RI}} \quad\text{(AUX path)} .$$

The approximation is **shared with the tensor implementation** and is basis independent. It is the
reason the two codes converge with supercell size in exactly the same way (Section 6).

---

## 2. RI-RS representation

A set of `N_l` real-space grid points `r_l` is placed per atom. The AO values on the grid are
Bloch-summed over the images of the AO centre,

$$\varphi_{\mu l}\;=\;\sum_{\mathbf R}\phi_\mu(\mathbf r_l+\mathbf R)
\qquad(\texttt{mat\_phi\_mu\_l}),$$

which for a molecule (`cell%perd = 0`) reduces to `phi_mu(r_l)`. The fit coefficients `Z_lP` solve a
Tikhonov-regularized least-squares problem against the three-centre integrals over the truncated
Coulomb metric,

$$(\mathbf D+\alpha\mathbf 1)\,\mathbf Z \;=\; \mathbf d,
\qquad D_{ll'}=\Bigl(\textstyle\sum_\mu \varphi_{\mu l}\varphi_{\mu l'}\Bigr)^{2},
\qquad d_{lP}=\sum_{\mathbf R_j\mathbf R_k}\bigl(\mu^{\mathbf R_j}\nu^{\mathbf R_k}\big|P^{\mathbf 0}\bigr),$$

with `alpha = TIKHONOV`. `Z_lP` therefore maps grid-space products onto the RI basis:
`(mu nu | P) ~ sum_l phi_mu(r_l) phi_nu(r_l) Z_lP`.

**Empirical status:** `Z_lP` and `d_lP` are validated. Swapping only the RI-RS `chi` into the tensor
code's `W` and `Sigma^c` reproduces the tensor result to **4 meV**, so the fit and its image sums are
correct.

---

## 3. Stage equations as implemented

### 3.1 Green's functions

$$G^{\mathrm{occ}}_{\mu\lambda}(i|\tau|)=\sum_n^{\mathrm{occ}} C_{\mu n}\,
e^{-|(\epsilon_n-\epsilon_F)\tau|}\,C_{\lambda n},
\qquad G^{\mathrm{vir}} \text{ analogously over virtuals,}$$

both at `k = 0` (AO x AO, dense but small). On the grid,

$$G_{ll'}=\sum_{\mu\nu}\varphi_{\mu l}\,G_{\mu\nu}\,\varphi_{\nu l'} .$$

### 3.2 Polarizability

$$\chi_{ll'}(i\tau)=G^{\mathrm{occ}}_{ll'}(i|\tau|)\;\circ\;G^{\mathrm{vir}}_{ll'}(i|\tau|),
\qquad
\boxed{\;\chi_{PQ}(i\tau,\mathbf k{=}0)=\sum_{ll'}Z_{lP}\,\chi_{ll'}(i\tau)\,Z_{l'Q}\;}$$

`∘` is the elementwise (Hadamard) product. The boxed projection back onto `P, Q` is the subject of
Section 5. Streamed over panels of contiguous grid row-blocks, so the `grid x grid` matrices are
never materialized.

### 3.3 Screened interaction (`get_W_MIC`, shared with the tensor code)

$$\chi_{PQ}(i\tau)\;\xrightarrow{\text{cos}}\;\chi_{PQ}(i\omega_j)
\;\xrightarrow{\text{one MIC phase}}\;\chi_{PQ}(i\omega_j,\mathbf k)$$
$$\varepsilon(i\omega_j,\mathbf k)=\mathbf 1-\mathbf V^{1/2}(\mathbf k)\mathbf M^{-1}(\mathbf k)\,
\chi(i\omega_j,\mathbf k)\,\mathbf M^{-1}(\mathbf k)\mathbf V^{1/2}(\mathbf k)$$
$$W(i\omega_j,\mathbf k)=\mathbf V^{1/2}\bigl(\varepsilon^{-1}-\mathbf 1\bigr)\mathbf V^{1/2},
\qquad
W^{\mathrm{MIC}}_{PQ}(i\omega_j)=\sum_{\mathbf k} w_{\mathbf k}\,
e^{-i\mathbf k\cdot\mathbf R_{\mathrm{MIC}}}W_{PQ}(i\omega_j,\mathbf k)$$

followed by a Fourier transform to `i tau` and

$$W^{\mathrm{MIC}}(i\tau)\;\leftarrow\;\mathbf M^{-1}(\mathbf k{=}0)\,W^{\mathrm{MIC}}(i\tau)\,
\mathbf M^{-1}(\mathbf k{=}0).$$

**`M^-1(k=0)` must not be regularized here.** `REGULARIZATION_RI` (default `1.0E-2` for periodic
systems) is applied to `M(k)` inside `compute_MinvVsqrt_Vsqrt`, ungated, by *both* codes -- that is
correct and stays. It was additionally applied to `M^-1(k=0)` in `compute_Minv_Gamma` behind
`IF (bs_env%do_gw_ri_rs)`, i.e. for RI-RS only, which made the two implementations differ in `W`
alone by up to 3.8 eV in `Sigma^c`. Removing that gate (commit `47123ddbdb`) is what brings the
periodic path to agreement. Do **not** test this by setting `REGULARIZATION_RI 0.0` in the input:
that also zeroes the shared `M(k)` and the comparison becomes meaningless.

### 3.4 Self-energies

$$\Sigma^x_{\lambda\sigma}=-\sum_{ll'}\varphi_{\lambda l}
\Bigl[D_{ll'}\circ\bigl(\mathbf Z\,\mathbf V^{\mathrm{tr}}\mathbf Z^{\!\top}\bigr)_{ll'}\Bigr]
\varphi_{\sigma l'},\qquad D_{ll'}=\sum_n^{\mathrm{occ}}\dots$$

$$\Sigma^{c,\mp}_{\lambda\sigma}(i\tau)=\mp\sum_{ll'}\varphi_{\lambda l}
\Bigl[G^{\mathrm{occ/vir}}_{ll'}(i|\tau|)\circ
\bigl(\mathbf Z\,W^{\mathrm{MIC}}(i\tau)\,\mathbf Z^{\!\top}\bigr)_{ll'}\Bigr]\varphi_{\sigma l'}$$

then analytic continuation to real frequency and the usual QP equation.

Note `Sigma^x` and `Sigma^c` have identical structure; they differ only in which operator is spread
onto the grid by `Z` -- the truncated Coulomb metric `V^tr` (short ranged) or `W^MIC` (long ranged).

### 3.5 Verified accuracy

Stage-swapping through the restart files (`read_chi` and `all_W_exist` key off file existence, so
renaming one code's `chi`/`W` files into the other's run directory substitutes that stage):

| configuration | d(gap) | median abs dSigma^c |
|---|---|---|
| RI-RS `chi` -> tensor `W` -> tensor `Sigma^c` | -0.001 eV | 4 meV |
| tensor `W` -> RI-RS `Sigma^c` | +0.007 eV | 10 meV |
| full RI-RS after the fix, Si n=2 (64 atoms) | **+0.000 eV** | **0.0 meV** (98 % within 10 meV) |
| full RI-RS after the fix, hBN n=5 (50 atoms) | +0.014 eV | 8.5 meV |

Medians are over states with `|Sigma^c| < 5 eV`; above that the Pade continuation dominates and a
mean over all states is not a meaningful metric. 62/62 regtests pass; the periodic H2O test now
reproduces the tensor value exactly (23.312 instead of 23.329).

---

## 4. Why `phi` is Bloch-summed and not minimum-imaged

Measured, not assumed: evaluating `phi` at the minimum image instead (`GW_RIRS_PHI_MIC=1`) degrades
agreement from 10 meV to 1436 meV once `W` is correct. The Bloch sum is right. Testing this on the
full pipeline before the `W` bug was fixed was useless, because a 3.8 eV error buries a 0.03 eV one.

---

## 5. Why project back to `PQ` instead of staying in `ll'`

The grid basis is used exactly where locality pays and abandoned exactly where it does not.

**The dielectric inversion is the reason.** `chi_ll'` is sparse -- that is the whole point of the
grid representation -- but `epsilon^-1` is **not**. Inverting destroys sparsity, because the Coulomb
kernel is long ranged: every grid point couples to every other through `V`. So the one stage that
can never be local is `(1 - V^(1/2) M^-1 chi M^-1 V^(1/2))^-1`, and it costs `O(n^3)` dense in
whatever basis it is performed.

Measured sizes on this branch:

| system | atoms | grid points `N_l` | `n_RI` | `N_l / n_RI` |
|---|---|---|---|---|
| Si, n=2 | 64 | 27 904 | 2 112 | 13.2 |
| hBN, n=5 | 50 | 16 800 | 1 775 | 9.5 |
| hBN, n=10 | 200 | 67 200 | 7 100 | 9.5 |
| hBN, n=14 | 392 | 131 712 | 13 916 | 9.5 |

A dense inversion in the grid basis therefore costs `(N_l/n_RI)^3 ~ 10^3` times the `PQ` inversion --
roughly **850x to 2300x** -- and it must be repeated for every `k` point and every frequency. On hBN
n=14 that is a dense `131712 x 131712` complex inversion per `(k, omega)`: about 280 GB for one
matrix, before counting the `n_k * n_omega` of them. It is not a constant-factor question.

Two further reasons:

1. **The `k`-point machinery is atom-pair indexed.** `V(k)` is built by a lattice sum over RI
   function pairs, the extrapolation weights `wkp_s_p` / `wkp_no_extra` are selected on the angular
   momentum `l_RI(P) + l_RI(Q)`, and `R_MIC` is defined per atom pair. Grid points carry no such
   atom-pair labelling, so the whole `get_W_MIC` apparatus -- shared with the tensor code and with
   RT-BSE -- would have to be rewritten.
2. **`W` is needed back on the grid anyway, but only as a product.** `Sigma^c` uses
   `(Z W^MIC Z^T)_ll'`, which is formed panel-by-panel and never stored whole. So the grid
   representation of `W` is recovered where it is cheap (a local contraction) rather than carried
   through the expensive step.

The division of labour is therefore: **build `chi` on the grid** (sparse, linear scaling, 27x faster
than the tensor contraction at 64 atoms), **invert in `PQ`** (small, dense, unavoidable), **contract
`Sigma` back on the grid** (sparse again). Staying in `ll'` throughout would trade a 27x win on one
stage for a ~10^3 loss on another.

---

## 6. Explicit `R`-dependence: what it buys and what it costs

### 6.1 Why the current scheme cannot converge in smaller cells

Both implementations compute the same quantity under the same MIC truncation, and they now agree to
meV at a given cell size. Their gap-versus-cell curves therefore coincide: **RI-RS cannot converge at
a smaller `MULTIPLE_UNIT_CELL` than the tensor code**, because it is computing the same number. The
locality of the grid basis buys cost, not cell-size convergence -- what limits convergence is the
one-image-per-pair reconstruction, which is basis independent.

To converge in smaller cells the method has to stop making that approximation.

### 6.2 The formulation

Keep the lattice index instead of collapsing it. With `phi_mu^R(r_l) = phi_mu(r_l + R)` kept
separate rather than summed,

$$G^{\mathbf R}_{ll'}=\sum_{\mu\nu}\varphi_{\mu l}\,G^{\mathbf R}_{\mu\nu}\,\varphi_{\nu l'},
\qquad
\chi^{\mathbf R}_{ll'}(i\tau)=G^{\mathbf R}_{ll'}\circ G^{-\mathbf R}_{l'l},
\qquad
\chi^{\mathbf R}_{PQ}=\sum_{ll'}Z_{lP}\,\chi^{\mathbf R}_{ll'}\,Z_{l'Q}$$

$$\chi_{PQ}(\mathbf k)=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf R}\chi^{\mathbf R}_{PQ}
\;\longrightarrow\;\varepsilon(\mathbf k)\;\longrightarrow\;W(\mathbf k)
\;\longrightarrow\;W^{\mathbf R}_{PQ}=\frac{1}{N_k}\sum_{\mathbf k}e^{-i\mathbf k\cdot\mathbf R}W_{PQ}(\mathbf k)$$

$$\Sigma^{c,\mathbf R}_{\lambda\sigma}(i\tau)=\sum_{ll'}\varphi_{\lambda l}
\Bigl[G^{\mathbf R}_{ll'}\circ\bigl(\mathbf Z W^{\mathbf R}\mathbf Z^{\!\top}\bigr)_{ll'}\Bigr]
\varphi_{\sigma l'},
\qquad \Sigma(\mathbf k)=\sum_{\mathbf R}e^{i\mathbf k\cdot\mathbf R}\Sigma^{\mathbf R}.$$

No single-image selection appears anywhere; `chi(k)` is now a genuine Bloch sum rather than one
matrix times one phase.

### 6.3 Cost

Let the physical volume needed for convergence be fixed. Write it either as one supercell of `M`
primitive cells, or as a primitive cell plus `n_R ~ M` images with `n_k ~ M` k-points.

| stage | supercell (current) | explicit `R` |
|---|---|---|
| `chi` build | `O(M)` (linear-scaling grid) | `O(M)` -- `n_R` images of an `M`-times-smaller cell |
| **dielectric inversion** | `O((M\,n_{RI}^{\mathrm{prim}})^3)=O(M^3)` | `O(n_k\,(n_{RI}^{\mathrm{prim}})^3)=O(M)` |
| `Sigma` contraction | `O(M)` | `O(M)` |
| memory, `W` | `(M\,n_{RI})^2` | `n_k\,(n_{RI})^2 = O(M)` |

**The inversion goes from `O(M^3)` to `O(M)`.** That is the real prize, and it is the same reason
k-point codes beat supercells generally. Concretely, the `n_RI^3` wall that stopped the hBN sweep --
`n_RI` 1 775 -> 13 916 from n=5 to n=14, a `10^3` growth in the inversion -- would flatten.

The costs added are: `n_R` times more `chi` and `Sigma` panel work (each on a much smaller cell, so
roughly neutral in total), `n_k` dielectric inversions instead of one (but of a matrix `M` times
smaller in each dimension), and storage of `n_R` real-space operator images instead of one.

### 6.4 This framework already exists in CP2K

`gw_tensor_small_cell_full_kp.F` is exactly this: explicit `R`-resolved operators, `V^tr(R)` via
`get_V_tr_R`, and `heuristic_RI_regularization` to condition `M(k)` across the k-mesh (which is why
`REGULARIZATION_RI` exists at all, and why it is applied there and not in the large-cell path).

So the proposal is **not** a new formulation -- it is to give RI-RS the small-cell full-k-point
treatment rather than the large-cell Gamma one. That de-risks it considerably: the k-point
bookkeeping, the `R`-vector lists, the regularization heuristic and the extrapolation weights are
all written and tested. What RI-RS would contribute is the linear-scaling, grid-based construction
of `chi^R` in place of the tensor contraction -- which is where the measured 27x already comes from.

The open questions before committing to it are the ones to settle on paper first:

1. whether `Z_lP` can be built once for the primitive cell and reused for every `R`, or whether the
   fit itself has to carry an `R` index (which would change the cost table);
2. how many `R` vectors are actually significant for `chi^R` at a given accuracy, since the grid
   locality bound differs from the AO-product bound the tensor code screens on;
3. whether the Hadamard product `G^R ∘ (Z W^R Z^T)` can reuse the existing panel machinery unchanged,
   given both operands now carry the same `R`.

---

## 7. Grid-resolved MIC (notes Sec. 7-8): measured feasibility

### 7.1 How crude is the atom-pair MIC?

The grids were reproduced exactly from `data/ri_rs_grid/*.ion` (336 pts/atom for B and N -> 16 800
for hBN n=5; 436 for Si -> 27 904 for Si n=2, both matching the runs). For sampled `(P,Q)` blocks the
combined lattice vector of Eqs.~(101)/(104),

$$\boldsymbol{\mathcal R}(\ell,\ell',P,Q)=\mathbf R_{\mathrm{MIC}}(\ell,\ell')
+\mathbf R_{\mathrm{MI}}(\ell',Q)-\mathbf R_{\mathrm{MI}}(\ell,P),$$

was histogrammed over grid pairs inside the `a_min/2` window:

| system | `a_min/2` | `cutoff_ri` | distinct `R` per `(P,Q)` block | weight in `R = 0` |
|---|---|---|---|---|
| hBN n=5 | 5.42 A | 11.21 A | median **8** (max 9) | **0.385** |
| Si n=2 | 5.43 A | 7.76 A | median **19** (max 26) | **0.230** |

The current code forces *all* of that weight onto a single phase per `(P,Q)` block. At these cells
only 23-38 % of the contribution actually belongs to `R = 0`, so 62-77 % is mis-phased at every
`k != 0`. That is a quantitative measure of how much accuracy the grid-resolved MIC can recover, and
it is large.

Consistency check: at n=5 the measured count is 8, which is also `(14/5)^2 = 7.8`, the number of
images needed to cover the n=14 volume. So `n_R ~ M` holds empirically.

### 7.2 Memory: large sparse cell vs small cell plus `R`

Per-rank at **8 MPI ranks on a 256 GB node (32 GB/rank)**, hBN, from measured components
(`Z_lP` solve peak and panel peak are printed by the code; `W(i tau)`, `eps`, `W(k)` are dense `cp_fm`,
block-cyclic distributed):

| cell | atoms | `n_RI` | `Z_lP` solve | `W(i tau)` set | `eps, W(k)` | panels | **total/rank** |
|---|---|---|---|---|---|---|---|
| n=5 | 50 | 1 775 | 2.3 GB | 0.03 | 0.01 | 0.5 | **2.8 GB** |
| n=8 | 128 | 4 544 | 14.8 GB | 0.21 | 0.08 | 1.2 | **16.2 GB** |
| n=10 | 200 | 7 100 | 18.8 GB | 0.50 | 0.20 | 1.8 | **21.3 GB** |
| n=12 | 288 | 10 224 | 18.8 GB | 1.05 | 0.42 | 2.6 | **22.9 GB** |
| n=14 | 392 | 13 916 | 18.8 GB | 1.94 | 0.77 | 7.0 | **28.5 GB** |

Same physical volume as n=14, reached instead with a small cell plus `n_R` images:

| cell | atoms | `n_RI` | `n_R` | `chi^R + W^R` | `Z_lP` solve | **total/rank** |
|---|---|---|---|---|---|---|
| n=5 | 50 | 1 775 | 8 | 0.28 GB | 2.3 GB | **2.5 GB** |
| n=7 | 98 | 3 479 | 4 | 0.53 GB | 8.7 GB | **9.2 GB** |
| n=10 | 200 | 7 100 | 2 | 1.11 GB | 18.8 GB | **19.9 GB** |
| n=14 | 392 | 13 916 | 1 | 2.13 GB | 18.8 GB | **21.0 GB** |

**Small cell plus `R` wins by about 11x** (2.5 GB vs 28.5 GB per rank for the same physics), and the
reason is specific:

* The dominant term is the **`Z_lP` Gram matrix** `D_local`, which is `N_l(cell)^2` and is *not*
  distributed (it is one atom's solve; the classifier only splits it over a ScaLAPACK subgroup).
  Shrinking the cell shrinks it **quadratically**: 48 524^2 -> 16 800^2, i.e. 18.8 GB -> 2.3 GB.
* The price of `R`-resolution is paid in **PQ space**, where the objects are small *and* distributed,
  and it grows only **linearly** in `n_R`: 0.28 GB per rank at n=5.
* `D_local` is bounded by the integration sphere (`cutoff_ri = r_c + r_RI = 11.2 A`, 48 524 points),
  so in the large-cell scheme it **saturates and never improves** -- it is 66 % of the n=14 budget and
  no larger cell will reduce it. The only other lever on it is lowering `CUTOFF_RADIUS_RI`.

Sparsity does not rescue the large cell at the sizes of interest: `G/chi/W` on the grid are **100 %
dense** at hBN n=5 and Si n=2, and only reach 34.9 % at n=14. The sparsity that the large-cell
argument relies on appears exactly where the memory has already saturated.

### 7.3 What is not yet established

The numbers above say the grid-resolved scheme is *affordable* and that the atom-pair MIC is *crude*.
They do **not** show that the grid-resolved MIC is converged at n=5: prerequisite Eq.~(74) still has
to hold for `G, chi, Sigma` at grid resolution, and with `r_AO = 9.6 A` against `a_min/2 = 5.4 A` it
may not. Testing that needs the implementation; the honest expectation is that it shifts the required
cell down by some factor, not that it removes the requirement.

One practical advantage over the notes as written: collecting the phases into
`chi_PQ^{\mathcal R}` and Bloch-summing in PQ space keeps every grid-space matrix **real**, which
sidesteps the concern raised in Sec.~8 that DBCSR has no complex support.
