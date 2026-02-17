# Integrating Alvarez-Melis et al. (2019) “OT with Global Invariances” into `manifoldalign`

Goal: implement the **Alvarez-Melis / Jegelka / Jaakkola (AISTATS 2019)** method
(*Towards Optimal Transport with Global Invariances*) in a way that is consistent
with the patterns, utilities, and backends already used in this package.

This note focuses on the **orthogonal invariance** specialization (the practical
word-embedding alignment case), and outlines an implementation path that can be
extended later.

---

## 1) Relevant existing patterns in this repo

### 1.1 Data conventions

- Throughout the package, domain matrices are **samples × features** (`n × d`).
- `hyperdesign` domain blocks store data as `block$x`.
- Many algorithms accept:
  - `hyperdesign` (multi-domain),
  - `list` of matrices (multi-domain), and/or
  - pairwise `(X_i, X_j)` via the `align_many()` adapter interface.

### 1.2 S3 “aligner adapter” framework

- `R/align_framework.R` defines:
  - `new_aligner(name, group, supports_multi)`
  - generics `fit_pair()`, `relative_transform()`, `pair_loss()`, `latent_dim()`
  - `align_many()` (in `R/align_many.R`) which:
    1) fits pairwise models on a dataset graph,
    2) extracts relative transforms, and
    3) synchronizes transforms by group (`"O"`, `"GL"`, `"perm"`).
- Existing adapters:
  - `grasp_aligner()` (group `"O"`)
  - `parrot_aligner()` (group `"perm"`)
  - etc.

**Takeaway:** the cleanest “package-consistent” way to add Alvarez-2019 is as a
new **pairwise aligner** with group `"O"` (orthogonal), so it works immediately
with `align_many()` + `rotation_sync()`.

### 1.3 OT utilities and solvers already present

Balanced OT (entropic Sinkhorn):
- `sinkhorn_unified()` (RcppArmadillo; `src/sinkhorn_unified.cpp`)
  - supports log-domain stabilization when `epsilon < 0.1`
  - expects marginals `a`, `b` with equal total mass
  - returns a **dense** plan `P` (transport matrix)

Exact OT (network simplex):
- `network_simplex_ot_cpp()` (Rcpp; `src/ot_solvers.cpp`) used as an optional
  backend in `fpgw` helpers.

Fast squared distances:
- `.Call("_manifoldalign_compute_squared_distances_cpp", X1, X2)` is wrapped by
  `compute_squared_distances_cpp_fallback()` in `R/zzz_cpp_fallbacks.R`.

Internal OT helper patterns (R code):
- `R/utils-ot.R` provides:
  - `validate_marginals()`, `prepare_ot_data()`
  - KNN helpers for out-of-sample mapping
  - cost builders like `compute_feature_cost()` (returns *Euclidean*, not squared)

**Takeaway:** Alvarez-2019 (orthogonal case) should rely on:
- `compute_squared_distances_cpp_fallback()` for `||x - y||^2` costs, and
- `sinkhorn_unified()` for entropic OT.

---

## 2) Alvarez-2019 method recap (orthogonal specialization)

For two domains `X ∈ R^{n×d}` and `Y ∈ R^{m×d}`, optimize over:
- coupling `π ∈ Π(a,b)` (balanced OT polytope), and
- orthogonal map `P ∈ O(d)`,

typically by alternating:

1) **Coupling update (OT):**
   - compute transformed points `Y_P = Y %*% P` (row convention)
   - cost matrix `C_{ij} = ||X_i - (Y_P)_j||^2`
   - solve entropic OT: `π ← Sinkhorn(C, a, b, epsilon)`

2) **Orthogonal update (Procrustes from coupling):**
   - compute cross-covariance-like matrix `M = t(X) %*% π %*% Y` (`d×d`)
   - SVD: `M = U Σ V^T`
   - set `P_col = U V^T` as the maximizer of `⟨P, M⟩` over `O(d)`
   - in row convention the applied operator is `P = t(P_col)` if you want `Y %*% P`
     (i.e., be explicit about row/column conventions; see §3.2 below)

Practical stabilization:
- entropy annealing schedule: start at large `epsilon0`, decay by factor `alpha`
  each outer iteration until `epsilon_min`.

---

## 3) Implementation design consistent with `manifoldalign`

### 3.1 Where it should live

Recommended minimal integration:
- Add a new adapter file, analogous to:
  - `R/adapters_grasp.R`
  - `R/adapters_parrot.R`

Proposed new file:
- `R/adapters_invar_ot.R`

This would export:
- `invar_ot_aligner()` (or `ot_procrustes_aligner()`; naming choice below)

And implement:
- `fit_pair.<aligner>()`
- `relative_transform.<fit_class>()`
- `pair_loss.<fit_class>()`
- `latent_dim.<fit_class>()`

This immediately unlocks:
- pairwise usage (`fit_pair()`),
- multi-domain synchronization (`align_many()` with group `"O"`),
- transformation application (`apply_transform()`).

Optional later (if desired):
- add a “direct” user function (non-adapter) like `invar_ot_pair(X, Y, ...)`
  that returns the same fit object but without the aligner wrapper.

### 3.2 Row/column convention: avoid transposition bugs

The package uses **row vectors** (`n×d`), i.e. data points are rows.

If we store an operator `op` intended for `apply_transform()` (which does
`X %*% op`), then `op` must be the **row-acting** map.

To keep the adapter API intuitive and consistent with existing `"O"` aligners:
- store `rotation` as the operator mapping domain **i → j**, so that
  `apply_transform(relative_transform(fit, "i","j"), X_i)` yields aligned data
  in domain j’s coordinate system.

Implementation detail:
- In Alvarez-2019 derivations, `P_col` is often written as a column-acting map.
  In row convention: `x_col ↦ P_col x_col` corresponds to `x_row ↦ x_row %*% t(P_col)`.
- Therefore, if you compute `P_col = U V^T` from SVD of `M = X^T π Y`,
  then the corresponding row-operator mapping `Y → X` is `t(P_col)`.
- For an adapter that returns `i → j`, return the transpose/inverse accordingly.

### 3.3 Backends to use (no new dependencies)

Cost computation:
- Use `compute_squared_distances_cpp_fallback(X, YP)` where `YP = Y %*% op`.

Sinkhorn step:
- Use `sinkhorn_unified(cost, a, b, epsilon, max_iter, tol, stabilized = TRUE)`.

Orthogonal projection:
- Use `svd(M)` and set `R = U %*% t(V)`.
- Optional `force_SO=TRUE`: flip last singular vector if `det(R) < 0` (same as
  `rotation_sync()` and `spectral_mnn_align_control(force_SO=TRUE)` patterns).

### 3.4 Proposed API surface (adapter)

Constructor:
- `invar_ot_aligner()` → `new_aligner("invar_ot", group="O", supports_multi=FALSE)`

Pair fit method (signature sketch):

```r
fit_pair.invar_ot_aligner <- function(
  algo, X_i, X_j, links = NULL,
  a = NULL, b = NULL,
  epsilon0 = 1, epsilon_min = 0.05, decay = 0.95,
  max_iter = 50, tol = 1e-6,
  sinkhorn_max_iter = 1000, sinkhorn_tol = 1e-9,
  stabilized = TRUE,
  init = c("identity", "random_orthogonal"),
  force_SO = FALSE,
  store_transport = TRUE,
  verbose = FALSE,
  ...
)
```

Notes:
- `links` can be ignored initially (Alvarez is unsupervised), but later could be
  used for warm-starts / priors.
- `a`, `b` default to uniform marginals (consistent with other OT code).
- `store_transport=FALSE` avoids huge fit objects when `n*m` is large; in that
  case store only summary diagnostics + final rotation.

Return object (class `invar_ot_pair_fit`) should likely contain:
- `rotation` (the `d×d` operator mapping i → j, row-acting; suitable for `apply_transform`)
- `transport` (the `n×m` plan; optional)
- `objective_trace` (vector; optional)
- `loss` (final `⟨π, C⟩` under final `epsilon_min` or final epsilon)
- `converged` (logical) + `iterations` (int)
- `epsilon_schedule` (vector of epsilons used)
- metadata: `n1`, `n2`, `d`

`relative_transform()` method:
- return `new_align_transform("O", rotation_or_transpose, from, to, k=d)`.

`pair_loss()` method:
- return `loss` (scalar), so `align_many()` can weight edges via `1/(1+loss)`.

`latent_dim()` method:
- return `d`.

### 3.5 Optional: exact OT backend (small problems)

For very small `n×m`, it may be useful to allow:
- `solver = c("sinkhorn","network_simplex")`
- if `network_simplex`, use `network_simplex_ot_cpp(cost, a, b, eps=1e-9)`

But note: Alvarez-2019’s optimization relies on entropic smoothing/annealing;
removing entropy can increase non-convexity and may worsen convergence without
careful initialization.

---

## 4) Testing plan consistent with the existing test suite

Add a new test file:
- `tests/testthat/test-invar-ot-aligner.R`

Suggested tests:
1) **Recovers known rotation on synthetic data**
   - sample `X` (n×d), generate random orthogonal `R_true`, set `Y = X %*% R_true`
   - run fit; check `||rotation - R_true||_F` is small (up to sign/perm symmetries),
     or check `rotation %*% t(R_true)` is close to identity.

2) **Transport plan has correct marginals**
   - if `store_transport=TRUE`, verify row/col sums match `a`, `b` (tolerance).

3) **Adapter works with align_many**
   - create 3 domains from one base `X` rotated by different random orthogonals
   - fit pairwise with `align_many(domains, invar_ot_aligner(), graph="complete")`
   - verify returned transforms are orthogonal and embeddings are aligned (e.g.,
     pairwise Frobenius distances reduced).

Use `skip_if_not_installed()` only when needed; for this method, all backends
are already compiled in-package, so skips should be minimal.

---

## 5) Naming recommendation

Two reasonable naming options:

1) **`invar_ot_*`**
   - Pros: matches paper framing (“invariances”).
   - Cons: could be confused with existing “translation-invariant” UOT code
     (`uot_ti_sinkhorn_kl`).

2) **`ot_procrustes_*`** (or `wasserstein_procrustes_*`)
   - Pros: immediately communicates “OT + orthogonal map”.
   - Cons: doesn’t emphasize the general Schatten framework, but that can be
     added later.

If the first implementation is orthogonal-only, `ot_procrustes_aligner()` is
arguably the clearest name; docs can cite Alvarez-Melis et al. (2019) and note
the generalization target.

