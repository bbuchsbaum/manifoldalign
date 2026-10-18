# manifoldalign NEWS

## Version 0.1.0.9000 (Development)

### Major Changes

* Added `sumcor_align()`, a native multi-domain alignment for observations of
  shared, keyed events. It learns explicit per-domain linear encoders that
  maximise reproducible cross-domain agreement of event-level means,
  normalised by repeat-level covariance, with an optional event-feature
  metric (`beta = 0` is the exact correspondence-only model). Correspondence
  is carried only by event identifiers. Training scores are produced by the
  stored encoder, so `oos_predict()` reproduces them exactly. New generic
  `calibrate_domain()` adds a domain against the frozen reference without
  changing existing encoders. `sumcor_aligner()` plugs the method into
  `fit_many()` and `align_many()`.
* `aligner_capabilities()` now also reports encoder type, prediction and
  calibration support, auxiliary-input needs, and invertibility for aligners
  that declare them; the two-field result is unchanged for other aligners.
* `global_geo_align()` and the multi-set UOT functions no longer assume
  neuroimaging data. `feature_correspondence$coords` replaces `voxel_coords`,
  coordinates may have any number of dimensions (the same in every domain),
  and `anchor_method = "grid"` builds a grid in that dimension.
  `uot_build_cost()`, `uot_fit_pair()` and `multiset_uot_align()` take
  `lambda_coord` in place of `lambda_anat`; documentation now speaks of
  datasets and coordinates. The old names still work with a deprecation
  warning, and supplying both an old and a new name is an error. The
  low-level `uot_build_sparse_cost_knn_cpp()` argument is renamed to
  `lambda_coord` without an alias.
* Added `rapid_ma_feature_scores()`, feature-side (dual) coordinates for
  RAPID-MA fits.

* Implemented Fused-Partial Gromov-Wasserstein (FPGW) distance for heterogeneous domain adaptation
* Added support for mass-constrained partial optimal transport with `rho` parameter

### Correctness fixes that change results

These fix defects found in a package review; several change fitted output.

* **KEMA**: out-of-sample prediction now applies the fitted kernel encoder
  (`oos_predict()` and `project_block()` compute k(x, reference rows) alpha);
  previously it used a linear surrogate that was wrong for nonlinear kernels
  and for the reduced solver. Fits have class `kema_fit`; generic `project()`
  now errors and `v` is `NA` for nonlinear kernels. New `seed` argument (stored
  REKEMA landmarks) and `backend_control$max_dense_bytes` (errors before an
  oversized dense allocation). `validate_out_of_sample_reconstruction()` and
  `validate_kema_eigenvalues()` are now defunct.
* **SSMA**: the dense eigen fallback no longer discards the smallest
  eigenvalues; returning fewer than `ncomp` components warns and is recorded,
  and trivial null-space vectors are never substituted. Label
  correspondences no longer pair a random row when a label has one member.
* **Spectral MNN**: loadings are an explicit ridge regression bridge to the
  aligned scores (`bridge_ridge`); heat-kernel signatures use one time grid
  shared by all domains; the Procrustes step applies its translation.
* **PARROT**: the proximal loop now iterates on the full objective (it
  previously reduced to a single entropic transport, so anchors and
  consistency terms had no effect); R and C++ backends solve the same
  problem; the default `center()` preprocessing is applied; `solver =
  "exact"` errors; new `max_outer`; non-convergence warns instead of being
  rebalanced away.
* **Sinkhorn**: `stabilized = TRUE` is honoured for every `epsilon`
  (it returned zero plans for large costs); convergence uses the marginal
  residual. Compiled kernels no longer fall back silently to R code on
  errors.
* **FPGW**: the mass-constrained variant starts and stays feasible and returns
  mass `rho`; the TV-penalised variant uses its own oracle and correct sign
  (larger `lambda` transports more mass); arguments are validated and LP
  failures are errors. FPGW and Gromov-Wasserstein objects are no longer
  classed `multiblock_biprojector`.
* **Gromov-Wasserstein**: entropy scaling is consistent with `epsilon`, the
  inner Sinkhorn is stabilised, and `predict()` always uses the fitted
  coupling (the raw-coordinate nearest-neighbour fallback is removed).
* **Network simplex**: exact at 2^-50 mass resolution with separate cost and
  mass scaling (small masses were rounded to zero), validated against the
  original marginals.
* **Multi-set UOT**: learned template features use consistent mass scaling;
  duplicate sparse edges are merged; iteration counts are reported correctly;
  `max_dist` is optional.
* **`align_many()`**: synchronisation, refinement and diagnostics share one
  documented transform contract (the inverse convention had let diagnostics
  certify wrong answers); permutation synchronisation works in all three
  modes; the permutation hub defaults to the largest domain; `links` is a
  formal argument; unsupported groups and consensus modes error.
  `relative_transform()` for GRASP pair fits returns the transpose of its
  previous operator.
* **`linear_sim_embed()`**: convergence no longer depends on `verbose` (the
  default run used to stop after one update); gradients are exact in R and
  C++; non-convergence is reported; `tol` reaches L-BFGS-B as `pgtol`.
* **Generalized Procrustes**: new `translation = "shared_tasks"` (default)
  centres partial-task designs on shared tasks; a stationarity
  `certificate` is returned; `final_diff` reports the last step.
* **`preproc = NULL`** now works in `gpca_align()`, `generalized_procrustes()`,
  `lowrank_align()` and `global_geo_align()`; `lowrank_align()` no longer builds
  a dense cross-domain mask.
* **RAPID-MA**: nested control lists are revalidated and unknown fields
  rejected; pattern relation matrices and sparse inputs are accepted without
  densifying; orthogonal maps stay orthogonal; the diffusion directedness rule
  changed, which alters default embeddings; outputs keep row names. Benchmark
  and validation runners no longer leak held-out labels and use shared metric
  definitions; the evidence documents carry notes on which numbers need
  regeneration.
* **`cv_alignment_rows()`**: targets are scored with the same projection as
  queries; folds are matched by name or position.
* **Pseudolabelling**: greedy selection is farthest-first; the large-n path
  no longer overflows or mislabels vertices.
* **`assess_ssma()`**: missing predictions count in the denominator and
  coverage and solver availability are reported.
* **Random numbers**: package functions no longer reset the caller's RNG
  stream; seeded steps restore it afterwards.
* **Dependencies**: PRIMME is optional (RSpectra fallback); `pryr` is no
  longer used; `src/Makevars` links LAPACK and BLAS; the vendored LEMON code
  no longer requires OpenMP and is credited in `inst/COPYRIGHTS`.

### Bug Fixes

* Fixed network simplex solver segfaults by improving error handling
* Completed the package-wide migration from the former `neighborweights`
  package name to the CRAN `adjoin` package. GlobalGeo anchor smoothness now
  uses the current spatial-adjacency contract, preserves one-way kNN edges by
  union symmetrization, and gives isolates an explicit zero-penalty policy.
* Withdrew `validate_out_of_sample_reconstruction()`, which previously returned
  `expected_error` plus random noise without fitting KEMA or evaluating held-out
  predictions. Historical values from this routine are invalid as accuracy
  evidence, and it is no longer included in `run_kema_validation_suite()`.
* Corrected KEMA component selection to solve on the numerical image of the
  kernel blocks and explicitly deflate graph-component null modes. KEMA now
  uses eigencore's certified full or constrained generalized LOBPCG solver,
  requires a strictly positive metric regularizer, rejects materially negative
  modes, and independently checks residuals and B-orthogonality at `1e-6`.
  Auto mode tries alternate backends and then stops if no backend passes;
  warning-only failed fits are no longer returned. Historical KEMA results with
  `fidelity$passed == FALSE`, or whose selected eigenvalues do not exceed their
  numerical-zero threshold, are invalid.
* Corrected TV penalty implementation to match paper's mathematical formulation (equation 11)

### Deprecations and Warnings

* `kema_orig()` is deprecated in favor of the single public `kema()` entry
  point and is scheduled for removal in manifoldalign 1.0.0. `solver =
  "regression"` is withdrawn because it never selected a distinct regression
  implementation. Formerly ignored extension arguments now error instead of
  misrepresenting the fitted objective.
* `validate_kema_eigenvalues()` is withdrawn because it did not reproduce a
  documented paper fixture and solved the wrong generalized right-hand side.
* **TV penalty (`lambda` parameter) is now marked as experimental**. The TV penalty `λ(|μ|² + |ν|² - 2|γ|²)` is a concave quadratic that biases optimization towards either zero or full mass, but cannot guarantee gradual sparsity. For predictable partial transport, use the mass-constrained formulation (`rho` parameter) instead.

### Internal Changes

* Removed incorrect augmentation method for TV penalty
* TV penalty now handled correctly in Frank-Wolfe step size computation
* Added runtime warning when `lambda > 0` to alert users about experimental status

### Documentation

* Updated FPGW documentation to clarify TV penalty limitations
* Examples now emphasize `rho` parameter for partial transport
* Added comprehensive technical documentation about TV penalty behavior

### Testing

* Modified TV penalty tests to check for experimental warnings rather than sparsity behavior
* All tests now pass with corrected expectations