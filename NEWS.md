# manifoldalign 0.1.0.9000

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
  The constraint now uses the degree matrix of the joint graph (within-domain
  degree plus incident correspondence weight), as in Wang & Mahadevan, so
  correspondence edges enter both sides of the generalized eigenproblem;
  this changes default embeddings. `ssma_align_control(degree_matrix =
  "within_domain")` reproduces the previous constraint, and fits record the
  choice in `degree_matrix`.
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
* **MMA (`mma_align_multiple()`)**: `assignment` maps each domain row to its
  reference row (it was inverted, and meaningless for non-square
  posteriors); EM posteriors and the log-likelihood include the full
  Gaussian normalisation; the M-step and the starting rotation were
  transposed, so EM results change; `converged` is reported correctly.
* **CONE-Align and GRASP**: score rows keep each domain's input order and the
  assignment is returned separately (`NA` when unmatched); `v` is a
  least-squares bridge to the aligned scores in CONE, CONE multiple, GRASP and
  MMA; the default `center()` preprocessing is fitted and stored, and
  preprocessors are concatenated by feature block, so out-of-sample projection
  works. GRASP's `lambda`, `max_iter` and `tol` now take effect, unequal node
  counts work, and `grasp_multiset()` handles four or more domains of unequal
  size. `solver = "auction"` now errors everywhere (it never worked).
* **GRASP graphs**: `grasp()`, `grasp_multiset()` and `grasp_aligner()` build
  Euclidean heat-kernel kNN graphs by default (`graph_metric = "euclidean"`)
  instead of graphs on row-standardised features. The old graphs were not
  invariant to rotations of the features and collapsed with two features, so
  GRASP matched a rotated copy of a domain at chance level; it now recovers
  rotated, permuted copies exactly. The bandwidth is the median k-th-neighbour
  distance under the graph's own metric (previously a random subsample of raw
  distances, which also consumed the caller's random numbers).
  `graph_metric = "correlation"` keeps row-standardised graphs and requires
  at least three features. GRASP results change; tournament rows for GRASP
  predate this.
* **Eigensolvers**: graph-alignment methods share one symmetric eigensolver
  that returns the requested end of the spectrum and checks convergence;
  coupled diagonalisation uses the smallest algebraic eigenvalues.
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
* **Out-of-sample preprocessing**: `oos_predict()` for fits that store one
  pre-processor per domain now applies the training statistics to new rows;
  it used to refit the pre-processor on the new rows (so, for example, a
  single new row was centred to zero). An unfitted stored pre-processor is an
  error.
* **multivarious**: pre-processors are fitted through `fit_transform()`, so
  fits no longer emit `prep()`/`init_transform()` deprecation warnings with
  multivarious >= 0.3.0; older versions fall back to the previous calls.
* **Dependencies**: PRIMME is optional (RSpectra fallback); `pryr` is no
  longer used or suggested; `src/Makevars` links LAPACK and BLAS; the vendored LEMON code
  no longer requires OpenMP and is credited in `inst/COPYRIGHTS`.

### Bug Fixes

* Fixed network simplex solver segfaults by improving error handling
* Withdrew FPGW's `lambda > 0` penalized route. It selected balanced transport
  vertices, so transported mass remained one and the parameter could not solve
  partial FPGW; its line-search term also disagreed with the stated objective.
  Historical penalized results are invalid as partial-transport evidence. Use
  the oracle-verified fixed-mass `rho` route instead.
* Fixed GlobalGeo spatial-anchor smoothness to use the declared `adjoin`
  dependency, count non-self neighbors correctly, preserve one-way kNN edges
  through union symmetrization, and expose a zero-penalty isolate policy.
* Corrected the Linear Similarity Embedding similarity-gradient sign and scale,
  the orthogonality-gradient scale in both backends, and the C++ L-BFGS-B
  stopping control. Nonzero optimizer statuses now use a verified
  exact-objective fallback or stop; they no longer return an apparently valid
  fit. Historical LSE fits carrying a nonzero convergence code are invalid.
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

### Deprecations and Warnings

* Reclassified GlobalGeo, CONE-Align, and CONE-Align Multiple as
  legacy/reference rather than recommended or experimental methods. GlobalGeo
  completed 8/21 intended frozen cases and its partial-overlap claim is
  withdrawn pending a disconnected-landmark policy. Pairwise CONE completed
  10/12; its two latent-dimension/Procrustes errors also gate the multiple-domain
  wrapper. These explicit failures are repair items, not silent wrong-result
  findings, so the methods are not deprecated solely for benchmark loss.
* `multiscale_manifold_align()` is deprecated and scheduled for removal in
  manifoldalign 1.0.0. It completed only 19 of 40 frozen validation runs and
  its scaling fixtures fail with a malformed solver formula. Existing code is
  retained temporarily for reproducibility, but new work should use
  `rapid_ma()` for supplied relations/positions/attributes, `ssma_align()` for
  semi-supervised feature manifolds, or `spectral_mnn_align()` for an
  unsupervised spectral baseline.
* FPGW's `lambda` parameter is retained only for compatibility and must be
  zero. Calls with `lambda > 0` now stop with a migration message directing
  users to the exact-mass `rho` formulation.
* `validate_out_of_sample_reconstruction()` is now a compatibility shim that
  warns and stops. Use `cv_alignment_rows()` for leakage-safe held-out alignment
  scoring until KEMA has a verified OOS reconstruction predictor and oracle.
* `kema_orig()` is deprecated in favor of the single public `kema()` entry
  point and is scheduled for removal in manifoldalign 1.0.0. `solver =
  "regression"` is withdrawn because it never selected a distinct regression
  implementation. Formerly ignored extension arguments (`use_laplacian`,
  `dweight`, `rweight`, `simfun`, `disfun`, and `centre_kernel`) now error
  instead of misrepresenting the fitted objective.
* `validate_kema_eigenvalues()` is withdrawn: it did not reproduce a documented
  paper fixture and solved the wrong generalized right-hand side. Historical
  paper-eigenvalue validation results from that function are invalid.

### Internal Changes

* The unused `partial_ot_tv()` augmentation now stops instead of exposing an
  objective that was not the package's documented FPGW penalty.
* Partial FPGW warm starts are rescaled into the exact-`rho` subcoupling
  feasible set before Frank-Wolfe iterations; infeasible internal starts now
  fail explicitly.
* Replaced the feature-cloud-only graph benchmark with attributed-graph
  protocol v2. Cases now carry supplied weighted relations, positions, mixed
  attributes, overlap/private-node truth, row identities, perturbation
  severity, train/dev/final roles, and content hashes. Twelve development cases
  and 42 zero-tuning final cases are frozen in installed manifests. The 30
  final accuracy/robustness cases include both partial-overlap and closed-world
  regimes so classic comparators are not judged only on an unsupported task.
* Added raw-attribute, position, topology-descriptor, anchor, row-identity,
  source-only and target-only classification, and chance controls. Benchmark rows now record
  typed outcomes and whether supplied topology was actually consumed; RAPID-MA
  is the supplied-graph comparator in the interactive dispatcher.
* Normalized package comparator adapters through the canonical method registry.
  Supplied-edge mutations now affect only RAPID-MA inputs; feature-derived and
  metric-measure methods are labeled accordingly, unsupported topology is
  rejected before execution, and prediction records preserve row identities
  plus anchor/non-anchor/private partitions. Frozen search-space hashes, equal
  tier wall-clock ceilings, convergence, timeout, and external peak-RSS status
  are carried in every result row. LSE is explicitly unsupported for
  cross-domain scoring, and LeMA remains an unavailable external reference
  until a reproducible adapter exists.
* Vectorized v2 edge canonicalization so large sparse scenario generation no
  longer performs one data-frame assembly per unique edge.
* Added governed real-data protocol 1.0 with checksum-pinned So2Sat-LCZ42,
  AUCS, and Sampson monastery tasks. Real network cases preserve observed
  direction, edge weights, natural partial overlap, unequal active-node sets,
  hidden scoring identities, and deterministic training anchors. Remote-
  sensing cases expose separate fit and scorer views so testing-east labels
  cannot enter fitting.
* Added resumable, license-gated acquisition, deterministic real-patch feature
  extraction, a small schema-only CI fixture, and a full frozen-split runner.
  Evidence rows distinguish unavailable or invalid datasets from algorithm
  failures and label CI accuracy as non-ranking diagnostics.
* Added the frozen package-wide alignment tournament: 588 typed rows over 42
  cases and 14 methods, with external runtime/peak-RSS measurement, resumable
  checkpoints, survivor-bias-safe summaries, paired effects, Pareto fronts,
  controls, and recursive checksums. RAPID-MA completed every case and is the
  only comparator that consumed supplied weighted relations, but its matching
  remains below the position-only control in every non-oracle regime. Its
  classification probe beats source-only transfer in every non-oracle regime;
  neither result is presented as a universal winner claim. See the checked
  [interpretation](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/INTERPRETATION.md)
  and [raw receipt](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/raw_results.csv).

### Documentation

* Updated FPGW documentation to distinguish metric-measure geometry from
  observed graph topology and to make `rho` the sole partial-transport API.

### Testing

* Added a dependency-independent FPGW oracle that exhaustively enumerates tiny
  weighted transport flows and checks compiled optimality, the independent
  four-index objective, exact mass and capacities, fusion endpoints,
  monotone objective traces, and target-permutation equivariance. The former
  balanced-mass `lambda` behavior is retained only as a rejection regression.
* Replaced KEMA completion-only and tautological solver tests with independent
  generalized-eigen oracles, hard failure-path checks, semi-supervised
  null-space coverage, and full/reduced/operator subspace agreement tests.
* Added graph-v2 generator-truth, common-random-severity, manifest-hash,
  no-final-tuning, negative-control, and observed-topology dispatch tests.
* Added real-data provenance, checksum, row/edge identity, weighted/directed
  topology, leakage, deterministic-control, and typed-failure tests. The
  enforced alignment matrix now names the governed real-data fixture and its
  CI-versus-release evidence scopes.
* Made `ALIGNMENT_TEST_MATRIX.yml` an enforced 22-method computational-test
  contract. PR, nightly, and release runners now report failures, warnings,
  skips, mandatory-oracle status, known gaps, and release availability by
  method and severity. Mandatory skips fail CI, stable promotion requires all
  computational families, and a dependency-independent coupled-diagonalization
  finite-difference oracle replaces optional-package-only evidence.
* Added tournament meta-tests for exact method-case coverage, denominator-correct
  missing predictions, semantic versus resource applicability, per-adapter tiny
  oracles, target-permutation invariance, typed errors, checkpoint resumption,
  and checksum-complete final evidence.
