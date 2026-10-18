# manifoldalign 0.1.0.9000

### Major Changes

* Implemented Fused-Partial Gromov-Wasserstein (FPGW) distance for heterogeneous domain adaptation
* Added support for mass-constrained partial optimal transport with `rho` parameter

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
