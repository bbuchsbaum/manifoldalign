# Alignment Method Test Assessment and Plan

This document summarizes our current confidence in each alignment method based on the existing test suite, and proposes concrete steps to bring every method to a **5/5** confidence level.

Confidence is about **correctness + robustness as demonstrated by tests**, not theoretical optimality.

- **5** – Very well validated by unit tests, edge cases, and representative benchmarks; failures would be surprising.
- **4** – Strong coverage on core behavior; some regimes or benchmarks are missing or only loosely checked.
- **3** – Basic structural and consistency tests exist; limited direct evaluation of alignment quality.

Methods covered here are those in the testing guide and matrix:

- `kema`
- `generalized_procrustes`
- `grasp`, `grasp_multiset`
- `cone_align`, `cone_align_multiple`
- `parrot`
- `gromov_wasserstein`
- `fpgw`
- `coupled_diagonalization`
- `gpca_align`
- `linear_sim_embed`
- `lowrank_align`

---

## Summary Table

| Method                             | Current Confidence | Notes |
|------------------------------------|--------------------|-------|
| KEMA                               | 4/5                | Strong numerical suite; OOS + paper-value checks not fully enforced |
| Generalized Procrustes             | 5/5                | Tight recovery tests, partial tasks, noise, diagnostics |
| GRASP                              | 4/5                | Identity/isomorphic graph tests + stress tests; some skips/fallbacks |
| grasp_multiset                     | 4/5                | Good structure + alignment-quality checks; no analytic baselines |
| CONE-Align                         | 4/5                | Identity, error reduction, permutations; less tied to canonical fixtures |
| CONE-Align Multiple                | 4/5                | Multi-domain alignment quality good; similar gaps as single CONE |
| PARROT                             | 4/5                | Comprehensive OT and anchor tests; thresholds intentionally loose |
| Gromov–Wasserstein                 | 4/5                | Strong OT structure tests; no analytic GW baselines |
| FPGW                               | 4/5                | Good coverage of variants; limit regimes not benchmarked vs references |
| Coupled Diagonalization            | 4/5                | Coupled/decoupled regimes + gradients; non-convex convergence accepted |
| GPCA Align                         | 5/5                | Extensive structural, stability, and cross-domain tests |
| Linear Similarity Embedding (LSE)  | 5/5                | Latent recovery, orthogonality, robust guards for edge cases |
| lowrank_align                      | 3/5                | Backend consistency only; little direct alignment-quality testing |

---

## Per-Method Assessment and Plan

### KEMA

**Current confidence: 4/5**

Evidence:
- Comprehensive numerical validation in `tests/testthat/test-kema-numerical-validation.R` and `test-kema.R`:
  - Synthetic two-domain spirals (paper-style) with eigenvalue ratio checks.
  - Exact vs regression vs REKEMA comparisons (subspace distance, variance ratios).
  - Auto-tuning / fail-soft logic (`choose_sigma`, solver retries).
- Edge-case handling: tiny datasets, rank-deficient cases, and benchmark gating.

Gaps:
- Out-of-sample reconstruction tests are approximate, not strict quantitative guarantees.
- Eigenvalue-ratio tests currently mostly log and sanity-check; they don’t enforce paper-like target ranges.
- KEMA is not yet a first-class citizen in the toy manifold benchmark suite (strict Procrustes-based evaluation).

**Plan to reach 5/5**

- [ ] Implement a proper `predict()` / OOS pipeline for KEMA:
  - Train/test split on the spiral and toy manifold datasets.
  - Procrustes-based reconstruction error and class-separation metrics on held-out data.
  - Hard thresholds (e.g., RMSE and correlation vs baseline).
- [ ] Tighten eigenvalue ratio tests:
  - For fixed synthetic spirals, assert eigenvalue ratios fall within tolerance bands derived from the paper or current stable behavior (store expected ranges in tests).
- [ ] Integrate KEMA into `test-toy-manifold-benchmarks.R`:
  - Evaluate KEMA on `isometric_curve` dataset with strict Procrustes metrics.
  - Add negative tests on `hard_nonisometric` to ensure it fails gracefully or is clearly worse than linear methods.
- [ ] Ensure full KEMA validation suite runs in CI:
  - Configure at least one CI job with `options(manifoldalign.run_benchmarks = TRUE)` and required deps (PRIMME, kernlab, multivarious).

---

### Generalized Procrustes

**Current confidence: 5/5**

Evidence:
- `tests/testthat/test-generalized_procrustes.R`:
  - Exact recovery of templates up to global rotation, including partial observation patterns.
  - Noise tolerance tests with explicit relative-error bounds.
  - Tightness certificate output and convergence checks.
  - Validation of error conditions (duplicate labels, single subject, mismatched dimensions, early stopping).

Gaps:
- Limited large-scale stress testing; performance and convergence at larger L and n not monitored.

**Plan to maintain 5/5**

- [ ] Add benchmark-sized test (gated by an option) with larger `L` and `n`:
  - Assert convergence and runtime within a documented budget.
  - Track typical iteration counts to catch performance regressions.
- [ ] Keep tightness-certificate tests in a mandatory CI job to catch numeric drifts.

---

### GRASP

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-grasp.R`:
  - Identity recovery on identical graphs with high accuracy requirements.
  - Isomorphic graphs with known permutations; compares accuracy against random baselines.
  - Multi-scenario, multi-parameter stress tests (tiny/small/medium/sparse/dense graphs).
  - Mathematical correctness: rotation orthogonality, descriptor normalization.

Gaps:
- Some GRASP tests use `tryCatch` + `skip` on error, which can mask regressions.
- Accuracy thresholds are relative to random baselines but not tied to canonical fixtures or tight targets.
- No explicit integration into the canonical fixtures from `ALIGNMENT_TEST_MATRIX.yml`.

**Plan to reach 5/5**

- [ ] Refactor GRASP tests to avoid `skip()` on computation errors:
  - Reduce problem sizes or tune parameters so tests are reliable enough to fail hard on regressions.
- [ ] Add tests using canonical fixtures:
  - Use `graph_cycle_mismatch` for GRASP alignment tests with explicit accuracy thresholds.
  - Assert assignment accuracy and top-k accuracy exceed fixed thresholds across multiple seeds.
- [ ] Add multi-seed stability checks:
  - Run GRASP on a fixed fixture over several seeds, assert median accuracy and variance stay within bounds.

---

### grasp_multiset

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-grasp_multiset.R`:
  - Hyperdesign and list-of-matrices interfaces validated.
  - Permutations are bijections; rotations are orthogonal.
  - Anchor choices and consensus-vs-direct mapping behavior tested.
  - Convergence behavior (max-iteration warnings, reasonable convergence).
  - Alignment quality vs permuted ground truth graphs.

Gaps:
- No connection yet to canonical fixtures (`high_dim_sparse`, etc.) beyond structural checks.
- Accuracy thresholds are formulated relatively, without analytic or cross-method baselines.

**Plan to reach 5/5**

- [ ] Add tests mapping `grasp_multiset` to fixtures in `ALIGNMENT_TEST_MATRIX.yml`:
  - E.g., `high_dim_sparse` for sparse graph metrics, with structural invariants and minimal alignment-quality expectations.
- [ ] Introduce cross-method benchmarks:
  - Compare `grasp_multiset` against GRASP on small multiset problems, ensuring comparable or better assignment accuracy.
- [ ] Add variance/stability tests across seeds on one multiset fixture.

---

### CONE-Align

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-cone_align.R`:
  - Basic return-structure checks; scores normalized, permutations valid.
  - Identity and permuted-graph recovery via `cone_align_iterate`.
  - Procrustes solver correctness (`solve_procrustes_cone`).
  - Parameter validation and solver-choice coverage.

Gaps:
- Tests are mostly internal/structural; limited use of canonical fixtures (`linear_maps`, `noisy_labels`, `anchor_partial_overlap`).
- Alignment quality is assessed via error reduction rather than explicit permutation recovery on fixtures with known ground truth.

**Plan to reach 5/5**

- [ ] Tie CONE tests to canonical fixtures:
  - Use `linear_maps` for linear alignment; assert Procrustes error vs latent is below a fixed threshold.
  - Use `noisy_labels` to check robustness to label imbalance (e.g., alignment-quality degradation within a controlled range).
  - Use `anchor_partial_overlap` where applicable to test alignment with partial anchors.
- [ ] Add explicit permutation-recovery metrics:
  - For synthetic permuted graphs, assert a minimum fraction of nodes are correctly matched.
- [ ] Integrate CONE into toy benchmark suite:
  - Evaluate on `linear_affine` and possibly `isometric_curve`, comparing Procrustes error and class separation to linear baselines.

---

### CONE-Align Multiple

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-cone_align_multiple.R`:
  - Multi-domain alignment to a reference space.
  - Rotations and permutations validated against pairwise `.align_two_embeddings`.
  - Node-count validations and error paths.

Gaps:
- Similar to single CONE: limited use of canonical multi-domain fixtures; no direct comparison to other multi-domain aligners (e.g., `mma_align_multiple`, `gpca_align`).

**Plan to reach 5/5**

- [ ] Add tests using multi-view fixtures (e.g., variations of `linear_maps` and multi-view synthetic datasets):
  - Check per-domain Procrustes error w.r.t. a reference domain and enforce thresholds.
- [ ] Compare to pairwise + consensus baselines:
  - For small 3-domain problems, compare CONE-multiple to repeated pairwise CONE or GPCA alignment and assert comparable quality.
- [ ] Integrate into toy benchmark suite as a multi-view graph alignment baseline.

---

### PARROT

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-parrot.R`:
  - Basic structural checks for `parrot` objects: `s`, `v`, `transport_plan`, `alignment_matrix`, anchors.
  - Parameter validation for all key arguments (sigma, lambda, tau, solver, etc.).
  - Anchor-constraint tests: higher mass on anchor pairs, robust to noise.
  - Transport-plan properties across tau values: marginals, non-negativity, entropy behavior.
  - Performance tests on easy/medium/hard aligned networks with relaxed accuracy thresholds.
  - Edge cases: very small networks, single anchor, high noise, different preprocessing options.

Gaps:
- Accuracy thresholds are intentionally conservative; they don’t reflect the desired “best-effort” performance.
- No cross-validation against a reference OT library or alternate implementation.
- No explicit evaluation on canonical fixtures (`anchor_partial_overlap`, `graph_cycle_mismatch`, `large_ot_cost`) for alignment accuracy, only for structural properties.

**Plan to reach 5/5**

- [ ] Strengthen performance thresholds:
  - On “easy” aligned networks, require significantly higher accuracy (e.g., >0.5 or >0.7) as implementation stabilizes.
  - Codify expected accuracy ranges per scenario based on current best results.
- [ ] Use canonical fixtures:
  - For `anchor_partial_overlap`, measure assignment or transport-based accuracy of anchors and non-anchors separately.
  - For `graph_cycle_mismatch` and `large_ot_cost`, assert that transport plans respect structural invariants (e.g., higher mass on consistent cycles, stable marginals across tau).
- [ ] Cross-check against a reference OT solver on tiny datasets:
  - Compare transport plans and objective values to an external implementation within tolerance.
- [ ] Add more explicit tests on how `lambda`, `lambda_p`, and `alpha` affect alignment quality and entropy.

---

### Gromov–Wasserstein

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-gromov_wasserstein.R`:
  - Basic hyperdesign input with different feature dimensions.
  - Doubly stochastic transport plans (within tolerance).
  - Identical-domain tests with small distances and diagonal-concentrated mass.
  - Different sample sizes and correct marginals.
  - Structural similarity between circle and spiral (valid transport and reasonable distances).
  - Input validation, convergence behavior, print method, and out-of-sample weights behaving as probabilities.

Gaps:
- No analytic GW baselines: we don’t test against known-optimal couplings on tiny problems.
- No multi-domain consistency checks (e.g., triangle-like inequalities on distances).

**Plan to reach 5/5**

- [ ] Construct tiny analytic test cases:
  - 1D grids or simple point sets where the optimal GW coupling is effectively diagonal/permutation-like.
  - Assert that the learned transport places ≥ some fraction of mass on the known matches, and that distances match expected ranges.
- [ ] Add multi-domain consistency tests:
  - For a simple three-domain setup, check that pairwise GW distances behave consistently (monotonicity, approximate triangle inequality).
- [ ] Add property-based tests:
  - Random small problems where we check marginals, convergence flags, and monotone behavior as `epsilon` changes.

---

### FPGW (Fused-Partial GW)

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-fpgw.R`:
  - Base hyperdesign functionality; transport-plan and distance matrix structure.
  - Mass-constrained variant (`rho`) and TV-penalized variant (`lambda`) with correct mass and metadata.
  - Mutual exclusivity of `rho` and `lambda`.
  - Behavior across `omega1` regimes (structural vs feature-dominant).
  - Convergence tests and print method.
  - Toy problem where transport mass is concentrated near the diagonal (expected structure).
  - Out-of-sample `predict` and `transform` tests with small shifts.

Gaps:
- No explicit connection to pure GW or pure OT baselines for edge-case parameter settings.
- No formal tests that FPGW’s distances respond smoothly and monotonically to `omega1` and `rho`.

**Plan to reach 5/5**

- [ ] Add limit-regime tests:
  - `omega1 → 0`: distances and transports comparable to a GW-only implementation.
  - `omega1 → 1`: behavior comparable to feature-only OT on the same distances.
  - `rho → 1` vs `rho < 1`: clear differences in total transported mass and alignment behavior.
- [ ] Add tiny analytic cases:
  - Compare FPGW plans and costs against a reference FGW implementation for small problems.
- [ ] Expand out-of-sample tests:
  - More challenging shifts and multiple `k` values, with hard bounds on RMSE and probability-vector behavior.

---

### Coupled Diagonalization

**Current confidence: 4/5**

Evidence:
- `tests/testthat/test-coupled_diagonalization.R` + `test-coupled_diag_mathematical.R`:
  - Basic hyperdesign input, different sample sizes, custom correspondence matrices.
  - Input validation and warnings.
  - Convergence behavior and diagnostic recording (cost history, eigenvalues, stiefel factors).
  - Recovery of known coupled latent structure with high coupling.
  - Decoupled vs high-coupling regimes clearly distinguished.
  - Partial correspondences, with per-sample embeddings aligned for corresponding samples.
  - Gradient checks via directional derivatives on the Stiefel manifold.

Gaps:
- Optimization is non-convex; some tests explicitly allow non-convergence as long as output is “reasonable”.
- Limited multi-run reproducibility checks (different initializations on the same problem).

**Plan to reach 5/5**

- [ ] Add multi-run reproducibility tests:
  - Run coupled diagonalization multiple times on the same synthetic coupled dataset with different seeds and measure:
    - Variation in final cost.
    - Alignment metrics (basis alignment, latent recovery).
  - Enforce that variation stays within acceptable bounds.
- [ ] Expand partial-correspondence tests:
  - Multiple correspondence densities and patterns, with explicit expectations on how alignment degrades.
- [ ] Add benchmark-sized tests (gated) to monitor runtime and convergence behavior for larger problems.

---

### GPCA Align

**Current confidence: 5/5**

Evidence:
- `tests/testthat/test-gpca_align.R`:
  - Basic structure and parameter handling (u, lambda, simfun, control).
  - Cross-domain alignment quality with centroid distances and correlations.
  - Edge cases: small n, no shared labels, single domain, sparse similarities, high-dimensional data.
  - Perfect correlation scenarios and custom preprocessing.
  - PSD correction behavior.

Gaps:
- Limited explicit involvement in the toy benchmark suite as a baseline.

**Plan to maintain 5/5**

- [ ] Ensure `genpca`/`PRIMME` are installed and GPCA tests are exercised in CI (no permanent skips).
- [ ] Add GPCA to toy benchmark suite:
  - Use it as the primary linear baseline on `linear_affine` and `isometric_curve` datasets with strict Procrustes metrics.

---

### Linear Similarity Embedding (LSE)

**Current confidence: 5/5**

Evidence:
- `tests/testthat/test-lse-recovery.R`, `test-lse-orthogonality.R`, `test-lse-robustness.R`:
  - Latent subspace recovery vs PCA on synthetic data, with Procrustes alignment and RMSE comparisons.
  - Near-orthogonality of learned weights under strong regularization.
  - Robust error/warning paths (zero masks, poorly scaled targets, predict dimension mismatch).

Gaps:
- Less coverage in extremely high-dimensional or p≫n regimes.
- Not yet integrated as a baseline in toy benchmarks.

**Plan to maintain 5/5**

- [ ] Add high-dimensional and p≫n tests:
  - Confirm convergence and orthogonality properties hold, and scores remain finite.
- [ ] Integrate LSE into toy benchmark suite:
  - Use as a linear baseline on `linear_affine` and compare to GPCA on alignment error.

---

### lowrank_align

**Current confidence: 3/5**

Evidence:
- `tests/testthat/test-lowrank_align.R`:
  - `createSimFun` preserves sparsity and symmetry of similarity matrices.
  - Label handling and NA behavior.
  - Agreement between explicit and operator solvers for `lowrank_align` (subspace overlap).

Gaps:
- No direct tests of alignment quality on canonical fixtures or synthetic datasets with known latent structure.
- No stress tests for noisy labels, partial correspondences, or high dimensionality.
- Not present in toy manifold benchmarks.

**Plan to reach 5/5**

- [ ] Add synthetic linear alignment tests:
  - Use a `linear_maps`-style generator with known latent space and labels.
  - Compare `lowrank_align` embeddings to the latent space via Procrustes, and assert error vs PCA/GPCA baselines.
- [ ] Use canonical fixtures:
  - Integrate `lowrank_align` with the `linear_maps` fixture from `ALIGNMENT_TEST_MATRIX.yml`, exercising `similarity_preservation` invariants.
- [ ] Add robustness tests:
  - Noisy/imbalanced labels, missing labels, and rank-deficient similarity matrices with appropriate warnings and stable outputs.
- [ ] Integrate into toy benchmark suite:
  - Evaluate on `linear_affine` as an additional linear method, with clear thresholds on alignment error.

---

## Cross-Cutting Plan

Beyond per-method work, several cross-cutting steps will raise overall confidence:

- [ ] **Toy benchmark integration**  
  Ensure each core method appears in `test-toy-manifold-benchmarks.R` (or an extended benchmark file) with strict Procrustes-based metrics:
  - Linear methods (GPCA, LSE, lowrank_align) on `linear_affine`.
  - Non-linear / manifold methods (KEMA, coupled_diagonalization, OT-based methods) on `isometric_curve`.
  - Failure-mode tests on `hard_nonisometric`.

- [ ] **CI configuration**  
  - At least one job with `options(manifoldalign.run_benchmarks = TRUE)` and required heavy dependencies installed.
  - Make sure benchmark skips are driven by explicit options, not silent missing packages.

- [ ] **Documentation alignment**  
  - Keep `tests/TESTING_GUIDE.md` and `ALIGNMENT_TEST_MATRIX.yml` updated as new fixtures/tests are added.
  - Document any changes in thresholds and rationale (especially when tightening performance expectations).

- [ ] **Stability across seeds**  
  - For stochastic methods (OT, GRASP, CONE, coupled_diagonalization), add multi-seed tests on small problems to monitor variance in key metrics (accuracy, Procrustes error, entropy).

