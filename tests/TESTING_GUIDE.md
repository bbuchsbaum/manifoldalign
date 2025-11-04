# Alignment Testing Guide

This document records the shared testing conventions for manifold alignment methods in this repository. It covers the canonical datasets, invariant assertions, method-specific checks, and contributor workflows that keep the alignment stack consistent.

---

## 1. Alignment Families & Capabilities

| Method | Type | Domains | Kernel | OT | Anchors | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| `gpca_align` | spectral linear | multiset | optional (pre-stage) | ✗ | optional | generalised PCA with within/between graphs |
| `grasp` / `grasp_multiset` | spectral + graph | multiset | ✗ | ✗ | optional | handles heterogeneous block metrics |
| `kema` | kernel spectral | multiset | ✓ | ✗ | optional | exact and REKEMA solvers |
| `cone_align` / `cone_align_multiple` | conic optimisation | multiset | ✗ | ✗ | optional | robust to label imbalance |
| `lowrank_align` | linear manifold | multiset | ✗ | ✗ | optional | scalable low-rank approximation |
| `parrot` | optimal transport | 2 domains | ✗ | ✓ | required | network alignment with entropic OT |
| `pseudolabel` | alignment + labelling | multiset | ✗ | ✗ | optional | leverages alignment for labels |
| `generalized_procrustes` | orthogonal | multiset | ✗ | ✗ | optional | task-based alignment with partial observations |
| `gromov_wasserstein` | optimal transport | multiset | ✗ | ✓ | optional | structural OT without features |
| `fpgw` | optimal transport | multiset | ✗ | ✓ | optional | fused/partial OT mixing feature + structure |
| `coupled_diagonalization` | spectral joint | multiset | ✗ | ✗ | optional | simultaneous multi-modal diagonalisation |
| `linear_sim_embed` | linear similarity | two set | ✗ | ✗ | optional | similarity-preserving linear maps |

Each method participates in the unified testing framework according to its capabilities (e.g. OT checks apply only to PARROT).

---

## 2. Canonical Test Datasets

Reusable fixtures live in `tests/testthat/helper-alignment.R` (to be extended) with these generators:

| Fixture | Description | Used for |
| --- | --- | --- |
| `spiral_two_domain(n, noise)` | Two interleaving spirals with rotation/scale | GPCA, KEMA (nonlinear separability) |
| `linear_maps(n, noise)` | Linear transforms of shared latent space | GPCA, cone_align, lowrank_align, linear_sim_embed |
| `anchor_partial_overlap()` | Two networks with partial anchors and missing labels | PARROT, cone_align, KEMA |
| `graph_cycle_mismatch()` | Graphs with cycle mismatches/edge perturbations | PARROT, GRASP |
| `high_dim_sparse()` | High-dimensional sparse similarities | GPCA PSD checks, grasp_multiset |
| `noisy_labels()` | Imbalanced/ noisy labels | cone_align balancing, pseudolabel |
| `singular_block()` | Duplicate/colinear rows causing rank deficiency | GPCA, KEMA (warning paths) |
| `large_ot_cost()` | Cost matrices spanning tau regimes | PARROT, Gromov-Wasserstein, FPGW |
| `task_partial_observation()` | Task matrices with missing entries | generalized_procrustes |
| `dual_graphs_structural()` | Pair of graphs with structural differences | gromov_wasserstein, fpgw |
| `multi_modal_views()` | Multiple modalities sharing latent space | coupled_diagonalization |

Datasets are tagged in a manifest (see Section 6) including size, modality, expected invariants, and any required args (anchors, kernel).

---

## 3. Core Invariants (All Methods)

Every aligner is expected to satisfy the following assertions where applicable:

1. **Structure & Metadata**
   - Returns an object inheriting `multiblock_biprojector` (or documented alternative).
   - Contains `v`, `s`, `sdev`, `block_indices`, `labels` (if labels supplied), and method-specific fields.

2. **Finite Outputs**
   - `s`, `v`, and `transport_plan` (if present) contain finite numeric values with no `NA`/`Inf`.

3. **Dimensional Consistency**
   - `nrow(s)` equals total samples post preprocessing.
   - Each `block_indices` span the correct feature ranges.

4. **PSD & Conditioning (spectral methods)**
   - Row-metric matrices are positive semi-definite; eigenvalues checked for non-negativity within tolerance.

5. **Stochastic Marginals (OT methods)**
   - For entropic OT, row/column sums of the transport plan approximate uniform marginals (`1/n`) within tolerance (default `5e-3`).

6. **Rank & Component Limits**
   - Requested components greater than available rank trigger informative warnings and reduce `ncomp` without failing tests.

7. **Solver Consistency**
   - Multiple solver paths (R vs C++, exact vs regression) produce numerically close outputs within scenario-specific tolerances.

Helper functions (`expect_alignment_ok()`, `expect_ot_plan_ok()`) should wrap these assertions to avoid duplication.

---

## 4. Method-Specific Checks

### 4.1 Spectral / Linear (GPCA, GRASP, Cone, Low-Rank)
- **Eigenvalue Ratios**: Compare dominant eigenvalues against baseline expectations (e.g. spiral dataset ratios ≈ paper values).
- **Subspace Alignment**: Principal angles between methods (exact vs regression) below threshold (e.g. sinθ ≤ 0.15).
- **Balancing**: For label balancing options, edge degree distributions before/after scaling validated.
- **Rank Deficiency**: Ensure warnings emitted when component request exceeds informative directions.

### 4.2 Kernel Methods (KEMA)
- **REKEMA Consistency**: Reduced sampling fractions converge toward full solution; Procrustes-aligned Frobenius error bounded.
- **Solver Retry**: Regression path triggers automatic exact retry on poor subspace distance; tests confirm both branches.
- **PSD Laplacians**: Within/between Laplacians remain PSD after adjustments.

### 4.3 Optimal Transport (PARROT)
- **Sinkhorn Equivalence**: R and C++ implementations match within tolerances across choices of `n`, `tau`.
- **Marginal Sums**: Transport plans satisfy uniform marginals within `5e-3`.
- **Cost Monotonicity**: Proximal loop reduces objective or keeps monotonic; history recorded.
- **Anchors**: Anchor priors integrate correctly; removing anchors should broaden transport support.

### 4.4 Multiset Alignment
- Tests should include: subset alignment accuracy, projector reconstruction, multi-block transformations, and sample weighting correctness.

### 4.5 Generalized Procrustes
- Verify task alignment handles partial observations: compare recovered task matrices against ground truth using Frobenius error.
- Ensure convergence metrics (objective reduction, iterations) remain within expected bounds.
- For orthogonal solutions, check `O` matrices are orthonormal and `det(O) = 1` up to tolerance.

### 4.6 Gromov-Wasserstein
- Validate transport plan marginals (uniform or supplied measures) within OT tolerance (`≤ 5e-3`).
- Compare GW objective values between R and C++ (if both paths exist).
- Use structured graph fixtures (e.g., `graph_cycle_mismatch`) to ensure structural similarities are preserved; compute Fréchet/GW distances for regression tracking.

### 4.7 Fused/Partial GW (FPGW)
- Confirm feature vs structure weighting behaves correctly: varying `rho` should shift mass between feature-only and structure-only alignments.
- Check transport plan mass constraint (partial transport) meets specification (mass ≤ `rho`).
- Compare R/C++ implementations if available.

### 4.8 Coupled Diagonalization
- Ensure coupled bases remain orthonormal and diagonal blocks match expected spectra.
- Verify convergence on Stiefel manifold (iterations < cap, gradient norms small).
- Use toy datasets where latent components are known to validate reconstruction error.

### 4.9 Linear Similarity Embedding
- Confirm similarity preservation: compare reconstructed similarity matrix against target (Frobenius error threshold).
- Ensure convergence status is success and embeddings are finite.
- For C++ path, assert numerical equivalence with R implementation.
---

## 5. Diagnostic Scenarios

| Scenario | Purpose | Expected Behaviour |
| --- | --- | --- |
| Zero-between similarities | GPCA should report "M_between has zero norm" or equivalent info | Method should not fail; within-dom alignment only |
| Negative similarities | Methods clip or regularise to maintain PSD | Message + corrected metric |
| OT marginal drift | Detect tolerance breaches, advise on tau scaling | Test ensures sums within threshold |
| Large k-NN `k` vs degree | Should warn but keep PSD | Graph degrees limited, controls override |
| Sparse anchors | PARROT should error with actionable message | Coverage ensures early detection |

---

## 6. Test Matrix & Manifest

Maintain `tests/testthat/ALIGNMENT_TEST_MATRIX.yml` with entries:
```yaml
fixtures:
  spiral_two_domain:
    generator: spiral_two_domain
    params: {n: 100, noise: 0.05}
    applies_to: [gpca_align, grasp, kema]
    invariants: [structure, finite, subspace, eigenvalues]
  anchor_partial_overlap:
    generator: anchor_partial_overlap
    applies_to: [parrot, cone_align, kema]
    invariants: [structure, finite, anchors, ot_marginals]
```

For each method provide required/optional flags. CI meta-tests should read this manifest to ensure coverage (future work).

---

## 7. Contributor Workflow

1. **Add Dataset**: Define generator in helper file; register in manifest with description + invariants; annotate the intended methods.
2. **Add Tests**: Use shared fixtures and helper expectations rather than bespoke data. Document new tests in this guide if they represent new invariants.
3. **Update Matrix**: Reflect changes (add rows/columns) to ensure the table mirrors actual coverage.
4. **Run Targeted Suites**: `devtools::test(filter = 'gpca_align')` etc; list recommended commands in each section.
5. **Document Diagnostics**: When addressing failures, note in Section 5 if new patterns emerge (e.g., tolerance updates, solver fallback warnings).

---

## 8. Next Steps

- [ ] Implement helper fixtures and expectations.
- [ ] Populate manifest with existing datasets.
- [ ] Migrate legacy tests to shared infrastructure.
- [ ] Add CI meta-test verifying matrix coverage.
- [ ] Periodically review tolerances (especially OT) as implementations evolve.

---

Keeping this guide current avoids drift between alignment methods and their test coverage, making it easier to add new methods or datasets with confidence.

