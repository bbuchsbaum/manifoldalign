# Alignment method sunset review

Status: accepted package review, 2026-07-22.

This review applies the executable lifecycle policy to all 22 canonical method
records. Correctness, completion, comparative accuracy, API duplication, and
task applicability are judged separately. The primary comparative evidence is
the frozen [tournament receipt](evidence/alignment-tournament-v2/raw_results.csv)
and its checked [interpretation](evidence/alignment-tournament-v2/INTERPRETATION.md).
Focused oracle paths are named below. No failed or missing run is removed from a
denominator.

## Decisions

| Canonical method | Decision | Evidence and completion | Capability retained or action required |
|---|---|---|---|
| RAPID-MA | Keep `experimental` | Tournament 42/42; invariant suites in [`test-rapid-ma.R`](../tests/testthat/test-rapid-ma.R), [`test-rapid-ma-uot.R`](../tests/testthat/test-rapid-ma-uot.R), and [`test-rapid-ma-relations.R`](../tests/testthat/test-rapid-ma-relations.R) | Only method here with native supplied sparse relations, weights, positions, mixed attributes, unequal sizes, and partial overlap. Matching loses to position-only control; blocked-inductive and million-node gates remain failed. |
| Token-OT Graph Align | Keep `experimental` | Tournament 18/18 applicable; sparse/dense oracle in [`test-token-ot-graph-align-contracts.R`](../tests/testthat/test-token-ot-graph-align-contracts.R) | Sparse candidate OT on feature-derived graphs. Low correspondence accuracy is a performance finding, not incoherence. |
| Spectral MNN Align | Keep `experimental` | Tournament 18/18 applicable; best exact-oracle matching at 0.750; invariants in [`test-spectral_mnn_align.R`](../tests/testthat/test-spectral_mnn_align.R) | Fast unsupervised/multiset spectral baseline; does not consume supplied topology. |
| SSMA | Keep `experimental` | Tournament 18/18 applicable; medium-closed classification OA 0.6635; operator/reduced checks in [`test-ssma-align.R`](../tests/testthat/test-ssma-align.R) | Semi-supervised multi-domain feature-manifold baseline. |
| GlobalGeo | Move to `legacy_reference` | Tournament 8/21: 6/6 oracle, 0/6 small closed, 2/9 small partial; geometry oracle in [`test-global_geo_align.R`](../tests/testthat/test-global_geo_align.R) | Preserve connected-geometry/OOS research use. Withdraw partial-overlap support. Restore experimental status only after a specified disconnected-component policy passes the frozen suite. |
| KEMA / REKEMA | Keep `experimental`; keep `kema_orig()` deprecated | Tournament 12/12 applicable; residual and B-orthogonality oracle in [`test-kema-numerical-validation.R`](../tests/testthat/test-kema-numerical-validation.R) | Repaired supervised kernel alignment remains coherent. Duplicate `kema_orig()` warns through 0.2.0 and is removed at 1.0.0 in favor of `kema()`. |
| GPCA Align | Keep `experimental` | Tournament 12/12 applicable; invariants in [`test-gpca_align.R`](../tests/testthat/test-gpca_align.R) | Supervised linear multiset alignment with native OOS transform. Independent exact oracle remains a promotion gap. |
| Low-rank Alignment | Keep `experimental` | Tournament 12/12 applicable; small-closed classification OA 0.7060; checks in [`test-lowrank_align.R`](../tests/testthat/test-lowrank_align.R) | Distinct low-rank supervised/OOS option; strong classification evidence prevents performance-only retirement. |
| Linear Similarity Embedding | Keep `experimental` as a utility, not a cross-domain aligner | Correctly `not_applicable` in 42/42 tournament cases; R/C++ objective and gradient oracle in [`test_linear_sim_embed_consistency.R`](../tests/testthat/test_linear_sim_embed_consistency.R) | Single-domain target-similarity estimator. Cross-domain claims remain forbidden. |
| Generalized Procrustes | Keep `experimental` | Not in tournament; template-recovery oracle in [`test-generalized_procrustes.R`](../tests/testthat/test-generalized_procrustes.R) | Unique orthogonal task alignment with partial task observation and OOS transforms. |
| GRASP | Keep `experimental` research baseline | Tournament 12/12 applicable; identity oracle in [`test-grasp.R`](../tests/testthat/test-grasp.R) | Functional-map graph matching is distinct; only two comparable tournament regimes and a unique mechanism, so the three-regime dominance retirement gate is not met. |
| GRASP Multiset | Keep `experimental` | Not in tournament; multi-domain checks in [`test-grasp_multiset.R`](../tests/testthat/test-grasp_multiset.R) | Multi-domain functional-map extension; independent oracle and comparative evidence remain gaps. |
| CONE-Align | Move to `legacy_reference` | Tournament 10/12; two equal-size closed cases stopped on mismatched Procrustes dimensions; rotation oracle in [`test-cone_align.R`](../tests/testthat/test-cone_align.R) | Preserve as a research baseline. It was second on oracle and small-closed matching, so performance deprecation is unjustified; repair latent-dimension reconciliation before restoration. |
| CONE-Align Multiple | Move to `legacy_reference` | Not directly in tournament; pairwise core failed 2/12 and wrapper checks are in [`test-cone_align_multiple.R`](../tests/testthat/test-cone_align_multiple.R) | Multi-domain capability is unique, but restoration inherits the pairwise dimension-reconciliation gate and must pass multi-domain regression cases. |
| PARROT | Keep `experimental` | Tournament 12/12; small-closed classification OA 0.8588, the regime leader; solver oracle in [`test-parrot-computational-contracts.R`](../tests/testthat/test-parrot-computational-contracts.R) | Anchor-guided dense OT/random-walk baseline has distinctive classification evidence despite weak correspondence accuracy. |
| Gromov-Wasserstein | Keep `experimental` | Tournament 12/12 applicable; transport invariants in [`test-gromov_wasserstein.R`](../tests/testthat/test-gromov_wasserstein.R) | Coherent balanced metric-measure OT. It neither consumes nor claims supplied graph topology; a dependency-independent exact oracle remains a promotion gap. |
| FPGW | Keep `experimental` | Tournament 21/21 applicable; exhaustive fixed-mass oracle in [`test-fpgw-oracle.R`](../tests/testthat/test-fpgw-oracle.R) | Coherent partial metric-measure OT. It is slow and inaccurate on this attributed-graph task, but exact partial mass and feature/geometry fusion are unique capabilities. |
| Coupled Diagonalization | Keep `experimental` | Not in tournament; finite-difference oracle in [`test-coupled-diagonalization-oracle.R`](../tests/testthat/test-coupled-diagonalization-oracle.R) | Joint spectral multi-domain alignment with native OOS transform; comparative accuracy is still missing. |
| Multiset Manifold Alignment | Keep `experimental` | Not in tournament; eigengate and multiset checks in [`test-mma-eigengate.R`](../tests/testthat/test-mma-eigengate.R) and [`test-mma-align-multiple.R`](../tests/testthat/test-mma-align-multiple.R) | Probabilistic registration for three or more feature-derived domains remains distinct. |
| Multiscale Manifold Align | Keep `deprecated` | Prior frozen validation 19/40; malformed solver failures begin at 550 total nodes; regression contract in [`test-multiscale-manifold-align-control.R`](../tests/testthat/test-multiscale-manifold-align-control.R) | Replace with RAPID-MA for supplied structure, SSMA for semi-supervised feature manifolds, or Spectral MNN for unsupervised spectral work. Warning spans 0.1.0--0.2.0; removal is 1.0.0; source/tests remain as the reproducibility archive. |
| Multiset UOT | Keep `experimental` | Not in tournament; dense/sparse oracle and stress checks in [`test-uot-ti-sinkhorn.R`](../tests/testthat/test-uot-ti-sinkhorn.R) and [`test-uot-stress.R`](../tests/testthat/test-uot-stress.R) | Unique sparse/dense unbalanced transport across pairwise or template supports. |
| OT-Procrustes | Keep `experimental` | Not in tournament; known-rotation oracle in [`test-ot-procrustes-aligner.R`](../tests/testthat/test-ot-procrustes-aligner.R) | Distinct alternating transport/orthogonal feature-space alignment with OOS transforms. |

## What is actually being retired

Only one canonical aligner is currently scheduled for removal:
`multiscale_manifold_align()`. Its implementation has an incoherent solver and
failed completion/scale contract, not merely poor accuracy. The duplicate
`kema_orig()` entry point is also scheduled for removal, but the repaired
canonical `kema()` method remains available.

GlobalGeo and CONE are demoted to legacy/reference, not deprecated. Their
measured failures are explicit errors rather than fabricated evidence or
silent wrong-result paths. Both also retain a plausible unique capability, and
CONE is not dominated in three applicable regimes. Runtime warnings and removal
versions would therefore overstate the evidence.

No other method meets the policy's performance-only retirement test of material
domination in at least three applicable regimes plus no unique supported
capability. Low accuracy alone is not a correctness finding.

## Repair and improvement gates

### RAPID-MA position-gated fusion

The next accuracy slice must be tuned only on development cases. It should add
an anchor-calibrated position-first candidate path, learn or gate view weights
using held-out training anchors, and guarantee a fallback when fusion does not
beat the position view. Acceptance requires a new versioned final receipt,
improved intent-to-run matching across all five non-oracle regimes, no
classification regression beyond a predeclared tolerance, and unchanged
supplied-topology/partial-overlap contracts.

### GlobalGeo connected-component policy

Restoration requires a documented mathematical policy for disconnected
landmark graphs, explicit behavior for cross-component correspondences, target
permutation equivariance, and at least 95% completion on the applicable frozen
cases. Partial-overlap support remains `unsupported` until all relevant cases
complete and score without using generator masks.

### CONE latent-dimension reconciliation

Restoration requires a deterministic shared-rank rule before Procrustes,
12/12 completion on the current applicable tournament cases, pairwise and
multiple-domain oracle/regression coverage, and no silent truncation that
changes the documented objective. Its current favorable small-case matching
must be retained within a predeclared tolerance.

## Governance result

The registry, generated comparison table, help pages, NEWS, and enforced test
matrix must agree with these states. A release check must fail if any
experimental, legacy/reference, deprecated, or quarantined method is described
as recommended. Future restoration or deprecation requires a new evidence
record rather than editing this decision around a failing result.
