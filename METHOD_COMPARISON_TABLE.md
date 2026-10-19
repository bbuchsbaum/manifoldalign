# Alignment method capabilities and evidence

Status: experimental package inventory, reviewed 2026-07-22.

This document is generated in part from `alignment_method_registry()`, the
machine-readable source of truth for public method contracts. The table between
the generated markers is checked in tests and must not be edited independently
of the registry.

## Current decision

No alignment method in the package is currently classified as stable or
recommended. The registry records what each implementation actually accepts and
what evidence exists; it does not assert a universal accuracy winner.

- RAPID-MA is currently the only implementation with a native API for supplied
  directed, nonnegative, weighted multirelation graphs together with mixed node
  attributes, missing attributes, and spatial positions. It remains
  experimental, does not support signed edge weights, and completed the frozen
  tournament without becoming a universal accuracy winner: its matching loses
  to the position-only control in every non-oracle regime.
- GRASP, CONE-Align, PARROT, Token-OT, Spectral MNN, SSMA, GlobalGeo, and
  related spectral methods construct topology from feature matrices. Their
  current APIs do not consume an arbitrary observed adjacency matrix.
- Gromov-Wasserstein, FPGW, and Multiset UOT operate on metric-measure geometry.
  Pairwise distances are not equivalent to native support for an observed,
  directed, attributed graph.
- KEMA, GPCA, Low-rank Alignment, Generalized Procrustes, and related methods
  remain candidates for supervised feature-manifold problems, subject to their
  stated label, anchor, and evidence requirements.
- `multiscale_manifold_align()` is deprecated because the implementation fails
  mandatory validation and scale cases. It is retained temporarily for
  reproducibility, not method selection.
- GlobalGeo, CONE-Align, and CONE-Align Multiple are legacy/reference methods
  pending the repair gates in the
  [sunset review](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/ALIGNMENT_METHOD_SUNSET_REVIEW.md).

## Reading the registry

Support values are contracts: `native` means the public method consumes the
input directly; `feature_matrix` means numeric columns are accepted as a
feature matrix; `derived_only` means the method derives the object internally;
`experimental` means an API exists without stable evidence; `unsupported`
means callers must not infer support; and `not_applicable` means the concept
does not belong to that input model.

C/A/S abbreviates correctness, accuracy, and scalability evidence. A stale
scale ceiling is historical evidence whose current runner or implementation no
longer matches completely. An unvalidated ceiling is deliberately not replaced
with a qualitative claim such as "large scale."

<!-- BEGIN GENERATED METHOD TABLE -->
| Method | Lifecycle/evidence | Input contract | Structure and node data | Supervision/shape | Output, cost, scale |
|---|---|---|---|---|---|
| RAPID-MA (`rapid_ma`) | experimental; C/A/S=invariant/comparative/comparative | supplied_graph | topology=native; weights=native; directed=native; multirelation=native; attributes=native/native; positions=native | labels=optional; anchors=optional; overlap=native; unequal=native; domains=2..Inf; OOS=native | embeddings, couplings, node_correspondences; O(E s t + N q I) for bounded relation width, sketch width s, prototypes q; storage O(E + N(s + q)) additional sparse storage; scale 22,600 total_nodes |
| Token-OT Graph Align (`token_ot_graph_align`) | experimental; C/A/S=invariant/comparative/comparative | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=optional; overlap=unsupported; unequal=native; domains=2..2; OOS=unsupported | embeddings, couplings, node_correspondences; O(E k + I C) with bounded candidate support C; storage O(E + C + N k) for sparse candidates and embeddings; scale 5,500 total_nodes |
| Spectral MNN Align (`spectral_mnn_align`) | experimental; C/A/S=invariant/comparative/comparative | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=optional; anchors=optional; overlap=unsupported; unequal=native; domains=2..Inf; OOS=experimental | embeddings, transforms; O(E k + N log N + N k^2); storage O(E + N k); scale 11,000 total_nodes |
| SSMA (`ssma_align`) | experimental; C/A/S=invariant/comparative/comparative | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=optional; anchors=optional; overlap=unsupported; unequal=native; domains=2..Inf; OOS=experimental | embeddings, transforms; O(E k + N r^2) for reduced rank r; storage O(E + N r); scale 11,000 total_nodes |
| GlobalGeo (`global_geo_align`) | legacy_reference; C/A/S=oracle/comparative/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=optional; anchors=optional; overlap=unsupported; unequal=native; domains=2..Inf; OOS=native | embeddings, transforms; O(E k + N r^2 + r^3); storage O(E + N r + r^2); scale unvalidated |
| KEMA / REKEMA (`kema`) | experimental; C/A/S=oracle/comparative/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=required; anchors=unsupported; overlap=unsupported; unequal=native; domains=2..Inf; OOS=experimental | embeddings, transforms; O(N^3) full; reduced/operator paths depend on landmark rank; storage O(N^2) full kernel/graph storage; reduced kernels remain stored; scale failed by 550 total_nodes |
| GPCA Align (`gpca_align`) | experimental; C/A/S=invariant/comparative/smoke | feature_manifold | topology=unsupported; weights=unsupported; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=required; anchors=unsupported; overlap=unsupported; unequal=native; domains=2..Inf; OOS=native | embeddings, transforms; O(N r^2 + r^3) after sparse similarity assembly; storage O(nnz(S) + N r + p r); scale unvalidated |
| Low-rank Alignment (`lowrank_align`) | experimental; C/A/S=invariant/comparative/smoke | feature_manifold | topology=unsupported; weights=unsupported; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=required; anchors=unsupported; overlap=unsupported; unequal=native; domains=2..Inf; OOS=native | embeddings, transforms; O(N r^2 + r^3) operator or explicit reduced solve; storage O(nnz(S) + N r + p r); scale unvalidated |
| Linear Similarity Embedding (`linear_sim_embed`) | experimental; C/A/S=oracle/single_fixture/smoke | feature_manifold | topology=unsupported; weights=unsupported; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=not_applicable; unequal=not_applicable; domains=1..1; OOS=native | embeddings, transforms; O(I N^2 k) for a dense target-similarity objective; storage O(N^2 + p k); scale unvalidated |
| Generalized Procrustes (`generalized_procrustes`) | experimental; C/A/S=oracle/single_fixture/smoke | feature_manifold | topology=unsupported; weights=unsupported; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=required; anchors=unsupported; overlap=native; unequal=native; domains=2..Inf; OOS=native | embeddings, transforms; O(I sum_d p k_d^2) plus SVD initialization; storage O(sum_d p k_d); scale unvalidated |
| GRASP (`grasp`) | experimental; C/A/S=invariant/comparative/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=unsupported; unequal=native; domains=2..2; OOS=experimental | embeddings, transforms, node_correspondences; O(E k + N k^2 + n_max^3) including dense assignment; storage O(E + N k + n_1 n_2); scale unvalidated |
| GRASP Multiset (`grasp_multiset`) | experimental; C/A/S=invariant/single_fixture/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=unsupported; unequal=native; domains=2..Inf; OOS=experimental | embeddings, transforms, node_correspondences; O(sum_d E_d k + I D k^3 + assignment); storage O(sum_d(E_d + n_d k) + pairwise assignment costs); scale unvalidated |
| CONE-Align (`cone_align`) | legacy_reference; C/A/S=invariant/comparative/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=optional; overlap=unsupported; unequal=native; domains=2..2; OOS=experimental | embeddings, transforms, node_correspondences; O(E k + I(n_1 n_2 k + assignment)); storage O(E + N k + n_1 n_2); scale unvalidated |
| CONE-Align Multiple (`cone_align_multiple`) | legacy_reference; C/A/S=invariant/single_fixture/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=unsupported; unequal=native; domains=3..Inf; OOS=experimental | embeddings, transforms, node_correspondences; O(sum_d E_d k + I sum_d n_ref n_d k); storage O(sum_d(E_d + n_d k + n_ref n_d)); scale unvalidated |
| PARROT (`parrot`) | experimental; C/A/S=invariant/comparative/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=required; overlap=unsupported; unequal=native; domains=2..2; OOS=experimental | embeddings, couplings, node_correspondences; O(I n_1 n_2 + random-walk feature construction); storage O(n_1 n_2 + n_1^2 + n_2^2); scale unvalidated |
| Gromov-Wasserstein (`gromov_wasserstein`) | experimental; C/A/S=invariant/comparative/smoke | metric_measure | topology=not_applicable; weights=not_applicable; directed=not_applicable; multirelation=not_applicable; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=unsupported; unequal=native; domains=2..Inf; OOS=experimental | couplings; O(I(n_1^2 n_2 + n_1 n_2^2)); storage O(n_1^2 + n_2^2 + n_1 n_2); scale unvalidated |
| Fused Partial Gromov-Wasserstein (`fpgw`) | experimental; C/A/S=oracle/comparative/smoke | metric_measure | topology=not_applicable; weights=not_applicable; directed=not_applicable; multirelation=not_applicable; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=native; unequal=native; domains=2..Inf; OOS=experimental | couplings; O(I(n_1^2 n_2 + n_1 n_2^2) + I LP(n_1 n_2)); storage O(n_1^2 + n_2^2 + n_1 n_2); scale unvalidated |
| Coupled Diagonalization (`coupled_diagonalization`) | experimental; C/A/S=oracle/single_fixture/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=optional; overlap=experimental; unequal=native; domains=2..Inf; OOS=native | embeddings, transforms; O(sum_d E_d k + I sum_d n_d k^2); storage O(sum_d(E_d + n_d k + p_d k)); scale unvalidated |
| Multiset Manifold Alignment (`mma_align_multiple`) | experimental; C/A/S=invariant/single_fixture/smoke | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=unsupported; unequal=native; domains=3..Inf; OOS=native | embeddings, transforms, node_correspondences; O(sum_d E_d k + I sum_d n_d K^2 + assignment); storage O(sum_d(E_d + n_d k + n_d K)); scale unvalidated |
| Multiscale Manifold Align (`multiscale_manifold_align`) | deprecated; C/A/S=smoke/comparative/comparative | feature_derived_graph | topology=unsupported; weights=derived_only; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=optional; anchors=optional; overlap=experimental; unequal=native; domains=2..Inf; OOS=experimental | embeddings, transforms; Declared sparse spectral hierarchy; current solver fails mandatory scale cases; storage Declared O(E + N r); not validated; scale failed by 550 total_nodes |
| Multiset UOT (`multiset_uot_align`) | experimental; C/A/S=oracle/single_fixture/smoke | metric_measure | topology=not_applicable; weights=not_applicable; directed=not_applicable; multirelation=not_applicable; attributes=native/unsupported; positions=native | labels=unsupported; anchors=optional; overlap=native; unequal=native; domains=2..Inf; OOS=unsupported | couplings, transforms; O(I nnz(C)) sparse or O(I n m) dense; storage O(nnz(C)) sparse or O(n m) dense; scale unvalidated |
| Event-keyed SUMCOR (`sumcor_align`) | experimental; C/A/S=oracle/single_fixture/none | feature_manifold | topology=unsupported; weights=unsupported; directed=unsupported; multirelation=unsupported; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=required; overlap=native; unequal=native; domains=2..Inf; OOS=native | embeddings, transforms; O(sum_s n_s p_s^2 + I m^2 E k p) for m domains, E events and basis width p; storage O(sum_s n_s p_s + m E p); scale unvalidated |
| OT-Procrustes (`ot_procrustes`) | experimental; C/A/S=invariant/single_fixture/smoke | feature_manifold | topology=not_applicable; weights=not_applicable; directed=not_applicable; multirelation=not_applicable; attributes=feature_matrix/unsupported; positions=unsupported | labels=unsupported; anchors=unsupported; overlap=unsupported; unequal=native; domains=2..2; OOS=native | couplings, transforms; O(I(n_1 n_2 d + d^3)); storage O(n_1 n_2 + d^2); scale unvalidated |
<!-- END GENERATED METHOD TABLE -->

## Benchmark implications

The table is a capability inventory, not a leaderboard. Comparative accuracy
must be reported separately by topology regime, supervision regime, overlap,
attribute type, graph size, and resource budget. Unsupported inputs and missing
optional dependencies must produce typed benchmark outcomes rather than numeric
scores. Stable or recommended status requires the lifecycle policy's
correctness, accuracy, and scalability gates, including current provenance and
at least one out-of-family fixture.

Protocol v2 supplies the frozen attributed-graph inputs, 30 final accuracy
cases, negative controls, provenance hashes, and no-final-tuning boundary. The
42-case, 14-method tournament is complete; its
[raw receipt](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/raw_results.csv) and checked
[interpretation](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/INTERPRETATION.md) report
regime-specific winners and failures without declaring one package-wide winner.
