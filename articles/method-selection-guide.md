# Method Selection Guide

This package provides more than a dozen alignment methods. Choosing the
right one depends on the structure of your problem, not on a single
accuracy leaderboard. This article walks through the decision space and
connects each branch to the method (or methods) that the current
benchmarks support.

## The four questions that matter

Before looking at any method, answer these about your data:

1.  **Correspondences** – Do you know which rows in domain A match which
    rows in domain B?

    - *Exact*: row IDs are shared across all domains (e.g. same
      subjects, same voxels).
    - *Partial anchors*: a subset of rows are matched; the rest are
      unmatched.
    - *Labels only*: you know class/condition labels but not row-level
      correspondences.
    - *None*: fully unsupervised.

2.  **Overlap** – Do all domains observe the same set of entities?

    - *Full*: every entity appears in every domain.
    - *Partial / open-set*: some entities are private to a single
      domain.

3.  **Number of domains** – Two, or three-or-more?

4.  **Output needed** – Do you need a shared latent embedding, a linear
    transform you can apply to new data, or an explicit sample-to-sample
    transport plan?

The table and flowchart below map combinations of these answers to
recommended methods.

## Method families at a glance

| Method family | Correspondences | Overlap | Domains | Output type | Feature dims | Vignette |
|:---|:---|:---|:---|:---|:---|:---|
| global_geo_align | Exact row IDs | Full | 2+ | Latent embedding | Any | – |
| spectral_mnn_align | Exact row IDs | Full or partial | 2+ | Latent embedding | Any | – |
| ssma_align | Exact row IDs | Full | 2+ | Latent embedding | Any | – |
| gpca_align | Labels | Full | 2+ | Latent embedding | Any | `gpca-align` |
| lowrank_align | Labels | Full | 2 | Latent embedding | Any | `lowrank-align` |
| kema | Labels | Full | 2 | Kernel embedding | Any | `kema-overview` |
| generalized_procrustes | Exact row IDs | Full | 2+ | Rotation + embedding | Must match | – |
| multiscale_manifold_align | Exact row IDs | Full | 2+ | Hierarchical embedding | Any | – |
| coupled_diagonalization | Exact row IDs (anchors) | Full (via anchors) | 2+ | Latent embedding | Any | `coupled-diagonalization` |
| mma_align_multiple | None (unsupervised) | Full | 3+ | Consensus embedding | Any | `mma-align` |
| grasp / grasp_multiset | None or partial | Full or partial | 2 (multiset: 3+) | Embedding + transport | Any | `grasp-align` |
| parrot | Anchors (partial) | Full or partial | 2 | Rotation + permutation | Any | `parrot-align` |
| fpgw | None (structural) | Full | 2 | Transport plan + embedding | Any | `fpgw_tutorial` |
| gromov_wasserstein | None (structural) | Full | 2 | Transport plan | Any | – |
| token_ot_graph_align | None or partial links | Full or partial | 2 | Transport plan + embedding | Any | – |
| cone_align / cone_align_multiple | Exact row IDs | Full | 2 (multiple: 3+) | Linear transform + embedding | Any | `cone-align` |
| multiset_uot_align | None (template geometry) | Full | 3+ | Transport-based mapping | Shared geometry dim | – |

## Decision flowchart

The flowchart below is a text rendering of the decision tree. Start at
the top and follow the branch that matches your problem.

    Do you have exact row-level correspondences across domains?
    │
    ├─ YES
    │  │
    │  ├─ Is overlap full (every entity in every domain)?
    │  │  │
    │  │  ├─ YES
    │  │  │  │
    │  │  │  ├─ Do feature dimensions match across domains?
    │  │  │  │  ├─ YES ──► generalized_procrustes (simple baseline)
    │  │  │  │  │         then try global_geo_align (benchmark leader)
    │  │  │  │  └─ NO ───► global_geo_align (benchmark leader)
    │  │  │  │             or ssma_align / spectral_mnn_align
    │  │  │  │
    │  │  │  └─ Do you want multiscale decomposition?
    │  │  │     ├─ YES ──► multiscale_manifold_align (research family)
    │  │  │     └─ NO ───► global_geo_align
    │  │  │
    │  │  └─ NO (partial / open-set overlap)
    │  │     └──────────► spectral_mnn_align (benchmark leader here)
    │  │
    │  └─ 3+ domains needing joint consensus?
    │     └──────────► coupled_diagonalization or cone_align_multiple
    │
    ├─ PARTIAL ANCHORS (some rows matched, others not)
    │  │
    │  ├─ Are anchors sparse / noisy?
    │  │  ├─ YES ──► parrot (robust to anchor noise)
    │  │  └─ NO ───► token_ot_graph_align or grasp
    │  │
    │  └─ 3+ domains?
    │     └──────────► grasp_multiset
    │
    ├─ LABELS ONLY (class/condition but no row IDs)
    │  │
    │  ├─ Labels clean?
    │  │  ├─ YES
    │  │  │  ├─ Want speed + simplicity? ──► gpca_align (safest baseline)
    │  │  │  ├─ Want latent geometry? ────► kema (linear kernel first)
    │  │  │  └─ Hard nonlinear regime? ──► lowrank_align (strongest here)
    │  │  │
    │  │  └─ NO (noisy / corrupted labels)
    │  │     └──────────► gpca_align (most robust to label noise)
    │  │
    │  └─ 2+ domains with labels?
    │     └──────────► gpca_align or lowrank_align via hyperdesign
    │
    └─ NO CORRESPONDENCES AT ALL (fully unsupervised)
       │
       ├─ Do domains share intrinsic geometry?
       │  ├─ YES ──► gromov_wasserstein or fpgw (structure matching)
       │  └─ NO ───► mma_align_multiple (3+ domains, consensus)
       │
       ├─ Do you have a shared geometry template?
       │  └─ YES ──► multiset_uot_align
       │
       └─ 3+ domains?
          └──────────► mma_align_multiple (consensus baseline)

## What the benchmarks say

The package tracks synthetic benchmark results across five fast-profile
scenarios that stress different problem axes. The `benchmark-guidance`
vignette has full tables; here is the executive summary. These are
frozen historical outputs. Current
[`alignment_method_registry()`](https://bbuchsbaum.github.io/manifoldalign/reference/alignment_method_registry.md)
lifecycle states take precedence over historical rankings.

### Exact correspondence, full overlap

| Scenario | Method | Rank | Latent cor | Centroid acc | Runtime (s) |
|:---|:---|---:|---:|---:|---:|
| hard_nonisometric | global_geo_rowid_exact | 1.500 | 0.960 | 1.000 | 0.008 |
| hard_nonisometric | ssma_exact | 2.333 | 0.984 | 0.960 | 0.016 |
| hard_nonisometric | spectral_mnn_exact | 2.833 | 0.791 | 1.000 | 0.015 |
| isometric_curve | global_geo_rowid_exact | 1.500 | 0.991 | 1.000 | 0.007 |
| isometric_curve | ssma_exact | 1.500 | 0.999 | 1.000 | 0.015 |
| isometric_curve | spectral_mnn_exact | 3.000 | 0.940 | 0.987 | 0.015 |
| linear_affine | global_geo_rowid_exact | 2.667 | 0.985 | 0.983 | 0.071 |
| linear_affine | ssma_exact | 2.667 | 0.972 | 0.983 | 0.063 |
| linear_affine | multiscale_spectral_energy | 4.000 | 0.837 | 0.917 | 0.053 |
| noisy_label_shift | global_geo_rowid_exact | 1.333 | 0.987 | 0.972 | 0.007 |
| noisy_label_shift | ssma_exact | 2.167 | 0.976 | 0.972 | 0.015 |
| noisy_label_shift | spectral_mnn_exact | 2.500 | 0.927 | 0.972 | 0.015 |
| partial_open_set | spectral_mnn_exact | 1.333 | 0.898 | 0.931 | 0.015 |
| partial_open_set | ssma_exact | 1.667 | 0.973 | 0.861 | 0.015 |
| partial_open_set | multiscale_spectral_hybrid | 3.833 | 0.770 | 0.764 | 0.018 |

**Takeaways:**

- `global_geo_align` leads on `linear_affine`, `isometric_curve`,
  `hard_nonisometric`, and `noisy_label_shift` with the fastest
  runtimes.
- `spectral_mnn_align` takes over on `partial_open_set` where it handles
  unmatched rows gracefully.
- `ssma_align` is a strong runner-up on curve and nonlinear regimes.

### Label-supervised

| Scenario          | Method            |  Rank | Latent cor | Centroid acc | Runtime (s) |
|:------------------|:------------------|------:|-----------:|-------------:|------------:|
| hard_nonisometric | lowrank_label     | 1.833 |      0.877 |        1.000 |       0.017 |
| hard_nonisometric | gpca_label        | 2.833 |      0.858 |        1.000 |       0.021 |
| hard_nonisometric | kema_linear_label | 2.833 |      0.995 |        0.973 |       0.022 |
| isometric_curve   | kema_linear_label | 1.667 |      0.992 |        1.000 |       0.020 |
| isometric_curve   | lowrank_label     | 2.000 |      0.890 |        1.000 |       0.018 |
| isometric_curve   | gpca_label        | 2.333 |      0.987 |        1.000 |       0.021 |
| linear_affine     | kema_linear_label | 1.667 |      0.966 |        0.983 |       0.024 |
| linear_affine     | gpca_label        | 2.167 |      0.961 |        1.000 |       0.062 |
| linear_affine     | lowrank_label     | 2.167 |      0.888 |        1.000 |       0.061 |
| noisy_label_shift | kema_linear_label | 1.667 |      0.989 |        0.889 |       0.021 |
| noisy_label_shift | gpca_label        | 2.000 |      0.918 |        0.972 |       0.022 |
| noisy_label_shift | lowrank_label     | 2.333 |      0.614 |        0.833 |       0.016 |

**Takeaways:**

- `gpca_align` is the safest label baseline: best on
  `noisy_label_shift`, competitive everywhere, and has the lowest
  seed-to-seed variance.
- `kema` with a linear kernel is fast and strong on latent correlation,
  but more seed-sensitive on curved data. Validate with multiple seeds.
- `lowrank_align` is the strongest on the hard nonlinear regime.
- `kema` with an RBF kernel is slow and highly variable on these small
  problems. Use it only when you have reason to believe the kernel is
  well-calibrated.

### Unsupervised multiset and transport

- `mma_align_multiple` is the current baseline for 3+ domain consensus
  alignment with no side information.
- `multiset_uot_align` occupies a separate regime: it requires a shared
  geometry template but achieves very high latent correlation when one
  is available.

### Unsupervised pairwise (no correspondences, 2 domains)

|  | Scenario | Method | Rank | Latent cor | Centroid acc | Runtime (s) |
|:---|:---|:---|---:|---:|---:|---:|
| 55 | hard_nonisometric | grasp_unsupervised | 1.000 | 0.406 | 1.000 | 0.015 |
| 36 | isometric_curve | cone_align_unsupervised | 1.000 | 0.766 | 1.000 | 0.010 |
| 37 | isometric_curve | grasp_unsupervised | 2.000 | 0.230 | 0.613 | 0.015 |
| 18 | linear_affine | cone_align_unsupervised | 1.000 | 0.899 | 1.000 | 0.029 |
| 19 | linear_affine | grasp_unsupervised | 2.000 | 0.381 | 0.850 | 0.055 |
| 85 | noisy_label_shift | cone_align_unsupervised | 1.000 | 0.894 | 0.986 | 0.011 |
| 86 | noisy_label_shift | grasp_unsupervised | 2.000 | 0.352 | 0.806 | 0.015 |
| 65 | partial_open_set | cone_align_unsupervised | 1.167 | 0.817 | 1.000 | 0.010 |
| 66 | partial_open_set | grasp_unsupervised | 1.833 | 0.520 | 1.000 | 0.014 |

These methods receive **no correspondences at all** — they discover
alignment purely from graph structure. `cone_align` uses iterative
Procrustes on spectral embeddings; `grasp` uses diffusion-based
descriptors. Both run on the first two views of each scenario.

### Partial anchor

|  | Scenario | Method | Rank | Latent cor | Centroid acc | Runtime (s) |
|:---|:---|:---|---:|---:|---:|---:|
| 38 | isometric_curve | parrot_partial_anchor | 1 | 0.508 | 0.827 | 0.114 |
| 20 | linear_affine | parrot_partial_anchor | 1 | 0.796 | 0.967 | 0.165 |
| 87 | noisy_label_shift | parrot_partial_anchor | 1 | 0.825 | 0.986 | 0.120 |
| 67 | partial_open_set | parrot_partial_anchor | 1 | 0.891 | 0.986 | 0.117 |

`parrot` receives 50% of shared row IDs as anchors (the rest are
hidden). This simulates the common real-world case where only a subset
of correspondences are known.

## Methods not yet in the synthetic benchmark

The benchmark now covers the major method families. The remaining gaps
are:

| Method | Why not benchmarked | Integration difficulty |
|:---|:---|:---|
| grasp_multiset | Multiset consensus; needs 3+ domain scenario integration | Medium – needs multiset scenario harness |
| fpgw | Pairwise OT; returns transport plan, needs embedding extraction | Medium – need to extract embedding from transport plan |
| gromov_wasserstein | Returns transport plan / barycentric weights, not latent embedding | High – output is coupling matrix, not embedding |
| token_ot_graph_align | Pairwise graph OT; needs link specification in benchmark harness | Medium – need link specification and embedding extraction |
| cone_align_multiple | Multiset; needs multi-domain scenario integration | Medium – needs multiset scenario harness |

**Recommended next steps for benchmark coverage:**

1.  **Transport methods**: `fpgw` and `token_ot_graph_align` both return
    `multiblock_biprojector` objects, so embedding extraction should
    work. The main work is constructing the right input format (cost
    matrices for FPGW, graph links for token OT).
2.  **Pure OT**: `gromov_wasserstein` returns couplings rather than
    embeddings. Benchmarking it requires either a transport-plan quality
    metric (e.g. fraction of correct correspondences recovered) or a
    barycentric projection step to produce comparable embeddings.
3.  **Multiset extensions**: `grasp_multiset` and `cone_align_multiple`
    need the full 3+ domain scenario infrastructure.

## Practical considerations

### Speed

From the benchmark timing data on small problems (n = 20–25 per domain):

- **Fastest** (\< 0.02s): `global_geo_align`, `spectral_mnn_align`,
  `ssma_align`, `kema` (linear), `mma_align_multiple`
- **Moderate** (0.02–0.2s): `gpca_align`, `lowrank_align`, multiscale
  variants
- **Slow** (\> 1s): `coupled_diagonalization` (~3–4s), `kema` (RBF, full
  landmarks, 10–50s on n=25)

For large problems (n \> 1000), methods with O(n^3) kernel operations
(`kema` with full landmarks, `gromov_wasserstein`) will become
bottlenecks. Use landmarking / reduced representations where available.

### Out-of-sample prediction

Not all methods support applying the learned alignment to new data.
Methods that return `multiblock_biprojector` or `projector` objects
support
[`oos_predict()`](https://bbuchsbaum.github.io/manifoldalign/reference/oos_predict.md):

- `spectral_mnn_align`, `ssma_align`, `global_geo_align`
- `gpca_align`, `lowrank_align`, `kema`
- `grasp`, `parrot`, `fpgw`, `cone_align`, `token_ot_graph_align`

Methods that produce point embeddings only (no projection operator):

- `generalized_procrustes` (rotation matrices can be applied manually)
- `mma_align_multiple` (consensus is fitted, no general OOS)
- `coupled_diagonalization` (fitted embeddings only)
- `gromov_wasserstein` (transport plan is sample-specific)
- `multiset_uot_align` (template-based mapping via `multiset_uot_map`)

### Interpretability

- **Linear transforms** (`generalized_procrustes`, `cone_align`,
  `lowrank_align`): the alignment is a matrix you can inspect.
- **Spectral embeddings** (`global_geo_align`, `spectral_mnn_align`,
  `ssma_align`, `gpca_align`): the latent dimensions are eigenvectors of
  a graph Laplacian or generalized eigenvalue problem.
- **Kernel methods** (`kema`): alignment lives in a kernel-induced
  feature space; harder to interpret directly.
- **Transport plans** (`fpgw`, `gromov_wasserstein`, `parrot`): you get
  an explicit sample-to-sample correspondence matrix.

## Quick-reference decision table

For the impatient, here is the shortest path to a reasonable method:

| Your situation | Start here | Then try |
|:---|:---|:---|
| Exact row IDs, full overlap | `global_geo_align` | `ssma_align` |
| Exact row IDs, partial overlap | `spectral_mnn_align` | `grasp` |
| Class labels, clean | `gpca_align` | `kema` (linear) |
| Class labels, noisy | `gpca_align` | `lowrank_align` |
| Hard nonlinear, labels | `lowrank_align` | `kema` (linear) |
| No supervision, 3+ domains | `mma_align_multiple` | – |
| No supervision, 2 domains | `fpgw` or `gromov_wasserstein` | `grasp` |
| Shared geometry template | `multiset_uot_align` | – |
| Need explicit correspondences | `parrot` or `fpgw` | `token_ot_graph_align` |
| Feature dims must match | `generalized_procrustes` | `cone_align` |

## Where to go next

- [`vignette("benchmark-guidance")`](https://bbuchsbaum.github.io/manifoldalign/articles/benchmark-guidance.md)
  for the full regime-specific benchmark tables and stability
  diagnostics.
- Individual method vignettes (see the “Vignette” column in the table
  above) for API details, parameter tuning, and worked examples.
- [`vignette("quickstart-hyperdesign")`](https://bbuchsbaum.github.io/manifoldalign/articles/quickstart-hyperdesign.md)
  for how to set up multi-domain data as a hyperdesign object, which
  most methods accept.
