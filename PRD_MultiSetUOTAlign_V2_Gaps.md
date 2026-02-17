# PRD v2: MultiSetUOTAlign (Gap-Closing + Production Hardening)

## Status (as of 2026-02-09)

**Already implemented in `manifoldalign` (WIP branch state)**

- `uot_build_cost()` supports `neighbor_mode = auto|dense|knn|radius|hybrid` with
  `radius`, `maxk`, `min_neighbors`, `dense_max_bytes`, and `ensure_cols`.
- Sparse cost builder returns **CSR+CSC** (`row_ptr/col_idx/cost` and
  `col_ptr/row_idx/cost_csc`) and supports distance filtering plus “extra edges”
  (used by `ensure_cols`).
- TI-Sinkhorn sparse solver supports CSR-only and CSR+CSC backends; the R wrapper
  dispatches to CSR+CSC when available.
- Hard constraints are supported via a general **group-equality** constraint
  (`group_X/group_Y`) implemented as **grouped neighbor search** (not post hoc
  filtering).
- `multiset_uot_align()` exposes neighborhood/constraint options and returns a
  `params` block for reproducible downstream mapping.
- Sparse time-series mapping is supported for **matrix signals** when the sparse
  cost includes CSC fields (C++ fast path).
- `uot_extract_coupling()` returns a reusable `dgCMatrix` map/coupling operator
  with optional pruning (`topk` / `threshold`), without materializing dense `π`.
- Synthetic tests cover: CSR↔CSC consistency, radius/hybrid behavior, grouped
  constraints, dense/sparse KKT identities, sparse `T×N` mapping, and operator reuse.

**Gaps remaining to meet PRD v2 milestones**

- Robust convergence upgrades: safeguarded Anderson acceleration, residual/objective
  stopping, warm-start across outer iterations, and optional ε-annealing.
- Additional coupling outputs: edge-list/CSR returns, row/col top‑k pruning modes,
  and diagnostics for isolate-heavy graphs.
- Sparse mapping ergonomics/performance: optional OpenMP and chunked mapping for very
  large `T` and/or `M`.
- Template regularization + structured diagnostics + benchmarks (as in Milestones 4+).

## 1) Summary

MultiSetUOTAlign v2 upgrades the current CPU MVP (KL translation-invariant
TI-Sinkhorn + star multiset template updates) into a production-ready functional
alignment engine for parcel and sparse-voxel regimes.

This PRD focuses on the **gaps** vs. PRD v1:

- Radius neighborhoods (in addition to kNN) using the existing
  `neighborweights` / `graphweights` “spatial_” machinery.
- Hard constraints (hemisphere, locality, block-diagonal masks).
- Robust convergence defaults (safeguarded Anderson acceleration; residual /
  objective-based stopping; warm-starts across outer loops).
- Explicit sparse coupling outputs and fast application to **time series** in
  sparse mode.
- Solver performance improvements (CSR+CSC, OpenMP-friendly kernels, memory
  caps).
- Diagnostics and reproducibility knobs that make failures debuggable.

**Target environment**

- Package: `manifoldalign`
- Compute: CPU only; heavy work in C++ (Rcpp/RcppArmadillo), OpenMP optional
- Typical scale:
  - Parcels: `N_k ~ 200–5000`, `M ~ 200–5000`, `K ~ 10–100`
  - Sparse voxel mode: `N_k ~ 20k–100k`, `M ~ 20k–100k`, `k ~ 16–128` edges/node

## 2) Goals and success metrics

### Goals

1. **Voxel-ready neighborhoods**
   - Support radius (`dthresh`) + `maxk` cap, not just fixed kNN.
   - Keep compute and memory near `O(N_k * k)` rather than `O(N_k * M)`.
2. **Plausible transport constraints**
   - Enforce “do not transport across hemispheres” and optional block masks.
   - Enforce hard locality gates (radius) even when feature distance is small.
3. **Reliable TI-Sinkhorn convergence**
   - Add robust acceleration and stopping criteria that correlate with solution
     quality (not just slow-moving iterates).
4. **Operational outputs**
   - Provide explicit sparse couplings (or pruned couplings) usable downstream.
   - Provide fast mapping of `T × N_k` time series into template space in sparse
     mode without densifying.
5. **Diagnostics**
   - Record per-subject convergence and mass-drop diagnostics; detect isolate
     nodes and empty neighborhoods early.

### Success metrics

- **Stability:** no NaNs/infs across standard parameter ranges; graceful failure
  messages for isolate-heavy graphs or zero-mass data.
- **Convergence quality:** fixed-point residual and marginal mismatch below
  thresholds; objective improvement plateaus consistently.
- **Performance:** on CPU, for sparse voxel mode (`N=100k`, `k=64`, `K=30`):
  - pair solve within minutes total (parallel over subjects),
  - memory stays within user-configured caps (no hidden dense allocations).
- **Downstream utility:** sparse couplings and mapping operators usable for:
  - activation maps, time series, low-dim embeddings, transported features.

## 3) Non-goals

- GPU/TPU support.
- Exact multi-marginal tensor UOT (exponential in K outside 1D structure).
- Full geodesic distances on surfaces (unless supplied as coordinates already).
- General smooth `φ` beyond KL in the initial v2 release (may be a follow-on).

## 4) Users and key use cases (unchanged, but stronger “voxel mode”)

- Researcher aligns:
  1) task contrasts (vector per node),
  2) connectivity embeddings (small `D`),
  3) time series (`T × N`).
- Needs to:
  - enforce hemisphere/locality constraints,
  - run in voxel-sparse mode, and
  - extract sparse mappings for downstream pipelines.

## 5) Functional requirements

### FR6 — Radius neighborhoods (and hybrid neighborhoods)

Add a new neighborhood mode for `uot_build_cost()` and the multiset loop:

- `neighbor_mode = "knn" | "radius" | "hybrid"`
- `radius` / `dthresh`: max anatomical distance allowed
- `maxk`: hard cap for radius search to prevent `O(n^2)` blowups
- `k_neighbors`: fallback kNN if radius returns too few neighbors

Implementation guidance

- Use the existing `neighborweights` / `graphweights` spatial stack:
  - In R, call `Rnanoflann::nn(data=Y, points=X, k=maxk, radius=radius)` to get
    `(indices, distances)`; this mirrors `graphweights::cross_spatial_adjacency`.
  - Feed those into a C++ cost builder that:
    - filters invalid indices / `dist > radius`,
    - computes `lambda_anat * dist^2 + lambda_feat * ||F_i - G_j||^2`,
    - emits CSR and CSC (see FR10).

### FR7 — Hard constraints and block masks

Support “allowed edges” constraints that can be applied before solving:

1. **Hemisphere constraint**
   - Inputs: `hemi_X` (length `N_k`), `hemi_Y` (length `M`)
   - Behavior: drop edges where `hemi_X[i] != hemi_Y[j]`
2. **Block-diagonal constraint**
   - Inputs: `block_X`, `block_Y` (e.g., network labels)
   - Behavior: only allow edges between matching blocks (or a user-provided
     compatibility matrix).
3. **Custom edge filter**
   - Allow user to provide an R callback (slow) for prototyping, and a
     fast-path C++ filtering mode for common constraints (hemi, blocks).

Design requirement

- Constraints must operate on sparse edge lists; they must never allocate dense
  `N×M` masks.

### FR8 — Robust TI-Sinkhorn acceleration + stopping

Upgrade the pairwise solver to provide robust defaults:

1. **Safeguarded Anderson acceleration (optional default-on)**
   - Parameters: `anderson = TRUE/FALSE`, `m = 3–6`, `reg = 1e-6`, `damping`
   - Safeguard rule:
     - compute candidate accelerated iterate `x_AA`,
     - accept only if fixed-point residual decreases sufficiently, else fallback
       to unaccelerated `T(x)`.
2. **Stopping criteria**
   - Replace/augment “iterate change” with:
     - fixed-point residual `||T(x) - x||_∞`,
     - marginal mismatch estimates `||pi1 - alpha||` and `||pi2 - beta||`
       (computed periodically from sparse edges),
     - objective change `ΔH_ε` (optional periodic check).
3. **Warm-start across outer iterations**
   - Store per-subject `(fbar,gbar)` and re-use as initialization when only
     `beta` and/or template features `G` have changed slightly.
4. **Epsilon annealing (optional)**
   - Schedule: start with larger `ε0`, decrease to `ε_final`.
   - Use warm-start across `ε` stages.

### FR9 — Explicit sparse coupling outputs

Expose coupling representations suitable for downstream use, without dense `π`:

- `uot_extract_coupling(...)` returns one of:
  1) CSR/CSC edge list with weights `w_ij`,
  2) `Matrix::sparseMatrix(i,j,x)` (dgCMatrix),
  3) pruned representation (top-k per row or per column).

Requirements

- Must support:
  - `pi2` and `pi1` retrieval,
  - coupling-weighted feature transport,
  - barycentric mapping operators.
- Provide pruning options:
  - `topk_row`, `topk_col`, or `threshold` on `w_ij`.

### FR10 — Sparse time-series mapping (C++ fast path)

Extend mapping to support sparse costs/couplings for matrix signals:

- Input: CSR coupling edges (or compute weights on the fly from potentials),
  signal matrix `S` (`T × N_k`), output `Ŝ` (`T × M`).
- Implement in C++:
  - streaming edge scan,
  - cache-friendly accumulation,
  - optional OpenMP over target columns or over edge blocks.

### FR11 — Template regularization (optional but important for fMRI plausibility)

Add optional regularizers on template weights and/or features:

1. **Spatial smoothness of `beta`**
   - Build template adjacency `A_Y` using `neighborweights::spatial_adjacency`
     with radius `dthresh` (already available).
   - Add an update step:
     - proximal smoothing, or
     - Laplacian/Tikhonov shrinkage: solve `(I + λ L) beta = beta_raw`.
2. **Feature smoothness of `G` (optional)**
   - Similar Laplacian smoothing on each feature dimension.

Requirement

- Regularizers must be configurable and must be off by default.

### FR12 — Diagnostics and reproducibility

Provide structured diagnostics:

- Per-subject:
  - iterations, residual, converged,
  - effective transported mass, mass dropped,
  - isolate counts (#rows with zero neighbors after constraints),
  - runtime (optional).
- Per-outer-iteration:
  - `||beta_new - beta||_∞`,
  - template mass normalization factor,
  - feature update norms.

Expose:

- `verbose` levels (`0,1,2`)
- `seed` for any randomized pruning or neighborhood sampling.

## 6) Non-functional requirements

### NFR4 — Memory caps

- Introduce explicit caps and early checks:
  - maximum nnz allowed per subject,
  - maximum `T×M` output size for mapping,
  - optional chunked mapping over timepoints.

### NFR5 — Threading model

- Prefer subject-parallelism in R (`mclapply`) for coarse parallelism.
- Within a subject, parallelize heavy kernels in C++ where safe:
  - CSR row-wise softmins: parallel over rows,
  - CSC col-wise softmins: parallel over columns.

## 7) API proposal (R level)

### Core builders

- `uot_build_cost(..., neighbor_mode = c("knn","radius","hybrid"), radius, maxk, ...)`
  - returns a sparse cost structure that includes both CSR and CSC.
- `uot_filter_edges(cost, hemi_X, hemi_Y, block_X, block_Y, ...)`
  - returns filtered sparse cost.

### Pairwise solver

- `uot_ti_sinkhorn_kl(cost, alpha, beta, epsilon, rho1, rho2, control=...)`
  - `control$anderson`, `control$stop`, `control$warm_start`

### Coupling + mapping

- `uot_extract_coupling(cost, alpha, beta, fbar, gbar, epsilon, prune=...)`
- `uot_apply_map(cost_or_coupling, signal, delta=..., mode=...)`

### Multiset driver

- `multiset_uot_align(..., control_pair=..., control_outer=..., constraints=..., regularize=...)`
  - warm-start per subject by default
  - optional annealing schedule

## 8) Implementation plan (C++ modules)

### C++ additions/changes

1. **Sparse cost structure with CSR+CSC**
   - New builder that emits both formats from `(idx,dists)`:
     - CSR for `Smin_beta_eps(C - g)` row-wise
     - CSC for `Smin_alpha_eps(C^T - f)` col-wise
2. **TI-Sinkhorn sparse solver refactor**
   - Accept CSR+CSC and compute row/col softmins without global reductions.
   - Add optional Anderson acceleration with minimal history allocations.
   - Add periodic checks for marginal mismatch and objective change.
3. **Coupling extraction**
   - Compute `w_ij` on edges and optionally output pruned edges.
4. **Sparse matrix mapping**
   - Add `uot_apply_map_sparse_mat_cpp` for `T×N` signals.

### R wrappers + glue

- Update `uot_build_cost()` to support `neighbor_mode="radius"` using
  `Rnanoflann::nn` (already the engine under `neighborweights` spatial tools).
- Expose constraint helpers and diagnostics.
- Keep the high-level multiset loop in R; push heavy scans to C++.

## 9) Testing and benchmarking

### Unit tests

- kNN vs radius mode sanity (small synthetic).
- Hemisphere/block constraint filtering correctness.
- Dense vs sparse equivalence on full bipartite graphs (already present; extend).
- Coupling extraction:
  - weights nonnegative, finite,
  - `pi1/pi2` match expected computed marginals within tolerance.
- Sparse time-series mapping:
  - compare sparse vs dense mapping on small problems.

### Benchmarks

- Parcel mode: `N=400, M=400, K=30`.
- Sparse voxel mode: `N=50k–100k, k=32–128, K=10–30`.
- Track:
  - time per pair solve,
  - outer loop time,
  - peak memory.

## 10) Risks and mitigations

- **Isolates after constraints (empty neighborhoods)**
  - Mitigation: detect early; fallback to kNN fill-in; warn with isolate counts.
- **Acceleration instability**
  - Mitigation: strict safeguards + damping; allow disabling globally.
- **Memory blowups in radius search**
  - Mitigation: require `maxk` cap; require `radius` sanity; enforce nnz caps.
- **Parameter sensitivity**
  - Mitigation: provide a recommended scaling pipeline and default schedules.

## 11) Milestones

1. Radius/hybrid neighborhoods + constraint filtering + CSR+CSC cost structure.
2. TI-Sinkhorn v2: Anderson + robust stopping + warm-start.
3. Coupling extraction + sparse `T×N` mapping (C++).
4. Template regularization + diagnostics + benchmarks.
