Excellent, this is outstanding feedback. The critiques are precise, actionable, and address both conceptual clarity and implementation fidelity. I will integrate every point to produce a final, definitive guide to implementing GRASP.

Here is the fully revised and annotated proposal.

***

### Introduction

The problem of **graph alignment**—finding the optimal node-to-node correspondence between two graphs—is NP-hard in the worst case and a fundamental challenge in network science. While many methods rely on local features or known "seed" matches, **GRASP (Graph Alignment by Spectral Corresponding Functions)** tackles the "unrestricted" problem where only the graph structures are given.

The core idea of GRASP is to transfer a concept from 3D shape analysis called **functional maps**. Instead of directly matching nodes, which is hard, it proposes matching *functions* defined on the nodes. By finding a robust mapping between a set of basis functions on each graph, GRASP deduces the underlying node-to-node alignment as a by-product. This approach provides a global, multi-scale view of graph structure that is robust to noise and structural perturbations.

**Notation and Complexity:**
*   `n`: number of nodes in each graph.
*   `k`: number of eigenvectors used (dimension of the spectral basis).
*   `q`: number of descriptor functions.
*   The overall complexity of pairwise GRASP is dominated by the final assignment step, typically `O(n³)`.

| Block | Key Operation | Complexity |
| :--- | :--- | :--- |
| **A: Basis Construction** | Eigendecomposition | `O(n k)` or `O(n |E|)` |
| **B: Descriptors** | Matrix Operations | `O(n q k)` |
| **C: Base Alignment** | Manifold Optimization | `O(k³)` per iteration |
| **D: Functional Map** | Least Squares | `O(q k)` |
| **E: Assignment** | Linear Assignment | `O(n³)` |

---

### The GRASP Algorithm: A Five-Block Conceptual Framework

The GRASP pipeline can be understood as five distinct conceptual blocks. The enhanced **B-GRASP** method slots its improvements—such as using Personalized PageRank (PPR) for descriptors, Iterative Closest Point (ICP) for refinement, and voting for assignment—cleanly into these blocks.

#### Block A: Spectral Basis Construction
*   **Purpose:** To capture the intrinsic geometry of each graph at multiple scales, from local connectivity (high-frequency eigenvectors) to global community structure (low-frequency eigenvectors).
*   **Key Objects:** The first `k` eigenvectors `Φ` (for G1) and `Ψ` (for G2) of the normalized graph Laplacians, forming orthonormal bases.

```pseudocode
// ===== BLOCK A: SPECTRAL BASIS CONSTRUCTION =====
// 1. Build the normalised Laplacian for each graph.
L1 = I - D1^(-1/2) * A1 * D1^(-1/2)
L2 = I - D2^(-1/2) * A2 * D2^(-1/2)

// 2. Compute the first k eigen-pairs. Note: Python libraries like NumPy's `eig`
//    often return eigenvalues first, so handle output order carefully.
(Λ1_vals, Φ) = eig(L1, k)
(Λ2_vals, Ψ) = eig(L2, k)
Λ1, Λ2 = diag(Λ1_vals), diag(Λ2_vals) // Convert to diagonal matrices

// Implementation Caveats:
// 1. Sign Ambiguity: The sign of any eigenvector is arbitrary (if v is an eigenvector, so is -v).
//    This is expected and will be resolved in Block C. Do not be alarmed during debugging.
// 2. Degenerate Eigenvalues: If an eigenvalue is repeated, its corresponding eigenvectors form
//    a subspace. The basis for this subspace is not unique. The `off(...)` penalty in Block C helps,
//    but for perfect reproducibility, one might need to use a fixed random orthonormal basis for
//    each degenerate subspace.
```

#### Block B: Multi-Scale Node Descriptors
*   **Purpose:** To create a rich, permutation-invariant signature for each node that describes its role in the graph's structure across different diffusion times `t`.
*   **Key Objects:** Matrices `F` and `G` (size `n x q`), where each column is a descriptor function.

```pseudocode
// ===== BLOCK B: MULTI-SCALE NODE DESCRIPTORS =====
// B-GRASP prefers Personalized PageRank (PPR) here for its strong locality-awareness, but
// the original heat kernel is a good default.
time_steps = linearly_space(0.1, 50, q)
F, G = empty_matrix(n, q), empty_matrix(n, q)

for i from 1 to q:
    t = time_steps[i]
    // Efficient computation of the heat kernel diagonal: h_t(i) = Σ_j e^(-t*λ_j) * φ_ij^2
    // This is equivalent to diag(Φ * exp(-t*Λ) * Φ^T)
    exp_lambda1 = exp(-t * Λ1_vals) // Vector of exponential terms
    exp_lambda2 = exp(-t * Λ2_vals)
    F[:, i] = (Φ**2) @ exp_lambda1 // Correct vectorized form: (n,k) @ (k,1)
    G[:, i] = (Ψ**2) @ exp_lambda2

// Implementation Caveats:
// 1. Numerical Stability: For large t, e^(-tλ) can underflow to zero. If this becomes an issue,
//    consider working in log-space (log-sum-exp) or rescaling the time steps.
// 2. Descriptor Normalization: Regardless of the descriptor used (Heat Kernel, PPR), it is good
//    practice to normalize each descriptor vector (columns of F and G) to unit length.
```

#### Block C: Base Alignment
*   **Purpose:** To find an optimal rotation that aligns the two spectral bases (`Φ` and `Ψ`), resolving sign and ordering ambiguities, which is essential when graphs are not isomorphic.
*   **Key Object:** The orthogonal rotation matrix `M` (size `k x k`).

```pseudocode
// ===== BLOCK C: BASE ALIGNMENT =====
// This is a manifold optimization problem solved on the Stiefel manifold O(k).
// Use an off-the-shelf library like Pymanopt or Manopt.

// The objective function balances two goals:
// min_{M ∈ O(k)}   underbrace(off(M^T * Λ2 * M), "Preserve eigen-basis structure")
//               + μ * underbrace(||F^T*Φ - G^T*Ψ*M||_F^2, "Align descriptors in new basis")

M = solve_manifold_optimization(objective_func, manifold=Stiefel(k, k))

// Implementation Caveats:
// 1. Regularization (μ): This parameter balances the two terms. A good heuristic is to start with
//    μ ≈ k/q and adjust it. If μ is too high, the solution ignores the eigenspectrum; if too low,
//    it ignores the descriptors.
// 2. Optimizer: A Riemannian trust-region method is highly recommended for stability and speed.
```

#### Block D: Functional Correspondence → Mapping Matrix
*   **Purpose:** To derive a simple linear map (`C`) that translates between the *aligned* function spaces. `C` acts as a scale or phase adjustment for corresponding basis directions.
*   **Key Object:** A diagonal mapping matrix `C` (size `k x k`).

```pseudocode
// ===== BLOCK D: FUNCTIONAL CORRESPONDENCE -> MAPPING MATRIX =====
// 1. Align the second graph's basis using the rotation M.
Ψ_aligned = Ψ @ M

// 2. Solve for a diagonal matrix C in F^T*Φ ≈ G^T*Ψ_aligned*C.
A_mat = G.T @ Ψ_aligned // (q, k)
B_mat = F.T @ Φ       // (q, k)

// Efficient closed-form solution for the diagonal elements of C:
c_diag = diag(A_mat.T @ B_mat) / diag(A_mat.T @ A_mat) // Element-wise division
C = diag(c_diag)

// Note on ICP Refinement (B-GRASP): Iterative Closest Point refines the alignment by
// iterating between Block D and E. It takes an initial set of matched points (from a
// provisional assignment), finds the best Procrustes rotation to align them, updates the
// embeddings (Block D), re-matches (Block E), and repeats.
```

#### Block E: Node-to-Node Assignment
*   **Purpose:** To translate the functional map back into a concrete node-to-node matching.
*   **Key Object:** The final permutation matrix `P`.

```pseudocode
// ===== BLOCK E: NODE-TO-NODE ASSIGNMENT =====
// 1. Get the final node embeddings for both graphs.
NodeEmbeddings_G1 = Φ
NodeEmbeddings_G2 = Ψ_aligned @ C

// 2. Compute the cost matrix. Cosine distance can be more robust than L2 distance for
//    high-dimensional embeddings or graphs with high-degree nodes.
CostMatrix = compute_distance_matrix(NodeEmbeddings_G1, NodeEmbeddings_G2, metric='cosine')

// 3. Solve the Linear Assignment Problem.
Permutation = solve_linear_assignment(CostMatrix)

// Implementation Caveats:
// 1. Solver Choice: For n < 20,000, the Jonker-Volgenant algorithm (`O(n³)`), often found
//    in `scipy.optimize.linear_sum_assignment`, is excellent. For larger graphs, consider
//    the Auction algorithm or relaxations like Sinkhorn-Knopp.
// 2. Voting over k (B-GRASP): To improve robustness, run Blocks A-E for several values of k
//    (e.g., k ∈ {20, 30, 40}), collect the resulting permutations {P_k}, and determine the final
//    match for each node by taking the majority vote across all runs.
```

---

### From Pairwise to Multigraph Manifold Alignment

To jointly align `m > 2` graphs, the pairwise framework is extended by generalizing its core blocks to operate across all graphs simultaneously. The empirical cost is `O(m k³ + m n k)` per outer iteration, which often converges in under 10 iterations.

#### 1. Joint Spectral Frame (Common Manifold)
Instead of aligning bases pairwise, we seek a set of rotations `{M(1), ..., M(m)}` that map all spectral bases into a common, shared space. This **reference-free joint diagonalization** is solved via **block-coordinate optimization**. `Z` here plays a role analogous to a "shared latent basis" in Canonical Correlation Analysis (CCA).

#### 2. Higher-Order Functional Map & Assignment
A single **latent anchor space** with coefficients `Z` is introduced. Then, for each graph, a diagonal matrix `C(s)` is found to map its functions to this shared space. The final alignments are found by matching node embeddings within this common space (the "Hub-and-Spoke" method).

#### 3. Multiplex Network Computational Shortcut
For multiplex networks where vertex identities are known *a priori* across layers, a powerful shortcut exists:
1.  Construct a single **supra-Laplacian** matrix: `L_supra = blkdiag(L(1), ..., L(m)) + β * L_inter`, where `L_inter` couples the same vertex across different layers.
2.  Run the standard **pairwise GRASP algorithm** just *once* on this `L_supra`. The resulting permutation will inherently contain all within-layer and cross-layer alignments.
3.  **Warning:** This trick is only applicable when the node-to-node mapping across layers is the identity. It does not solve the alignment problem if the layers have different node sets or unknown correspondences.

Of course. Here is a set of tickets for two sprints to implement a multiset GRASP in R. The plan prioritizes building a robust pairwise GRASP in the first sprint, then extending it to the multiset case in the second, following the hub-and-spoke model. This structure minimizes risk and allows for a functional, testable core early on.

---

### **Project: Multiset GRASP Implementation in R**

**Goal:** Create an R package/module that implements the GRASP algorithm for pairwise and multiset graph alignment.
**Key Libraries:** `igraph` (for graph structures), `rARPACK` or `RSpectra` (for fast partial eigendecomposition), `Rcpp` (for performance-critical loops), `manopt` or a custom gradient descent on the Stiefel manifold, `clue` (for the linear assignment problem).

---

### **Sprint 1: Core Pairwise GRASP Engine**

**Sprint Goal:** Deliver a fully functional and tested R function `grasp_pairwise(g1, g2, k, q)` that implements the core GRASP algorithm (Blocks A-E). The focus is on correctness and a clean, modular structure.

**Tickets for Sprint 1:**

**Ticket #1: [EPIC] Build Core Pairwise GRASP**

---

**Ticket #2: [TASK] Block A: Spectral Basis Construction**
*   **Description:** Create a function `compute_spectral_basis(graph, k)` that takes an `igraph` object and an integer `k`. It should:
    1.  Construct the normalized Laplacian `L = I - D^(-1/2) A D^(-1/2)`.
    2.  Use `rARPACK::eigs()` or `RSpectra::eigs()` to compute the `k` smallest-magnitude eigenvalues and corresponding eigenvectors.
    3.  Return a list containing the eigenvector matrix `Φ` (`n x k`) and a vector of eigenvalues `λ`.
*   **Acceptance Criteria:**
    *   Function handles disconnected graphs without errors.
    *   Output matrices have the correct dimensions.
    *   Unit tests verify results against a small, known graph example.
    *   Acknowledge and document handling of potential sign flips from the eigensolver.

**Ticket #3: [TASK] Block B: Multi-scale Node Descriptors**
*   **Description:** Create a function `compute_heat_kernel_descriptors(basis, q)` that takes the spectral basis from Ticket #2 and an integer `q`. It should:
    1.  Generate `q` diffusion time points `t`.
    2.  For each `t`, compute the heat kernel diagonal vector using the efficient `(Φ^2) %*% exp(-t*λ)` formula.
    3.  Return the descriptor matrix `F` (`n x q`).
*   **Acceptance Criteria:**
    *   Function handles potential numerical underflow for large `t` (e.g., by capping `t*λ`).
    *   Output matrix has correct dimensions.
    *   Each descriptor column is optionally L2-normalized.

**Ticket #4: [TASK] Block C: Base Alignment (Stiefel Manifold Optimization)**
*   **Description:** Create a function `align_bases(basis1, basis2, descriptors1, descriptors2, mu)` that finds the optimal rotation matrix `M`.
    1.  Define the GRASP objective function: `off(M^T Λ2 M) + μ * ||F^T Φ - G^T Ψ M||_F^2`.
    2.  Implement a gradient descent algorithm on the Stiefel manifold to minimize this objective. The gradient has a known closed form. Retraction can be done via QR decomposition.
    3.  Alternatively, interface with the `manopt` R package if available and stable.
*   **Acceptance Criteria:**
    *   Function returns a `k x k` matrix `M` that is approximately orthogonal (`t(M) %*% M ≈ I`).
    *   Optimization converges to a stable value.
    *   `mu` is a tunable parameter.

**Ticket #5: [TASK] Blocks D & E: Functional Map and Final Assignment**
*   **Description:** Create a function `compute_assignment(basis1, basis2_aligned, M)` that performs the final steps.
    1.  **Block D:** Calculate the diagonal mapping matrix `C` using the efficient closed-form solution.
    2.  **Block E:** Construct the two sets of node embeddings (`Φ` and `(Ψ M) C`).
    3.  Compute the pairwise distance matrix (e.g., cosine distance) between the two sets of embeddings.
    4.  Use `clue::solve_LSAP()` to solve the linear assignment problem.
*   **Acceptance Criteria:**
    *   Function returns a permutation vector indicating the final node-to-node mapping.
    *   The cost matrix calculation is correct.
    *   Handles the case where `clue` is not installed gracefully.

**Ticket #6: [TASK] Wrapper & Documentation: `grasp_pairwise()`**
*   **Description:** Combine the functions from tickets #2-5 into a single, user-friendly wrapper function `grasp_pairwise(g1, g2, k=30, q=100, mu=0.1)`.
    1.  Add clear Roxygen2 documentation explaining all parameters, the return value, and a simple usage example.
    2.  Include parameter validation (e.g., graphs must have the same number of nodes, `k` must be less than `n`).
*   **Acceptance Criteria:**
    *   Function can be called with just two `igraph` objects.
    *   Default parameters are sensible.
    *   `?grasp_pairwise` displays helpful documentation.

---

### **Sprint 2: Multiset Extension and Enhancements**

**Sprint Goal:** Extend the core engine to handle an arbitrary number of graphs. Implement the "hub-and-spoke" model for multiset alignment. Add key B-GRASP enhancements (PPR descriptors and voting) to improve robustness.

**Tickets for Sprint 2:**

**Ticket #7: [EPIC] Extend GRASP to Multiset Alignment**

---

**Ticket #8: [TASK] Enhancement: Implement PPR Descriptors**
*   **Description:** Create an alternative descriptor function `compute_ppr_descriptors(graph, q)`.
    1.  Use `igraph::page_rank()` to compute Personalized PageRank.
    2.  Define `q` "source" node sets (e.g., top `n/q` nodes by degree, etc.).
    3.  For each source set, compute the PPR vector and use it as a descriptor.
    4.  Modify `grasp_pairwise()` to accept a `descriptor_method` argument (`"heat"` or `"ppr"`).
*   **Acceptance Criteria:**
    *   PPR function is reasonably efficient for sparse graphs.
    *   The main `grasp_pairwise` function can seamlessly switch between descriptor types.

**Ticket #9: [TASK] Multiset: Joint Base Alignment (Hub-and-Spoke)**
*   **Description:** Create a function `align_bases_multiset(basis_list, descriptor_list, mu)` that finds rotations `{M(s)}` mapping all bases to a single latent space.
    1.  Implement the block-coordinate descent algorithm described in the guide.
    2.  In each inner step, fix all but one `M(s)` and solve a simpler optimization problem to update it. The objective for updating `M(s)` will involve aligning its descriptors to the *mean* of all other aligned descriptors.
    3.  The core optimization logic from Ticket #4 can be reused here.
*   **Acceptance Criteria:**
    *   Function takes a list of bases and a list of descriptor matrices as input.
    *   Returns a list of `k x k` orthogonal rotation matrices `{M(s)}`.
    *   Algorithm converges reliably.

**Ticket #10: [TASK] Multiset: Hub-and-Spoke Assignment**
*   **Description:** Create a function `compute_multiset_assignment(basis_list, M_list, C_list)` that generates all pairwise alignments.
    1.  For each graph `s`, compute its node embeddings in the shared latent space: `Embeddings(s) = (Φ(s) M(s)) C(s)`.
    2.  For any pair of graphs `(s, t)`, compute the distance matrix between `Embeddings(s)` and `Embeddings(t)`.
    3.  Solve the linear assignment problem to get the permutation `P_st`.
    4.  Return a list of all pairwise permutation matrices/vectors.
*   **Acceptance Criteria:**
    *   Function correctly computes alignments between any two specified graphs in the set.
    *   The resulting alignments are consistent (i.e., `P_st ≈ P_su %*% P_ut`).

**Ticket #11: [TASK] Wrapper & Documentation: `grasp_multiset()`**
*   **Description:** Combine the multiset components into a top-level function `grasp_multiset(graph_list, k, q, ...)`.
    1.  The function should first call `compute_spectral_basis` and `compute_..._descriptors` for all graphs.
    2.  It then calls `align_bases_multiset` to get the joint rotations.
    3.  It then computes the functional maps `C(s)` for each graph relative to the latent space.
    4.  Finally, it orchestrates the hub-and-spoke assignment.
    5.  Add clear Roxygen2 documentation.
*   **Acceptance Criteria:**
    *   The function takes a list of `igraph` objects and returns a list of pairwise alignments.
    *   The internal flow is logical and follows the hub-and-spoke model.

**Ticket #12: [TASK] Enhancement: Voting Over `k`**
*   **Description:** Add a voting mechanism to the `grasp_pairwise` and `grasp_multiset` functions.
    1.  Modify the wrappers to accept a vector of `k` values (e.g., `k_vals = c(20, 30, 40)`).
    2.  Loop over `k_vals`, performing a full alignment for each.
    3.  Collect the resulting permutations.
    4.  Implement a voting scheme: for each node, choose the match that appeared most frequently across all runs. Use a greedy assignment to resolve ties and ensure a valid one-to-one mapping.
*   **Acceptance Criteria:**
    *   The function can be called with a single `k` (no voting) or multiple `k`s (voting).
    *   The final output is a single, consensus permutation.

# GRASP Implementation Tickets - Revised for manifoldalign Package Consistency

## Overview

Based on analysis of existing manifoldalign package patterns (KEMA, GPCA, lowrank_align, generalized_procrustes), the GRASP implementation has been revised to ensure consistency with established conventions.

## Key Consistency Requirements

### 1. Parameter Naming Conventions
- `ncomp` (not `k`) for number of components
- `preproc` for preprocessing function (default: `center()`)
- `sigma` for diffusion/kernel parameters
- `lambda` for regularization parameters
- `use_laplacian` for Laplacian normalization choice
- `solver` for algorithm variants

### 2. Function Structure Patterns
- S3 methods with generic functions
- `hyperdesign` and `multidesign` method variants
- Consistent parameter validation using `chk` package
- Error handling with `safe_compute` wrapper pattern

### 3. Matrix Operations
- Extensive use of `Matrix` package for sparse operations
- `PRIMME::eigs_sym` for eigenvalue computations
- Pattern: `Matrix::bdiag()` for block diagonal construction
- Consistent sparse matrix validation and conversion

### 4. Return Values
- `multiblock_biprojector` objects for alignment methods
- Include `preproc`, `block_indices`, `labels` metadata
- Standard structure: `v` (loadings), `s` (scores), `sdev` (standard deviations)

### 5. Graph Construction
- Use `neighborweights::graph_weights` and `neighborweights::adjacency`
- Pattern from KEMA: separate `compute_local_similarity` functions
- Sparse graph representations throughout

### 6. Dependencies and Imports
- `igraph` required for graph operations (already in DESCRIPTION)
- Pattern: `RSpectra` (not `rARPACK`) for eigenvalue computations
- `clue` for assignment problems
- Follow `requireNamespace()` patterns for suggested packages

---

## Updated Sprint 1: Core Pairwise GRASP Engine

**Sprint Goal:** Deliver `grasp.hyperdesign()` following manifoldalign conventions.

### Ticket #1: [EPIC] Core GRASP Implementation

**Epic Goal:** Implement GRASP (Graph Alignment by Spectral Corresponding Functions) following manifoldalign package patterns and conventions.

---

### Ticket #2: [TASK] Block A: Spectral Basis Construction (Consistent with KEMA patterns)

**Description:** Create function `compute_grasp_basis(strata, ncomp, use_laplacian=TRUE)` following KEMA's eigenvalue computation patterns.

**Implementation Requirements:**
1. Use `Matrix` package for sparse Laplacian construction
2. Follow KEMA's `normalize_laplacian()` pattern for consistency
3. Use `PRIMME::eigs_sym()` with same parameters as other package methods
4. Handle isolated nodes gracefully (following KEMA pattern)
5. Return sparse matrices consistently

```r
compute_grasp_basis <- function(strata, ncomp, use_laplacian = TRUE) {
  # Follow KEMA pattern: validate input
  if (!is.list(strata) || length(strata) == 0) {
    stop("strata must be a non-empty list", call. = FALSE)
  }
  
  # Construct adjacency matrices using neighborweights pattern
  graphs <- lapply(strata, function(stratum) {
    # Use igraph for graph construction, then convert to Matrix
    adj_matrix <- # ... construct from stratum$x
    if (use_laplacian) {
      normalize_laplacian(adj_matrix)  # Reuse KEMA function
    } else {
      adj_matrix
    }
  })
  
  # Eigenvalue computation using PRIMME (consistent with package)
  bases <- lapply(graphs, function(L) {
    decomp <- safe_compute(
      PRIMME::eigs_sym(L, NEig = ncomp + 1, which = "SA"),
      "Eigenvalue computation failed for GRASP basis"
    )
    
    # Follow KEMA pattern: filter trivial eigenvectors
    list(vectors = decomp$vectors, values = decomp$values)
  })
  
  bases
}
```

**Acceptance Criteria:**
- Function follows `safe_compute` error handling pattern
- Uses `Matrix` package sparse operations throughout
- Eigenvalue computation matches KEMA methodology
- Handles edge cases (disconnected components) like other package functions

---

### Ticket #3: [TASK] Block B: Multi-scale Descriptors (Following KEMA descriptor patterns)

**Description:** Create `compute_grasp_descriptors(bases, q_descriptors, sigma=0.73)` following KEMA's diffusion time parameter patterns.

**Implementation Requirements:**
1. Follow KEMA's `sigma` parameter naming convention
2. Use vectorized Matrix operations for efficiency
3. Pattern after KEMA's heat kernel computations
4. Include numerical stability checks from KEMA

```r
compute_grasp_descriptors <- function(bases, q_descriptors, sigma = 0.73) {
  # Follow KEMA validation patterns
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  
  # Generate time steps (following KEMA diffusion patterns)
  time_steps <- seq(0.1, 50, length.out = q_descriptors) * sigma
  
  descriptors <- lapply(bases, function(basis) {
    phi <- basis$vectors
    lambda_vals <- basis$values
    
    # Vectorized descriptor computation (KEMA-style efficiency)
    desc_matrix <- Matrix(0, nrow = nrow(phi), ncol = q_descriptors, sparse = TRUE)
    
    for (i in seq_len(q_descriptors)) {
      t <- time_steps[i]
      # Follow KEMA numerical stability patterns
      exp_vals <- exp(-t * pmax(lambda_vals, 1e-12))
      desc_matrix[, i] <- Matrix::rowSums(phi^2 %*% Matrix::Diagonal(x = exp_vals))
    }
    
    desc_matrix
  })
  
  descriptors
}
```

**Acceptance Criteria:**
- Parameter naming matches KEMA conventions (`sigma`)
- Uses `chk` package for validation (package pattern)
- Sparse matrix operations throughout
- Numerical stability handling matches package standards

---

### Ticket #4: [TASK] Block C: Base Alignment (Following Procrustes optimization patterns)

**Description:** Create `align_grasp_bases(basis1, basis2, desc1, desc2, lambda=0.1)` following generalized_procrustes optimization patterns.

**Implementation Requirements:**
1. Use `lambda` parameter naming (consistent with package)
2. Follow Procrustes convergence patterns and tolerance settings
3. Matrix operations using `Matrix` package
4. Error handling following package patterns

```r
align_grasp_bases <- function(basis1, basis2, desc1, desc2, lambda = 0.1, 
                             max_iter = 100, tol = 1e-6) {
  # Validation following package patterns
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_number(max_iter)
  chk::chk_number(tol)
  
  # Extract matrices (following KEMA matrix extraction patterns)
  phi1 <- as(basis1$vectors, "dgCMatrix")
  phi2 <- as(basis2$vectors, "dgCMatrix")
  lambda1 <- Matrix::Diagonal(x = basis1$values)
  lambda2 <- Matrix::Diagonal(x = basis2$values)
  
  # Stiefel manifold optimization (following Procrustes patterns)
  ncomp <- ncol(phi1)
  M <- Matrix::Diagonal(ncomp)  # Initialize as identity
  
  for (iter in seq_len(max_iter)) {
    # Objective function (following package optimization patterns)
    obj1 <- Matrix::norm(Matrix::crossprod(M, lambda2 %*% M) - lambda1, "F")
    obj2 <- Matrix::norm(Matrix::crossprod(desc1, phi1) - 
                        Matrix::crossprod(desc2, phi2 %*% M), "F")
    
    # Combined objective
    objective <- obj1 + lambda * obj2
    
    # Gradient computation and Stiefel projection
    # ... (following Procrustes optimization patterns)
    
    # Convergence check (package pattern)
    if (iter > 1 && abs(prev_obj - objective) < tol) {
      break
    }
    prev_obj <- objective
  }
  
  list(rotation = M, iterations = iter, converged = iter < max_iter)
}
```

**Acceptance Criteria:**
- Uses package convergence patterns (`max_iter`, `tol`)
- Matrix operations consistent with package style
- Error handling follows `safe_compute` patterns
- Return structure matches package conventions

---

### Ticket #5: [TASK] Blocks D & E: Functional Map and Assignment (Following KEMA final assignment patterns)

**Description:** Create `compute_grasp_assignment(basis1, basis2, M, method="cosine")` following KEMA's assignment methodology.

**Implementation Requirements:**
1. Distance computation following package patterns
2. Use similar assignment algorithms as other methods
3. Sparse matrix efficiency throughout
4. Error handling consistency

```r
compute_grasp_assignment <- function(basis1, basis2, M, method = "cosine") {
  # Validation (package pattern)
  if (!method %in% c("cosine", "euclidean")) {
    stop("method must be 'cosine' or 'euclidean'", call. = FALSE)
  }
  
  # Block D: Functional correspondence (following KEMA coefficient computation)
  phi1 <- basis1$vectors
  phi2_aligned <- basis2$vectors %*% M
  
  # Compute diagonal mapping matrix C (efficient closed form)
  C_diag <- Matrix::diag(Matrix::crossprod(phi2_aligned, phi1)) / 
           Matrix::diag(Matrix::crossprod(phi2_aligned))
  C <- Matrix::Diagonal(x = C_diag)
  
  # Block E: Final embeddings and assignment
  embed1 <- phi1
  embed2 <- phi2_aligned %*% C
  
  # Distance computation (following package patterns)
  if (method == "cosine") {
    # Use coop package (already imported in package)
    cost_matrix <- 1 - coop::cosine(Matrix::t(embed1), Matrix::t(embed2))
  } else {
    cost_matrix <- as.matrix(dist(rbind(as.matrix(embed1), as.matrix(embed2))))
    cost_matrix <- cost_matrix[1:nrow(embed1), (nrow(embed1)+1):ncol(cost_matrix)]
  }
  
  # Linear assignment (following package dependency patterns)
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue package is required for assignment computation. ",
         "Please install it with: install.packages('clue')", call. = FALSE)
  }
  
  assignment <- clue::solve_LSAP(cost_matrix)
  
  list(assignment = as.integer(assignment), 
       cost_matrix = cost_matrix,
       mapping_matrix = C)
}
```

**Acceptance Criteria:**
- Distance computation uses package-imported functions (`coop`)
- Assignment follows package dependency patterns
- Return structure consistent with package conventions
- Error messages match package style

---

### Ticket #6: [TASK] Main Interface: `grasp.hyperdesign()` (Following package S3 method patterns)

**Description:** Create main S3 method following exact patterns from KEMA, GPCA, and lowrank_align.

**Implementation Requirements:**
1. Exact parameter naming conventions from package
2. S3 method structure matching existing methods
3. Return `multiblock_biprojector` object like other methods
4. Documentation following package Roxygen2 patterns

```r
#' Graph Alignment by Spectral Corresponding Functions (GRASP)
#'
#' Performs GRASP alignment on hyperdesign data structures. 
#' Projects graph data from multiple domains into aligned spectral spaces.
#'
#' @param data A hyperdesign object containing multiple graph domains
#' @param preproc Preprocessing function to apply to the data (default: center())
#' @param ncomp Number of spectral components to extract (default: 30)
#' @param q_descriptors Number of diffusion-time descriptors (default: 100)
#' @param sigma Diffusion parameter for descriptor computation (default: 0.73)
#' @param lambda Regularization parameter for base alignment (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Assignment method: "linear" for exact assignment (default)
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the GRASP alignment
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign graph data
#' library(multidesign)
#' 
#' # Create synthetic graph domains
#' set.seed(123)
#' domain1 <- list(
#'   x = matrix(rnorm(100), 50, 2),  # Node features
#'   design = data.frame(node_id = 1:50)
#' )
#' domain2 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(node_id = 1:50)  
#' )
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(domain1 = domain1, domain2 = domain2))
#' 
#' # Run GRASP alignment
#' result <- grasp(hd, ncomp = 20, q_descriptors = 50)
#' }
#'
#' @export
#' @importFrom multivarious init_transform
grasp.hyperdesign <- function(data, 
                             preproc = center(), 
                             ncomp = 30,
                             q_descriptors = 100,
                             sigma = 0.73,
                             lambda = 0.1,
                             use_laplacian = TRUE,
                             solver = "linear",
                             ...) {
  
  # Input validation (following package patterns)
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_logical(use_laplacian)
  
  # Validate solver parameter (following KEMA pattern)
  if (!solver %in% c("linear", "hungarian")) {
    stop("solver must be either 'linear' or 'hungarian'", call. = FALSE)
  }
  
  # Validate input data (following package pattern)
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty list of hyperdesign objects", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("GRASP currently supports exactly 2 domains. For multiset alignment, ", 
         "use grasp_multiset() [coming in Sprint 2]", call. = FALSE)
  }
  
  # Preprocess data (following exact package pattern)
  pdata <- safe_compute(
    multivarious::init_transform(data, preproc),
    "Data preprocessing failed. Check preprocessing function or data format."
  )
  
  # Extract preprocessing information (package pattern)
  proclist <- attr(pdata, "preproc")
  names(proclist) <- names(pdata)
  
  # Block indices computation (following package pattern)
  block_indices <- block_indices(pdata)
  proc <- multivarious::concat_pre_processors(proclist, 
                                             split(block_indices, row(block_indices)))
  names(block_indices) <- names(pdata)
  
  # Call GRASP fitting function
  grasp_fit(pdata, proc, ncomp, q_descriptors, sigma, lambda, 
           use_laplacian, solver, block_indices)
}

#' @keywords internal
grasp_fit <- function(strata, proc, ncomp, q_descriptors, sigma, lambda, 
                     use_laplacian, solver, block_indices) {
  
  # Block A: Spectral basis construction
  bases <- compute_grasp_basis(strata, ncomp, use_laplacian)
  
  # Block B: Multi-scale descriptors  
  descriptors <- compute_grasp_descriptors(bases, q_descriptors, sigma)
  
  # Block C: Base alignment
  alignment_result <- align_grasp_bases(bases[[1]], bases[[2]], 
                                       descriptors[[1]], descriptors[[2]], lambda)
  
  # Blocks D & E: Functional map and assignment
  assignment_result <- compute_grasp_assignment(bases[[1]], bases[[2]], 
                                               alignment_result$rotation, "cosine")
  
  # Compute scores (following package pattern)
  embed1 <- bases[[1]]$vectors
  embed2 <- bases[[2]]$vectors %*% alignment_result$rotation %*% assignment_result$mapping_matrix
  scores <- rbind(embed1, embed2)
  
  # Primal vectors computation (following KEMA pattern)
  v <- do.call(rbind, lapply(seq_along(strata), function(i) {
    xi <- strata[[i]]$x
    alpha_i <- bases[[i]]$vectors
    Matrix::crossprod(xi, alpha_i)
  }))
  
  # Return multiblock_biprojector (exact package pattern)
  multivarious::multiblock_biprojector(
    v = v,
    s = scores,
    sdev = apply(scores, 2, sd),
    preproc = proc,
    block_indices = block_indices,
    assignment = assignment_result$assignment,
    rotation = alignment_result$rotation,
    mapping_matrix = assignment_result$mapping_matrix,
    classes = "grasp"
  )
}

#' @export
grasp <- function(data, ...) {
  UseMethod("grasp")
}
```

**Acceptance Criteria:**
- Exact parameter naming matches package conventions
- S3 method structure identical to existing methods
- Returns `multiblock_biprojector` with standard components
- Documentation follows package Roxygen2 patterns
- Error handling and validation matches package style

---

## Updated Sprint 2: Multiset Extension

### Ticket #7: [EPIC] Multiset GRASP Extension

**Epic Goal:** Extend GRASP to handle arbitrary number of domains following hub-and-spoke pattern, consistent with package multiset handling approaches.

### Ticket #8: [TASK] Enhanced Descriptors (Following neighborweights integration patterns)

**Description:** Add PPR descriptors using `neighborweights` functions already imported by the package.

**Implementation Requirements:**
1. Use existing `neighborweights::graph_weights` for consistency
2. Follow package's graph construction patterns
3. Add `descriptor_method` parameter to main function

```r
compute_ppr_descriptors <- function(strata, q_descriptors, sigma = 0.73) {
  # Use neighborweights pattern from KEMA
  descriptors <- lapply(strata, function(stratum) {
    # Construct graph using package patterns
    graph_weights <- neighborweights::graph_weights(
      stratum$x,
      weight_mode = "normalized",
      neighbor_mode = "knn", 
      k = 5,  # Default from KEMA
      type = "normal",
      sigma = sigma
    )
    
    # Extract adjacency matrix
    adj_matrix <- neighborweights::adjacency(graph_weights)
    
    # Compute PPR descriptors using igraph
    # ... (PPR computation following package efficiency patterns)
  })
  
  descriptors
}
```

### Ticket #9: [TASK] Multiset Joint Alignment (Following generalized_procrustes patterns)

**Description:** Extend base alignment to multiple domains using block coordinate descent, following generalized_procrustes optimization patterns.

### Ticket #10: [TASK] Hub-and-Spoke Assignment (Following package multiset patterns)

**Description:** Implement multiset assignment following the hub-and-spoke model used in other package methods.

### Ticket #11: [TASK] Main Interface: `grasp_multiset()` (Following package conventions)

**Description:** Create multiset wrapper following exact patterns from package multiset methods.

### Ticket #12: [TASK] Voting Enhancement (Following package robustness patterns)

**Description:** Add k-voting mechanism following package robustness enhancement patterns.

---

## Integration Notes

1. **DESCRIPTION Updates:** Add `clue` to Suggests field
2. **NAMESPACE Updates:** Add appropriate imports following package patterns
3. **Documentation:** Follow exact Roxygen2 patterns from existing methods
4. **Testing:** Follow package testing patterns with synthetic data
5. **Error Messages:** Use exact error message format from package

This revised implementation ensures GRASP integrates seamlessly with manifoldalign while maintaining its algorithmic integrity and following all established package conventions.