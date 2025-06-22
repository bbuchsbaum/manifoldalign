# Method Comparison Table: manifoldalign Package

## Algorithm Overview

| Method | Type | Key Innovation | Use Case | Complexity |
|--------|------|----------------|----------|------------|
| **KEMA** | Kernel/Manifold | Kernel embeddings with manifold preservation | Semi-supervised multi-domain alignment | O(n²) kernel, O(n³) eigen |
| **Generalized Procrustes** | Orthogonal | Handles partial task observations | Task-based alignment with missing data | O(n²L) per iteration |
| **GRASP** | Graph/Spectral | Functional maps with multi-scale descriptors | Near-isomorphic graph matching | O(n³) + O(n³) assignment |
| **CONE-Align** | Graph/Iterative | Alternating optimization on embeddings | General graph alignment | O(n² × iterations) |
| **PARROT** | Graph/Transport | Position-aware features with optimal transport | Semi-supervised network alignment | O(n² × iterations) |
| **Gromov-Wasserstein** | Optimal Transport | Structural optimal transport | Domain adaptation without correspondence | O(n² × iterations) |
| **FPGW** | Optimal Transport | Combines features and structure | Partial transport with mixed costs | O(n² × iterations) |
| **Coupled Diagonalization** | Spectral | Joint diagonalization with coupling | Multi-modal manifold analysis | O(n³) + Stiefel optimization |
| **GPCA Align** | Linear/PCA | Generalized PCA with constraints | Linear multi-domain alignment | O(n³) eigendecomposition |
| **Linear Similarity Embedding** | Linear/Similarity | Similarity-preserving linear maps | Fast alignment with known similarities | O(n³) or iterative |

## Common Parameters

| Parameter | KEMA | Procrustes | GRASP | CONE | PARROT | GW | FPGW | Coupled Diag | GPCA | LSE |
|-----------|------|------------|-------|------|--------|----|----|--------------|------|-----|
| `ncomp` | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ |
| `preproc` | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ | ✓ | ✗ |
| `max_iter` | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| `tol` | ✗ | ✓ | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✓ |
| `sigma` | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ | ✗ | ✓ | ✗ | ✗ |
| `lambda` | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ | ✓ | ✗ | ✗ | ✓ |
| `solver` | ✓ | ✗ | ✓ | ✓ | ✓ | ✗ | ✗ | ✗ | ✗ | ✓ |
| `verbose` | ✓ | ✓ | ✗ | ✗ | ✓ | ✓ | ✓ | ✗ | ✗ | ✓ |

## Input/Output Specifications

| Method | Input Type | Primary Output | Additional Outputs |
|--------|------------|----------------|-------------------|
| **KEMA** | hyperdesign/multidesign | multiblock_biprojector | alpha coefficients, kernel info |
| **Generalized Procrustes** | hyperdesign/list | list | O_mats, A_est, convergence info |
| **GRASP** | hyperdesign | multiblock_biprojector | assignment, rotation, mapping |
| **CONE-Align** | hyperdesign | multiblock_biprojector | assignment, rotation |
| **PARROT** | hyperdesign | multiblock_biprojector | transport_plan, alignment_matrix |
| **Gromov-Wasserstein** | hyperdesign | gromov_wasserstein | transport_plans, distances |
| **FPGW** | hyperdesign | fpgw | transport_plans, distances |
| **Coupled Diagonalization** | hyperdesign | multiblock_biprojector | coupled_bases |
| **GPCA Align** | hyperdesign | gpca_aligned | loadings, scores |
| **Linear Similarity Embedding** | matrices + similarity | list | W, objective_values |

## Method Selection Guide

### By Data Type
- **Labeled data**: KEMA (semi-supervised), Generalized Procrustes (task-based)
- **Graph data**: GRASP, CONE-Align, PARROT
- **General domains**: Gromov-Wasserstein, FPGW
- **Multi-modal**: Coupled Diagonalization
- **Linear assumptions**: GPCA Align, Linear Similarity Embedding

### By Problem Type
- **Known correspondences**: KEMA, GPCA Align, Linear Similarity Embedding
- **Partial correspondences**: PARROT (anchors), Generalized Procrustes (tasks)
- **No correspondences**: GRASP, CONE-Align, Gromov-Wasserstein
- **Partial transport**: FPGW with rho parameter
- **Semi-supervised**: KEMA, PARROT

### By Scale
- **Small (n < 1000)**: All methods feasible
- **Medium (1000 < n < 10000)**: Prefer KEMA with REKEMA, iterative methods
- **Large (n > 10000)**: Linear methods, REKEMA, careful parameter tuning

## Performance Characteristics

| Method | Time Complexity | Space Complexity | Scalability | Parallelizable |
|--------|----------------|------------------|-------------|----------------|
| **KEMA** | O(n²k + n³) | O(n²) | Medium (REKEMA helps) | Partially |
| **Generalized Procrustes** | O(I × n²L) | O(nL) sparse | High | Yes |
| **GRASP** | O(n³) | O(n²) | Low | Eigendecomp only |
| **CONE-Align** | O(I × n²) | O(n²) | Medium | No |
| **PARROT** | O(I × n²) | O(n²) | Medium | C++ accelerated |
| **Gromov-Wasserstein** | O(I × n²) | O(n²) | Medium | Partially |
| **FPGW** | O(I × n²) | O(n²) | Medium | Partially |
| **Coupled Diagonalization** | O(n³ + I×n²k) | O(nk) | Low | Partially |
| **GPCA Align** | O(n³) | O(n²) | Low | SVD only |
| **Linear Similarity Embedding** | O(n³) or O(I×n²) | O(n²) | Medium | Solver dependent |

## Robustness Features

| Method | Handles Missing Data | Noise Robust | Initialization | Convergence Guarantee |
|--------|---------------------|--------------|----------------|---------------------|
| **KEMA** | ✓ (NA labels) | ✓ | Auto-tuned | Local (global for convex) |
| **Generalized Procrustes** | ✓ (partial tasks) | ✓ | SVD/random fallback | Local minimum |
| **GRASP** | ✗ | ✓ | Spectral | Local minimum |
| **CONE-Align** | ✗ | Moderate | Spectral | Local minimum |
| **PARROT** | ✓ (partial anchors) | ✓ | Warm-start | Local minimum |
| **Gromov-Wasserstein** | ✗ | Moderate | Uniform | Local minimum |
| **FPGW** | ✓ (partial transport) | ✓ | Entropic warm-start | Local minimum |
| **Coupled Diagonalization** | ✓ (partial corresp.) | ✓ | Spectral | Local minimum |
| **GPCA Align** | ✗ | Moderate | PCA-based | Global (convex) |
| **Linear Similarity Embedding** | ✗ | Depends on similarity | Identity/learned | Depends on solver |

## Unique Features

| Method | Unique Capability |
|--------|-------------------|
| **KEMA** | Automatic kernel selection, regression solver for speed |
| **Generalized Procrustes** | Handles arbitrary partial observation patterns |
| **GRASP** | Multi-scale heat diffusion descriptors |
| **CONE-Align** | Simple and effective for general graphs |
| **PARROT** | Position-aware features via Random Walk with Restart |
| **Gromov-Wasserstein** | Pure structural alignment without features |
| **FPGW** | Flexible feature/structure trade-off, mass constraints |
| **Coupled Diagonalization** | Simultaneous multi-modal analysis |
| **GPCA Align** | Interpretable linear transformations |
| **Linear Similarity Embedding** | Direct optimization of similarity preservation |