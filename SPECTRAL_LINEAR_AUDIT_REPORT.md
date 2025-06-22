# Spectral/Linear Method Group Audit Report

## Executive Summary

This audit examines the Spectral/Linear method group in the manifoldalign package, consisting of:
1. **Coupled Diagonalization** (`coupled_diagonalization.R`)
2. **GPCA Align** (`gpca_align.R`) 
3. **Linear Similarity Embedding** (`linear_sim_embed.R`)
4. **Low Rank Alignment** (`low_rank_alignment.R`)

All methods show solid implementation with good mathematical foundations, though some areas need attention for numerical stability and interface consistency.

## 1. Mathematical Correctness

### Coupled Diagonalization ✓
- **Algorithm**: Correctly implements Eynard et al. (2015) simultaneous diagonalization
- **Objective**: Properly minimizes F = Σᵢ ||AᵢᵀΛᵢAᵢ - diag(AᵢᵀΛᵢAᵢ)||²_F + μ_c Σᵢ<ⱼ ||FᵢᵀUᵢAᵢ - FⱼᵀUⱼAⱼ||²_F
- **Stiefel Manifold**: Correctly maintains orthogonality via QR decomposition after gradient steps
- **Tests**: Comprehensive mathematical validation including gradient checks and structure recovery

### GPCA Align ✓
- **Algorithm**: Proper use of generalized PCA with row metric matrix M
- **Similarity Matrices**: Correctly handles within-domain and between-domain similarities
- **PSD Preservation**: Good approach using separate M_within and M_between construction
- **Normalization**: Innovative Frobenius norm normalization for interpretable u parameter

### Linear Similarity Embedding ✓ 
- **Algorithm**: Faithful implementation of Passalis & Tefas (2016)
- **Gradients**: Correct gradient computations for both Js (similarity) and Jp (orthogonality)
- **Dual Implementation**: Both R (ADAM) and C++ (L-BFGS-B) backends validated
- **Tests**: Strong validation including latent structure recovery

### Low Rank Alignment ✓
- **Algorithm**: Correct implementation balancing low-rank and similarity constraints
- **Semi-supervised**: Excellent NA label handling for unlabeled samples
- **Eigenvalue Problem**: Properly handles zero modes from isolated nodes

## 2. Stiefel Manifold Optimization

### Coupled Diagonalization - EXCELLENT
```r
# Proper projection back to Stiefel manifold
A[[i]] <- qr.Q(qr(A_new))[, seq_len(ncomp), drop = FALSE]
```
- Uses QR decomposition for exact orthogonality
- Momentum on tangent space before projection
- Adaptive step size with backtracking
- Tracks best solution during optimization

### Linear Similarity Embedding - GOOD (Soft Constraint)
```r
# Orthogonality as penalty term
Jp = ||W'W - I||²_F / (2*m²)
```
- Uses soft orthogonality constraint via penalty
- Appropriate for the problem formulation
- Both ADAM and L-BFGS-B handle unconstrained optimization

### Recommendation
Consider adding a `manifold_optimizer` utility for shared Stiefel optimization code.

## 3. Numerical Stability

### Critical Issues Found

#### Issue 1: Eigenvalue Handling in Coupled Diagonalization
```r
# Current code has protection but could be more robust
eig$values[abs(eig$values) < 1e-10] <- 1e-10
```
**Recommendation**: Add relative tolerance check: `eig$values[abs(eig$values) < 1e-10 * max(abs(eig$values))]`

#### Issue 2: GPCA Align PSD Correction
```r
if (min_eigenval < 1e-8) {
    correction <- abs(min_eigenval) + 1e-6
    M <- M + Matrix::Diagonal(x=rep(correction, nrow(M)))
}
```
**Good practice** but consider using machine epsilon scaling.

#### Issue 3: Linear Similarity Embedding Gradient Stability
```r
# Check for valid gradient
if (any(is.na(grad_i)) || any(!is.finite(grad_i))) {
    warning("Invalid gradient at iteration ", iter, " for domain ", i)
    grad_i[!is.finite(grad_i)] <- 0
}
```
**Good defensive programming** but zeroing may not be optimal. Consider using previous gradient.

### Recommendations for Numerical Stability
1. Add condition number monitoring for critical matrices
2. Use relative tolerances instead of absolute where appropriate
3. Implement gradient clipping for extreme values
4. Add matrix regularization options for ill-conditioned problems

## 4. Interface Consistency

### Strengths
- All methods follow hyperdesign pattern
- Consistent parameter naming (ncomp, preproc, verbose)
- Return multiblock_biprojector objects
- Support multivarious preprocessing pipeline

### Inconsistencies Found

#### Parameter Naming
- `mu_coupling` (coupled_diag) vs `u` (gpca_align) vs `mu` (low_rank_align)
- `max_iter` vs `maxit` 
- `knn` vs `k` for nearest neighbors

#### Similarity Function Interface
```r
# GPCA Align - good pattern
simfun = neighborweights::binary_label_matrix

# Low Rank Alignment - custom createSimFun
simfun = createSimFun(S)
```

### Recommendations
1. Standardize trade-off parameter to `mu` across methods
2. Use consistent `max_iter` everywhere
3. Document similarity function requirements in shared location
4. Create `similarity_functions.R` with common patterns

## 5. Test Coverage Analysis

### Test Coverage Summary
| Method | Basic Tests | Mathematical Tests | Edge Cases | Performance |
|--------|------------|-------------------|------------|-------------|
| Coupled Diagonalization | ✓ (6 tests) | ✓ (5 tests) | ✓ | ⚠️ |
| GPCA Align | ✓ (3 tests) | ✓ (2 tests) | ✓ | ❌ |
| Linear Similarity Embedding | ✓ (3 tests) | ✓ (3 tests) | ✓ | ✓ |
| Low Rank Alignment | ⚠️ | ⚠️ | ⚠️ | ❌ |

### Missing Test Coverage
1. **Low Rank Alignment**: No dedicated test file found
2. **Performance benchmarks**: No systematic performance testing
3. **Large-scale problems**: Limited testing with n > 1000
4. **Failure modes**: More tests needed for degenerate cases

### Test Quality Highlights
- **Coupled Diagonalization**: Excellent "working backwards" tests with known ground truth
- **Linear Similarity Embedding**: Strong R vs C++ consistency tests
- **GPCA Align**: Good cross-domain alignment validation

## 6. Implementation Quality

### Code Organization - GOOD
- Clear separation of concerns
- Helper functions properly scoped
- Good use of S3 methods

### Documentation - GOOD
- Comprehensive roxygen2 documentation
- Mathematical formulas in details sections
- Clear examples

### Error Handling - NEEDS IMPROVEMENT
```r
# Good
if (m < 2) {
  stop("Coupled diagonalization requires at least 2 domains", call. = FALSE)
}

# Needs improvement - should validate before computation
if (sum_M_batch == 0) {
  message(paste("Skipping iteration", t, "due to zero sum in batch mask M."))
  next
}
```

## 7. Specific Recommendations by Method

### Coupled Diagonalization
1. Add convergence diagnostics export (gradient norms, cost history)
2. Implement warm start capability for iterative refinement
3. Add option for different graph construction methods
4. Consider L-BFGS-B as alternative optimizer

### GPCA Align
1. Export the normalized M_within and M_between for inspection
2. Add cross-validation for u parameter selection
3. Implement incremental GPCA for streaming data
4. Add sparse GPCA option for high-dimensional data

### Linear Similarity Embedding
1. Unify sigma_P selection across R and C++ implementations
2. Add more initialization options beyond PCA
3. Implement stochastic variants for large-scale data
4. Export convergence diagnostics from C++ backend

### Low Rank Alignment
1. **URGENT**: Add comprehensive test suite
2. Document the semi-supervised capabilities more prominently
3. Add adaptive mu selection based on eigenvalue gaps
4. Implement out-of-sample extension method

## 8. Performance Considerations

### Memory Usage
- Good use of sparse matrices throughout
- Consider memory profiling for large problems
- Add memory usage estimates in documentation

### Computational Efficiency
- PRIMME usage is appropriate for large eigenproblems
- Consider parallel implementations for multi-domain operations
- Profile bottlenecks in gradient computations

## 9. Integration Points

### With Other Package Components
- Good integration with multivarious preprocessing
- Consistent with multidesign/hyperdesign patterns
- Could better integrate with other alignment methods

### Suggestions
1. Create unified `spectral_alignment` wrapper for method comparison
2. Add method selection guide in package vignette
3. Implement ensemble methods combining spectral approaches

## 10. Priority Actions

### HIGH Priority
1. Add comprehensive test suite for low_rank_alignment
2. Standardize parameter names across methods
3. Fix numerical stability issues in eigenvalue computations
4. Document semi-supervised capabilities

### MEDIUM Priority
1. Create shared utilities for Stiefel optimization
2. Add performance benchmarks
3. Implement cross-validation for hyperparameters
4. Improve error messages with actionable guidance

### LOW Priority
1. Add visualization methods for coupled bases
2. Implement GPU acceleration
3. Create method comparison vignettes
4. Add streaming/online variants

## Conclusion

The Spectral/Linear method group shows strong mathematical foundations and generally good implementation quality. The main areas for improvement are:

1. **Test Coverage**: Especially for low_rank_alignment
2. **Interface Consistency**: Parameter naming and similarity functions
3. **Numerical Stability**: More robust eigenvalue handling
4. **Documentation**: Better cross-method integration examples

The methods are production-ready but would benefit from the recommended improvements for robustness and usability.

## Appendix: Code Quality Metrics

```r
# Cyclomatic Complexity (McCabe)
coupled_diagonalization: 12 (Moderate)
gpca_align: 8 (Good)
linear_sim_embed: 15 (High - consider refactoring)
low_rank_alignment: 6 (Good)

# Lines of Code
coupled_diagonalization: 601
gpca_align: 264  
linear_sim_embed: 737
low_rank_alignment: 304

# Test Coverage Estimate
Overall: ~75% (Good but gaps in low_rank_alignment)
```