# Pseudolabeling Module Improvement TODO

Based on comprehensive technical audit of `assign_pseudolabels()` and related functions.

## 1. FUNCTIONAL CORRECTNESS & NUMERICAL ISSUES

### 🟢 HIGH FEASIBILITY (Easy fixes) - **PHASE 1 COMPLETED** ✅

- [x] **Quantile validation in `.threshold_similarity()`** ✅
  - **Issue**: Silent failure when `q` outside [0,1]
  - **Fix**: Add `if (q <= 0 || q >= 1) stop("q must be in (0,1)")`
  - **Effort**: 5 minutes
  - **Impact**: Prevents silent NA threshold bugs

- [x] **Use `components$no` instead of `max(cluster_assignments)`** ✅
  - **Issue**: Incorrect cluster count if component IDs non-contiguous
  - **Fix**: Replace `max(cluster_assignments, na.rm = TRUE)` with `components$no`
  - **Effort**: 2 minutes
  - **Impact**: Accurate reporting

- [x] **Fix edge buffer overflow guidance** ✅
  - **Issue**: Contradictory advice to "increase k" when buffer grows linearly with k
  - **Fix**: Better error message or dynamic reallocation
  - **Effort**: 10 minutes
  - **Impact**: Better UX

- [x] **Complete namespace qualification** ✅
  - **Issue**: Functions like `sample()`, `sprintf()`, `which.max()` not namespaced
  - **Fix**: Add `base::` prefix for base R functions
  - **Effort**: 15 minutes
  - **Impact**: Prevents conflicts, cleaner imports

- [x] **Optimized sparse thresholding** ✅
  - **Issue**: `sim_matrix >= threshold` creates 8x larger lgCMatrix
  - **Fix**: In-place modification of `@x` slot
  - **Effort**: 30 minutes
  - **Impact**: 8x memory reduction during thresholding

- [x] **Efficient symmetry enforcement with `forceSymmetric()`** ✅
  - **Issue**: `(sim_matrix + t(sim_matrix))/2` doubles memory transiently
  - **Fix**: Use `Matrix::forceSymmetric(sim_matrix, uplo = "U")` (Matrix ≥ 1.6)
  - **Effort**: 15 minutes + testing
  - **Impact**: 50% memory reduction for large sparse matrices

- [x] **Use `sum()` for sparse diagonal operations** ✅
  - **Issue**: `diag(matrix) <- 0` traverses all zeros 
  - **Fix**: Use summary() to access only non-zero entries
  - **Effort**: 10 minutes
  - **Impact**: Faster diagonal operations

- [x] **Improved input validation with clear messages** ✅
  - **Issue**: Cryptic error messages for invalid inputs
  - **Fix**: Add descriptive errors with suggested remedies
  - **Effort**: 20 minutes
  - **Impact**: Better user experience

### 🟡 MEDIUM FEASIBILITY (Moderate effort) - **PHASE 2 COMPLETED** ✅

- [x] **Improve average similarity computation** ✅
  - **Issue**: `sum(cluster_sims)/(m*(m-1))` treats missing edges as 0
  - **Fix**: Expose both "observed mean" and "density-adjusted mean"
  - **Effort**: 45 minutes
  - **Impact**: Better cluster quality assessment
  - **Enhancement**: Added `observed_similarity` and `edge_density` columns to cluster_info

- [x] **Annoy version compatibility check** ✅
  - **Issue**: Distance formula assumes specific Annoy version behavior
  - **Fix**: Add version check or distance validation
  - **Effort**: 30 minutes
  - **Impact**: Prevents false positives/negatives
  - **Enhancement**: Added automatic distance validation with fallback formulas

- [x] **Return S3 class with methods** ✅
  - **Issue**: Raw list return provides poor user experience
  - **Fix**: `structure(..., class = "pseudolabels")` with `print()`, `summary()`, `as.data.frame()`
  - **Effort**: 2-3 hours
  - **Impact**: Professional interactive experience
  - **Enhancement**: Full S3 class with comprehensive methods

- [x] **Continuous diversity weighting** ✅
  - **Issue**: Binary switch `diversity_weight > 0` vs `== 0` 
  - **Fix**: Accept α ∈ [0,1] mixing diversity and confidence
  - **Effort**: 1-2 hours
  - **Impact**: More nuanced representative selection
  - **Enhancement**: Weighted combination of diversity and cluster confidence scores

### 🔴 LOW FEASIBILITY (Complex changes)

- [ ] **Split oversized connected components**
  - **Issue**: No mechanism to handle `max_cluster_size` violations
  - **Fix**: Recursive re-clustering (spectral split) or sampling
  - **Effort**: 2-4 hours
  - **Impact**: Better handling of "blob" anchors
  - **Complexity**: Requires spectral clustering implementation

## 2. ALGORITHMIC & SCALABILITY IMPROVEMENTS

### 🟡 MEDIUM FEASIBILITY

- [ ] **Replace igraph with union-find for giant graphs (n ≥ 50k)**
  - **Issue**: igraph incurs >2x memory overhead for large graphs
  - **Fix**: Direct union-find on CSR arrays or use tidygraph
  - **Effort**: 2-3 hours
  - **Impact**: Major memory reduction for large problems
  - **Implementation**: Either Rcpp union-find or tidygraph integration

- [x] **Improve diversity selection algorithm** ✅
  - **Issue**: O(k²) complexity in candidate selection
  - **Fix**: MaxMin sampling with priority queue or farthest-first traversal
  - **Effort**: 3-4 hours
  - **Impact**: Faster diversity selection
  - **Enhancement**: Implemented advanced farthest-first traversal (FFS) algorithm with automatic algorithm selection based on problem size (≥50 candidates). Uses efficient distance tracking and weighted combination of diversity + confidence scores.

- [x] **Continuous diversity weighting** ✅ *(Moved to Phase 2 - Already Completed)*
  - **Issue**: Binary switch `diversity_weight > 0` vs `== 0` 
  - **Fix**: Accept α ∈ [0,1] mixing diversity and confidence
  - **Effort**: 1-2 hours
  - **Impact**: More nuanced representative selection
  - **Enhancement**: Completed in Phase 2 with weighted combination formula

### 🔴 LOW FEASIBILITY (Major changes)

- [ ] **Adaptive thresholding improvements**
  - **Issue**: Single global quantile inadequate for heterogeneous data
  - **Fix**: Per-sample k-NN radius or per-domain quantiles
  - **Effort**: 4-6 hours
  - **Impact**: Better anchor quality for mixed data
  - **Note**: We have `neighborweights` which offers mutual k-NN

- [ ] **Parallelization of major loops**
  - **Issue**: Row-wise operations not parallelized
  - **Fix**: RcppParallel or future.apply for n ≥ 5k
  - **Effort**: 3-5 hours
  - **Impact**: Significant speed improvements

- [ ] **GPU/torch integration**
  - **Issue**: CPU-only implementation limits scale
  - **Fix**: Optional Faiss integration for embeddings
  - **Effort**: 8-12 hours
  - **Impact**: 10-50x speedup for 100k+ samples
  - **Complexity**: Major dependency and code restructuring

## 3. CODE QUALITY & API IMPROVEMENTS

### 🟢 HIGH FEASIBILITY

- [x] **Complete namespace qualification** ✅
  - **Fix**: Always use `Matrix::`, `stats::`, etc.
  - **Effort**: 30 minutes
  - **Impact**: Better reproducibility

- [x] **Enhanced documentation clarity** ✅
  - **Fix**: Clarify that `sim_matrix` must contain all relevant similarities
  - **Effort**: 20 minutes
  - **Impact**: Prevent user confusion

### 🟡 MEDIUM FEASIBILITY

- [ ] **Return S3 class with methods**
  - **Fix**: `structure(..., class = "pseudolabels")` with `print()`, `plot()`, `as.data.frame()`
  - **Effort**: 2-3 hours
  - **Impact**: Better interactive experience

- [ ] **Replace verbose with cli logging**
  - **Fix**: Use `cli::cli_inform()` for colored, suppressible progress
  - **Effort**: 1 hour
  - **Impact**: Professional logging experience

- [ ] **Comprehensive unit tests**
  - **Fix**: Add testthat specs for edge cases
  - **Effort**: 3-4 hours
  - **Impact**: Robust edge case handling
  - **Tests needed**: 
    - Zero-edge matrices
    - Diagonal-only matrices
    - Huge sparse matrices
    - Symmetric-but-unsorted CSR

## 4. MICRO-OPTIMIZATIONS (Drop-in fixes)

### 🟢 HIGH FEASIBILITY

- [x] **Faster in-place thresholding** ✅
  - **Fix**: Implemented efficient sparse thresholding without intermediate matrices
  - **Effort**: 15 minutes
  - **Impact**: No memory duplication

- [x] **Pre-allocation optimizations** ✅
  - **Fix**: Better initial size estimates for edge vectors with auto-expansion
  - **Effort**: 30 minutes
  - **Impact**: Reduced memory allocations

## IMPLEMENTATION PRIORITY

### Phase 1: Quick Wins (1-2 days) - **COMPLETED** ✅
1. ✅ Quantile validation
2. ✅ Component count fix  
3. ✅ Namespace qualification
4. ✅ Optimized sparse thresholding
5. ✅ Better error messages
6. ✅ Matrix::forceSymmetric() for symmetry
7. ✅ Sparse diagonal operations optimization
8. ✅ Edge buffer auto-expansion

### Phase 2: Medium Impact (1 week) - **COMPLETED** ✅
1. ✅ forceSymmetric() integration (moved to Phase 1)
2. ✅ Improved similarity computation with observed/density metrics
3. ✅ S3 class with comprehensive methods (print, summary, as.data.frame)
4. ✅ Enhanced documentation
5. ✅ Continuous diversity weighting implementation

### Phase 3: Major Improvements (2-3 weeks) - **PARTIALLY COMPLETED** ✅
1. Union-find for large graphs
2. ✅ Diversity selection improvements
3. ✅ Continuous diversity weighting (moved to Phase 2)
4. Parallelization
5. Comprehensive test suite

### Phase 4: Advanced Features (1+ months)
1. Adaptive thresholding with neighborweights
2. GPU integration
3. Advanced clustering algorithms
4. Performance benchmarking suite

## DEPENDENCIES TO CONSIDER

- **Matrix package version**: forceSymmetric() requires Matrix ≥ 1.6
- **neighborweights integration**: Leverage existing mutual k-NN capabilities
- **cli package**: For enhanced logging
- **RcppParallel**: For parallelization
- **testthat**: For comprehensive testing

## ESTIMATED TOTAL EFFORT

- **Phase 1**: ~2 developer days (immediate fixes)
- **Phase 2**: ~5 developer days (medium improvements) 
- **Phase 3**: ~15 developer days (major algorithmic improvements)
- **Phase 4**: ~30+ developer days (advanced features)

**Recommendation**: Start with Phase 1 quick wins for immediate impact, then evaluate Phase 2 based on user feedback and performance needs. 