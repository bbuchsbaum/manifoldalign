# Audit Action Items: manifoldalign Package

## Critical Issues (Implement Immediately)

### 1. Interface Fixes
- [ ] Add `kema.default()` method with helpful error message
- [ ] Add `generalized_procrustes.default()` method
- [ ] Fix generalized_procrustes interface inconsistency
- [ ] Add missing default methods for all generics

### 2. Documentation Fixes
- [ ] Update KEMA examples to use proper hyperdesign constructors
- [ ] Document TV penalty limitations prominently in FPGW
- [ ] Complete documentation for all C++ functions in RcppExports.R
- [ ] Fix C++ style comment in genprocrustes.R line 164

### 3. Testing Gaps
- [ ] Add comprehensive test suite for low_rank_alignment
- [ ] Implement skipped hyperdesign tests in test-generalized_procrustes.R
- [ ] Add error condition tests for KEMA
- [ ] Fix expect_no_warning syntax in test-fpgw-validation.R

## High Priority (Week 1-2)

### 1. Parameter Standardization
- [ ] Standardize coupling parameter: mu_coupling → mu across all methods
- [ ] Standardize max iterations: maxit → max_iter
- [ ] Standardize components: ncomp across all methods
- [ ] Create parameter naming guide

### 2. Code Organization
- [ ] Split all_generic.R into smaller focused files
- [ ] Extract long functions in kema_fit() into sub-functions
- [ ] Extract helper functions in generalized_procrustes
- [ ] Create validation_utils.R with common patterns

### 3. Common Utilities
- [ ] Create `validate_hyperdesign()` function
- [ ] Create `extract_data_from_hyperdesign()` utility
- [ ] Add progress reporting utilities
- [ ] Standardize matrix validation patterns

## Medium Priority (Month 1)

### 1. Algorithm Improvements
- [ ] Make KEMA lambda scaling optional (add parameter)
- [ ] Make KEMA regression quality thresholds configurable
- [ ] Improve cone_align_multiple to use joint optimization
- [ ] Add more solver options to GRASP

### 2. Documentation Enhancement
- [ ] Add `@family` tags to group related methods
- [ ] Create unified parameter tuning guide
- [ ] Document semi-supervised capabilities better
- [ ] Add method selection flowchart

### 3. Testing Enhancement
- [ ] Tighten numerical tolerances in tests
- [ ] Add performance benchmarks
- [ ] Add cross-method consistency tests
- [ ] Implement out-of-sample reconstruction tests

## Low Priority (Quarter 1)

### 1. Code Quality
- [ ] Standardize comment styles
- [ ] Remove all TODO/FIXME comments
- [ ] Apply consistent code formatting
- [ ] Remove legacy mode code to vignette

### 2. Advanced Features
- [ ] Add GPU acceleration support
- [ ] Implement automatic parameter selection
- [ ] Add cross-validation utilities
- [ ] Create Shiny app for method comparison

### 3. Documentation
- [ ] Create comprehensive developer guide
- [ ] Add troubleshooting vignette
- [ ] Create performance tuning guide
- [ ] Document internal design patterns

## File-Specific Actions

### kema.R
- [ ] Add kema.default() method
- [ ] Extract graph construction logic
- [ ] Make lambda scaling optional
- [ ] Update examples

### genprocrustes.R
- [ ] Fix interface inconsistency
- [ ] Fix C++ style comment
- [ ] Extract helper functions
- [ ] Add memory checks

### grasp.R / cone_align.R / parrot.R
- [ ] Standardize error handling
- [ ] Add progress indicators
- [ ] Document parameter interactions
- [ ] Improve multi-graph versions

### fpgw.R
- [ ] Document TV penalty as non-functional
- [ ] Consider deprecating TV penalty
- [ ] Add marginals parameter
- [ ] Improve line search

### coupled_diagonalization.R
- [ ] Standardize parameter names
- [ ] Add numerical stability checks
- [ ] Improve error messages
- [ ] Add cross-validation

### utils.R
- [ ] Add validation utilities
- [ ] Add hyperdesign utilities
- [ ] Add progress utilities
- [ ] Document all utilities

## Testing Checklist

### For Each Method
- [ ] Basic functionality test
- [ ] Parameter validation test
- [ ] Error condition test
- [ ] Edge case test
- [ ] Numerical validation test
- [ ] Performance benchmark
- [ ] Cross-method consistency

## Documentation Checklist

### For Each Exported Function
- [ ] Complete @param descriptions
- [ ] Mathematical formulation in @details
- [ ] Working @examples
- [ ] @seealso cross-references
- [ ] @references to papers
- [ ] @family groupings

## Code Quality Checklist

- [ ] No magic numbers
- [ ] Clear variable names
- [ ] Consistent indentation
- [ ] No dead code
- [ ] Proper error messages
- [ ] Input validation
- [ ] Memory efficiency
- [ ] Numerical stability