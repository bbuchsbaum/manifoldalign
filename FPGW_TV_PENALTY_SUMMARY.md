# FPGW TV Penalty Implementation Summary

## Overview
This document summarizes the work done to fix the TV penalty implementation in the Fused-Partial Gromov-Wasserstein (FPGW) algorithm.

## Problem Statement
The TV-penalized partial optimal transport was showing inverted behavior - higher λ values increased transported mass instead of reducing it, which was opposite to the expected sparsity-inducing behavior.

## Investigation Process

### 1. Initial Consultation (Gemini & O3)
Both AI models identified that the augmentation method was incorrect:
- **Gemini**: Suggested λ should be added to transport costs with zero escape costs
- **O3**: Identified misplaced dummy penalty in the augmented formulation
- **Consensus**: The implementation was penalizing NOT transporting mass instead of penalizing transported mass

### 2. Mathematical Analysis
The user provided a crucial insight: the TV penalty in the paper is actually:
```
λ * (|μ|² + |ν|² - 2|γ|²)
```
This is a **smooth quadratic term**, not an L1 penalty. The augmentation method is fundamentally incompatible with this formulation.

### 3. Implementation Fix
Modified `fw_fpgw()` to:
- Use standard OT oracle for both classical and TV-penalized cases
- Add TV penalty contributions to step size computation
- Include quadratic TV penalty in objective function

Key code changes:
```r
# Step size modification
if (lambda > 0) {
  mass_P <- sum(P)
  mass_dP <- sum(dP)
  a <- a + 2 * lambda * mass_dP^2      # quadratic term
  b <- b - 4 * lambda * mass_P * mass_dP  # linear term
}

# Objective computation
if (lambda > 0) {
  current_obj <- current_obj + lambda * (sum(p)^2 + sum(q)^2 - 2 * sum(P)^2)
}
```

## Results

### ✅ What Works
- Implementation is now mathematically correct according to the paper
- No more augmentation-related errors
- Step size computation properly includes TV penalty terms
- Objective function correctly computes the quadratic penalty

### ❌ What Doesn't Work
- Tests still fail - the quadratic penalty doesn't encourage sparsity
- The penalty creates a valley at mass=1, causing convergence to full transport
- The formulation `-2λ|γ|²` is minimized at either |γ|=0 or |γ|=1

## Root Cause
The TV penalty as formulated in the paper is not actually a sparsity-inducing penalty. The quadratic term creates two minima (no transport or full transport), and typical FGW objectives push the solution toward full transport.

## Actions Taken
1. **Documentation Updated**: Added warnings that TV penalty is experimental
2. **Code Comments**: Clarified the mathematical formulation and limitations
3. **Test Status**: Tests remain failing but now we understand why
4. **User Warning**: Added runtime warning when λ > 0

## Recommendations

### For Users
- Use `rho` parameter for mass-constrained transport instead of TV penalty
- TV penalty with λ > 0 is experimental and may not work as expected

### For Developers
1. Contact paper authors about the TV penalty behavior
2. Consider implementing alternative sparsity-inducing penalties:
   - L1 penalty: `λ * ||γ||₁`
   - Entropy regularization
   - KL divergence penalty
3. Update or remove failing TV penalty tests

## Files Modified
- `/Users/bbuchsbaum/code/manifoldalign/R/fpgw.R` - Main implementation
- Created multiple documentation files explaining the issue

## Conclusion
The TV penalty implementation is now mathematically correct but reveals a fundamental limitation in the penalty formulation itself. The quadratic penalty does not reliably encourage sparse transport plans, which appears to be a theoretical issue rather than an implementation bug.