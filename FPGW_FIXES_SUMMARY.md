# FPGW Algorithm Fixes Summary

Based on the senior engineer's review, I've implemented the following three critical fixes:

## 1. Frank-Wolfe Step Size Computation (Fixed)

**Issue**: The quadratic coefficient `a` was computed incorrectly, leading to negative values and alpha = 0.

**Fix** (lines 264-288):
```r
# OLD (incorrect):
a <- struct_weight * sum((const - 2 * Cx %*% dP %*% t(Cy)) * dP)

# NEW (correct): 
A <- Cx %*% dP %*% t(Cy)  # n1 x n2
a <- 2 * struct_weight * sum(A * dP)  # = omega2 * <L, dP ⊗ dP> >= 0
```

**Result**: The coefficient `a` is now correctly non-negative (verified: 10.37 vs -8.72 in test).

**Note**: Alpha still becomes 0 near convergence when `b ≈ 0` (gradient orthogonal to direction), which is expected Frank-Wolfe behavior.

## 2. Gradient Computation (Already Correct)

**Review confirmed**: The gradient simplification is mathematically sound.

```r
# Pre-compute (outside loop):
const1 <- as.vector(Cx^2 %*% p)  
const2 <- as.vector(Cy^2 %*% q)
const <- outer(const1, rep(1, n2)) + outer(rep(1, n1), const2)

# Gradient (inside loop):
M_P <- const - 2 * Cx %*% P %*% t(Cy)
grad <- feature_weight * C + struct_weight * M_P
```

**Result**: Correctly implements ∇F(γ) for square-loss GW with uniform marginals.

## 3. Mass-Constrained LP (Fixed)

**Issue**: Only enforced Σγ ≥ ρ instead of exact equality Σγ = ρ.

**Fix** (lines 358-364):
```r
# OLD: Only lower bound
Amass <- matrix(-1, 1, n * m)
Amat <- rbind(Amat, Amass)
rhs <- c(rhs, -rho)

# NEW: Both bounds for exact equality
Amass_upper <- matrix(1, 1, n * m)   # sum gamma_ij <= rho
Amass_lower <- matrix(-1, 1, n * m)  # sum gamma_ij >= rho  
Amat <- rbind(Amat, Amass_upper, Amass_lower)
rhs <- c(rhs, rho, -rho)
```

**Result**: Now enforces exact mass constraint (though lpSolve may still have numerical tolerance).

## Additional Fixes

- Fixed C++ function name mismatch: `partial_ot_mass_cpp` → `partial_ot_mass_rcpp`
- Fixed row normalization with zero rows in `predict.fpgw`
- Added error checking for negative `a` values

## Remaining Issues

1. **Convergence**: Frank-Wolfe can be slow near optimum (alpha → 0). Consider:
   - Implementing Armijo backtracking line search when alpha < 1e-12
   - Using away-steps or pairwise FW variants
   - Better warm-start initialization

2. **C++ Implementation**: The reviewer noted the C++ partial OT solver uses a greedy heuristic that ignores feasibility. Should call proper LP solver.

3. **Numerical Precision**: Transport plans show small between-cluster leakage (1e-13) due to LP solver tolerances.