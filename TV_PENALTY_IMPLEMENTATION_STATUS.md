# TV Penalty Implementation Status

## Summary
The TV penalty in the FPGW implementation has been mathematically corrected according to the paper's equation (11), but fundamental issues remain with the penalty formulation itself.

## Mathematical Formulation
The paper defines the TV penalty as:
```
λ * (|μ|² + |ν|² - 2|γ|²)
```
This is a **smooth quadratic penalty**, not an L1 penalty as the name "Total Variation" might suggest.

## Implementation Status

### ✅ Completed
1. **Removed incorrect augmentation method** - The augmentation approach was fundamentally wrong for this quadratic formulation
2. **Implemented correct step size computation** - Added TV penalty contributions to the Frank-Wolfe step size:
   ```r
   if (lambda > 0) {
     mass_P <- sum(P)
     mass_dP <- sum(dP)
     a <- a + 2 * lambda * mass_dP^2      # quadratic term
     b <- b - 4 * lambda * mass_P * mass_dP  # linear term
   }
   ```
3. **Updated objective function** - Includes the quadratic TV penalty term
4. **Uses standard OT oracle** - Both classical and TV-penalized cases use the same oracle

### ❌ Known Issues
1. **Penalty doesn't encourage sparsity** - The quadratic formulation creates a valley at mass=1, causing the algorithm to converge to full transport
2. **Test failures** - Tests expecting monotonic mass decrease with increasing λ fail
3. **Initialization dependency** - Results depend heavily on warm-start initialization

## Root Cause Analysis

The fundamental issue is that the TV penalty `-2λ|γ|²` is minimized when either:
- |γ| = 0 (no transport)
- |γ| = 1 (full transport)

For typical problems where the FGW objective favors some transport, the algorithm converges to |γ| = 1 regardless of λ value.

## Recommendations

### Short Term
1. **Mark TV penalty as experimental** in documentation
2. **Add warning** when users set λ > 0
3. **Update tests** to reflect actual behavior or skip TV tests

### Long Term
1. **Contact paper authors** to clarify if this behavior is intended
2. **Consider alternative formulations**:
   - True L1 penalty: `λ * ||γ||₁`
   - Entropy regularization: `-λ * H(γ)`
   - KL divergence penalty
3. **Implement proximal methods** if pursuing L1 penalties

## Code Locations
- Main implementation: `/Users/bbuchsbaum/code/manifoldalign/R/fpgw.R` (lines 310-320)
- TV oracle (deprecated): `/Users/bbuchsbaum/code/manifoldalign/R/fpgw.R` (lines 548-616)
- Tests: `/Users/bbuchsbaum/code/manifoldalign/tests/testthat/test-fpgw-validation.R` (lines 268-305)

## Technical Details
The current implementation is mathematically correct according to the paper but doesn't produce the expected sparsity-inducing behavior. This appears to be a limitation of the penalty formulation rather than an implementation bug.