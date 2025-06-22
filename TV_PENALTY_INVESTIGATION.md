# TV Penalty Investigation Report

## Summary
After extensive investigation and consultation with Gemini Pro, O3, and Pro models, the TV penalty implementation remains problematic. Despite implementing the consensus fix, the behavior is still inverted.

## Consensus from AI Models
All three models agreed on the core issue:
- **Root Cause**: The current implementation penalizes NOT transporting mass instead of penalizing transported mass
- **Solution**: Add lambda to transport costs (grad + lambda) with zero-cost escape routes

## Implementation Attempted
```r
partial_ot_tv <- function(cost, p, q, lambda) {
  n <- length(p)
  m <- length(q)
  
  # Create augmented cost matrix with zero escape costs
  C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
  
  # Add lambda to transport costs to penalize transported mass
  C_aug[1:n, 1:m] <- cost + lambda
  
  # Standard augmentation for marginals
  p_aug <- c(p, sum(q))
  q_aug <- c(q, sum(p))
  
  # Solve augmented OT problem
  gamma_aug <- classical_ot_lp(C_aug, p_aug, q_aug)
  
  # Extract partial transport plan
  gamma <- gamma_aug[1:n, 1:m]
  return(gamma)
}
```

## Test Results
### Simple Test Cases
- With uniform cost matrices, the implementation correctly routes all mass through dummy nodes when lambda > 0
- This produces zero transport in the top-left block (correct behavior)

### FPGW Integration
- lambda=0: transported mass = 0.626 (partial due to outliers)
- lambda=1: transported mass = 0.626 (same as lambda=0)
- lambda=5: transported mass = 1.0 (full transport!)

This is still backwards - higher lambda should reduce transport.

## Key Findings

1. **The augmentation method works correctly in isolation**
   - When tested with simple cost matrices, mass correctly escapes through dummy nodes
   - The issue appears only in the Frank-Wolfe context

2. **Gradient values are all positive**
   - The FPGW gradient ranges from 0.44 to 7.55 (no negative values)
   - Even with lambda=5, all transport costs (grad + lambda) are positive
   - Yet the solver still prefers full transport

3. **Possible deeper issue with Frank-Wolfe**
   - The Frank-Wolfe algorithm may have additional logic that interferes
   - The step size computation or convergence criteria might be affected
   - The warm-start initialization might bias toward full transport

## Remaining Questions

1. **Is the TV penalty formulation correct for Frank-Wolfe?**
   - The paper mentions TV penalty but doesn't provide implementation details
   - The augmentation method might not be compatible with Frank-Wolfe's assumptions

2. **Alternative formulations to consider:**
   - Direct penalty in objective without augmentation
   - Proximal gradient methods for the TV term
   - Different regularization that's more compatible with Frank-Wolfe

## Current Status
- All other FPGW tests pass (52/54)
- Only TV penalty tests fail (2/54)
- The implementation is functional for classical FGW and mass-constrained partial OT
- TV penalty remains an open issue requiring further theoretical investigation

## Recommendation
Given the complexity of getting TV penalty to work correctly with Frank-Wolfe, and that all other functionality is working, I recommend:
1. Documenting the TV penalty as experimental/not fully supported
2. Focusing on the working features (classical FGW, mass-constrained partial OT)
3. Potentially reaching out to the paper authors for clarification on the TV penalty implementation