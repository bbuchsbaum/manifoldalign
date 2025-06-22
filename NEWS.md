# manifoldalign NEWS

## Version 0.1.0.9000 (Development)

### Major Changes

* Implemented Fused-Partial Gromov-Wasserstein (FPGW) distance for heterogeneous domain adaptation
* Added support for mass-constrained partial optimal transport with `rho` parameter

### Bug Fixes

* Fixed network simplex solver segfaults by improving error handling
* Corrected TV penalty implementation to match paper's mathematical formulation (equation 11)

### Deprecations and Warnings

* **TV penalty (`lambda` parameter) is now marked as experimental**. The TV penalty `λ(|μ|² + |ν|² - 2|γ|²)` is a concave quadratic that biases optimization towards either zero or full mass, but cannot guarantee gradual sparsity. For predictable partial transport, use the mass-constrained formulation (`rho` parameter) instead.

### Internal Changes

* Removed incorrect augmentation method for TV penalty
* TV penalty now handled correctly in Frank-Wolfe step size computation
* Added runtime warning when `lambda > 0` to alert users about experimental status

### Documentation

* Updated FPGW documentation to clarify TV penalty limitations
* Examples now emphasize `rho` parameter for partial transport
* Added comprehensive technical documentation about TV penalty behavior

### Testing

* Modified TV penalty tests to check for experimental warnings rather than sparsity behavior
* All tests now pass with corrected expectations