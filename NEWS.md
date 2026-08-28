# manifoldalign NEWS

## Version 0.1.0.9000 (Development)

### Major Changes

* Implemented Fused-Partial Gromov-Wasserstein (FPGW) distance for heterogeneous domain adaptation
* Added support for mass-constrained partial optimal transport with `rho` parameter

### Bug Fixes

* Fixed network simplex solver segfaults by improving error handling
* Completed the package-wide migration from the former `neighborweights`
  package name to the CRAN `adjoin` package. GlobalGeo anchor smoothness now
  uses the current spatial-adjacency contract, preserves one-way kNN edges by
  union symmetrization, and gives isolates an explicit zero-penalty policy.
* Withdrew `validate_out_of_sample_reconstruction()`, which previously returned
  `expected_error` plus random noise without fitting KEMA or evaluating held-out
  predictions. Historical values from this routine are invalid as accuracy
  evidence, and it is no longer included in `run_kema_validation_suite()`.
* Corrected KEMA component selection to solve on the numerical image of the
  kernel blocks and explicitly deflate graph-component null modes. KEMA now
  uses eigencore's certified full or constrained generalized LOBPCG solver,
  requires a strictly positive metric regularizer, rejects materially negative
  modes, and independently checks residuals and B-orthogonality at `1e-6`.
  Auto mode tries alternate backends and then stops if no backend passes;
  warning-only failed fits are no longer returned. Historical KEMA results with
  `fidelity$passed == FALSE`, or whose selected eigenvalues do not exceed their
  numerical-zero threshold, are invalid.
* Corrected TV penalty implementation to match paper's mathematical formulation (equation 11)

### Deprecations and Warnings

* `kema_orig()` is deprecated in favor of the single public `kema()` entry
  point and is scheduled for removal in manifoldalign 1.0.0. `solver =
  "regression"` is withdrawn because it never selected a distinct regression
  implementation. Formerly ignored extension arguments now error instead of
  misrepresenting the fitted objective.
* `validate_kema_eigenvalues()` is withdrawn because it did not reproduce a
  documented paper fixture and solved the wrong generalized right-hand side.
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