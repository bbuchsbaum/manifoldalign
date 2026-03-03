# Utility functions for manifoldalign package

This file contains shared utility functions used across multiple
modules.

## Usage

``` r
safe_compute(expr, error_msg, warning_fn = NULL)
```

## Arguments

- expr:

  Expression to evaluate

- error_msg:

  Custom error message to display on failure

- warning_fn:

  Optional function to handle warnings (default: standard warning)

## Value

Result of expression evaluation

## Usage of safe_compute

The `safe_compute` function should be used for operations that:

- May fail due to external dependencies or numerical issues

- Should stop execution with a clear error message (not provide
  fallbacks)

- Benefit from consistent error reporting across the package

Do NOT use `safe_compute` for:

- Operations that need custom fallback behavior

- Cases where you want to continue execution after an error

- Simple parameter validation (use chk:: functions instead)

Examples of good candidates:

- Eigenvalue computations:
  `safe_compute(PRIMME::eigs_sym(...), "Eigenvalue computation failed")`

- Graph construction:
  `safe_compute(neighborweights::graph_weights(...), "Graph construction failed")`

- Matrix operations: `safe_compute(solve(A), "Matrix solve failed")`

Safe computation wrapper with enhanced error handling

Wraps expressions in tryCatch with informative error messages and
optional warning handling. Provides consistent error reporting across
the package.
