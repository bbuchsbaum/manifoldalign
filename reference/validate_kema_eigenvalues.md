# Withdrawn KEMA Paper-Eigenvalue Validator

This compatibility shim is retained so historical callers receive an
explicit failure. The former routine did not reproduce a documented
paper fixture and solved a different generalized problem whose
right-hand side was only \`lambda \* I\`; its claimed comparison with
paper eigenvalues was therefore not valid accuracy evidence.

## Usage

``` r
validate_kema_eigenvalues(
  expected_eigenvals = c(0.82, 0.41),
  tolerance = 0.1,
  n_per_domain = 100
)
```

## Arguments

- expected_eigenvals:

  Deprecated and ignored.

- tolerance:

  Deprecated and ignored.

- n_per_domain:

  Deprecated and ignored.

## Value

This function emits a deprecation warning and then errors.

## See also

\[run_kema_validation_suite()\]

## Examples

``` r
if (FALSE) { # \dontrun{
validate_kema_eigenvalues()
} # }
```
