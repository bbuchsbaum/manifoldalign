# Run Enforced KEMA Numerical Validation

Fits deterministic KEMA fixtures and records the same numerical
contracts required of every returned fit: generalized-eigen backward
residual, B-orthogonality, numerical-zero exclusion, and agreement
between matrix and operator backends up to subspace indeterminacy.

## Usage

``` r
run_kema_validation_suite(verbose = TRUE)
```

## Arguments

- verbose:

  Whether to print a concise report.

## Value

A list containing \`fidelity_validation\`, \`backend_agreement\`,
\`withdrawn_validations\`, and \`overall_success\`.

## Details

Withdrawn paper-eigenvalue, regression-solver, and fabricated
reconstruction checks are reported for auditability but do not
contribute to \`overall_success\`.

## Examples

``` r
# \donttest{
results <- run_kema_validation_suite(verbose = FALSE)
results$overall_success
#> [1] TRUE
# }
```
