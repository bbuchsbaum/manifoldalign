# Validate KEMA Eigenvalues Against Paper Specifications

Tests whether KEMA produces eigenvalues matching those reported in
Figure 2 of Tuia & Camps-Valls (2016).

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

  Expected eigenvalue ratios from paper (default: c(0.82, 0.41))

- tolerance:

  Numerical tolerance for comparison

- n_per_domain:

  Number of samples per domain for test

## Value

List with validation results

## Examples

``` r
# \donttest{
result <- validate_kema_eigenvalues()
#> Warning: Eigenvalue extraction failed: upper value must be greater than lower value
print(result$success)
#> [1] FALSE
# }
```
