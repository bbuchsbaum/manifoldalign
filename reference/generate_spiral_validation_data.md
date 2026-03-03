# Generate Synthetic Two-Domain Spiral Data

Creates the synthetic spiral dataset used in KEMA paper Figure 2 for
numerical validation.

## Usage

``` r
generate_spiral_validation_data(
  n_per_domain = 100,
  noise_level = 0.1,
  seed = 42
)
```

## Arguments

- n_per_domain:

  Number of samples per domain

- noise_level:

  Gaussian noise standard deviation

- seed:

  Random seed for reproducibility

## Value

List with domain data and labels

## Examples

``` r
data <- generate_spiral_validation_data(n_per_domain = 50, seed = 123)
str(data)
#> List of 4
#>  $ strata :List of 2
#>   ..$ :List of 2
#>   .. ..$ x     : num [1:50, 1:2] -0.056 0.225 0.603 0.56 0.545 ...
#>   .. ..$ design:'data.frame':    50 obs. of  1 variable:
#>   .. .. ..$ labels: chr [1:50] "early" "early" "early" "early" ...
#>   ..$ :List of 2
#>   .. ..$ x     : num [1:50, 1:2] -0.071 0.1292 0.0858 -0.0249 -0.2906 ...
#>   .. ..$ design:'data.frame':    50 obs. of  1 variable:
#>   .. .. ..$ labels: chr [1:50] "early" "early" "early" "early" ...
#>  $ labels : chr [1:100] "early" "early" "early" "early" ...
#>  $ domain1:List of 2
#>   ..$ x     : num [1:50, 1:2] -0.056 0.225 0.603 0.56 0.545 ...
#>   ..$ design:'data.frame':   50 obs. of  1 variable:
#>   .. ..$ labels: chr [1:50] "early" "early" "early" "early" ...
#>  $ domain2:List of 2
#>   ..$ x     : num [1:50, 1:2] -0.071 0.1292 0.0858 -0.0249 -0.2906 ...
#>   ..$ design:'data.frame':   50 obs. of  1 variable:
#>   .. ..$ labels: chr [1:50] "early" "early" "early" "early" ...
```
