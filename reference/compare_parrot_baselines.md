# Compare PARROT with Baseline Methods

Benchmarks PARROT against simpler alignment methods

## Usage

``` r
compare_parrot_baselines(val_data)
```

## Arguments

- val_data:

  Validation data from generate_parrot_validation_data

## Value

Data frame comparing different methods

## Examples

``` r
# \donttest{
vdata <- generate_parrot_validation_data(n_nodes = 30, n_anchors = 5)
comparison <- compare_parrot_baselines(vdata)
#> Running PARROT...
#> Running feature-only baseline...
#> Running anchor propagation baseline...
#> 
#> METHOD COMPARISON
#> =================
#>                          method      runtime top1_accuracy top5_accuracy
#> parrot                   PARROT 0.2009155750     0.5000000     0.9666667
#> feature_only       Feature-Only 0.0002756119     0.7000000            NA
#> anchor_prop  Anchor-Propagation 0.0044350624     0.1666667            NA
#>                    mrr
#> parrot       0.7166667
#> feature_only        NA
#> anchor_prop         NA
# }
```
