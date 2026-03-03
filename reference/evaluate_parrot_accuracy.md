# Evaluate PARROT Alignment Accuracy

Computes various accuracy metrics for PARROT alignment results

## Usage

``` r
evaluate_parrot_accuracy(parrot_result, ground_truth, k = c(1, 5, 10))
```

## Arguments

- parrot_result:

  Result from parrot() function

- ground_truth:

  Ground truth from generate_parrot_validation_data

- k:

  Top-k accuracy levels to compute (default: c(1, 5, 10))

## Value

A list with alignment accuracy metrics including hit rate at various k
values

## Examples

``` r
# \donttest{
vdata <- generate_parrot_validation_data(n_nodes = 20, n_anchors = 5)
# evaluate_parrot_accuracy(result, vdata$ground_truth)
# }
```
