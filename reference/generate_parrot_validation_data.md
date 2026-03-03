# Generate Synthetic Network Alignment Data

Creates two networks with known correspondence for validating PARROT

## Usage

``` r
generate_parrot_validation_data(
  n_nodes = 100,
  n_anchors = 10,
  noise_level = 0.1,
  structure = c("ring", "grid", "random", "community"),
  permute_fraction = 0.3
)
```

## Arguments

- n_nodes:

  Number of nodes per network

- n_anchors:

  Number of anchor correspondences

- noise_level:

  Noise level for second network (0-1)

- structure:

  Type of network structure: "ring", "grid", "random", "community"

- permute_fraction:

  Fraction of nodes to permute in second network

## Value

List with two networks and ground truth alignment

## Examples

``` r
# \donttest{
vdata <- generate_parrot_validation_data(n_nodes = 20, n_anchors = 5)
str(vdata$ground_truth)
#> List of 3
#>  $ permutation   : int [1:20] 1 2 3 4 5 6 7 8 20 10 ...
#>  $ inverse_perm  : int [1:20] 1 2 3 4 5 6 7 8 20 10 ...
#>  $ anchor_indices:List of 2
#>   ..$ net1: int [1:5] 10 12 13 14 18
#>   ..$ net2: int [1:5] 10 12 13 14 18
# }
```
