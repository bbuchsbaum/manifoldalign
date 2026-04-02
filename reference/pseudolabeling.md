# Pseudolabeling for Unsupervised Domain Adaptation

Provides pseudolabeling system for unsupervised domain adaptation with
KEMA. Identifies high-confidence anchor samples to guide domain
alignment.

## Value

NULL (documentation page only).

## Details

The pseudolabeling system addresses unsupervised domain adaptation by
identifying reliable correspondences between samples from different
domains. The approach uses similarity-based clustering, diversity-aware
selection, adaptive thresholding, and quality control filtering.

Main functions:

- [`assign_pseudolabels()`](https://bbuchsbaum.github.io/manifoldalign/reference/assign_pseudolabels.md):
  General-purpose pseudolabeling from sparse similarity matrices

- [`high_sim_pseudolabels()`](https://bbuchsbaum.github.io/manifoldalign/reference/high_sim_pseudolabels.md):
  Specialized function for multi-domain data using cosine similarity

- [`create_synthetic_similarity_matrix()`](https://bbuchsbaum.github.io/manifoldalign/reference/create_synthetic_similarity_matrix.md):
  Generate synthetic data for testing

- [`evaluate_pseudolabeling()`](https://bbuchsbaum.github.io/manifoldalign/reference/evaluate_pseudolabeling.md):
  Evaluate pseudolabeling performance against ground truth

Integration with KEMA:

    # Generate pseudolabels
    plabs <- assign_pseudolabels(similarity_matrix, min_clusters = 20)

    # Use with KEMA
    fit <- kema.hyperdesign(
      data = strata,
      y = plabs$labels,
      u = 0.8,           # Trust geometry over pseudolabels
      dweight = 0.2,     # Mild class separation
      simfun = function(lab) binary_label_matrix(lab, type = "s"),
      disfun = function(lab) binary_label_matrix(lab, type = "d")
    )

Key parameters:

- **sim_threshold**: Controls which similarities are considered "high".
  Can be adaptive.

- **diversity_weight**: Balances cluster coherence vs. representative
  diversity

- **min_clusters/max_clusters**: Controls the number of anchor points

- **min_cluster_size**: Ensures clusters are large enough to be reliable

## Examples

``` r
# \donttest{
library(Matrix)
# Create synthetic similarity matrix
n <- 100
sim_matrix <- Matrix::rsparsematrix(n, n, density = 0.1, rand.x = runif)
sim_matrix <- (sim_matrix + Matrix::t(sim_matrix)) / 2
Matrix::diag(sim_matrix) <- 1

# Assign pseudolabels
result <- assign_pseudolabels(sim_matrix, min_clusters = 5)
#> Warning: Final number of representatives (1) is below min_clusters (5). Consider lowering sim_threshold or min_cluster_size.
table(result$labels, useNA = "always")
#> 
#> anchor_001       <NA> 
#>         98          2 
# }
```
