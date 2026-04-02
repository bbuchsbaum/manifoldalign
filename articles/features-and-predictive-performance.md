# Features, Correspondences, and Predictive Performance

You have multiple data domains, and each row also has an external
feature vector. Those features might come from text embeddings, anatomy,
metadata, or another modality entirely. The alignment problem is:

- use the observed matrices to learn a shared manifold,
- use feature similarity to define cross-domain supervision,
- use held-out feature agreement to define predictive performance.

The important point is that predictive performance is **not**
“reconstruct the feature vector.” It is “does a held-out row land near
rows in other domains whose features are similar?”

## Setup

``` r
library(manifoldalign)
library(multidesign)
library(ggplot2)
```

## Simulate two aligned domains

We will generate two observed matrices from the same latent variables
and keep separate feature vectors for scoring.

``` r
set.seed(20260420)
n <- 60
latent_dim <- 3

Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
A1 <- matrix(rnorm(latent_dim * 8), latent_dim, 8)
A2 <- matrix(rnorm(latent_dim * 6), latent_dim, 6)

X1 <- Z %*% A1 + matrix(rnorm(n * 8, sd = 0.05), n, 8)
X2 <- Z %*% A2 + matrix(rnorm(n * 6, sd = 0.05), n, 6)

features_true <- list(
  a = Z + matrix(rnorm(n * latent_dim, sd = 0.02), n, latent_dim),
  b = Z + matrix(rnorm(n * latent_dim, sd = 0.02), n, latent_dim)
)
```

The observed data go into a `hyperdesign`. The feature matrices stay
outside the fit and will be used for correspondences and held-out
evaluation.

``` r
hd <- hyperdesign(list(
  a = multidesign(X1, data.frame(row_id = seq_len(n))),
  b = multidesign(X2, data.frame(row_id = seq_len(n)))
))

hd
#> 
#> === Hyperdesign Object ===
#> 
#> Number of blocks:  2 
#> 
#> +- Block  1  (a)  -----------------
#> |  Dimensions: 60 x 8 
#> |  Design Variables: row_id 
#> |  Design Structure: 
#> |   * row_id: 60 levels (1, 2, 3...59, 60)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> +- Block  2  (b)  -----------------
#> |  Dimensions: 60 x 6 
#> |  Design Variables: row_id 
#> |  Design Structure: 
#> |   * row_id: 60 levels (1, 2, 3...59, 60)
#> |  Column Design: Present
#> |   Variables:  .index 
#> 
#> =======================
#> 
```

## Turn feature similarity into correspondences

Here we use mutual top-1 cosine matches between the feature vectors. In
a real problem you might use a threshold, mutual top-`k`, or weighted
anchors.

``` r
full_corr <- feature_correspondences(features_true$a, features_true$b)

nrow(full_corr)
#> [1] 57
head(full_corr)
#>   domain_i index_i domain_j index_j
#> 1        1       1        2       1
#> 2        1       2        2       2
#> 3        1       3        2       3
#> 4        1       4        2       4
#> 5        1       5        2       5
#> 6        1       6        2       6

stopifnot(nrow(full_corr) >= 0.8 * n)
```

That correspondence table is only for the full data. During
cross-validation we must rebuild it on each training split so held-out
rows do not leak into the fit.

## Score held-out rows

[`cv_alignment_rows()`](https://bbuchsbaum.github.io/manifoldalign/reference/cv_alignment_rows.md)
does row-wise cross-validation. We hold out rows from domain `a`, fit on
the remaining data, project the held-out rows, and score them against
the training rows in domain `b`.

``` r
rows <- list(
  list(a = 1:10),
  list(a = 11:20),
  list(a = 21:30)
)

fit_fn <- function(analysis) {
  ids_a <- analysis$a$design$row_id
  ids_b <- analysis$b$design$row_id

  corr <- feature_correspondences(
    features_true$a[ids_a, , drop = FALSE],
    features_true$b[ids_b, , drop = FALSE]
  )

  manifoldalign::ssma_align(
    analysis,
    correspondences = corr,
    preproc = multivarious::center(),
    ncomp = 2,
    control = manifoldalign::ssma_align_control(
      knn = 8,
      rank_per_domain = 12,
      verbose = FALSE
    )
  )
}

cv_signal <- manifoldalign::cv_alignment_rows(
  hd,
  rows = rows,
  fit_fn = fit_fn,
  features = features_true,
  k = 3,
  target_pool = "analysis"
)

cv_signal$scores
#> # A tibble: 3 × 9
#>   mean_top1_similarity mean_topk_similarity oracle_top1_similarity
#>                  <dbl>                <dbl>                  <dbl>
#> 1                0.675                0.680                  0.999
#> 2                0.398                0.473                  1.000
#> 3                0.815                0.780                  1.000
#> # ℹ 6 more variables: oracle_topk_similarity <dbl>, top1_gap <dbl>,
#> #   topk_gap <dbl>, n_queries <dbl>, n_pairs <dbl>, .fold <int>

stopifnot(all(is.finite(cv_signal$scores$mean_top1_similarity)))
stopifnot(all(is.finite(cv_signal$scores$mean_topk_similarity)))
stopifnot(mean(cv_signal$scores$mean_top1_similarity) > 0.5)
```

The key score is `mean_top1_similarity`: for each held-out query row,
project it into the fitted space, retrieve its nearest neighbor in the
target block, and measure the cosine similarity of their external
feature vectors.

## Does the score detect broken semantic structure?

A good held-out metric should drop when the feature-defined target
structure is destroyed. We can test that by shuffling the target feature
rows while keeping the fitted model fixed.

``` r
set.seed(20260421)
features_shuffled <- list(
  a = features_true$a,
  b = features_true$b[sample(seq_len(n)), , drop = FALSE]
)

cv_shuffled <- manifoldalign::cv_alignment_rows(
  hd,
  rows = rows,
  fit_fn = fit_fn,
  features = features_shuffled,
  k = 3,
  target_pool = "analysis"
)

comparison <- rbind(
  aggregate_cv(cv_signal, "Aligned features"),
  aggregate_cv(cv_shuffled, "Shuffled target features")
)

comparison
#>                   scenario metric     value
#> 1         Aligned features  top-1 0.6292889
#> 2         Aligned features  top-k 0.6444415
#> 3 Shuffled target features  top-1 0.1063780
#> 4 Shuffled target features  top-k 0.0671158

stopifnot(
  mean(cv_signal$scores$mean_top1_similarity) >
    mean(cv_shuffled$scores$mean_top1_similarity) + 0.2
)
```

![Held-out feature agreement drops when the target feature rows are
shuffled, which is what we want from a predictive scoring
rule.](features-and-predictive-performance_files/figure-html/plot-comparison-1.png)

Held-out feature agreement drops when the target feature rows are
shuffled, which is what we want from a predictive scoring rule.

This is the predictive story:

- the alignment is fit on the observed matrices,
- feature similarity defines the semantic notion of a good cross-domain
  match,
- held-out performance asks whether new rows land near semantically
  similar rows in the fitted space.

## Design takeaway

If your features are only there to define supervision and evaluation,
keep them **outside** the fit and use them to:

1.  build correspondence edges on the training rows,
2.  score held-out rows with
    [`cv_alignment_rows()`](https://bbuchsbaum.github.io/manifoldalign/reference/cv_alignment_rows.md).

If the features are a real modality you want represented in the
embedding, you can also add them as another domain in the alignment
itself. That is a different model, and cross-validation then has to hold
out the matching feature rows too.
