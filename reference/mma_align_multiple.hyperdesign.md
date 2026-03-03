# Multiset Manifold Alignment (MMA)

Aligns 3+ domains by: (i) building spectral embeddings (commute-time or
Laplacian), (ii) aligning eigenbases via eigensignatures (histograms +
Hungarian), and (iii) refining per-domain alignments with an EM-based
probabilistic point-registration on the orthogonal group.

## Usage

``` r
# S3 method for class 'hyperdesign'
mma_align_multiple(
  data,
  ref_idx = 1L,
  preproc = center(),
  ncomp = 10,
  sigma = 0.73,
  knn = NULL,
  embedding = c("ctd", "lap"),
  normalize = c("hypersphere", "none"),
  hist_bins = NULL,
  hist_max_bins = 128L,
  hist_similarity = c("cor", "cosine"),
  signature_method = c("hist", "w2", "hybrid"),
  eig_penalty = 0,
  signature_gate_multiplicity = FALSE,
  signature_gate_tau = NULL,
  signature_retry_if_weak = FALSE,
  em_max_iter = 50L,
  em_tol = 1e-05,
  em_sigma0 = NULL,
  em_outlier = 0.01,
  force_SO = FALSE,
  target_var = 0.95,
  match_to = c("reference", "consensus"),
  consensus_centers = c("ref", "min"),
  n_centers = NULL,
  consensus_init = c("ref", "kmeans"),
  final_assignment = c("none", "hungarian", "argmax"),
  ...
)
```

## Arguments

- data:

  A hyperdesign object containing three or more domains

- ref_idx:

  Integer; which domain acts as reference (default 1)

- preproc:

  Preprocessing function (default multivarious::center())

- ncomp:

  Number of components; set NULL to choose automatically

- sigma:

  Gaussian affinity bandwidth (in data units)

- knn:

  Number of nearest neighbors for graph construction (default adaptive)

- embedding:

  One of "ctd" (commute-time; default) or "lap"

- normalize:

  One of "hypersphere" (default) or "none" for embeddings

- hist_bins:

  Integer number of bins for eigensignature histograms; NULL for auto

- hist_max_bins:

  Clamp for auto-bins (default 128)

- hist_similarity:

  One of "cor" or "cosine" (default "cor")

- signature_method:

  One of "hist", "w2", or "hybrid" for eigenvector matching

- eig_penalty:

  Eigenvalue penalty weight for signature alignment (default 0)

- signature_gate_multiplicity:

  Logical; gate matches by eigenvalue clusters

- signature_gate_tau:

  Threshold for eigenvalue clustering (NULL for auto)

- signature_retry_if_weak:

  Logical; retry with alternative method if weak match

- em_max_iter:

  Integer EM iterations (default 50)

- em_tol:

  Numeric EM tolerance on relative change (default 1e-5)

- em_sigma0:

  Optional initial stddev in embedding units; NULL for auto

- em_outlier:

  Prior mass for uniform outlier component (default 0.01)

- force_SO:

  Logical; if TRUE, force det(R)=+1 (default FALSE)

- target_var:

  Lower bound criterion for auto K when ncomp=NULL (default 0.95)

- match_to:

  One of "reference" (default) or "consensus"

- consensus_centers:

  Strategy to set number of consensus centers when match_to="consensus":
  "ref" (use ref domain size) or "min" (use minimum size)

- n_centers:

  Optional integer overriding consensus_centers; must be \<= min sizes

- consensus_init:

  Initialization for consensus centers: "ref" or "kmeans"

- final_assignment:

  Optional post-processing: "none" (default), "hungarian", or "argmax".
  Affects only extras\$assignment; core EM remains soft.

- ...:

  Unused.

## Value

A multiblock alignment object with stacked scores, block loadings, and
extras (rotations, posteriors, histogram alignment info, and optionally
consensus and final assignment).

## Details

Two matching modes are supported: - match_to = "reference": pairwise EM
to a chosen reference domain. - match_to = "consensus": joint EM that
learns a shared template (consensus centers) along with per-domain
orthogonal transforms; supports unequal node counts naturally.

The implementation reuses shared helpers in this package: -
neighborweights::graph_weights + adjacency for graph building -
graph_laplacian (internal) for Laplacian construction - safe_compute,
compute_block_indices, feature_block_indices, new_alignment_result

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(30), 10, 3)
X2 <- matrix(rnorm(30), 10, 3)
X3 <- matrix(rnorm(30), 10, 3)
data_list <- list(X1, X2, X3)
result <- mma_align_multiple(data_list, ncomp = 2)
dim(result$scores)
#> NULL
# }
```
