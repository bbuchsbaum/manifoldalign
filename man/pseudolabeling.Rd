% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudolabel.R
\name{pseudolabeling}
\alias{pseudolabeling}
\title{Pseudolabeling for Unsupervised Domain Adaptation}
\description{
Provides pseudolabeling system for unsupervised domain adaptation with 
KEMA. Identifies high-confidence anchor samples to guide domain alignment.
}
\details{
The pseudolabeling system addresses unsupervised domain adaptation by 
identifying reliable correspondences between samples from different domains. 
The approach uses similarity-based clustering, diversity-aware selection, 
adaptive thresholding, and quality control filtering.

Main functions:
\itemize{
  \item \code{assign_pseudolabels()}: General-purpose pseudolabeling from 
    sparse similarity matrices
  \item \code{high_sim_pseudolabels()}: Specialized function for 
    multi-domain data using cosine similarity
  \item \code{create_synthetic_similarity_matrix()}: Generate synthetic 
    data for testing
  \item \code{evaluate_pseudolabeling()}: Evaluate pseudolabeling 
    performance against ground truth
}

Integration with KEMA:
\preformatted{
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
}

Key parameters:
\itemize{
  \item \strong{sim_threshold}: Controls which similarities are considered 
    "high". Can be adaptive.
  \item \strong{diversity_weight}: Balances cluster coherence vs. 
    representative diversity
  \item \strong{min_clusters/max_clusters}: Controls the number of anchor 
    points
  \item \strong{min_cluster_size}: Ensures clusters are large enough to be 
    reliable
}
}
