% Generated manually for package development
\name{cv_alignment_rows}
\alias{cv_alignment_rows}
\title{Cross-Validated Row-Wise Alignment Scoring}
\usage{
cv_alignment_rows(
  data,
  rows,
  fit_fn,
  features,
  k = 5L,
  feature_similarity = c("cosine", "correlation"),
  target_pool = c("analysis", "assessment", "both"),
  prediction_mode = c("auto", "embedding", "weights")
)
}
\arguments{
\item{data}{A \code{hyperdesign}-compatible object. If \code{data} is not
already a \pkg{multidesign} hyperdesign, it is coerced internally.}

\item{rows}{Explicit held-out row specification forwarded to
\code{multidesign::cv_rows()}. For hyperdesigns, each fold should be a named
list mapping block names or positions to held-out row indices.}

\item{fit_fn}{Function taking the analysis split for a fold and returning a
fitted alignment object.}

\item{features}{External feature matrices aligned to the original rows of each
block. Supply either a named list of matrices, one per block, or a
hyperdesign-compatible object whose \code{$x} matrices are treated as
features.}

\item{k}{Positive integer number of latent nearest neighbours used when scoring
held-out rows.}

\item{feature_similarity}{Similarity used for the external feature space.
Either \code{"cosine"} (default) or \code{"correlation"}.}

\item{target_pool}{Which rows are available as retrieval targets for each
held-out query block: \code{"analysis"} (default) uses training rows from the
fitted fold, \code{"assessment"} uses only other held-out rows, and
\code{"both"} concatenates the two pools when latent projection is available.}

\item{prediction_mode}{How neighbour rankings are produced. \code{"auto"}
(default) uses latent \code{\link{oos_predict}} projection for
embedding-oriented fits and prefers \code{predict(..., type = "weights")} for
transport-style fits when \code{target_pool = "analysis"}. \code{"embedding"}
forces latent neighbour search and \code{"weights"} forces weight-based
ranking against training targets.}
}
\value{
A \code{multidesign::cv_result} with one row per fold. The score table includes
\code{mean_top1_similarity}, \code{mean_topk_similarity},
\code{oracle_top1_similarity}, \code{oracle_topk_similarity},
\code{top1_gap}, \code{topk_gap}, \code{n_queries}, and \code{n_pairs}.
}
\description{
Runs explicit row-index cross-validation for multi-domain alignment models by
delegating fold construction and execution to \code{multidesign::cv_rows()} and
\code{multidesign::cross_validate()}. Held-out rows are scored by how similar
their cross-block neighbours are in an external feature space, using either
latent \code{\link{oos_predict}} projections or method-specific
\code{predict(..., type = "weights")} support when available.
}
\details{
The external feature matrices are used only for held-out evaluation. If your
fitting procedure uses correspondence tables keyed by original rows, include a
stable row identifier in each block's design before calling this function so
the training correspondences can be rebuilt inside \code{fit_fn}.
}
