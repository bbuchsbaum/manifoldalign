#' Benchmark Data for Alignment Methods
#'
#' A small, deterministic multi-domain dataset used across alignment examples and
#' tests. It contains three domains that share two latent classes but exhibit
#' different linear transformations and noise levels. The object stores the raw
#' matrices/design frames and the class labels so alignment methods can generate
#' comparable hyperdesigns on demand.
#'
#' @name alignment_benchmark
#' @docType data
#' @format A list with the following components:
#' \describe{
#'   \item{domains}{List of three elements named `domain1`, `domain2`, and
#'   `domain3`. Each element contains a centered numeric matrix `x` (samples
#'   x features) and a `design` data frame with `id`, `condition`, and
#'   `domain` columns.}
#'   \item{labels}{Factor of length 80 giving the latent class for each sample.}
#' }
#'
#' @usage data(alignment_benchmark)
#' @keywords datasets
#'
#' @examples
#' data <- manifoldalign::alignment_benchmark
#' hd <- multidesign::hyperdesign(lapply(data$domains, function(dom) {
#'   multidesign::multidesign(dom$x, dom$design)
#' }))
#' str(hd, max.level = 1)
#'
#' @export
alignment_benchmark <- local({
  set.seed(20240219)

  n_per_class <- 40L
  n_domains <- 3L
  classes <- factor(rep(c("class_A", "class_B"), each = n_per_class))
  n_total <- length(classes)

  # Latent 4D structure shared across domains
  latent <- matrix(rnorm(n_total * 4L), nrow = n_total, ncol = 4L)

  class_shift <- matrix(rep(c(1, -1), each = n_per_class), nrow = n_total, ncol = 4L)
  class_shift[, 3:4] <- class_shift[, 3:4] / 1.5

make_domain <- function(transform, noise_sd, suffix) {
  X <- latent %*% transform + class_shift + matrix(rnorm(n_total * ncol(transform), sd = noise_sd),
                                                   nrow = n_total)
  X <- scale(X, center = TRUE, scale = FALSE)
  design <- data.frame(
    id = seq_len(n_total),
    condition = classes,
    domain = suffix,
    stringsAsFactors = FALSE
  )
  list(x = X, design = design)
}

  transforms <- list(
    diag(4),
    matrix(c(0.9, 0.2, 0, 0,
             -0.2, 1.1, 0, 0,
             0, 0, 0.8, 0.3,
             0, 0, -0.3, 0.9), nrow = 4, byrow = TRUE),
    matrix(c(1, 0, 0.4, 0,
             0, 0.7, 0, 0.2,
             -0.4, 0, 1, 0,
             0, -0.2, 0, 0.85), nrow = 4, byrow = TRUE)
  )

  noises <- c(0.2, 0.25, 0.3)
  domain_names <- paste0("domain", seq_len(n_domains))

  domain_list <- lapply(seq_along(transforms), function(idx) {
    make_domain(transforms[[idx]], noises[idx], suffix = domain_names[idx])
  })
  names(domain_list) <- domain_names

  list(
    domains = domain_list,
    labels = classes
  )
})

# Ensure object is visible for data() calls
utils::globalVariables("alignment_benchmark")
