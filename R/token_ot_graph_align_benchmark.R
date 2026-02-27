#' Benchmark token_ot_graph_align scalability on synthetic permuted graphs
#'
#' Generates paired domains where the target is a permuted/noisy copy of the
#' source and reports runtime and basic alignment accuracy as problem size grows.
#'
#' Benchmarks are intended for interactive profiling rather than automated CI.
#'
#' @param sizes Integer vector of node counts to benchmark.
#' @param d Integer feature dimension for synthetic node features.
#' @param noise_sd Numeric noise standard deviation added to target features.
#' @param candidate_k Integer number of candidates per source node.
#' @param n_levels Integer number of multilevel hierarchy levels.
#' @param prior_mode Prior lifting mode for multilevel (`"none"`, `"hard"`, `"soft"`).
#' @param token_mode Tokenization mode passed to [token_ot_graph_align_control()].
#' @param views Views passed to [token_ot_graph_align_control()].
#' @param n_reps Integer number of replications per size.
#' @param seed Integer seed for reproducibility.
#' @param verbose Logical; print per-run progress.
#'
#' @return A data.frame with per-run runtime and summary metrics.
#' @export
benchmark_token_ot_graph_align_scalability <- function(
  sizes = c(50L, 100L, 200L),
  d = 4L,
  noise_sd = 0.01,
  candidate_k = 80L,
  n_levels = 1L,
  prior_mode = c("none", "hard", "soft"),
  token_mode = c("view_only", "view_plus_neighbors"),
  views = "raw",
  n_reps = 3L,
  seed = 1L,
  verbose = FALSE
) {
  prior_mode <- match.arg(prior_mode)
  token_mode <- match.arg(token_mode)

  chk::chk_number(d)
  chk::chk_true(d >= 1)
  chk::chk_number(noise_sd)
  chk::chk_true(noise_sd >= 0)
  chk::chk_number(candidate_k)
  chk::chk_true(candidate_k >= 1)
  chk::chk_number(n_levels)
  chk::chk_true(n_levels >= 1)
  chk::chk_number(n_reps)
  chk::chk_true(n_reps >= 1)
  chk::chk_number(seed)
  chk::chk_logical(verbose)

  sizes <- as.integer(sizes)
  sizes <- sizes[is.finite(sizes) & sizes >= 3L]
  if (!length(sizes)) stop("`sizes` must contain integers >= 3.", call. = FALSE)

  rows <- vector("list", length(sizes) * as.integer(n_reps))
  ptr <- 0L

  for (s in seq_along(sizes)) {
    n <- sizes[[s]]
    for (r in seq_len(as.integer(n_reps))) {
      set.seed(as.integer(seed) + s * 1009L + r * 17L)
      X1 <- matrix(rnorm(n * as.integer(d)), n, as.integer(d))
      perm21 <- sample.int(n)
      map12 <- match(seq_len(n), perm21)
      X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * as.integer(d), sd = noise_sd), n, as.integer(d))

      hd <- as_hyperdesign(list(source = X1, target = X2))
      ctrl <- token_ot_graph_align_control(
        n_levels = n_levels,
        prior_mode = prior_mode,
        views = views,
        candidate_k = candidate_k,
        token_mode = token_mode,
        verbose = FALSE
      )

      if (isTRUE(verbose)) {
        message("token_ot_graph_align benchmark: n=", n, " rep=", r)
      }

      t0 <- proc.time()[[3]]
      fit <- token_ot_graph_align(hd, ncomp = min(10L, as.integer(d)), control = ctrl)
      runtime <- proc.time()[[3]] - t0

      assignment <- fit$assignment
      acc <- NA_real_
      if (!is.null(assignment) && length(assignment) == n) {
        ok <- !is.na(assignment)
        if (any(ok)) {
          acc <- mean(assignment[ok] == map12[ok])
        }
      }

      nnz <- if (inherits(fit$transport_plan, "Matrix")) length(fit$transport_plan@x) else NA_integer_
      ptr <- ptr + 1L
      rows[[ptr]] <- data.frame(
        n = n,
        rep = r,
        runtime_sec = runtime,
        accuracy = acc,
        nnz = nnz,
        n_levels = as.integer(n_levels),
        prior_mode = prior_mode,
        candidate_k = as.integer(candidate_k),
        token_mode = token_mode,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows[seq_len(ptr)])
}

