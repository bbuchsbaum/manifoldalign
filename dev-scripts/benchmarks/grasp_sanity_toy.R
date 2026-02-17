#!/usr/bin/env Rscript

# Tiny deterministic sanity check for core GRASP assignment quality.
# This script uses the internal pairwise GRASP pipeline directly because it
# exposes the exact assignment vector for a known permutation.

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION") &&
      grepl("^Package:\\s*manifoldalign", paste(readLines("DESCRIPTION", n = 5L), collapse = "\n"))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(manifoldalign)
})

build_toy_case <- function(seed = 15L, n = 8L, noise_sd = 0.01) {
  set.seed(seed)
  angles <- 2 * pi * (seq_len(n) / n)

  # Slightly asymmetric circular features to avoid perfect rotational symmetry.
  X1 <- cbind(
    cos(angles),
    sin(angles),
    0.2 * (seq_len(n) / n) + stats::rnorm(n, sd = 0.01)
  )

  # Non-trivial, fixed permutation.
  perm <- n:1
  X2 <- X1[perm, , drop = FALSE] + matrix(stats::rnorm(n * ncol(X1), sd = noise_sd), nrow = n)

  hd <- list(
    domain1 = list(x = X1, design = data.frame(node_id = seq_len(n))),
    domain2 = list(x = X2, design = data.frame(node_id = seq_len(n)))
  )
  class(hd) <- c("hyperdesign", "list")

  list(hyperdesign = hd, perm = perm, X1 = X1, X2 = X2)
}

run_grasp_core <- function(toy, ncomp = 5L, q_descriptors = 30L, sigma = 1.0, lambda = 0.05) {
  bases <- manifoldalign:::compute_grasp_basis(toy$hyperdesign, ncomp = ncomp, use_laplacian = TRUE)
  descriptors <- manifoldalign:::compute_grasp_descriptors(bases, q_descriptors = q_descriptors, sigma = sigma)
  alignment <- manifoldalign:::align_grasp_bases(
    bases[[1]], bases[[2]],
    descriptors[[1]], descriptors[[2]],
    lambda = lambda
  )
  manifoldalign:::compute_grasp_assignment(
    bases[[1]], bases[[2]],
    descriptors[[1]], descriptors[[2]],
    alignment$rotation,
    distance_method = "cosine",
    solver_method = "linear",
    alpha = 0.7
  )
}

run_raw_baseline <- function(toy) {
  n <- nrow(toy$X1)
  d12 <- as.matrix(stats::dist(rbind(toy$X1, toy$X2)))[seq_len(n), n + seq_len(n), drop = FALSE]
  manifoldalign:::solve_lsap_with_padding(d12)
}

accuracy <- function(pred, target) {
  mean(as.integer(pred) == as.integer(target))
}

toy <- build_toy_case()
grasp_out <- run_grasp_core(toy)
raw_pred <- run_raw_baseline(toy)

acc_grasp <- accuracy(grasp_out$assignment, toy$perm)
acc_raw <- accuracy(raw_pred, toy$perm)
acc_random_expectation <- 1 / length(toy$perm)

cat("GRASP Toy Sanity Check\n")
cat("----------------------\n")
cat(sprintf("n nodes: %d\n", length(toy$perm)))
cat(sprintf("target permutation: %s\n", paste(toy$perm, collapse = ",")))
cat(sprintf("predicted (grasp):  %s\n", paste(as.integer(grasp_out$assignment), collapse = ",")))
cat(sprintf("predicted (raw):    %s\n", paste(as.integer(raw_pred), collapse = ",")))
cat(sprintf("accuracy (grasp):   %.3f\n", acc_grasp))
cat(sprintf("accuracy (raw):     %.3f\n", acc_raw))
cat(sprintf("random expectation: %.3f\n", acc_random_expectation))

# Sanity threshold: should be comfortably above chance on this tiny deterministic case.
if (acc_grasp < 0.90) {
  stop(sprintf("Sanity check failed: grasp accuracy %.3f < 0.90", acc_grasp), call. = FALSE)
}

cat("Status: PASS\n")
