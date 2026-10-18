# GRASP base alignment (lambda, max_iter, tol are live), solver truthfulness,
# and multiset robustness to unequal domain sizes.

grasp_fixture <- function(n = 30, seed = 2, ncomp = 8, q = 15) {
  set.seed(seed)
  perm <- sample(n)
  X1 <- matrix(rnorm(n * 4), n, 4)
  X1[, 1] <- X1[, 1] + 2 * cos(2 * pi * seq_len(n) / n)
  X1[, 2] <- X1[, 2] + 2 * sin(2 * pi * seq_len(n) / n)
  X2 <- X1[perm, ] + matrix(rnorm(n * 4, sd = 0.02), n, 4)
  strata <- list(list(x = X1), list(x = X2))
  bases <- manifoldalign:::compute_grasp_basis(strata, ncomp, TRUE)
  desc <- manifoldalign:::compute_grasp_descriptors(bases, q, 0.73)
  list(bases = bases, desc = desc)
}

base_objective <- function(M, fx, lambda) {
  M <- as.matrix(M)
  lam2 <- fx$bases[[2]]$values
  Fm <- as.matrix(Matrix::crossprod(fx$desc[[1]], fx$bases[[1]]$vectors))
  Gm <- as.matrix(Matrix::crossprod(fx$desc[[2]], fx$bases[[2]]$vectors))
  S <- crossprod(M, lam2 * M)
  diag(S) <- 0
  sum(S^2) + lambda * sum((Fm - Gm %*% M)^2)
}

test_that("base alignment minimises the documented objective and reports it truthfully", {
  fx <- grasp_fixture()
  start <- manifoldalign:::align_grasp_bases(fx$bases[[1]], fx$bases[[2]],
                                             fx$desc[[1]], fx$desc[[2]],
                                             lambda = 0.5, max_iter = 0)
  expect_identical(start$iterations, 0L)
  expect_false(start$converged)
  M0 <- as.matrix(start$rotation)
  expect_true(all(M0 %in% c(-1, 0, 1)))                 # signed permutation start
  expect_equal(start$objective, base_objective(M0, fx, 0.5))

  fit <- manifoldalign:::align_grasp_bases(fx$bases[[1]], fx$bases[[2]],
                                           fx$desc[[1]], fx$desc[[2]],
                                           lambda = 0.5, max_iter = 500, tol = 1e-10)
  M <- as.matrix(fit$rotation)
  expect_equal(crossprod(M), diag(ncol(M)), tolerance = 1e-10)
  expect_equal(fit$objective, base_objective(M, fx, 0.5), tolerance = 1e-10)
  expect_equal(fit$objective_init, start$objective)
  expect_lt(fit$objective, fit$objective_init)
  expect_gte(fit$iterations, 1L)
  expect_lte(fit$iterations, 500L)

  # lambda changes the solution (it is no longer ignored)
  fit_small <- manifoldalign:::align_grasp_bases(fx$bases[[1]], fx$bases[[2]],
                                                 fx$desc[[1]], fx$desc[[2]],
                                                 lambda = 0.01, max_iter = 500, tol = 1e-10)
  expect_gt(max(abs(as.matrix(fit_small$rotation) - M)), 1e-6)

  # a tighter iteration budget is reported as not converged
  short <- manifoldalign:::align_grasp_bases(fx$bases[[1]], fx$bases[[2]],
                                             fx$desc[[1]], fx$desc[[2]],
                                             lambda = 0.5, max_iter = 1, tol = 1e-14)
  expect_identical(short$iterations, 1L)
  expect_false(short$converged)
})

test_that("base alignment with lambda = 0 stops at the signed permutation", {
  fx <- grasp_fixture(seed = 5)
  fit <- manifoldalign:::align_grasp_bases(fx$bases[[1]], fx$bases[[2]],
                                           fx$desc[[1]], fx$desc[[2]], lambda = 0)
  expect_true(fit$converged)
  expect_equal(fit$objective, 0)
  expect_true(all(as.matrix(fit$rotation) %in% c(-1, 0, 1)))
})

test_that("grasp() exposes the base-alignment diagnostics", {
  set.seed(6)
  X1 <- matrix(rnorm(25 * 3), 25, 3)
  X2 <- X1[sample(25), ] + rnorm(75, sd = 0.01)
  hd <- create_alignment_hyperdesign(X1, X2)
  res <- suppressWarnings(grasp(hd, ncomp = 6, q_descriptors = 10, lambda = 0.2,
                                max_iter = 50))
  ba <- res$base_alignment
  expect_true(is.logical(ba$converged))
  expect_lte(ba$iterations, 50L)
  expect_lte(ba$objective, ba$objective_init)
  expect_error(grasp(hd, ncomp = 6, q_descriptors = 10, solver = "auction"), "not available")
  expect_s3_class(suppressWarnings(grasp(hd, ncomp = 6, q_descriptors = 10,
                                         solver = "hungarian")), "grasp")
})

test_that("LSAP helper accepts only exact solvers", {
  C <- matrix(c(3, 1, 2, 2, 3, 1, 1, 2, 3), 3, 3)
  expect_identical(manifoldalign:::solve_lsap_with_padding(C), c(3L, 1L, 2L))
  expect_identical(manifoldalign:::solve_lsap_with_padding(C, method = "dense"), c(3L, 1L, 2L))
  expect_error(manifoldalign:::solve_lsap_with_padding(C, method = "auction"), "not available")
  expect_error(manifoldalign:::compute_grasp_assignment(NULL, NULL, NULL, NULL, NULL,
                                                        solver_method = "auction"),
               "not available")
})

test_that("grasp_multiset votes across four domains of unequal size", {
  set.seed(3)
  ms <- list(matrix(rnorm(20 * 3), 20, 3), matrix(rnorm(18 * 3), 18, 3),
             matrix(rnorm(22 * 3), 22, 3), matrix(rnorm(19 * 3), 19, 3))
  res <- suppressWarnings(grasp_multiset(ms, ncomp = 4, q_descriptors = 8, max_iter = 5))
  n_anchor <- 20
  for (s in seq_along(ms)) {
    p <- res$permutations[[s]]
    expect_length(p, n_anchor)
    ok <- p[!is.na(p)]
    expect_true(all(ok >= 1 & ok <= nrow(ms[[s]])))
    expect_false(anyDuplicated(ok) > 0)
    expect_equal(sum(!is.na(p)), min(n_anchor, nrow(ms[[s]])))
    # anchor-ordered embeddings: matched rows copied, unmatched rows NA
    E <- res$embeddings[[s]]
    expect_equal(nrow(E), n_anchor)
    expect_equal(E[!is.na(p), ], res$embeddings_raw[[s]][ok, ], ignore_attr = TRUE)
    expect_true(all(is.na(E[is.na(p), ])))
  }
})

test_that("grasp_multiset pairwise fallback accepts lambda = 0", {
  set.seed(4)
  ms <- lapply(1:3, function(i) matrix(rnorm(15 * 3), 15, 3))
  res <- suppressWarnings(grasp_multiset(ms, ncomp = 4, q_descriptors = 8,
                                         max_iter = 1, lambda = 0))
  expect_false(res$converged)
  expect_identical(res$rotation_source, "pairwise_fallback")
})

test_that("grasp_multiset is quiet unless verbose", {
  set.seed(5)
  ms <- lapply(1:2, function(i) matrix(rnorm(15 * 3), 15, 3))
  expect_no_message(suppressWarnings(grasp_multiset(ms, ncomp = 3, q_descriptors = 6,
                                                    max_iter = 5)))
  expect_message(suppressWarnings(grasp_multiset(ms, ncomp = 3, q_descriptors = 6,
                                                 max_iter = 5, verbose = TRUE)),
                 "Only 2 domains")
  expect_error(grasp_multiset(ms, ncomp = 3, q_descriptors = 6, solver = "auction"),
               "not available")
})
