test_that("domain convergence is invariant to cost scale and feasible-set constants", {
  set.seed(742)
  M <- matrix(rnorm(144), 12, 12); A <- crossprod(M)
  reduced <- function(H) list(A_r = H, ranks = rep(4L, 3),
                             basis = rep(list(list(T = diag(4))), 3))
  fit <- function(H) .ssma_solve_domain_scores(reduced(H), 2, 1e-8, 30000)
  original_cost <- function(z) sum(z$alpha * (A %*% z$alpha))
  base <- fit(A)
  for (scale in c(1e-12, 1e12)) {
    z <- fit(A * scale)
    expect_equal(original_cost(z), original_cost(base), tolerance = 1e-6)
    expect_gt(z$domain_diagnostics$iterations, 1)
    expect_equal(z$domain_diagnostics$objective / scale, original_cost(z), tolerance = 1e-7)
  }
  shifts <- rep(c(1e8, 2e8, 3e8), each = 4)
  z <- fit(A + diag(shifts))
  expect_equal(original_cost(z), original_cost(base), tolerance = 1e-6)
  expect_equal(z$domain_diagnostics$objective - 2 * sum(c(1e8, 2e8, 3e8)),
               original_cost(z), tolerance = 1e-6)
})

test_that("public SSMA preserves its solution when all linear costs are scaled", {
  set.seed(742)
  x <- list(a = diag(4), b = diag(4), c = diag(4))
  h <- lapply(x, function(z) matrix(rnorm(40), 10, 4))
  links <- function(s) list(linear_correspondence(c("a", "b"), h[1:2], strength = s, name = "ab"),
                           linear_correspondence(c("a", "c"), h[c(1, 3)], strength = s, name = "ac"))
  control <- ssma_align_control(compression = "none", normalization = "domain_scores",
                               auto_correspondences = FALSE, domain_tol = 1e-8, domain_maxiter = 30000)
  fit <- function(s) ssma_align(x, preproc = NULL, ncomp = 2, control = control,
                               constraints = links(s),
                               within_graphs = lapply(x, function(z) Matrix::Matrix(0, 4, 4, sparse = TRUE)))
  cost <- function(f) sum(vapply(f$constraint_diagnostics, function(z) z$residual, numeric(1)))
  a <- fit(1); b <- fit(1e-12)
  expect_equal(cost(a), cost(b), tolerance = 1e-6)
  expect_gt(b$domain_diagnostics$iterations, 1)
})

test_that("full-rank domain alignment attains the exact two-domain Procrustes oracle", {
  for (seed in 1:20) {
    set.seed(seed)
    M <- matrix(rnorm(64), 8, 8); A <- crossprod(M)
    z <- .ssma_solve_domain_scores(list(A_r = A, ranks = c(4L, 4L),
                  basis = list(list(T = diag(4)), list(T = diag(4)))), 4, 1e-8, 30000)
    oracle <- sum(diag(A)) - 2 * sum(svd(A[1:4, 5:8])$d)
    expect_equal(z$domain_diagnostics$objective, oracle, tolerance = 1e-7, info = paste("seed", seed))
    expect_equal(sum(z$alpha * (A %*% z$alpha)), oracle, tolerance = 1e-7)
    expect_true(all(diff(z$domain_diagnostics$objective_history) <= 1e-10))
  }
  C <- diag(c(3, 2))
  H <- rbind(cbind(matrix(0, 2, 2), C), cbind(C, matrix(0, 2, 2)))
  Y <- rbind(diag(2), diag(c(1, -1)))
  z <- .ssma_update_square_blocks(Y, H, list(1:2, 3:4))
  expect_gt(z$updates, 0)
  expect_equal(sum(z$Y * (H %*% z$Y)), -10, tolerance = 1e-10)
})

test_that("constant domain objectives retain correct costs and normalized scores", {
  z <- .ssma_solve_domain_scores(list(A_r = diag(rep(c(7, 5), each = 3)), ranks = c(3L, 3L),
                  basis = list(list(T = diag(3)), list(T = diag(3)))), 2, 1e-8, 100)
  expect_equal(z$domain_diagnostics$objective, 24)
  expect_equal(crossprod(z$alpha[1:3, ]), diag(2), tolerance = 1e-10)
  expect_equal(crossprod(z$alpha[4:6, ]), diag(2), tolerance = 1e-10)
})

test_that("global SSMA retains exact nonconstant zero-cost correspondence modes", {
  set.seed(804)
  X <- scale(matrix(rnorm(96), 16, 6), center = TRUE, scale = FALSE)
  x <- list(a = X, b = X)
  corr <- data.frame(domain_i = 1L, index_i = 1:16, domain_j = 2L, index_j = 1:16)
  for (solver in c("reduced", "operator")) {
    f <- ssma_align(x, preproc = NULL, correspondences = corr,
          within_graphs = lapply(x, function(z) Matrix::Matrix(0, 16, 16, sparse = TRUE)),
          ncomp = 2, control = ssma_align_control(compression = "none",
          auto_correspondences = FALSE, solver = solver, operator_tol = 1e-10))
    a <- f$s[1:16, , drop = FALSE]; b <- f$s[17:32, , drop = FALSE]
    expect_lt(sum((a - b)^2) / (sum(a^2) + sum(b^2)), 1e-12)
    expect_equal(ncol(f$s), 2)
    expect_gt(sum(a^2), 0)
  }
})

test_that("singleton labels retain their true sample identity through the public API", {
  set.seed(777)
  x <- list(a = matrix(rnorm(32), 8, 4), b = matrix(rnorm(32), 8, 4))
  labs <- list(letters[1:8], letters[1:8])
  hd <- as_hyperdesign(x, labels = labs, domain_names = names(x))
  f <- ssma_align(hd, y = label, preproc = NULL, ncomp = 2,
      control = ssma_align_control(compression = "none", normalization = "domain_scores",
                                  domain_tol = 1e-7, domain_maxiter = 30000, seed = 1),
      within_graphs = lapply(x, function(z) Matrix::Matrix(0, 8, 8, sparse = TRUE)))
  cc <- f$correspondences
  expect_equal(nrow(cc), 8)
  expect_identical(labs[[1]][cc$index_i], labs[[2]][cc$index_j])
})

test_that("invalid explicit correspondence identity is rejected before automatic matching", {
  set.seed(32)
  x <- list(a = matrix(rnorm(32), 8, 4), b = matrix(rnorm(32), 8, 4))
  good <- data.frame(domain_i = 1, index_i = 3, domain_j = 2, index_j = 3)
  for (field in names(good)) {
    for (value in c(NA, Inf, 0, 1.5, 1e20)) {
      bad <- good; bad[[field]] <- value
      expect_error(ssma_align(x, correspondences = bad, preproc = NULL, ncomp = 2),
                   "Correspondence")
    }
  }
  for (field in c("index_i", "index_j", "domain_i", "domain_j")) {
    bad <- good; bad[[field]] <- 100
    expect_error(ssma_align(x, correspondences = bad, preproc = NULL, ncomp = 2), "in range")
  }
  bad <- good; bad$domain_j <- 1
  expect_error(ssma_align(x, correspondences = bad, preproc = NULL, ncomp = 2), "different domains")
})

test_that("zero-mode numerical validation uses the full operator scale", {
  # A partial eigensolve may return only near-zero modes, while the omitted
  # spectrum is large. Its tolerance must not be inferred from returned modes.
  expect_identical(.ssma_eigenpair_indices(c(-1e-7, 0, 1e-8), 2, 1e-8, 1e6), 1:2)
  expect_error(.ssma_eigenpair_indices(c(-1, 0), 1, 1e-8, 1), "negative")
})
