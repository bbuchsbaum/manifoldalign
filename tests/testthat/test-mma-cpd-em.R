# MMA hard assignments, CPD posteriors and EM rotation convention.

posterior_from_map <- function(ref_to_g, n_g, hit = 0.9, base = 0.01) {
  A <- matrix(base, length(ref_to_g), n_g)
  for (i in seq_along(ref_to_g)) {
    if (!is.na(ref_to_g[i])) A[i, ref_to_g[i]] <- hit
  }
  A
}

test_that("hungarian and argmax both map domain rows to reference rows (square)", {
  ref_to_g <- c(3L, 1L, 5L, 2L, 4L)
  A <- posterior_from_map(ref_to_g, 5L)
  g_to_ref <- match(seq_len(5), ref_to_g)            # 2 4 1 5 3

  hung <- manifoldalign:::hungarian_from_posterior(A)
  amax <- manifoldalign:::build_final_assignment("argmax", list(A), "consensus")[[1]]
  expect_identical(hung, g_to_ref)
  expect_identical(amax, g_to_ref)
})

test_that("hungarian handles fewer domain rows than reference rows", {
  # 6 reference rows, 4 domain rows; domain row j belongs to reference row truth[j]
  truth <- c(5L, 2L, 6L, 1L)
  A <- matrix(0.01, 6, 4)
  A[cbind(truth, seq_along(truth))] <- 0.9
  hung <- manifoldalign:::hungarian_from_posterior(A)
  expect_identical(hung, truth)
  expect_length(hung, ncol(A))
  expect_false(anyDuplicated(hung) > 0)
})

test_that("hungarian handles more domain rows than reference rows", {
  # 3 reference rows, 5 domain rows; rows 2 and 4 of the domain have no partner
  truth <- c(3L, NA, 1L, NA, 2L)
  A <- matrix(0.01, 3, 5)
  ok <- !is.na(truth)
  A[cbind(truth[ok], which(ok))] <- 0.9
  hung <- manifoldalign:::hungarian_from_posterior(A)
  expect_identical(hung, truth)
})

test_that("reference domain assignment is the identity", {
  set.seed(4)
  mats <- lapply(1:3, function(i) matrix(rnorm(15 * 3), 15, 3))
  res <- suppressWarnings(mma_align_multiple(mats, ncomp = 3, em_max_iter = 5,
                                             final_assignment = "hungarian"))
  expect_identical(res$assignment[[1]], seq_len(15))
  for (g in 2:3) {
    expect_length(res$assignment[[g]], 15)
    expect_true(all(res$assignment[[g]] %in% seq_len(15)))
  }
})

test_that("CPD E-step matches a hand-computed Gaussian mixture with uniform outliers", {
  set.seed(21)
  K <- 3; M <- 4; N <- 6
  Y <- matrix(rnorm(M * K), M, K)     # centroids
  X <- matrix(rnorm(N * K), N, K)     # data
  D2 <- outer(rowSums(Y^2), rowSums(X^2), "+") - 2 * tcrossprod(Y, X)
  for (sigma2 in c(0.05, 0.7, 4)) {
    for (w in c(0, 0.1, 0.4)) {
      dens <- matrix(0, M, N)
      for (m in seq_len(M)) for (n in seq_len(N)) {
        dens[m, n] <- (1 - w) / M * prod(dnorm(X[n, ], Y[m, ], sqrt(sigma2)))
      }
      total <- colSums(dens) + w / N
      P_ref <- sweep(dens, 2, total, "/")
      ll_ref <- sum(log(total))

      es <- manifoldalign:::.mma_cpd_estep(D2, sigma2, K = K, w = w)
      expect_equal(es$alpha, P_ref, tolerance = 1e-10)
      expect_equal(es$loglik, ll_ref, tolerance = 1e-10)
    }
  }
})

test_that("CPD E-step is stable when sigma2 is tiny", {
  D2 <- matrix(c(0, 4, 4, 0.01), 2, 2)
  es <- manifoldalign:::.mma_cpd_estep(D2, sigma2 = 1e-6, K = 2, w = 0.01)
  expect_true(all(is.finite(es$alpha)))
  expect_true(is.finite(es$loglik))
  expect_true(all(colSums(es$alpha) <= 1 + 1e-12))
})

test_that("EM rotation step is consistent with its own E-step model", {
  # X_g = X_ref Q  <=>  X_g ~ X_ref t(R) with R = t(Q). Started at the truth the
  # EM must stay there instead of flipping to the transpose.
  set.seed(7)
  n <- 30; K <- 3
  X_ref <- matrix(rnorm(n * K), n, K)
  Q <- qr.Q(qr(matrix(rnorm(K * K), K, K)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  expect_gt(max(abs(Q - t(Q))), 0.1)                 # genuinely non-symmetric
  X_g <- X_ref %*% Q

  em <- manifoldalign:::em_mma_pairwise_R(X_ref, X_g, R_init = t(Q), sigma2_0 = 0.01,
                                          outlier = 0.01, max_iter = 20, tol = 1e-12)
  expect_equal(em$R, t(Q), tolerance = 1e-6)
  expect_identical(apply(em$alpha, 2, which.max), seq_len(n))

  # A rotation that is close to the truth converges to it.
  E <- matrix(c(0, -0.08, 0, 0.08, 0, 0, 0, 0, 0), K, K)
  sv <- svd(diag(K) + E)
  R_near <- t(Q) %*% sv$u %*% t(sv$v)
  em2 <- manifoldalign:::em_mma_pairwise_R(X_ref, X_g, R_init = R_near, sigma2_0 = 0.05,
                                           outlier = 0.01, max_iter = 200, tol = 1e-12)
  expect_equal(em2$R, t(Q), tolerance = 1e-4)
})

test_that("Procrustes step returns argmax over O(K) and respects force_SO", {
  set.seed(3)
  K <- 4
  M <- matrix(rnorm(K * K), K, K)
  R <- manifoldalign:::.mma_procrustes_rotation(M)
  expect_equal(crossprod(R), diag(K), tolerance = 1e-10)
  # tr(R M) is maximal: no random orthogonal matrix does better
  best <- sum(diag(R %*% M))
  others <- replicate(200, {
    O <- qr.Q(qr(matrix(rnorm(K * K), K, K)))
    sum(diag(O %*% M))
  })
  expect_true(all(others <= best + 1e-10))
  R_so <- manifoldalign:::.mma_procrustes_rotation(M, force_SO = TRUE)
  expect_equal(det(R_so), 1, tolerance = 1e-10)
})

test_that("EM convergence flag is truthful at the final iteration", {
  set.seed(9)
  n <- 25; K <- 2
  X_ref <- matrix(rnorm(n * K), n, K)
  X_g <- X_ref + matrix(rnorm(n * K, sd = 0.05), n, K)
  run <- function(max_iter) {
    manifoldalign:::em_mma_pairwise_R(X_ref, X_g, R_init = diag(K), sigma2_0 = 0.2,
                                      outlier = 0.05, max_iter = max_iter, tol = 1e-6)
  }
  free <- run(500)
  expect_true(free$converged)
  it <- free$iterations
  expect_gte(it, 2)
  expect_true(run(it)$converged)          # converges exactly on the last allowed iteration
  expect_false(run(it - 1L)$converged)    # budget exhausted before the criterion
})

test_that("eigensignature initial rotation maps the reference basis onto the domain basis", {
  set.seed(11)
  n <- 400
  Uref <- cbind(rexp(n) / 10, runif(n, -0.2, 0.05), rbeta(n, 5, 1) / 10, rnorm(n, 0.03, 0.01))
  perm <- c(3L, 1L, 4L, 2L)
  s <- c(-1, 1, 1, -1)
  Ug <- matrix(0, n, 4)
  for (i in 1:4) Ug[, perm[i]] <- s[i] * Uref[, i]
  for (meth in c("hist", "w2")) {
    ag <- manifoldalign:::eigensignature_align(Uref, Ug, method = meth)
    expect_identical(ag$perm, perm)
    R0 <- manifoldalign:::.mma_signature_rotation(ag)
    # MMA convention: E_g ~ E_ref t(R), aligned scores E_g R
    expect_equal(Uref %*% t(R0), Ug, tolerance = 1e-12)
    expect_equal(Ug %*% R0, Uref, tolerance = 1e-12)
  }
})
