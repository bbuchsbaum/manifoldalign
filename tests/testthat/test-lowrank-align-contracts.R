# Contracts of lowrank_align(): preprocessing lists, cross-domain edge
# filtering, sv_thresh semantics, eigensolver fallback and RNG hygiene.

lowrank_fixture <- function(seed = 7) {
  set.seed(seed)
  X1 <- matrix(rnorm(6 * 3), 6, 3)
  X2 <- matrix(rnorm(5 * 4), 5, 4)
  labels1 <- c("L1", "L2", "L3", "L1", NA, "L2")
  labels2 <- c("L2", "L1", "L3", NA, "L1")
  pool <- c("L1", "L2", "L3")
  S <- diag(3)
  S[1, 2] <- S[2, 1] <- 0.4
  dimnames(S) <- list(pool, pool)
  list(
    X1 = X1, X2 = X2,
    simfun = createSimFun(S),
    hd = function(X1, X2) {
      multidesign::hyperdesign(list(
        b1 = multidesign::multidesign(X1, data.frame(y = labels1)),
        b2 = multidesign::multidesign(X2, data.frame(y = labels2))
      ))
    }
  )
}

# Eigenvector columns agree up to sign.
expect_same_columns_up_to_sign <- function(A, B, tolerance = 1e-6) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  expect_equal(dim(A), dim(B))
  for (j in seq_len(ncol(A))) {
    s <- sign(sum(A[, j] * B[, j]))
    expect_equal(A[, j], s * B[, j], tolerance = tolerance)
  }
}

fit_lra <- function(hd, simfun, ...) {
  suppressMessages(lowrank_align(hd, y, simfun = simfun, ncomp = 2, lambda = 0.01, ...))
}

test_that("lowrank_align accepts preproc = NULL and NULL list entries", {
  skip_if_not_installed("multidesign")
  fx <- lowrank_fixture()
  X1c <- sweep(fx$X1, 2, colMeans(fx$X1))
  X2c <- sweep(fx$X2, 2, colMeans(fx$X2))
  hd_c <- fx$hd(X1c, X2c)

  fit_center <- fit_lra(hd_c, fx$simfun, preproc = multivarious::center())
  # Regression: `proclist[[i]] <- NULL` shrank the list and `names<-` failed.
  fit_null <- fit_lra(hd_c, fx$simfun, preproc = NULL)
  fit_mixed <- fit_lra(hd_c, fx$simfun, preproc = list(multivarious::center(), NULL))

  # Already-centred data: centring is a no-op, so all three fits coincide.
  expect_same_columns_up_to_sign(fit_null$s, fit_center$s)
  expect_same_columns_up_to_sign(fit_mixed$s, fit_center$s)
  expect_equal(nrow(fit_null$v), ncol(fx$X1) + ncol(fx$X2))

  # The stored preprocessor has one identity block per unprocessed domain.
  Xnew <- cbind(fx$X1[1:2, ], fx$X2[1:2, ])
  expect_equal(unname(multivarious::transform(fit_null$preproc, Xnew)), unname(Xnew))
  mixed_tf <- multivarious::transform(fit_mixed$preproc, Xnew)
  expect_equal(unname(mixed_tf[, 4:7]), unname(Xnew[, 4:7]))
  expect_equal(unname(mixed_tf[, 1:3]),
               unname(sweep(Xnew[, 1:3], 2, colMeans(X1c))), tolerance = 1e-12)
})

test_that("cross-domain filtering of sparse edges matches the dense mask", {
  set.seed(11)
  sizes <- c(4L, 3L, 5L)
  n <- sum(sizes)
  C <- abs(Matrix::rsparsematrix(n, n, density = 0.4, symmetric = TRUE))
  Matrix::diag(C) <- 0
  block_id <- rep.int(seq_along(sizes), sizes)

  # Oracle: the previous construction (cross-domain mask by sub-assignment).
  offsets <- cumsum(c(0L, sizes))
  C_mask <- Matrix::Matrix(0, n, n, sparse = TRUE)
  for (a in seq_along(sizes)) {
    ia <- (offsets[a] + 1L):offsets[a + 1L]
    for (b in seq_along(sizes)) {
      if (a == b) next
      jb <- (offsets[b] + 1L):offsets[b + 1L]
      C_mask[ia, jb] <- 1
    }
  }
  expected <- as.matrix(C * C_mask)

  got <- manifoldalign:::.lowrank_cross_domain_edges(C, block_id)
  expect_s4_class(got, "sparseMatrix")
  expect_equal(as.matrix(got), expected, ignore_attr = TRUE)
  expect_true(Matrix::isSymmetric(got))
  expect_lte(Matrix::nnzero(got), Matrix::nnzero(C))

  # A general (non-symmetric-class) input gives the same answer.
  got_general <- manifoldalign:::.lowrank_cross_domain_edges(methods::as(C, "generalMatrix"), block_id)
  expect_equal(as.matrix(got_general), expected, ignore_attr = TRUE)

  # A single domain has no cross-domain edges.
  one <- manifoldalign:::.lowrank_cross_domain_edges(C, rep(1L, n))
  expect_equal(Matrix::nnzero(one), 0)
})

test_that("sv_thresh below 1 behaves exactly like sv_thresh = 1; invalid values error", {
  skip_if_not_installed("multidesign")
  fx <- lowrank_fixture()
  hd <- fx$hd(fx$X1, fx$X2)
  f1 <- fit_lra(hd, fx$simfun, sv_thresh = 1)
  f05 <- fit_lra(hd, fx$simfun, sv_thresh = 0.5)
  expect_equal(f05$s, f1$s)

  # The R factor keeps exactly the singular values above max(sv_thresh, 1).
  X <- diag(c(3, 1.5, 0.5))
  fac <- manifoldalign:::.compute_R_factors(X, sv_thresh = 0.2)
  expect_equal(sort(fac$a), sort(1 - 1 / c(3, 1.5)^2))
  fac2 <- manifoldalign:::.compute_R_factors(X, sv_thresh = 2)
  expect_equal(fac2$a, 1 - 1 / 9)

  expect_error(lowrank_align(hd, y, simfun = fx$simfun, sv_thresh = -1), "sv_thresh")
  expect_error(lowrank_align(hd, y, simfun = fx$simfun, sv_thresh = NA_real_), "sv_thresh")
})

test_that("explicit solver gives the same embedding without PRIMME", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("PRIMME")
  fx <- lowrank_fixture()
  hd <- fx$hd(fx$X1, fx$X2)
  with_primme <- fit_lra(hd, fx$simfun, scale_M = TRUE)
  local_mocked_bindings(.primme_available = function() FALSE)
  without <- fit_lra(hd, fx$simfun, scale_M = TRUE)
  expect_same_columns_up_to_sign(without$s, with_primme$s)

  # The fallback returns the smallest/largest eigenpairs in the right order.
  set.seed(3)
  A <- crossprod(matrix(rnorm(64), 8, 8))
  ev <- eigen(A, symmetric = TRUE)$values
  sa <- manifoldalign:::.lowrank_eigs_sym(A, 3, which = "SA")
  la <- manifoldalign:::.lowrank_eigs_sym(A, 2, which = "LA")
  expect_equal(sa$values, rev(ev)[1:3])
  expect_equal(la$values, ev[1:2])
  expect_equal(A %*% sa$vectors, sa$vectors %*% diag(sa$values), tolerance = 1e-8)
})

test_that("operator solver leaves the caller's RNG stream untouched", {
  skip_if_not_installed("multidesign")
  fx <- lowrank_fixture()
  hd <- fx$hd(fx$X1, fx$X2)
  set.seed(99)
  expected <- runif(3)
  set.seed(99)
  fit1 <- fit_lra(hd, fx$simfun, solver = "operator")
  expect_equal(runif(3), expected)
  fit2 <- fit_lra(hd, fx$simfun, solver = "operator")
  expect_equal(fit1$s, fit2$s)
})
