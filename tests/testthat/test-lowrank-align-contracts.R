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

# ---------------------------------------------------------------------------
# Trivial (per-domain constant) directions are excluded from the embedding.
# ---------------------------------------------------------------------------

# Independent dense oracle: centre each domain, build R_k = U diag(1 - 1/s^2) U'
# from the left singular vectors with s > 1, M = (I - R)'(I - R), the
# cross-domain Laplacian L, Z = (1 - mu) M + 2 mu L, and take the smallest
# eigenpairs of Z restricted to the orthogonal complement of the per-domain
# constant vectors.
lra_dense_oracle <- function(Xs, labels, simfun, mu, ncomp) {
  Xs <- lapply(Xs, function(X) sweep(X, 2, colMeans(X)))
  sizes <- vapply(Xs, nrow, integer(1))
  n <- sum(sizes)
  dom <- rep(seq_along(sizes), sizes)
  Rk <- lapply(Xs, function(X) {
    s <- svd(X)
    keep <- s$d > 1
    U <- s$u[, keep, drop = FALSE]
    U %*% diag(1 - 1 / s$d[keep]^2, sum(keep)) %*% t(U)
  })
  R <- matrix(0, n, n)
  for (k in seq_along(sizes)) R[dom == k, dom == k] <- Rk[[k]]
  M <- crossprod(diag(n) - R)
  C <- as.matrix(simfun(labels))
  C[C < 0] <- 0
  C <- (C + t(C)) / 2
  C[outer(dom, dom, "==")] <- 0
  L <- diag(rowSums(C)) - C
  Z <- (1 - mu) * M + 2 * mu * L
  E <- vapply(seq_along(sizes), function(k) as.numeric(dom == k), numeric(n))
  N <- qr.Q(qr(E), complete = TRUE)[, -seq_along(sizes), drop = FALSE]
  e <- eigen(crossprod(N, Z %*% N), symmetric = TRUE)
  idx <- order(e$values)[seq_len(ncomp)]
  list(values = e$values[idx], vectors = N %*% e$vectors[, idx, drop = FALSE],
       all_values = sort(e$values), Z = Z)
}

# Per-domain column sums of the scores (zero when orthogonal to every e_k).
domain_sums <- function(S, sizes) {
  rowsum(as.matrix(S), rep(seq_along(sizes), sizes))
}

test_that("lowrank_align matches a dense oracle restricted to the non-constant complement", {
  skip_if_not_installed("multidesign")
  fx <- lowrank_fixture()
  hd <- fx$hd(fx$X1, fx$X2)
  labels <- c(c("L1", "L2", "L3", "L1", NA, "L2"), c("L2", "L1", "L3", NA, "L1"))
  sizes <- c(6L, 5L)

  for (mu in c(0.2, 0.5, 0.8)) {
    orc <- lra_dense_oracle(list(fx$X1, fx$X2), labels, fx$simfun, mu = mu, ncomp = 2)
    # The selected oracle eigenvalues are separated from their neighbours.
    expect_gt(min(diff(orc$all_values[1:3])), 1e-3)

    fit_exp <- fit_lra(hd, fx$simfun, mu = mu, solver = "explicit")
    fit_op <- fit_lra(hd, fx$simfun, mu = mu, solver = "operator")
    expect_same_columns_up_to_sign(fit_exp$s, orc$vectors, tolerance = 1e-5)
    expect_same_columns_up_to_sign(fit_op$s, orc$vectors, tolerance = 1e-6)
    expect_same_columns_up_to_sign(fit_op$s, fit_exp$s, tolerance = 1e-5)

    # Without PRIMME the dense fallback solves the same restricted problem.
    local({
      local_mocked_bindings(.primme_available = function() FALSE)
      fit_dense <- fit_lra(hd, fx$simfun, mu = mu, solver = "explicit")
      expect_same_columns_up_to_sign(fit_dense$s, orc$vectors, tolerance = 1e-8)
    })

    # Every component has zero mean within each domain.
    expect_lt(max(abs(domain_sums(fit_exp$s, sizes))), 1e-8)
    expect_lt(max(abs(domain_sums(fit_op$s, sizes))), 1e-8)
  }
})

test_that("restricted eigen helper agrees across PRIMME, dense and RSpectra paths", {
  set.seed(5)
  sizes <- c(9L, 14L, 7L)
  n <- sum(sizes)
  A <- crossprod(matrix(rnorm(n * n), n)) / n
  A <- Matrix::Matrix(A, sparse = TRUE)
  dom <- rep(seq_along(sizes), sizes)
  Q <- manifoldalign:::.lowrank_constant_basis(dom)
  expect_s4_class(Q, "sparseMatrix")
  expect_equal(as.matrix(Matrix::crossprod(Q)), diag(3), ignore_attr = TRUE)
  # The projector removes per-domain means for vectors and matrices.
  X <- matrix(rnorm(n * 2), n, 2)
  expect_equal(manifoldalign:::.lowrank_project_complement(X, dom),
               X - as.matrix(Q %*% Matrix::crossprod(Q, X)), ignore_attr = TRUE)
  expect_equal(manifoldalign:::.lowrank_project_complement(X[, 1], dom),
               as.numeric(X[, 1] - Q %*% Matrix::crossprod(Q, X[, 1])))

  E <- vapply(seq_along(sizes), function(k) as.numeric(dom == k), numeric(n))
  N <- qr.Q(qr(E), complete = TRUE)[, -(1:3)]
  ref <- eigen(crossprod(N, as.matrix(A) %*% N), symmetric = TRUE)
  ref_vals <- rev(ref$values)[1:4]
  ref_vecs <- N %*% ref$vectors[, rev(seq_len(n - 3))[1:4]]

  sigma <- max(Matrix::rowSums(abs(A))) + 1
  check <- function(res, tol) {
    expect_equal(res$values, ref_vals, tolerance = tol)
    expect_lt(max(abs(as.matrix(Matrix::crossprod(Q, res$vectors)))), 1e-10)
    expect_same_columns_up_to_sign(res$vectors, ref_vecs, tolerance = tol)
  }
  local({
    local_mocked_bindings(.primme_available = function() FALSE)
    check(manifoldalign:::.lowrank_eigs_restricted(A, 4, dom, sigma), 1e-10)
    skip_if_not_installed("RSpectra")
    check(manifoldalign:::.lowrank_eigs_restricted(A, 4, dom, sigma, dense_max = 0L), 1e-7)
  })
  skip_if_not_installed("PRIMME")
  check(manifoldalign:::.lowrank_eigs_restricted(A, 4, dom, sigma), 1e-5)
})

test_that("lowrank_align never returns the constant eigenvector (toy linear fixture)", {
  skip_if_not_installed("multidesign")
  # Regression: on this two-class, two-feature fixture the label term
  # dominates, the constant vector (an exact eigenvector of Z with eigenvalue
  # 1 - mu) was the second smallest eigenpair, and component 2 was constant at
  # -1/sqrt(120) with sdev ~ 1e-15.
  ld <- manifoldalign:::synthetic_alignment_scenario(
    "linear_affine", profile = "full", seed = 20240502, n = 60, noise = 0.01
  )
  labels <- as.character(ld$true_labels)
  hd <- multidesign::hyperdesign(list(
    view1 = multidesign::multidesign(ld$view1, data.frame(y = labels)),
    view2 = multidesign::multidesign(ld$view2, data.frame(y = labels))
  ))
  simfun <- function(lbl) {
    lbl <- as.character(lbl)
    S <- outer(lbl, lbl, "==")
    S[is.na(S)] <- FALSE
    Matrix::Matrix(S * 1, sparse = TRUE)
  }
  sizes <- c(60L, 60L)
  mu <- 0.5

  # The mechanism: with centred features 1 is an exact eigenvector of Z and is
  # among its two smallest eigenpairs.
  all_labels <- c(labels, labels)
  orc <- lra_dense_oracle(list(ld$view1, ld$view2), all_labels, simfun, mu = mu, ncomp = 2)
  one <- rep(1, 120)
  expect_equal(as.numeric(orc$Z %*% one), (1 - mu) * one, tolerance = 1e-10)
  z_vals <- sort(eigen(orc$Z, symmetric = TRUE, only.values = TRUE)$values)
  expect_lt(1 - mu, z_vals[3])

  for (solver in c("explicit", "operator")) {
    fit <- suppressMessages(lowrank_align(hd, y, simfun = simfun, ncomp = 2, mu = mu,
                                          solver = solver, sv_thresh = 0.5, lambda = 0.01))
    S <- as.matrix(fit$s)
    expect_equal(dim(S), c(120L, 2L))
    # No zero-variance component, and none in the span of the domain constants.
    expect_true(all(fit$sdev > 0.05))
    expect_lt(max(abs(domain_sums(S, sizes))), 1e-8)
    expect_equal(unname(crossprod(S)), diag(2), tolerance = 1e-8)
    expect_same_columns_up_to_sign(S, orc$vectors, tolerance = 1e-5)
  }
})

test_that("lowrank_align rejects mu outside [0, 1]", {
  skip_if_not_installed("multidesign")
  fx <- lowrank_fixture()
  hd <- fx$hd(fx$X1, fx$X2)
  expect_error(fit_lra(hd, fx$simfun, mu = -0.1), "mu")
  expect_error(fit_lra(hd, fx$simfun, mu = 1.5), "mu")
  expect_error(fit_lra(hd, fx$simfun, mu = NA_real_), "mu")
  expect_error(fit_lra(hd, fx$simfun, mu = c(0.2, 0.3)), "mu")
})
