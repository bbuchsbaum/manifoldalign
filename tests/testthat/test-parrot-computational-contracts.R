library(testthat)
library(manifoldalign)

inverse_perm_local <- function(perm) {
  inv <- integer(length(perm))
  inv[perm] <- seq_along(perm)
  inv
}

make_parrot_hyperdesign <- function(X1, X2, map12, anchor_idx) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  stopifnot(length(map12) == n1)

  anchors1 <- rep(NA_integer_, n1)
  anchors2 <- rep(NA_integer_, n2)

  for (k in seq_along(anchor_idx)) {
    i <- anchor_idx[k]
    j <- map12[i]
    anchors1[i] <- k
    anchors2[j] <- k
  }

  hd <- list(
    domain1 = list(
      x = X1,
      design = data.frame(node_id = seq_len(n1), anchors = anchors1)
    ),
    domain2 = list(
      x = X2,
      design = data.frame(node_id = seq_len(n2), anchors = anchors2)
    )
  )
  class(hd) <- "hyperdesign"
  hd
}

extract_parrot_state <- function(hd) {
  strata <- list(
    list(x = hd$domain1$x, design = hd$domain1$design),
    list(x = hd$domain2$x, design = hd$domain2$design)
  )

  anchor_vec1 <- hd$domain1$design$anchors
  anchor_vec2 <- hd$domain2$design$anchors

  anchor_info <- list(
    vec1 = anchor_vec1,
    vec2 = anchor_vec2,
    idx1 = which(!is.na(anchor_vec1)),
    idx2 = which(!is.na(anchor_vec2)),
    n1 = nrow(hd$domain1$x),
    n2 = nrow(hd$domain2$x)
  )

  list(
    networks = manifoldalign:::extract_parrot_networks(strata),
    anchor_info = anchor_info
  )
}

test_that("PARROT Sylvester solver matches analytic oracle on a small system", {
  set.seed(20260219)

  n1 <- 4
  n2 <- 3
  W1_raw <- matrix(runif(n1 * n1), n1, n1)
  W2_raw <- matrix(runif(n2 * n2), n2, n2)
  W1 <- W1_raw / rowSums(W1_raw)
  W2 <- W2_raw / rowSums(W2_raw)
  Cnode <- matrix(runif(n1 * n2, min = 0.1, max = 1.7), n1, n2)

  beta <- 0.2
  gamma <- 0.15
  kappa <- (1 - beta) * gamma

  X_iter <- manifoldalign:::solve_sylvester_rwr(
    W1 = W1,
    W2T = t(W2),
    Cnode = Cnode,
    beta = beta,
    gamma = gamma,
    tol = 1e-12,
    max_iter = 2000,
    use_cpp = FALSE
  )

  A <- diag(n1 * n2) - kappa * kronecker(W2, W1)
  b <- as.vector((1 + beta) * Cnode)
  X_exact <- matrix(solve(A, b), nrow = n1, ncol = n2)

  expect_equal(X_iter, X_exact, tolerance = 1e-7)
  expect_true(all(is.finite(X_iter)))
})

test_that("PARROT transport is permutation equivariant for domain-2 relabeling", {
  set.seed(444)

  n <- 24
  theta <- 2 * pi * seq_len(n) / n
  X1 <- cbind(
    cos(theta),
    sin(theta),
    sin(2 * theta),
    seq_len(n) / n
  )
  X1 <- X1 + matrix(rnorm(n * ncol(X1), sd = 0.015), n, ncol(X1))

  perm21 <- sample(n)
  map12 <- match(seq_len(n), perm21)
  X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * ncol(X1), sd = 0.01), n, ncol(X1))

  anchor_idx <- seq(1, n, by = 4)
  hd_base <- make_parrot_hyperdesign(X1, X2, map12, anchor_idx)

  fit_base <- parrot(
    hd_base,
    anchors = anchors,
    sigma = 0.2,
    lambda = 0.1,
    tau = 0.01,
    alpha = 0.25,
    max_iter = 60,
    use_cpp = FALSE
  )
  S_base <- as.matrix(fit_base$transport_plan)

  q <- sample(n)
  q_inv <- inverse_perm_local(q)
  X2_perm <- X2[q, , drop = FALSE]
  map12_perm <- q_inv[map12]
  hd_perm <- make_parrot_hyperdesign(X1, X2_perm, map12_perm, anchor_idx)

  fit_perm <- parrot(
    hd_perm,
    anchors = anchors,
    sigma = 0.2,
    lambda = 0.1,
    tau = 0.01,
    alpha = 0.25,
    max_iter = 60,
    use_cpp = FALSE
  )
  S_perm <- as.matrix(fit_perm$transport_plan)

  expect_equal(S_perm, S_base[, q, drop = FALSE], tolerance = 2e-2)
})

test_that("PARROT internal cost matrix remains finite and strictly positive", {
  set.seed(2024)

  n <- 12
  d <- 3
  X1 <- matrix(rnorm(n * d), n, d)
  perm21 <- sample(n)
  map12 <- match(seq_len(n), perm21)
  X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * d, sd = 0.05), n, d)

  hd <- make_parrot_hyperdesign(X1, X2, map12, anchor_idx = c(1, 4, 7, 10))
  state <- extract_parrot_state(hd)

  rwr <- manifoldalign:::compute_parrot_rwr(
    state$networks,
    state$anchor_info,
    sigma = 0.2,
    max_iter = 80,
    tol = 1e-8,
    use_cpp = FALSE
  )

  C <- manifoldalign:::compute_parrot_cost(
    state$networks,
    rwr,
    state$anchor_info,
    alpha = 0.3,
    sigma = 0.2,
    gamma = 0.1,
    use_cpp = FALSE
  )

  expect_equal(dim(C), c(n, n))
  expect_true(all(is.finite(C)))
  expect_true(all(C > 0))
})

test_that("PARROT preserves transport invariants across randomized stress fixtures", {
  seeds <- c(11, 17, 29)

  for (seed in seeds) {
    set.seed(seed)

    n <- 12
    d <- 3

    base <- matrix(rnorm(n * d), n, d)
    scales <- 10^runif(d, min = -5, max = 5)
    X1 <- sweep(base, 2, scales, "*")
    X1[1, ] <- 0
    X1[2, ] <- X1[3, ] + 1e-10

    perm21 <- sample(n)
    map12 <- match(seq_len(n), perm21)
    X2 <- X1[perm21, , drop = FALSE]
    X2 <- X2 + matrix(rnorm(n * d, sd = 1e-4), n, d)
    X2 <- sweep(X2, 2, c(1e5, -5e4, 2e4), "+")

    hd <- make_parrot_hyperdesign(X1, X2, map12, anchor_idx = c(1, 4, 7, 10))
    fit <- parrot(
      hd,
      anchors = anchors,
      sigma = 0.2,
      lambda = 0.05,
      tau = 0.02,
      alpha = 0.25,
      max_iter = 40,
      use_cpp = FALSE
    )

    S <- as.matrix(fit$transport_plan)
    target <- rep(1 / n, n)

    expect_true(all(is.finite(S)))
    expect_true(all(S >= 0))
    expect_equal(rowSums(S), target, tolerance = 8e-3)
    expect_equal(colSums(S), target, tolerance = 8e-3)
  }
})
