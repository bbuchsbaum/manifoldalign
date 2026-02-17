test_that("uot_build_sparse_cost_knn_cpp returns consistent CSR and CSC", {
  n <- 3L
  m <- 4L
  k <- 2L

  nn_idx <- matrix(c(
    1L, 3L,
    2L, NA_integer_,
    4L, 1L
  ), nrow = n, byrow = TRUE)
  nn_dists <- matrix(c(
    0.1, 0.9,
    0.3, NA_real_,
    2.5, 0.2
  ), nrow = n, byrow = TRUE)

  F <- matrix(c(
    0, 0,
    1, 0,
    0, 1
  ), nrow = n, byrow = TRUE)
  G <- matrix(c(
    0, 0,
    1, 0,
    0, 1,
    1, 1
  ), nrow = m, byrow = TRUE)

  extra_i <- c(1L, 2L)
  extra_j <- c(4L, 3L)
  extra_d <- c(0.5, 2.0)

  out <- uot_build_sparse_cost_knn_cpp(
    nn_idx = nn_idx,
    nn_dists = nn_dists,
    n_cols = m,
    lambda_anat = 2.0,
    F = F,
    G = G,
    lambda_feat = 0.25,
    max_dist = 1.0,
    extra_i = extra_i,
    extra_j = extra_j,
    extra_dists = extra_d
  )

  expect_equal(out$n_rows, n)
  expect_equal(out$n_cols, m)
  expect_equal(out$nnz, out$row_ptr[n + 1L])
  expect_equal(out$nnz, out$col_ptr[m + 1L])

  csr_triplets <- function(cost) {
    rp <- cost$row_ptr
    out <- do.call(rbind, lapply(seq_len(cost$n_rows), function(i) {
      start <- rp[i] + 1L
      end <- rp[i + 1L]
      if (end < start) return(NULL)
      data.frame(
        i = rep.int(i, end - start + 1L),
        j = cost$col_idx[start:end],
        c = cost$cost[start:end]
      )
    }))
    if (is.null(out)) out <- data.frame(i = integer(), j = integer(), c = numeric())
    out
  }

  csc_triplets <- function(cost) {
    cp <- cost$col_ptr
    out <- do.call(rbind, lapply(seq_len(cost$n_cols), function(j) {
      start <- cp[j] + 1L
      end <- cp[j + 1L]
      if (end < start) return(NULL)
      data.frame(
        i = cost$row_idx[start:end],
        j = rep.int(j, end - start + 1L),
        c = cost$cost_csc[start:end]
      )
    }))
    if (is.null(out)) out <- data.frame(i = integer(), j = integer(), c = numeric())
    out
  }

  A <- csr_triplets(out)
  B <- csc_triplets(out)

  A <- A[order(A$i, A$j, A$c), , drop = FALSE]
  B <- B[order(B$i, B$j, B$c), , drop = FALSE]

  expect_equal(A$i, B$i)
  expect_equal(A$j, B$j)
  expect_equal(A$c, B$c, tolerance = 1e-12)

  # Spot-check one computed cost: edge (1,3) uses dist=0.9 and squared feature diff.
  i <- 1L
  j <- 3L
  d <- 0.9
  feat <- sum((F[i, ] - G[j, ])^2)
  expected <- 2.0 * d^2 + 0.25 * feat
  got <- A$c[A$i == i & A$j == j][1]
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("uot_build_sparse_cost_knn_cpp adds a soft prior_map penalty when provided", {
  n <- 2L
  m <- 3L
  k <- 3L

  nn_idx <- matrix(c(
    1L, 2L, 3L,
    1L, 2L, 3L
  ), nrow = n, byrow = TRUE)
  nn_dists <- matrix(0, nrow = n, ncol = k)

  # Template positions in 1D.
  Y <- matrix(c(0, 1, 3), ncol = 1)

  # Row 1 prefers column 1, row 2 prefers column 3.
  prior_map <- c(1L, 3L)

  out <- uot_build_sparse_cost_knn_cpp(
    nn_idx = nn_idx,
    nn_dists = nn_dists,
    n_cols = m,
    lambda_anat = 0,
    lambda_feat = 0,
    max_dist = Inf,
    Y_template = Y,
    prior_map = prior_map,
    lambda_prior = 2.0,
    prior_sigma = 1.0
  )

  rp <- out$row_ptr
  get_cost <- function(i, j) {
    start <- rp[i] + 1L
    end <- rp[i + 1L]
    jj <- out$col_idx[start:end]
    cc <- out$cost[start:end]
    cc[which(jj == j)][1]
  }

  # Row 1: pref = 1 (Y=0). Penalty is 2 * (Yj - 0)^2.
  expect_equal(get_cost(1L, 1L), 0, tolerance = 1e-12)
  expect_equal(get_cost(1L, 2L), 2 * (1^2), tolerance = 1e-12)
  expect_equal(get_cost(1L, 3L), 2 * (3^2), tolerance = 1e-12)

  # Row 2: pref = 3 (Y=3). Penalty is 2 * (Yj - 3)^2.
  expect_equal(get_cost(2L, 1L), 2 * (3^2), tolerance = 1e-12)
  expect_equal(get_cost(2L, 2L), 2 * (2^2), tolerance = 1e-12)
  expect_equal(get_cost(2L, 3L), 0, tolerance = 1e-12)
})

test_that("uot_build_cost dense applies prior_map penalty and biases TI-Sinkhorn", {
  n <- 5L
  m <- 5L
  X <- matrix(0, n, 1)
  Y <- matrix(seq_len(m), ncol = 1)

  prior_map <- seq_len(n)
  C <- uot_build_cost(
    X = X, Y = Y,
    lambda_anat = 0,
    lambda_feat = 0,
    neighbor_mode = "dense",
    prior_map = prior_map,
    lambda_prior = 10,
    prior_sigma = 1
  )

  expect_equal(dim(C), c(n, m))
  expect_equal(apply(C, 1, which.min), prior_map)

  alpha <- rep(1, n)
  beta <- rep(1, m)
  eps <- 0.6

  fit <- uot_ti_sinkhorn_kl(
    cost = C,
    alpha = alpha,
    beta = beta,
    epsilon = eps,
    rho1 = 1,
    rho2 = 1,
    max_iter = 2000,
    tol = 1e-12
  )

  fbar <- as.numeric(fit$fbar)
  gbar <- as.numeric(fit$gbar)
  pi <- outer(alpha, beta) * exp((outer(fbar, rep(1, m)) + outer(rep(1, n), gbar) - C) / eps)

  diag_mean <- mean(diag(pi))
  off_mean <- mean(pi[row(pi) != col(pi)])
  expect_gt(diag_mean, off_mean)
})

test_that("uot_build_cost respects group constraints via grouped neighbor search", {
  X <- matrix(c(9.9, 20.1), ncol = 1)
  Y <- matrix(c(0, 10, 20, 30), ncol = 1)
  group_X <- c("A", "B")
  group_Y <- c("A", "B", "A", "B")

  cost <- uot_build_cost(
    X = X, Y = Y,
    lambda_anat = 1,
    neighbor_mode = "knn",
    k_neighbors = 1L,
    group_X = group_X,
    group_Y = group_Y,
    ensure_cols = FALSE
  )

  rp <- cost$row_ptr
  nbr1 <- cost$col_idx[(rp[1] + 1L):rp[2]]
  nbr2 <- cost$col_idx[(rp[2] + 1L):rp[3]]
  expect_equal(nbr1, 1L) # forced to group A (Y[1]=0), not the global nearest (Y[2]=10)
  expect_equal(nbr2, 4L) # forced to group B (Y[4]=30), not the global nearest (Y[3]=20)
})

test_that("uot_build_cost radius mode enforces radius and errors on isolates", {
  skip_if_not_installed("Rnanoflann")

  X_ok <- matrix(c(0, 1), ncol = 1)
  Y <- matrix(c(0, 1.1, 2.05), ncol = 1)
  r <- 0.2

  cost_ok <- uot_build_cost(
    X = X_ok, Y = Y,
    lambda_anat = 1,
    neighbor_mode = "radius",
    radius = r,
    maxk = 8L
  )

  rp <- cost_ok$row_ptr
  edges <- do.call(rbind, lapply(seq_len(nrow(X_ok)), function(i) {
    start <- rp[i] + 1L
    end <- rp[i + 1L]
    data.frame(i = rep.int(i, end - start + 1L), j = cost_ok$col_idx[start:end])
  }))
  dists <- sqrt(rowSums((X_ok[edges$i, , drop = FALSE] - Y[edges$j, , drop = FALSE])^2))
  expect_true(all(dists <= r + 1e-12))

  X_bad <- matrix(c(0, 1, 10), ncol = 1)
  expect_error(
    uot_build_cost(
      X = X_bad, Y = Y,
      lambda_anat = 1,
      neighbor_mode = "radius",
      radius = r,
      maxk = 8L
    ),
    "zero edges"
  )
})

test_that("uot_build_cost hybrid falls back to kNN when radius too small", {
  skip_if_not_installed("Rnanoflann")

  X <- matrix(c(0, 1, 10), ncol = 1)
  Y <- matrix(c(0, 1.1, 2.05), ncol = 1)
  r <- 0.2

  cost <- uot_build_cost(
    X = X, Y = Y,
    lambda_anat = 1,
    neighbor_mode = "hybrid",
    k_neighbors = 2L,
    radius = r,
    maxk = 8L,
    min_neighbors = 1L
  )

  rp <- cost$row_ptr
  # Row 3 is far from all Y under radius r, so it should get kNN fallback edges.
  row3_start <- rp[3] + 1L
  row3_end <- rp[4]
  expect_equal(row3_end - row3_start + 1L, 2L)
  j3 <- cost$col_idx[row3_start:row3_end]
  d3 <- sqrt(rowSums(sweep(Y[j3, , drop = FALSE], 2, X[3, , drop = FALSE], "-")^2))
  expect_true(all(d3 > r))
})

test_that("TI-Sinkhorn dense solution satisfies KL KKT marginal identities", {
  set.seed(11)
  n <- 7
  m <- 6
  eps <- 0.5
  rho1 <- 1.3
  rho2 <- 0.9

  C <- matrix(runif(n * m), n, m)
  alpha <- runif(n)
  beta <- runif(m)

  fit <- uot_ti_sinkhorn_kl(
    cost = C,
    alpha = alpha,
    beta = beta,
    epsilon = eps,
    rho1 = rho1,
    rho2 = rho2,
    max_iter = 2000,
    tol = 1e-12
  )

  fbar <- as.numeric(fit$fbar)
  gbar <- as.numeric(fit$gbar)
  f <- as.numeric(fit$f)
  g <- as.numeric(fit$g)

  # Check lambda* equals the closed form and enforces mass match.
  logsumexp <- function(x) {
    mx <- max(x)
    if (!is.finite(mx)) return(mx)
    mx + log(sum(exp(x - mx)))
  }
  loga <- ifelse(alpha > 0, log(alpha), -Inf)
  logb <- ifelse(beta > 0, log(beta), -Inf)
  La <- logsumexp(loga - fbar / rho1)
  Lb <- logsumexp(logb - gbar / rho2)
  lambda_expected <- (rho1 * rho2) / (rho1 + rho2) * (La - Lb)
  expect_equal(as.numeric(fit$lambda), lambda_expected, tolerance = 1e-10)

  m1 <- sum(alpha * exp(-f / rho1))
  m2 <- sum(beta * exp(-g / rho2))
  expect_equal(m1, m2, tolerance = 1e-10)

  # KKT: pi = alpha ⊗ beta ⊙ exp((fbar ⊕ gbar - C)/eps)
  pi <- outer(alpha, beta) * exp((outer(fbar, rep(1, m)) + outer(rep(1, n), gbar) - C) / eps)
  pi1 <- rowSums(pi)
  pi2 <- colSums(pi)
  expect_equal(pi1, alpha * exp(-f / rho1), tolerance = 1e-8)
  expect_equal(pi2, beta * exp(-g / rho2), tolerance = 1e-8)
})

test_that("TI-Sinkhorn sparse CSR+CSC backend matches dense identities", {
  set.seed(12)
  n <- 6
  m <- 5
  eps <- 0.7
  rho1 <- 1.1
  rho2 <- 0.8

  C <- matrix(runif(n * m), n, m)
  alpha <- runif(n)
  beta <- runif(m)

  # Full bipartite graph with both CSR and CSC representations.
  row_ptr <- as.integer(seq(0, n * m, by = m))
  col_idx <- as.integer(rep(seq_len(m), times = n))
  cost_csr <- as.numeric(t(C))

  col_ptr <- as.integer(seq(0, n * m, by = n))
  row_idx <- as.integer(rep(seq_len(n), times = m))
  cost_csc <- as.numeric(C)

  cost_sparse <- list(
    row_ptr = row_ptr,
    col_idx = col_idx,
    cost = cost_csr,
    col_ptr = col_ptr,
    row_idx = row_idx,
    cost_csc = cost_csc,
    n_rows = n,
    n_cols = m
  )

  fit <- uot_ti_sinkhorn_kl(
    cost = cost_sparse,
    alpha = alpha,
    beta = beta,
    epsilon = eps,
    rho1 = rho1,
    rho2 = rho2,
    max_iter = 2000,
    tol = 1e-12
  )

  expect_equal(fit$backend, "sparse_csr_csc")

  fbar <- as.numeric(fit$fbar)
  gbar <- as.numeric(fit$gbar)
  f <- as.numeric(fit$f)
  g <- as.numeric(fit$g)

  pi <- outer(alpha, beta) * exp((outer(fbar, rep(1, m)) + outer(rep(1, n), gbar) - C) / eps)
  pi1 <- rowSums(pi)
  pi2 <- colSums(pi)
  expect_equal(pi1, alpha * exp(-f / rho1), tolerance = 1e-8)
  expect_equal(pi2, beta * exp(-g / rho2), tolerance = 1e-8)
})

test_that("uot_apply_map matches explicit coupling (dense and sparse vector)", {
  set.seed(13)
  n <- 5
  m <- 4
  T <- 3
  eps <- 0.6
  rho1 <- 1.2
  rho2 <- 0.7

  C <- matrix(runif(n * m), n, m)
  alpha <- runif(n)
  beta <- runif(m)

  fit <- uot_ti_sinkhorn_kl(
    cost = C,
    alpha = alpha,
    beta = beta,
    epsilon = eps,
    rho1 = rho1,
    rho2 = rho2,
    max_iter = 2000,
    tol = 1e-12
  )

  fbar <- as.numeric(fit$fbar)
  gbar <- as.numeric(fit$gbar)
  pi <- outer(alpha, beta) * exp((outer(fbar, rep(1, m)) + outer(rep(1, n), gbar) - C) / eps)
  pi2 <- colSums(pi)

  s_vec <- rnorm(n)
  s_mat <- matrix(rnorm(T * n), T, n)
  delta <- 1e-8

  y_exp_vec <- as.numeric(crossprod(pi, s_vec) / (pi2 + delta))
  y_hat_vec <- uot_apply_map(C, alpha, beta, fit$fbar, fit$gbar, epsilon = eps, signal = s_vec, delta = delta)
  expect_equal(y_hat_vec, y_exp_vec, tolerance = 1e-8)

  y_exp_mat <- s_mat %*% sweep(pi, 2, pi2 + delta, "/")
  y_hat_mat <- uot_apply_map(C, alpha, beta, fit$fbar, fit$gbar, epsilon = eps, signal = s_mat, delta = delta)
  expect_equal(y_hat_mat, y_exp_mat, tolerance = 1e-8)

  # Same potentials applied over an equivalent full-bipartite sparse cost list.
  row_ptr <- as.integer(seq(0, n * m, by = m))
  col_idx <- as.integer(rep(seq_len(m), times = n))
  cost_csr <- as.numeric(t(C))
  C_sparse <- list(
    row_ptr = row_ptr,
    col_idx = col_idx,
    cost = cost_csr,
    n_rows = n,
    n_cols = m
  )
  y_hat_sparse <- uot_apply_map(C_sparse, alpha, beta, fit$fbar, fit$gbar, epsilon = eps, signal = s_vec, delta = delta)
  expect_equal(y_hat_sparse, y_exp_vec, tolerance = 1e-8)

  # Sparse matrix mapping (requires CSC fields).
  col_ptr <- as.integer(seq(0, n * m, by = n))
  row_idx <- as.integer(rep(seq_len(n), times = m))
  cost_csc <- as.numeric(C)
  C_sparse_csc <- list(
    row_ptr = row_ptr,
    col_idx = col_idx,
    cost = cost_csr,
    col_ptr = col_ptr,
    row_idx = row_idx,
    cost_csc = cost_csc,
    n_rows = n,
    n_cols = m
  )
  y_hat_sparse_mat <- uot_apply_map(C_sparse_csc, alpha, beta, fit$fbar, fit$gbar,
                                    epsilon = eps, signal = s_mat, delta = delta)
  expect_equal(y_hat_sparse_mat, y_exp_mat, tolerance = 1e-8)
})

test_that("multiset_uot_align wires constraint_fields into sparse neighborhoods", {
  set.seed(14)

  X <- matrix(c(9.9, 20.1), ncol = 1)
  Y <- matrix(c(0, 10, 20, 30), ncol = 1)
  hemi_X <- c("L", "R")
  hemi_Y <- c("L", "R", "L", "R")

  datasets <- list(list(X = X, alpha = rep(1, nrow(X)), hemi = hemi_X))
  template <- list(Y = Y, beta = rep(1, nrow(Y)), hemi = hemi_Y)

  fit <- multiset_uot_align(
    datasets = datasets,
    template = template,
    epsilon = 0.5,
    rho1 = 1,
    rho2 = 1,
    lambda_anat = 1,
    neighbor_mode = "knn",
    k_neighbors = 1L,
    constraint_fields = "hemi",
    max_outer = 1,
    max_inner = 50,
    tol_inner = 1e-10,
    tol_outer = 1e-10,
    parallel = FALSE
  )

  cost <- fit$fits[[1]]$cost
  rp <- cost$row_ptr
  edges <- do.call(rbind, lapply(seq_len(nrow(X)), function(i) {
    start <- rp[i] + 1L
    end <- rp[i + 1L]
    data.frame(i = rep.int(i, end - start + 1L), j = cost$col_idx[start:end])
  }))

  expect_true(all(hemi_X[edges$i] == hemi_Y[edges$j]))
})

test_that("uot_extract_coupling produces a reusable map operator", {
  set.seed(15)
  n <- 6
  m <- 5
  T <- 4
  eps <- 0.55
  rho1 <- 1.0
  rho2 <- 0.9
  delta <- 1e-8

  C <- matrix(runif(n * m), n, m)
  alpha <- runif(n)
  beta <- runif(m)
  fit <- uot_ti_sinkhorn_kl(C, alpha, beta, epsilon = eps, rho1 = rho1, rho2 = rho2,
                            max_iter = 2000, tol = 1e-12)

  row_ptr <- as.integer(seq(0, n * m, by = m))
  col_idx <- as.integer(rep(seq_len(m), times = n))
  cost_csr <- as.numeric(t(C))
  col_ptr <- as.integer(seq(0, n * m, by = n))
  row_idx <- as.integer(rep(seq_len(n), times = m))
  cost_csc <- as.numeric(C)
  cost_sparse <- list(
    row_ptr = row_ptr,
    col_idx = col_idx,
    cost = cost_csr,
    col_ptr = col_ptr,
    row_idx = row_idx,
    cost_csc = cost_csc,
    n_rows = n,
    n_cols = m
  )

  op <- uot_extract_coupling(cost_sparse, alpha, beta, fit$fbar, fit$gbar,
                             epsilon = eps, weights = "map", delta = delta)
  Q <- op$coupling
  expect_s4_class(Q, "dgCMatrix")
  expect_equal(dim(Q), c(n, m))

  s_vec <- rnorm(n)
  s_mat <- matrix(rnorm(T * n), T, n)

  y1 <- uot_apply_map(C, alpha, beta, fit$fbar, fit$gbar, epsilon = eps, signal = s_vec, delta = delta)
  y2 <- as.numeric(crossprod(Q, s_vec))
  expect_equal(y2, y1, tolerance = 1e-8)

  Y1 <- uot_apply_map(C, alpha, beta, fit$fbar, fit$gbar, epsilon = eps, signal = s_mat, delta = delta)
  Y2 <- as.matrix(s_mat %*% Q)
  expect_equal(Y2, Y1, tolerance = 1e-8)

  op1 <- uot_extract_coupling(cost_sparse, alpha, beta, fit$fbar, fit$gbar,
                              epsilon = eps, weights = "map", delta = delta,
                              prune_topk_col = 1L)
  Q1 <- op1$coupling
  expect_true(length(Q1@x) <= m)
})

test_that("uot_fit_pair + uot_map provide a simple pairwise workflow", {
  set.seed(16)
  n <- 7
  m <- 6
  eps <- 0.5
  rho1 <- 1.2
  rho2 <- 0.8
  X <- matrix(rnorm(n * 3), n, 3)
  Y <- matrix(rnorm(m * 3), m, 3)
  alpha <- runif(n)
  beta <- runif(m)

  fit <- uot_fit_pair(
    X = X, Y = Y,
    alpha = alpha, beta = beta,
    lambda_anat = 1,
    neighbor_mode = "dense",
    epsilon = eps, rho1 = rho1, rho2 = rho2,
    max_iter = 500, tol = 1e-10
  )

  expect_true(inherits(fit, "uot_pair_fit"))

  s <- rnorm(n)
  y1 <- uot_map(fit, s)
  y2 <- uot_apply_map(fit$cost, fit$alpha, fit$beta, fit$sol$fbar, fit$sol$gbar,
                      epsilon = fit$epsilon, signal = s)
  expect_equal(y1, y2, tolerance = 1e-10)

  Q <- uot_operator(fit, weights = "map")
  y3 <- as.numeric(crossprod(Q, s))
  expect_equal(y3, y2, tolerance = 1e-8)
})

test_that("multiset_uot_map maps signals using stored params", {
  set.seed(17)
  K <- 2
  n <- 10
  m <- 12

  datasets <- lapply(seq_len(K), function(k) {
    X <- matrix(rnorm(n * 3), n, 3)
    list(X = X, alpha = rep(1, n))
  })
  template <- list(Y = matrix(rnorm(m * 3), m, 3), beta = rep(1, m))

  eps <- 0.6
  fit <- multiset_uot_align(
    datasets = datasets,
    template = template,
    epsilon = eps, rho1 = 1, rho2 = 1,
    lambda_anat = 1, lambda_feat = 0,
    neighbor_mode = "knn",
    k_neighbors = 4L,
    max_outer = 1,
    max_inner = 200,
    tol_inner = 1e-10,
    tol_outer = 1e-10,
    parallel = FALSE
  )

  signals <- lapply(seq_len(K), function(k) rnorm(n))
  mapped <- multiset_uot_map(fit, datasets = datasets, signals = signals, delta = 1e-8)
  expect_equal(length(mapped), K)
  expect_equal(length(mapped[[1]]), m)

  manual <- lapply(seq_len(K), function(k) {
    uot_apply_map(
      cost = fit$fits[[k]]$cost,
      alpha = datasets[[k]]$alpha,
      beta = fit$template$beta,
      fbar = fit$fits[[k]]$sol$fbar,
      gbar = fit$fits[[k]]$sol$gbar,
      epsilon = eps,
      signal = signals[[k]]
    )
  })

  expect_equal(mapped, manual, tolerance = 1e-10)
})
