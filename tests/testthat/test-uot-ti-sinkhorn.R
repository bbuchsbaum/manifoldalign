test_that("TI-Sinkhorn dense and sparse agree on pi2 up to tolerance", {
  set.seed(1)

  n <- 6
  m <- 5
  eps <- 0.7
  rho1 <- 1.1
  rho2 <- 0.9

  C <- matrix(runif(n * m), n, m)
  alpha <- runif(n)
  beta <- runif(m)

  fit_dense <- uot_ti_sinkhorn_kl(
    C, alpha, beta,
    epsilon = eps, rho1 = rho1, rho2 = rho2,
    max_iter = 500, tol = 1e-10
  )

  # Full bipartite graph as sparse CSR (row-major edge order)
  row_ptr <- as.integer(seq(0, n * m, by = m))
  col_idx <- as.integer(rep(seq_len(m), times = n))
  cost_vec <- as.numeric(t(C))
  C_sparse <- list(
    row_ptr = row_ptr,
    col_idx = col_idx,
    cost = cost_vec,
    n_rows = n,
    n_cols = m
  )

  fit_sparse <- uot_ti_sinkhorn_kl(
    C_sparse, alpha, beta,
    epsilon = eps, rho1 = rho1, rho2 = rho2,
    max_iter = 500, tol = 1e-10
  )

  logsumexp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(x - m)))
  }

  pi2_from_potentials <- function(fbar, gbar) {
    loga <- ifelse(alpha > 0, log(alpha), -Inf)
    logb <- ifelse(beta > 0, log(beta), -Inf)
    log_pi2 <- vapply(seq_len(m), function(j) {
      z <- loga + logb[j] + (fbar + gbar[j] - C[, j]) / eps
      logsumexp(z)
    }, numeric(1))
    # Compare as a probability vector to remove any global scaling drift
    p <- exp(log_pi2 - logsumexp(log_pi2))
    p[!is.finite(p)] <- 0
    p / sum(p)
  }

  pi2_dense <- pi2_from_potentials(fit_dense$fbar, fit_dense$gbar)
  pi2_sparse <- pi2_from_potentials(fit_sparse$fbar, fit_sparse$gbar)

  expect_equal(pi2_dense, pi2_sparse, tolerance = 1e-6)
})

test_that("multiset_uot_align runs and returns finite template updates", {
  set.seed(2)

  K <- 3
  datasets <- lapply(seq_len(K), function(k) {
    X <- matrix(rnorm(30), 10, 3)
    F <- matrix(rnorm(40), 10, 4)
    list(X = X, alpha = rep(1, 10), F = F)
  })
  template <- list(
    Y = matrix(rnorm(36), 12, 3),
    beta = rep(1, 12),
    G = matrix(0, 12, 4)
  )

  fit <- multiset_uot_align(
    datasets, template,
    epsilon = 0.5, rho1 = 1, rho2 = 1,
    lambda_anat = 1, lambda_feat = 0.1,
    k_neighbors = 6,
    max_outer = 2,
    max_inner = 200,
    tol_inner = 1e-8,
    tol_outer = 1e-6,
    learn_template_features = TRUE,
    parallel = FALSE
  )

  expect_equal(length(fit$template$beta), 12)
  expect_true(all(is.finite(fit$template$beta)))
  expect_equal(sum(fit$template$beta), 10, tolerance = 1e-6)

  expect_equal(dim(fit$template$G), c(12, 4))
  expect_true(all(is.finite(fit$template$G)))

  expect_equal(length(fit$fits), K)
  expect_true(all(vapply(fit$fits, function(x) is.list(x$sol), logical(1))))
})

