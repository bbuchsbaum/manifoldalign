test_that("TI-Sinkhorn remains stable across a stiff parameter sweep (dense)", {
  set.seed(101)

  logsumexp <- function(x) {
    mx <- max(x)
    if (!is.finite(mx)) return(mx)
    mx + log(sum(exp(x - mx)))
  }

  check_kkt_logs <- function(C, alpha, beta, fit, eps, rho1, rho2, tol_log = 5e-5) {
    fbar <- as.numeric(fit$fbar)
    gbar <- as.numeric(fit$gbar)
    f <- as.numeric(fit$f)
    g <- as.numeric(fit$g)

    loga <- ifelse(alpha > 0, log(alpha), -Inf)
    logb <- ifelse(beta > 0, log(beta), -Inf)

    n <- nrow(C)
    m <- ncol(C)

    log_pi1 <- vapply(seq_len(n), function(i) {
      z <- logb + (fbar[i] + gbar - C[i, ]) / eps
      loga[i] + logsumexp(z)
    }, numeric(1))

    log_pi2 <- vapply(seq_len(m), function(j) {
      z <- loga + (fbar + gbar[j] - C[, j]) / eps
      logb[j] + logsumexp(z)
    }, numeric(1))

    rhs1 <- loga - f / rho1
    rhs2 <- logb - g / rho2

    idx1 <- which(is.finite(log_pi1) & is.finite(rhs1))
    idx2 <- which(is.finite(log_pi2) & is.finite(rhs2))
    expect_true(length(idx1) > 0)
    expect_true(length(idx2) > 0)
    expect_lt(max(abs(log_pi1[idx1] - rhs1[idx1])), tol_log)
    expect_lt(max(abs(log_pi2[idx2] - rhs2[idx2])), tol_log)

    m1 <- sum(alpha * exp(-f / rho1))
    m2 <- sum(beta * exp(-g / rho2))
    expect_true(is.finite(m1) && is.finite(m2))
    expect_equal(m1, m2, tolerance = 1e-8)
  }

  combos <- list(
    list(eps = 0.02, rho1 = 0.1, rho2 = 0.1),
    list(eps = 0.02, rho1 = 0.1, rho2 = 10),
    list(eps = 0.02, rho1 = 10, rho2 = 0.1),
    list(eps = 0.05, rho1 = 1, rho2 = 1),
    list(eps = 0.2, rho1 = 0.1, rho2 = 1),
    list(eps = 0.2, rho1 = 1, rho2 = 10),
    list(eps = 1.0, rho1 = 10, rho2 = 1)
  )

  for (trial in seq_len(3)) {
    n <- 9
    m <- 8
    C <- matrix(rexp(n * m, rate = 0.5), n, m) * 10

    alpha <- runif(n)
    beta <- runif(m)
    alpha[sample.int(n, size = 2)] <- 0
    beta[sample.int(m, size = 2)] <- 0
    alpha <- alpha / sum(alpha)
    beta <- beta / sum(beta)

    for (p in combos) {
      fit <- uot_ti_sinkhorn_kl(
        cost = C,
        alpha = alpha,
        beta = beta,
        epsilon = p$eps,
        rho1 = p$rho1,
        rho2 = p$rho2,
        max_iter = 5000,
        tol = 1e-12
      )
      expect_true(all(is.finite(as.numeric(fit$lambda))))
      expect_true(all(is.finite(as.numeric(fit$fbar))))
      expect_true(all(is.finite(as.numeric(fit$gbar))))
      expect_true(is.finite(fit$residual))

      check_kkt_logs(C, alpha, beta, fit, eps = p$eps, rho1 = p$rho1, rho2 = p$rho2)
    }
  }
})

test_that("TI-Sinkhorn remains stable across a stiff parameter sweep (sparse CSR+CSC)", {
  set.seed(102)

  logsumexp <- function(x) {
    mx <- max(x)
    if (!is.finite(mx)) return(mx)
    mx + log(sum(exp(x - mx)))
  }

  combos <- list(
    list(eps = 0.02, rho1 = 0.2, rho2 = 5),
    list(eps = 0.05, rho1 = 1, rho2 = 1),
    list(eps = 0.2, rho1 = 10, rho2 = 0.5)
  )

  for (trial in seq_len(2)) {
    n <- 8
    m <- 7
    C <- matrix(rexp(n * m, rate = 0.8), n, m) * 8

    alpha <- runif(n)
    beta <- runif(m)
    alpha[sample.int(n, size = 1)] <- 0
    beta[sample.int(m, size = 1)] <- 0
    alpha <- alpha / sum(alpha)
    beta <- beta / sum(beta)

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

    for (p in combos) {
      fit <- uot_ti_sinkhorn_kl(
        cost = cost_sparse,
        alpha = alpha,
        beta = beta,
        epsilon = p$eps,
        rho1 = p$rho1,
        rho2 = p$rho2,
        max_iter = 5000,
        tol = 1e-12
      )

      fbar <- as.numeric(fit$fbar)
      gbar <- as.numeric(fit$gbar)
      f <- as.numeric(fit$f)
      g <- as.numeric(fit$g)

      loga <- ifelse(alpha > 0, log(alpha), -Inf)
      logb <- ifelse(beta > 0, log(beta), -Inf)

      log_pi1 <- vapply(seq_len(n), function(i) {
        z <- logb + (fbar[i] + gbar - C[i, ]) / p$eps
        loga[i] + logsumexp(z)
      }, numeric(1))

      log_pi2 <- vapply(seq_len(m), function(j) {
        z <- loga + (fbar + gbar[j] - C[, j]) / p$eps
        logb[j] + logsumexp(z)
      }, numeric(1))

      rhs1 <- loga - f / p$rho1
      rhs2 <- logb - g / p$rho2

      idx1 <- which(is.finite(log_pi1) & is.finite(rhs1))
      idx2 <- which(is.finite(log_pi2) & is.finite(rhs2))
      expect_true(length(idx1) > 0)
      expect_true(length(idx2) > 0)
      expect_lt(max(abs(log_pi1[idx1] - rhs1[idx1])), 5e-5)
      expect_lt(max(abs(log_pi2[idx2] - rhs2[idx2])), 5e-5)
    }
  }
})

