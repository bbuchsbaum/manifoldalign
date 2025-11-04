test_that("alignment_quality works for reference mode with assignments", {
  skip_if_missing_cd_deps()
  set.seed(12)
  n <- 30; d <- 5
  t <- runif(n, 0, 2*pi)
  Xb <- cbind(cos(t), sin(t), matrix(rnorm(n*(d-2), sd = 0.03), n, d-2))
  orth <- function(p) { M <- matrix(rnorm(p*p), p, p); sv <- svd(M); sv$u %*% t(sv$v) }
  X1 <- Xb
  X2 <- Xb %*% orth(d)
  X3 <- Xb %*% orth(d)

  fit <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = 5,
                            sigma = 0.6,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            signature_method = "w2",
                            final_assignment = "hungarian",
                            em_max_iter = 15)
  aq <- alignment_quality(fit)
  expect_equal(aq$mode, "reference")
  expect_true(is.data.frame(aq$per_domain))
  # For easy synthetic setup, errors should be very small
  pd <- aq$per_domain
  expect_true(all(pd$mse[-1] < 1e-2, na.rm = TRUE))
  expect_true(all(pd$mse_soft[-1] < 1e-2, na.rm = TRUE))
  expect_true(is.list(aq$global))
  expect_true(aq$global$mse_weighted < 1e-2)
})

test_that("alignment_quality works for consensus mode", {
  skip_if_missing_cd_deps()
  set.seed(13)
  n1 <- 40; n2 <- 35; n3 <- 30; d <- 4
  X <- matrix(rnorm(1000), 1000, d)
  X1 <- X[1:n1,]
  X2 <- X[(n1+1):(n1+n2),]
  X3 <- X[(n1+n2+1):(n1+n2+n3),]

  fit <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = 5,
                            sigma = 0.7,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            match_to = "consensus",
                            consensus_centers = "min",
                            final_assignment = "hungarian",
                            em_max_iter = 10)
  aq <- alignment_quality(fit)
  expect_equal(aq$mode, "consensus")
  expect_true(is.data.frame(aq$per_domain))
  pd <- aq$per_domain
  expect_equal(nrow(pd), 3L)
  expect_true(all(pd$n == c(n1, n2, n3)))
  # Soft errors should be defined for at least one domain and be non-negative when present
  expect_true(any(is.finite(pd$mse_soft)))
  expect_true(all(pd$mse_soft[is.finite(pd$mse_soft)] >= 0))
  expect_true(is.list(aq$global))
  expect_true(is.finite(aq$global$mse_soft_weighted))
})
