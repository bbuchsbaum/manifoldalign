test_that("mma_align_multiple achieves high correspondence on rotated domains (reference mode)", {
  skip_if_missing_cd_deps()
  set.seed(123)

  n <- 40; d <- 5
  # Base data: points on a noisy circle embedded in d-dim
  t <- runif(n, 0, 2*pi)
  Xbase2 <- cbind(cos(t), sin(t))
  Xbase <- cbind(Xbase2, matrix(rnorm(n * (d - 2), sd = 0.05), n, d - 2))

  # Random orthogonal transforms for domains 2 and 3
  rand_orth <- function(p) {
    M <- matrix(rnorm(p*p), p, p)
    sv <- svd(M)
    sv$u %*% t(sv$v)
  }
  Q2 <- rand_orth(d)
  Q3 <- rand_orth(d)

  X1 <- Xbase
  X2 <- Xbase %*% Q2
  X3 <- Xbase %*% Q3

  res <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = 5,
                            sigma = 0.7,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            em_max_iter = 25,
                            em_tol = 1e-6,
                            final_assignment = "hungarian")

  expect_true(inherits(res, "mma_align_multiple"))

  # EM diagnostics present; iterations should be positive for non-ref domains
  expect_true(is.list(res$em_stats))
  iters <- vapply(res$em_stats[-1], function(z) as.integer(z$iterations), integer(1))
  expect_true(all(iters > 0))

  # Assignment accuracy against identity should be high
  assign <- res$assignment
  expect_equal(length(assign), 3L)
  acc2 <- mean(assign[[2]] == seq_len(n))
  acc3 <- mean(assign[[3]] == seq_len(n))
  expect_gte(acc2, 0.9)
  expect_gte(acc3, 0.9)

  # Reconstruct aligned embeddings per domain and check MSE after applying assignment
  ns <- c(nrow(res$posteriors[[1]]), ncol(res$posteriors[[2]]), ncol(res$posteriors[[3]]))
  idx2 <- (ns[1] + 1):(ns[1] + ns[2])
  idx3 <- (ns[1] + ns[2] + 1):(ns[1] + ns[2] + ns[3])
  E1 <- res$s[seq_len(ns[1]), , drop = FALSE]
  E2 <- res$s[idx2, , drop = FALSE][assign[[2]], , drop = FALSE]
  E3 <- res$s[idx3, , drop = FALSE][assign[[3]], , drop = FALSE]

  mse12 <- mean((E1 - E2)^2)
  mse13 <- mean((E1 - E3)^2)
  expect_lte(mse12, 1e-2)
  expect_lte(mse13, 1e-2)
})
