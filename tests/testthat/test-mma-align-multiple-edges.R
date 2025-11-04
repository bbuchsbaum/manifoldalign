test_that("mma_align_multiple handles disconnected graphs (auto-K)", {
  skip_if_missing_cd_deps()
  set.seed(42)

  # Two separated clusters per domain to induce disconnected kNN graphs with small k
  mk_domain <- function(n1 = 20, n2 = 20, d = 3) {
    X1 <- matrix(rnorm(n1 * d, mean = -3, sd = 0.3), n1, d)
    X2 <- matrix(rnorm(n2 * d, mean = +3, sd = 0.3), n2, d)
    rbind(X1, X2)
  }

  X1 <- mk_domain()
  X2 <- mk_domain()
  X3 <- mk_domain()

  res <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = NULL,            # auto-K selection
                            target_var = 0.9,
                            sigma = 0.5,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            em_max_iter = 10)

  expect_true(inherits(res, "mma_align_multiple"))
  expect_true(is.numeric(res$K) && res$K >= 2)
  expect_equal(ncol(res$s), res$K)
})

test_that("mma_align_multiple remains numerically stable with high outlier prior", {
  skip_if_missing_cd_deps()
  set.seed(100)
  n <- 30; d <- 4
  X <- matrix(rnorm(n * d), n, d)
  # Make domain 3 noisy/outlier-heavy by adding large noise
  X1 <- X
  X2 <- X + matrix(rnorm(n * d, sd = 0.05), n, d)
  X3 <- matrix(rnorm(n * d, sd = 2.0), n, d)

  res <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = 4,
                            sigma = 0.7,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            em_max_iter = 15,
                            em_tol = 1e-5,
                            em_outlier = 0.3)

  expect_true(inherits(res, "mma_align_multiple"))
  # EM stats present; sigma2 finite
  expect_true(is.list(res$em_stats))
  s2 <- vapply(res$em_stats, function(z) if (is.list(z)) z$sigma2 else NA_real_, numeric(1))
  expect_true(all(is.finite(s2[!is.na(s2)])))
})

