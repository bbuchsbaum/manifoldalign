test_that("cone_align_multiple handles disconnected graphs", {

  set.seed(123)
  n  <- 60
  m  <- 3

  # Two equal-size disconnected blobs per graph ------------------------
  make_blob <- function() rbind(
      matrix(rnorm(n/2 * 3,  3), n/2, 3),   # cluster A
      matrix(rnorm(n/2 * 3, -3), n/2, 3))   # cluster B

  mats <- replicate(m, make_blob(), simplify = FALSE)

  # No ground‑truth permutation here; we only require convergence
  # Suppress the valid warning about ncomp > informative eigenvectors
  res <- suppressWarnings(
    cone_align_multiple(
      mats, ncomp = 6, max_iter = 30, tol = 1e-2
    )
  )

  # Check that the result is a valid object despite the warning
  expect_s3_class(res, "cone_align_multiple")
  expect_type(res$assignment, "list")
  expect_length(res$assignment, m)
  expect_true(all(vapply(res$assignment, is.integer, logical(1))))
}) 