test_that(".with_local_seed fixes the seed and restores the caller's stream", {
  with_local_seed <- get(".with_local_seed", envir = asNamespace("manifoldalign"))
  set.seed(42)
  expected_next <- runif(1)
  set.seed(42)
  a <- with_local_seed(7, runif(3))
  b <- with_local_seed(7, runif(3))
  expect_identical(a, b)
  expect_identical(runif(1), expected_next)
  set.seed(42)
  expect_error(with_local_seed(1, stop("boom")), "boom")
  expect_identical(runif(1), expected_next)
})

test_that(".irlba_det is reproducible and leaves the caller's RNG alone", {
  skip_if_not_installed("irlba")
  set.seed(3)
  A <- matrix(rnorm(60 * 20), 60, 20)
  set.seed(100)
  before <- .Random.seed
  s1 <- manifoldalign:::.irlba_det(A, nv = 3, nu = 0)
  expect_identical(.Random.seed, before)
  set.seed(999)
  s2 <- manifoldalign:::.irlba_det(A, nv = 3, nu = 0)
  expect_identical(s1$v, s2$v)
  expect_equal(s1$d, svd(A)$d[1:3], tolerance = 1e-6)
})

test_that("SSMA and multiscale fits do not depend on the caller's RNG state", {
  skip_if_not_installed("irlba")
  fx <- manifoldalign::rapid_ma_benchmark_fixture(seed = 3005)
  run <- function(seed) {
    set.seed(seed)
    suppressWarnings(manifoldalign::benchmark_rapid_ma(
      fx, methods = c("ssma", "multiscale")
    ))$results$oa
  }
  expect_identical(run(1), run(2))
})
