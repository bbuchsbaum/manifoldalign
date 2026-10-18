test_that(".local_seed restores the caller's RNG stream on exit and error", {
  local_seed <- get(".local_seed", envir = asNamespace("manifoldalign"))
  f <- function(s) {
    local_seed(s)
    a <- runif(2)
    local_seed(s + 1)
    c(a, runif(1))
  }
  set.seed(42)
  expected_next <- runif(1)
  set.seed(42)
  first <- f(5)
  expect_identical(runif(1), expected_next)
  set.seed(42)
  expect_identical(f(5), first)
  g <- function() {
    local_seed(3)
    stop("boom")
  }
  set.seed(42)
  expect_error(g(), "boom")
  expect_identical(runif(1), expected_next)
})

test_that("benchmark registries no longer reset the caller's RNG", {
  case <- get("synthetic_graph_alignment_case", envir = asNamespace("manifoldalign"))
  set.seed(11)
  expected_next <- runif(1)
  set.seed(11)
  first <- case(n = 20L, seed = 3L)
  expect_identical(runif(1), expected_next)
  second <- case(n = 20L, seed = 3L)
  expect_identical(first, second)
})
