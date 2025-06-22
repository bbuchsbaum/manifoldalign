library(testthat)
library(manifoldalign)

skip_if_not_installed("lpSolve")

test_that("lpSolve fallback matches LEMON solver on 3x3 problem", {
  if (!manifoldalign:::has_lemon_network_simplex()) {
    skip("LEMON not available")
  }
  set.seed(1)
  cost <- matrix(runif(9), 3, 3)
  a <- rep(1/3, 3)
  b <- rep(1/3, 3)
  sol_lemon <- manifoldalign:::network_simplex_ot_cpp(cost, a, b)
  sol_lp <- manifoldalign:::fallback_network_simplex_ot_cpp(cost, a, b)
  expect_equal(sol_lp, sol_lemon, tolerance = 1e-6)
})

test_that("lpSolve fallback matches on rectangular problem", {
  if (!manifoldalign:::has_lemon_network_simplex()) {
    skip("LEMON not available")
  }
  set.seed(2)
  cost <- matrix(runif(12), 3, 4)
  a <- rep(1/3, 3)
  b <- rep(1/4, 4)
  sol_lemon <- manifoldalign:::network_simplex_ot_cpp(cost, a, b)
  sol_lp <- manifoldalign:::fallback_network_simplex_ot_cpp(cost, a, b)
  expect_equal(sol_lp, sol_lemon, tolerance = 1e-6)
})
