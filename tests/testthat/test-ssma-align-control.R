library(testthat)

test_that("ssma_serial_control validates and stores fields", {
  sc <- ssma_serial_control(
    enabled = TRUE,
    row_whiten = "ar1",
    lag_mode = "hard",
    lag_exclusion = 2,
    lag_scale = 3
  )

  expect_s3_class(sc, "ssma_serial_control")
  expect_true(isTRUE(sc$enabled))
  expect_equal(sc$row_whiten, "ar1")
  expect_equal(sc$lag_mode, "hard")
  expect_equal(sc$lag_exclusion, 2L)

  expect_error(ssma_serial_control(ar1_shrink = 1), "ar1_shrink")
  expect_error(ssma_serial_control(lag_exclusion = -1), "lag_exclusion")
  expect_error(ssma_serial_control(lag_scale = 0), "lag_scale")
})

test_that("resolve_ssma_align_control merges serial config and rejects unknown fields", {
  ctrl <- manifoldalign:::resolve_ssma_align_control(list(
    knn = 9,
    rank_per_domain = 64,
    solver = "operator",
    operator_tol = 1e-5,
    serial_control = list(enabled = TRUE, row_whiten = "ar1", lag_mode = "hard", lag_exclusion = 1)
  ))

  expect_s3_class(ctrl, "ssma_align_control")
  expect_equal(ctrl$knn, 9L)
  expect_equal(ctrl$rank_per_domain, 64L)
  expect_equal(ctrl$solver, "operator")
  expect_equal(ctrl$operator_tol, 1e-5)
  expect_s3_class(ctrl$serial_control, "ssma_serial_control")
  expect_true(isTRUE(ctrl$serial_control$enabled))
  expect_equal(ctrl$serial_control$lag_exclusion, 1L)

  expect_error(
    manifoldalign:::resolve_ssma_align_control(list(not_a_field = 1)),
    "Unknown control argument"
  )
  expect_error(
    manifoldalign:::resolve_ssma_align_control(list(serial_control = list(not_here = 1))),
    "Unknown serial control argument"
  )
  expect_error(ssma_align_control(operator_tol = 0), "operator_tol")
  expect_error(ssma_align_control(operator_maxiter = 0), "operator_maxiter")
})
