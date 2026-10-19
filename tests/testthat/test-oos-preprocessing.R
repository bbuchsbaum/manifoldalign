# oos_predict() must preprocess new rows with the statistics fitted on the
# training rows, never refit the preprocessor on the new rows.

test_that("oos_predict applies per-domain fitted preprocessors without refitting", {
  set.seed(5)
  X1 <- matrix(rnorm(20 * 3, mean = 4), 20, 3)
  X2 <- matrix(rnorm(15 * 2, mean = -2), 15, 2)
  p1 <- manifoldalign:::.preproc_fit_transform(multivarious::center(), X1)$preproc
  p2 <- manifoldalign:::.preproc_fit_transform(multivarious::center(), X2)$preproc
  V <- matrix(rnorm(5 * 2), 5, 2)
  fit <- structure(
    list(v = V, preproc = list(p1, p2), block_indices = list(1:3, 4:5)),
    class = c("multiblock_biprojector", "bi_projector", "projector")
  )

  newX <- matrix(rnorm(4 * 3, mean = 10), 4, 3)
  expected <- sweep(newX, 2, colMeans(X1)) %*% V[1:3, ]
  expect_equal(oos_predict(fit, newX, side = 1), expected)

  # A single new row is a valid query (refitting would centre it to zero).
  expect_equal(oos_predict(fit, newX[1, , drop = FALSE], side = 1),
               expected[1, , drop = FALSE])

  unfitted <- fit
  unfitted$preproc <- list(multivarious::center(), multivarious::center())
  expect_error(oos_predict(unfitted, newX, side = 1), "not fitted")
})

test_that("the preprocessing helper fits on its input and uses no deprecated API", {
  X <- matrix(c(1, 2, 3, 10, 20, 30), 3, 2)
  expect_no_warning(ft <- manifoldalign:::.preproc_fit_transform(multivarious::center(), X))
  expect_equal(ft$transformed, sweep(X, 2, colMeans(X)), ignore_attr = TRUE)
  expect_s3_class(ft$preproc, "pre_processor")
  expect_equal(multivarious::transform(ft$preproc, X + 1), sweep(X, 2, colMeans(X)) + 1,
               ignore_attr = TRUE)
})
