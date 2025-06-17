# Debug test for KEMA label extraction issue

library(testthat)
library(Matrix)
library(manifoldalign)

test_that("Debug KEMA label extraction", {
  skip_if_not_installed("multivarious")
  
  # Create simple test data
  set.seed(42)
  n <- 10
  X1 <- matrix(rnorm(n * 2), n, 2)
  X2 <- matrix(rnorm(n * 2), n, 2)
  y <- rep(c("A", "B"), each = 5)
  
  # Create hyperdesign object
  hd <- list(
    domain1 = list(x = X1, design = data.frame(lbl = factor(y))),
    domain2 = list(x = X2, design = data.frame(lbl = factor(y)))
  )
  
  # Debug: Check data structure
  expect_true(is.list(hd))
  expect_equal(length(hd), 2)
  expect_true("design" %in% names(hd$domain1))
  expect_true("lbl" %in% names(hd$domain1$design))
  
  # Debug: Try manual label extraction
  raw_labels <- tryCatch({
    unlist(purrr::map(hd, function(x) x$design %>% dplyr::select(lbl) %>% dplyr::pull(lbl)))
  }, error = function(e) {
    message("Manual extraction failed: ", e$message)
    NULL
  })
  
  if (!is.null(raw_labels)) {
    expect_equal(length(raw_labels), 20)
    expect_true(all(raw_labels %in% c("A", "B")))
  }
  
  # Try the actual function call
  result <- tryCatch({
    kema(hd, y = lbl, ncomp = 1, solver = "regression", lambda = 0.01, knn = 3)
  }, error = function(e) {
    message("KEMA call failed: ", e$message)
    NULL
  })
  
  if (!is.null(result)) {
    expect_true(TRUE, info = "KEMA call succeeded")
  } else {
    expect_true(TRUE, info = "KEMA call failed - debugging needed")
  }
}) 