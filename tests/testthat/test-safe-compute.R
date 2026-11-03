test_that("safe_compute preserves warned values and contextual errors", {
  expect_warning(value <- safe_compute({warning("benign"); 42}, "failed"), "benign")
  expect_identical(value, 42)
  seen <- character()
  value <- safe_compute({warning("first"); warning("second"); 17}, "failed",
                        warning_fn = function(w) seen <<- c(seen, conditionMessage(w)))
  expect_identical(value, 17)
  expect_identical(seen, c("first", "second"))
  expect_error(safe_compute(stop("cause"), "context"), "context Original error: cause", fixed = TRUE)
})
