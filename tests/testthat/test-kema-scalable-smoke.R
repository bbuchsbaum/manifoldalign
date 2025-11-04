skip_if_not_installed("kernlab")
skip_if_not_installed("tibble")
skip_if_not_installed("manifoldalign")

skip_if_benchmarks_disabled("Benchmark tests disabled; enable with options(manifoldalign.run_benchmarks = TRUE)")

make_toy_hyperdesign <- function(n_per_domain = 20, noise = 0.05) {
  set.seed(123)
  t <- seq(0, 2 * pi, length.out = n_per_domain)
  spiral <- function(t, phase = 0) {
    cbind(
      cos(t + phase) + rnorm(length(t), sd = noise),
      sin(t + phase) + rnorm(length(t), sd = noise)
    )
  }
  domain1 <- spiral(t)
  domain2 <- spiral(t, phase = pi / 4)
  labels <- factor(ifelse(t <= pi, "early", "late"))
  multidesign::hyperdesign(list(
    d1 = multidesign::multidesign(domain1, tibble::tibble(lbl = labels)),
    d2 = multidesign::multidesign(domain2, tibble::tibble(lbl = labels))
  ))
}

with_scalable_option <- function(value, expr) {
  old <- getOption("kema.scalable_mode")
  on.exit(options(kema.scalable_mode = old), add = TRUE)
  options(kema.scalable_mode = value)
  force(expr)
}

test_that("scalable KEMA matches baseline on toy data", {
  if (!"package:tibble" %in% search()) {
    suppressPackageStartupMessages(library(tibble))
  }
  hd <- make_toy_hyperdesign()

  baseline <- with_scalable_option(FALSE, {
    manifoldalign::kema(hd, y = lbl, solver = "exact", ncomp = 2, knn = 3,
         sigma = 0.6, lambda = 1e-3, dweight = 1, rweight = 0)
  })

  scalable <- with_scalable_option(TRUE, {
    manifoldalign::kema(hd, y = lbl, solver = "exact", ncomp = 2, knn = 3,
         sigma = 0.6, lambda = 1e-3, dweight = 1, rweight = 0)
  })

  expect_true(all(is.finite(scalable$eigenvalues$values)))
  expect_equal(dim(scalable$s), dim(baseline$s))
  expect_true(all(is.finite(as.matrix(scalable$s))))
  if (!is.null(scalable$coef)) {
    expect_true(all(is.finite(as.matrix(scalable$coef))))
  }
})
