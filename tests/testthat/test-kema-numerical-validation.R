library(testthat)
library(Matrix)
library(multidesign)
library(manifoldalign)

generate_two_domain_spiral <- function(n_per_domain = 40, noise_level = 0.05) {
  t <- seq(0, 4 * pi, length.out = n_per_domain)
  X1 <- cbind(
    t * cos(t) + rnorm(n_per_domain, 0, noise_level),
    t * sin(t) + rnorm(n_per_domain, 0, noise_level)
  )
  X2 <- cbind(
    0.8 * t * cos(t + pi / 4) + rnorm(n_per_domain, 0, noise_level),
    0.8 * t * sin(t + pi / 4) + rnorm(n_per_domain, 0, noise_level)
  )
  labels <- ifelse(t < 2 * pi, "early", "late")

  list(
    domain1 = list(x = X1, labels = labels),
    domain2 = list(x = X2, labels = labels),
    all_labels = c(labels, labels)
  )
}

create_kema_hyperdesign <- function(data) {
  blocks <- lapply(data[c("domain1", "domain2")], function(domain) {
    multidesign::multidesign(
      domain$x,
      data.frame(lbl = factor(domain$labels))
    )
  })
  multidesign::hyperdesign(blocks)
}

expect_fidelity_gated_kema <- function(fit, ncomp) {
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(ncol(fit$s), ncomp)
  expect_equal(ncol(fit$alpha), ncomp)
  expect_true(all(is.finite(fit$s)))
  expect_true(all(is.finite(fit$alpha)))
  expect_true(fit$fidelity$passed)
  expect_lte(fit$fidelity$max_rel_residual, 1e-6)
  expect_lte(fit$fidelity$max_b_orth_offdiag, 1e-6)
  expect_true(all(
    fit$eigenvalues$values > fit$fidelity$eigenvalue_zero_tol
  ))
  expect_identical(fit$eigenvalues$solver, "exact")
}

test_that("spiral fixture and hyperdesign construction are deterministic", {
  set.seed(42)
  data <- generate_two_domain_spiral(n_per_domain = 12)
  hd <- create_kema_hyperdesign(data)

  expect_equal(length(hd), 2)
  expect_equal(
    unname(vapply(hd, function(x) nrow(x$x), integer(1))),
    c(12L, 12L)
  )
  expect_equal(
    unname(vapply(hd, function(x) ncol(x$x), integer(1))),
    c(2L, 2L)
  )
  expect_true(all(data$all_labels %in% c("early", "late")))
  expect_true(all(vapply(hd, function(x) is.factor(x$design$lbl), logical(1))))
})

test_that("public KEMA returns only fidelity-gated spiral components", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(42)
  data <- generate_two_domain_spiral(n_per_domain = 30)
  hd <- create_kema_hyperdesign(data)

  fit <- kema(
    hd,
    y = lbl,
    ncomp = 3,
    kernel = kernlab::rbfdot(sigma = 0.5),
    lambda = 1e-3,
    knn = 5,
    u = 0.5
  )

  expect_fidelity_gated_kema(fit, ncomp = 3)
  expect_equal(nrow(fit$s), length(data$all_labels))
  expect_gt(stats::var(as.numeric(fit$s[, 1L])), 1e-10)
})

test_that("semi-supervised KEMA resolves a multi-dimensional null space", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(42)
  X1 <- matrix(rnorm(40), 20, 2)
  X2 <- matrix(rnorm(40), 20, 2)
  labels <- c(rep(c("A", "B"), each = 5), rep(NA_character_, 10))
  hd <- create_kema_hyperdesign(list(
    domain1 = list(x = X1, labels = labels),
    domain2 = list(x = X2, labels = labels)
  ))

  fit <- kema(
    hd,
    y = lbl,
    ncomp = 2,
    lambda = 1e-2,
    knn = 3,
    u = 0.5
  )

  expect_fidelity_gated_kema(fit, ncomp = 2)
  expect_identical(fit$fidelity$eigensolver, "eigencore_dense_full")
  expect_true(fit$fidelity$eigencore_certificate_passed)
  expect_gte(fit$fidelity$nullity_deflated, 1L)
  expect_gt(fit$fidelity$n_eigenpairs_examined, 3)
  expect_gt(fit$fidelity$eigenvalue_zero_tol, 0)
})

test_that("KEMA skips quotient pairs that fail the original-pencil residual", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(1)
  X1 <- matrix(rnorm(60), 30, 2)
  X2 <- matrix(rnorm(60), 30, 2)
  labels <- sample(c("A", "B"), 30, TRUE)
  hd <- multidesign::hyperdesign(list(
    d1 = multidesign::multidesign(X1, data.frame(labels = labels)),
    d2 = multidesign::multidesign(X2, data.frame(labels = labels))
  ))

  fit <- fit_many(kema_aligner(), hd, y = labels, ncomp = 2)

  expect_true(fit$fidelity$passed)
  expect_lte(fit$fidelity$max_rel_residual, 1e-6)
  expect_gte(fit$fidelity$n_residual_rejected, 1L)
  expect_true(all(fit$eigenvalues$values > fit$fidelity$eigenvalue_zero_tol))
})

test_that("ill-conditioned KEMA either meets the gate or is rejected", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(42)
  X1 <- matrix(c(rep(1, 10), rnorm(10, 0, 1e-3)), 10, 2)
  X2 <- matrix(c(rep(2, 10), rnorm(10, 0, 1e-3)), 10, 2)
  labels <- rep(c("A", "B"), each = 5)
  hd <- create_kema_hyperdesign(list(
    domain1 = list(x = X1, labels = labels),
    domain2 = list(x = X2, labels = labels)
  ))

  outcome <- tryCatch(
    kema(
      hd,
      y = lbl,
      ncomp = 1,
      kernel = kernlab::rbfdot(sigma = 0.73),
      lambda = 1e-2,
      knn = 3,
      u = 0.5
    ),
    error = identity
  )

  if (inherits(outcome, "error")) {
    expect_match(
      conditionMessage(outcome),
      "No KEMA fit was returned|Failed to extract.*non-trivial"
    )
  } else {
    expect_fidelity_gated_kema(outcome, ncomp = 1)
  }
})

test_that("withdrawn KEMA controls are rejected rather than ignored", {
  set.seed(11)
  data <- generate_two_domain_spiral(n_per_domain = 12)
  hd <- create_kema_hyperdesign(data)

  expect_error(
    kema(hd, y = lbl, solver = "regression"),
    "never selected a distinct regression implementation"
  )
  expect_error(
    kema(hd, y = lbl, use_laplacian = FALSE),
    "Unsupported KEMA argument.*use_laplacian"
  )
  expect_error(
    kema(hd, y = lbl, rweight = 0.1),
    "Unsupported KEMA argument.*rweight"
  )
  expect_error(
    kema(hd, y = lbl, centre_kernel = TRUE),
    "Unsupported KEMA argument.*centre_kernel"
  )
})

test_that("withdrawn OOS validator cannot fabricate a passing result", {
  expect_true(
    "validate_out_of_sample_reconstruction" %in%
      getNamespaceExports("manifoldalign")
  )

  set.seed(4102)
  rng_before <- .Random.seed

  expect_warning(
    result <- try(validate_out_of_sample_reconstruction(), silent = TRUE),
    regexp = "withdrawn|deprecated"
  )
  expect_s3_class(result, "try-error")
  expect_match(as.character(result), "withdrawn")
  expect_identical(.Random.seed, rng_before)
})

test_that("withdrawn paper-eigenvalue validator cannot claim accuracy", {
  expect_warning(
    result <- try(validate_kema_eigenvalues(), silent = TRUE),
    regexp = "withdrawn|deprecated"
  )
  expect_s3_class(result, "try-error")
  expect_match(as.character(result), "invalid|withdrawn")
})

test_that("KEMA validation suite runs enforced numerical checks only", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  suite_body <- paste(deparse(body(run_kema_validation_suite)), collapse = "\n")
  expect_false(grepl(
    "validate_out_of_sample_reconstruction",
    suite_body,
    fixed = TRUE
  ))
  expect_false(grepl("reconstruction_validation", suite_body, fixed = TRUE))
  expect_true(grepl("withdrawn_validations", suite_body, fixed = TRUE))
  expect_true(grepl("expected_error plus random noise", suite_body, fixed = TRUE))

  suite <- run_kema_validation_suite(verbose = FALSE)
  expect_true(suite$overall_success)
  expect_true(suite$fidelity_validation$success)
  expect_lte(suite$fidelity_validation$max_rel_residual, 1e-6)
  expect_lte(suite$fidelity_validation$max_b_orth_offdiag, 1e-6)
  expect_true(suite$backend_agreement$success)
  expect_gt(min(suite$backend_agreement$canonical_correlations), 1 - 1e-6)
  expect_identical(
    suite$withdrawn_validations$paper_eigenvalue_ratios$status,
    "withdrawn"
  )
})

test_that("choose_sigma returns finite positive fallbacks", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  expect_gt(choose_sigma(X), 0)
  expect_equal(choose_sigma(matrix(c(1, 2), 1, 2)), 1)
  expect_equal(
    choose_sigma(matrix(rep(c(1, 2), 5), 5, 2, byrow = TRUE)),
    1
  )
})
