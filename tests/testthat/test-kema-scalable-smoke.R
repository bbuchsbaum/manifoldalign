make_operator_kema_fixture <- function(n_per_domain = 18, n_features = 4) {
  set.seed(91)
  X1 <- matrix(rnorm(n_per_domain * n_features), n_per_domain, n_features)
  X2 <- matrix(
    rnorm(n_per_domain * n_features) + 0.25,
    n_per_domain,
    n_features
  )
  labels <- factor(rep(1:3, length.out = n_per_domain))
  multidesign::hyperdesign(list(
    d1 = multidesign::multidesign(
      X1,
      data.frame(lbl = labels)
    ),
    d2 = multidesign::multidesign(
      X2,
      data.frame(lbl = labels)
    )
  ))
}

test_that("operator KEMA matches the full matrix backend on a toy oracle", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  hd <- make_operator_kema_fixture()
  common <- list(
    data = hd,
    ncomp = 2,
    knn = 4,
    sigma = 0.5,
    lambda = 1e-3
  )

  full <- rlang::inject(manifoldalign::kema(
    !!!common,
    y = lbl,
    backend = "full_exact"
  ))
  operator <- rlang::inject(manifoldalign::kema(
    !!!common,
    y = lbl,
    backend = "operator_exact"
  ))

  expect_true(full$fidelity$passed)
  expect_true(operator$fidelity$passed)
  expect_lte(full$fidelity$max_rel_residual, 1e-6)
  expect_lte(operator$fidelity$max_rel_residual, 1e-6)
  expect_equal(
    operator$eigenvalues$values,
    full$eigenvalues$values,
    tolerance = 1e-7
  )

  Q_full <- qr.Q(qr(as.matrix(full$s)))
  Q_operator <- qr.Q(qr(as.matrix(operator$s)))
  correlations <- svd(
    crossprod(Q_full, Q_operator),
    nu = 0,
    nv = 0
  )$d
  expect_gt(min(correlations), 1 - 1e-7)
})
