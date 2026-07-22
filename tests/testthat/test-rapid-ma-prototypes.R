rapid_initialize_prototypes <- function(...) {
  manifoldalign:::.rapid_initialize_prototypes(...)
}

rapid_update_prototypes <- function(...) {
  manifoldalign:::.rapid_update_prototypes(...)
}

expect_valid_prototype_bank <- function(bank) {
  K <- nrow(bank$embedding)
  expect_s3_class(bank, "rapid_ma_prototypes")
  expect_gt(K, 0L)
  expect_identical(nrow(bank$structure), K)
  expect_identical(nrow(bank$attribute), K)
  expect_identical(nrow(bank$class_prob), K)
  expect_length(bank$prior, K)
  expect_length(bank$active, K)
  expect_true(all(is.finite(bank$embedding)))
  expect_true(all(is.finite(bank$structure)))
  expect_true(all(is.finite(bank$attribute)))
  expect_true(all(is.finite(bank$class_prob)))
  expect_true(all(is.finite(bank$prior)))
  expect_equal(rowSums(bank$class_prob), rep(1, K), tolerance = 1e-12)
  expect_equal(sum(bank$prior), 1, tolerance = 1e-12)
}

test_that("rare and one-node classes retain balanced valid prototypes", {
  set.seed(1201)
  X <- matrix(rnorm(13L * 5L), nrow = 13L)
  labels <- c(rep("common", 12L), "rare")
  bank <- rapid_initialize_prototypes(
    X,
    labels,
    control = list(prototypes_per_class = 4L, unknown_prototypes = 0L)
  )

  expect_valid_prototype_bank(bank)
  expect_identical(as.vector(table(bank$prototype_class)), c(4L, 1L))
  expect_equal(
    as.vector(bank$diagnostics$class_prior_mass[c("common", "rare")]),
    c(0.5, 0.5),
    tolerance = 1e-12
  )
  rare <- which(bank$prototype_class == "rare")
  expect_length(rare, 1L)
  expect_equal(bank$embedding[rare, ], X[13L, ])
  expect_true(bank$fixed_class[[rare]])
  expect_equal(unname(bank$class_prob[rare, "rare"]), 1)
})

test_that("prototype initialization is reproducible and domain-order invariant", {
  set.seed(1202)
  X1 <- matrix(rnorm(21L * 4L), nrow = 21L)
  X2 <- matrix(rnorm(17L * 4L), nrow = 17L)
  S1 <- cbind(seq_len(21L) / 21L, sin(seq_len(21L)))
  S2 <- cbind(seq_len(17L) / 17L, sin(seq_len(17L)))
  A1 <- matrix(rnorm(21L * 3L), nrow = 21L)
  A2 <- matrix(rnorm(17L * 3L), nrow = 17L)
  y1 <- c(rep(c("forest", "water"), 10L), NA)
  y2 <- c(rep(c("water", "forest"), 8L), NA)
  control <- list(prototypes_per_class = 3L, unknown_prototypes = 2L, seed = 81L)

  forward <- rapid_initialize_prototypes(
    list(first = X1, second = X2),
    list(y1, y2),
    structures = list(S1, S2),
    attributes = list(A1, A2),
    control = control
  )
  repeated <- rapid_initialize_prototypes(
    list(first = X1, second = X2),
    list(y1, y2),
    structures = list(S1, S2),
    attributes = list(A1, A2),
    control = control
  )
  reversed <- rapid_initialize_prototypes(
    list(second = X2, first = X1),
    list(y2, y1),
    structures = list(S2, S1),
    attributes = list(A2, A1),
    control = control
  )

  expect_identical(forward$embedding, repeated$embedding)
  expect_identical(forward$prototype_class, repeated$prototype_class)
  expect_equal(reversed$embedding, forward$embedding, tolerance = 1e-12)
  expect_equal(reversed$structure, forward$structure, tolerance = 1e-12)
  expect_equal(reversed$attribute, forward$attribute, tolerance = 1e-12)
  expect_identical(reversed$prototype_class, forward$prototype_class)
  expect_equal(reversed$prior, forward$prior)
})

test_that("unknown-only and fully supervised modes are explicit", {
  set.seed(1203)
  X <- matrix(rnorm(30L * 3L), nrow = 30L)

  unknown <- rapid_initialize_prototypes(
    X,
    labels = rep(NA, nrow(X)),
    control = list(unknown_prototypes = 3L)
  )
  expect_valid_prototype_bank(unknown)
  expect_identical(nrow(unknown$embedding), 3L)
  expect_true(unknown$diagnostics$unknown_only)
  expect_true(all(!unknown$fixed_class))
  expect_true(all(unknown$prototype_class == ".unknown"))
  expect_identical(unknown$class_levels, ".unknown")

  supervised_labels <- rep(c("a", "b", "c"), each = 10L)
  supervised <- rapid_initialize_prototypes(
    X,
    supervised_labels,
    control = list(prototypes_per_class = 2L, unknown_prototypes = 5L)
  )
  expect_valid_prototype_bank(supervised)
  expect_true(supervised$diagnostics$fully_supervised)
  expect_true(all(supervised$fixed_class))
  expect_false(any(supervised$prototype_class == ".unknown"))
  expect_identical(nrow(supervised$embedding), 6L)

  no_unknown <- rapid_initialize_prototypes(
    X,
    c(supervised_labels[-30L], NA),
    control = list(prototypes_per_class = 1L, unknown_prototypes = 0L)
  )
  expect_false(any(no_unknown$prototype_class == ".unknown"))
})

test_that("empty prototypes recover without NaN or class loss", {
  set.seed(1204)
  X <- rbind(
    matrix(rnorm(36L, -2, 0.2), ncol = 3L),
    matrix(rnorm(36L, 2, 0.2), ncol = 3L)
  )
  labels <- rep(c("left", "right"), each = 12L)
  bank <- rapid_initialize_prototypes(
    X, labels, control = list(prototypes_per_class = 3L)
  )
  K <- nrow(bank$embedding)
  coupling <- Matrix::sparseMatrix(
    i = seq_len(nrow(X)), j = rep(1L, nrow(X)), x = 1,
    dims = c(nrow(X), K)
  )

  updated <- rapid_update_prototypes(bank, coupling, X, labels)

  expect_valid_prototype_bank(updated)
  expect_identical(updated$prototype_class, bank$prototype_class)
  expect_identical(updated$diagnostics$empty_before_recovery, 2:K)
  expect_identical(updated$diagnostics$n_recovered, K - 1L)
  expect_identical(updated$diagnostics$n_retained, 0L)
  expect_true(all(updated$active))
  expect_true(all(updated$recovered[2:K]))
  expect_true(all(updated$mass > 0))
  for (k in which(updated$fixed_class)) {
    expect_equal(
      unname(updated$class_prob[k, updated$prototype_class[[k]]]), 1
    )
  }
})

test_that("empty recovery is bounded and retaining old state stays finite", {
  set.seed(1205)
  X <- matrix(rnorm(40L * 4L), nrow = 40L)
  labels <- rep(c("a", "b"), each = 20L)
  bank <- rapid_initialize_prototypes(
    X,
    labels,
    control = list(prototypes_per_class = 3L, max_recoveries = 1L)
  )
  K <- nrow(bank$embedding)
  coupling <- matrix(0, nrow(X), K)
  coupling[, 1L] <- 1
  updated <- rapid_update_prototypes(bank, coupling, X, labels)

  expect_valid_prototype_bank(updated)
  expect_identical(updated$diagnostics$n_recovered, 1L)
  expect_identical(updated$diagnostics$n_retained, K - 2L)
  expect_true(all(updated$active))
  expect_true(all(updated$mass > 0))
})

test_that("unknown prototypes learn soft classes while fixed classes stay fixed", {
  set.seed(1206)
  X <- rbind(
    matrix(rnorm(30L, -1, 0.1), ncol = 3L),
    matrix(rnorm(30L, 1, 0.1), ncol = 3L),
    matrix(rnorm(6L, 0, 0.1), ncol = 3L)
  )
  labels <- c(rep("a", 10L), rep("b", 10L), NA, NA)
  bank <- rapid_initialize_prototypes(
    X,
    labels,
    control = list(prototypes_per_class = 1L, unknown_prototypes = 1L)
  )
  unknown <- which(!bank$fixed_class)
  coupling <- matrix(0, nrow(X), nrow(bank$embedding))
  coupling[, unknown] <- 1

  updated <- rapid_update_prototypes(bank, coupling, X, labels)

  expect_gt(updated$class_prob[unknown, "a"], 0.45)
  expect_gt(updated$class_prob[unknown, "b"], 0.45)
  expect_lt(updated$class_prob[unknown, ".unknown"], 0.01)
  for (k in which(updated$fixed_class)) {
    expect_equal(
      unname(updated$class_prob[k, updated$prototype_class[[k]]]), 1
    )
    expect_equal(sum(updated$class_prob[k, ]), 1)
  }
})

test_that("embedding, structure, and attribute summaries follow coupling mass", {
  X <- rbind(c(0, 1), c(2, 3), c(10, 11), c(12, 13))
  S <- matrix(c(0, 2, 8, 10), ncol = 1L)
  A <- cbind(c(1, 1, 5, 5), c(2, 4, 6, 8))
  labels <- c("a", "a", "b", "b")
  bank <- rapid_initialize_prototypes(
    X,
    labels,
    structures = S,
    attributes = A,
    control = list(prototypes_per_class = 1L)
  )
  Q <- rbind(c(1, 0), c(1, 0), c(0, 2), c(0, 1))
  updated <- rapid_update_prototypes(
    bank, Q, X, labels, structures = S, attributes = A
  )

  expect_equal(updated$embedding[1L, ], colMeans(X[1:2, , drop = FALSE]))
  expect_equal(updated$embedding[2L, ], (2 * X[3L, ] + X[4L, ]) / 3)
  expect_equal(updated$structure[2L, ], (2 * S[3L, ] + S[4L, ]) / 3)
  expect_equal(updated$attribute[2L, ], (2 * A[3L, ] + A[4L, ]) / 3)
  expect_equal(updated$diagnostics$total_coupling_mass, 5)
})

test_that("prototype validation rejects silent dimensional and class drift", {
  X <- matrix(seq_len(60L), nrow = 20L)
  labels <- rep(c("a", "b"), each = 10L)
  bank <- rapid_initialize_prototypes(
    X, labels, control = list(prototypes_per_class = 1L)
  )
  Q <- matrix(0.5, nrow(X), nrow(bank$embedding))

  expect_error(
    rapid_update_prototypes(bank, Q, cbind(X, 1), labels),
    "dimensions"
  )
  expect_error(
    rapid_update_prototypes(bank, Q, X, replace(labels, 1L, "new-class")),
    "absent from the prototype bank"
  )
  expect_error(
    rapid_initialize_prototypes(X, labels, control = list(unknown_level = "a")),
    "collides"
  )
})
