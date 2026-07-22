rapid_uot_build_cost <- function(...) {
  manifoldalign:::.rapid_uot_build_cost(...)
}

rapid_uot_solve <- function(...) {
  manifoldalign:::.rapid_uot_solve(...)
}

rapid_sparse_cost_matrix <- function(cost_state) {
  cost <- cost_state$cost
  out <- matrix(Inf, cost$n_rows, cost$n_cols)
  rows <- rep(seq_len(cost$n_rows), diff(cost$row_ptr))
  out[cbind(rows, cost$col_idx)] <- cost$cost
  out
}

test_that("composite costs combine all requested evidence exactly", {
  reference <- rbind(c(0, 0), c(2, 0), c(0, 2), c(2, 2))
  labels <- c("a", "a", "b", "b")
  reference_structure <- matrix(c(0, 2, 8, 10), ncol = 1L)
  reference_attribute <- matrix(c(1, 3, 7, 9), ncol = 1L)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    reference,
    labels,
    structures = reference_structure,
    attributes = reference_attribute,
    control = list(prototypes_per_class = 1L)
  )
  X <- rbind(c(0.5, 0), c(1.5, 2))
  S <- matrix(c(1, 9), ncol = 1L)
  A <- matrix(c(2, 8), ncol = 1L)
  position <- rbind(c(0, 0), c(2, 2))
  prototype_position <- rbind(c(0, 1), c(2, 1))
  node_labels <- c("a", NA)
  confidence <- c(0.9, 0)
  anchors <- c(NA, 1L)
  control <- list(
    latent_weight = 2,
    structure_weight = 3,
    attribute_weight = 4,
    position_weight = 5,
    label_weight = 6,
    anchor_weight = 7,
    label_confidence_threshold = 0.8
  )

  built <- rapid_uot_build_cost(
    X,
    prototypes,
    structure = S,
    attributes = A,
    positions = position,
    prototype_positions = prototype_position,
    labels = node_labels,
    label_confidence = confidence,
    anchors = anchors,
    control = control,
    mode = "dense"
  )
  expected <- matrix(0, nrow(X), nrow(prototypes$embedding))
  for (i in seq_len(nrow(X))) {
    for (j in seq_len(nrow(prototypes$embedding))) {
      expected[i, j] <-
        2 * mean((X[i, ] - prototypes$embedding[j, ])^2) +
        3 * mean((S[i, ] - prototypes$structure[j, ])^2) +
        4 * mean((A[i, ] - prototypes$attribute[j, ])^2) +
        5 * mean((position[i, ] - prototype_position[j, ])^2)
      if (i == 1L) {
        expected[i, j] <- expected[i, j] +
          6 * confidence[[i]] * (1 - prototypes$class_prob[j, "a"])
      }
      if (i == 2L && j != anchors[[i]]) {
        expected[i, j] <- expected[i, j] + 7
      }
    }
  }

  expect_equal(built$cost, expected, tolerance = 1e-12)
  expect_setequal(
    built$diagnostics$active_views,
    c("latent", "structure", "attribute", "position")
  )
})

test_that("full-support sparse and dense UOT solutions agree", {
  set.seed(1302)
  reference <- rbind(
    matrix(rnorm(24L, -1, 0.2), ncol = 3L),
    matrix(rnorm(24L, 1, 0.2), ncol = 3L)
  )
  reference_labels <- rep(c("a", "b"), each = 8L)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    reference,
    reference_labels,
    control = list(prototypes_per_class = 2L)
  )
  X <- matrix(rnorm(33L), nrow = 11L)
  control <- list(
    q = nrow(prototypes$embedding),
    epsilon = 0.35,
    rho_source = 0.9,
    rho_target = 1.2,
    max_iter = 3000L,
    tol = 1e-10
  )

  dense <- rapid_uot_solve(X, prototypes, control = control, mode = "dense")
  sparse <- rapid_uot_solve(X, prototypes, control = control, mode = "sparse")

  expect_equal(
    as.matrix(sparse$coupling),
    as.matrix(dense$coupling),
    tolerance = 1e-7
  )
  expect_equal(sparse$row_mass, dense$row_mass, tolerance = 1e-7)
  expect_equal(
    sparse$prototype_mass, dense$prototype_mass, tolerance = 1e-7
  )
  expect_identical(sparse$diagnostics$backend, "sparse_csr_csc")
  expect_true(sparse$diagnostics$converged)
})

test_that("CSR and CSC sparse costs retain the same bounded support", {
  set.seed(1303)
  reference <- matrix(rnorm(60L), nrow = 20L)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    reference,
    rep(c("a", "b", "c", "d"), each = 5L),
    control = list(prototypes_per_class = 2L)
  )
  X <- matrix(rnorm(39L), nrow = 13L)
  built <- rapid_uot_build_cost(
    X, prototypes, control = list(q = 3L), mode = "sparse"
  )
  cost <- built$cost
  csr <- rapid_sparse_cost_matrix(built)
  csc <- matrix(Inf, cost$n_rows, cost$n_cols)
  columns <- rep(seq_len(cost$n_cols), diff(cost$col_ptr))
  csc[cbind(cost$row_idx, columns)] <- cost$cost_csc

  expect_equal(csc, csr)
  expect_true(all(diff(cost$row_ptr) >= 1L))
  expect_true(all(diff(cost$col_ptr) >= 1L))
  expect_lte(length(cost$cost), nrow(X) * 3L + nrow(prototypes$embedding))
  expect_false(built$diagnostics$dense_node_prototype_retained)
})

test_that("outliers drop transport mass under partial overlap", {
  set.seed(1304)
  reference <- rbind(
    matrix(rnorm(60L, -1, 0.15), ncol = 3L),
    matrix(rnorm(60L, 1, 0.15), ncol = 3L)
  )
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    reference,
    rep(c("a", "b"), each = 20L),
    control = list(prototypes_per_class = 1L)
  )
  X <- rbind(
    matrix(rnorm(90L, -1, 0.2), ncol = 3L),
    matrix(rnorm(90L, 1, 0.2), ncol = 3L),
    matrix(rnorm(30L, 12, 0.2), ncol = 3L)
  )
  fit <- rapid_uot_solve(
    X,
    prototypes,
    control = list(
      q = 2L,
      epsilon = 0.1,
      rho_source = 0.08,
      rho_target = 1,
      max_iter = 3000L,
      tol = 1e-9
    )
  )
  inlier_mass <- median(fit$row_mass[1:60])
  outlier_mass <- median(fit$row_mass[61:70])

  expect_lt(outlier_mass, inlier_mass * 1e-4)
  expect_gt(fit$diagnostics$dropped_fraction, 0.05)
  expect_true(all(is.finite(fit$coupling@x)))
})

test_that("prototype coverage is added or declared inactive explicitly", {
  set.seed(1305)
  reference <- matrix(rnorm(60L), nrow = 20L)
  reference_labels <- rep(c("a", "b"), each = 10L)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    reference,
    reference_labels,
    control = list(prototypes_per_class = 4L)
  )
  X <- matrix(rnorm(15L), nrow = 5L)
  covered <- rapid_uot_build_cost(
    X, prototypes, control = list(q = 1L), mode = "sparse"
  )

  expect_true(all(covered$diagnostics$prototype_edge_counts > 0L))
  expect_gt(covered$diagnostics$coverage_edges_added, 0L)
  expect_length(covered$inactive_prototypes, 0L)

  constrained <- rapid_uot_build_cost(
    X,
    prototypes,
    labels = rep("a", nrow(X)),
    control = list(q = 1L, hard_labels = TRUE),
    mode = "sparse"
  )
  inactive_classes <- prototypes$prototype_class[constrained$inactive_prototypes]
  expect_true(all(inactive_classes == "b"))
  expect_true(all(
    prototypes$prototype_class[constrained$prototype_index] == "a"
  ))
  expect_true(all(constrained$diagnostics$prototype_edge_counts > 0L))
})

test_that("hard anchors constrain support and invalid anchors fail", {
  X <- matrix(seq_len(36L) / 10, nrow = 12L)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    X,
    rep(c("a", "b"), each = 6L),
    control = list(prototypes_per_class = 2L)
  )
  anchors <- rep(NA_integer_, nrow(X))
  anchors[1:3] <- c(1L, 2L, 3L)
  built <- rapid_uot_build_cost(
    X,
    prototypes,
    anchors = anchors,
    control = list(q = 4L, hard_anchors = TRUE),
    mode = "sparse"
  )
  support <- rapid_sparse_cost_matrix(built)

  for (i in 1:3) {
    global_columns <- built$prototype_index[which(is.finite(support[i, ]))]
    expect_identical(global_columns, anchors[[i]])
  }
  expect_error(
    rapid_uot_build_cost(X, prototypes, anchors = rep(99L, nrow(X))),
    "valid integer prototype"
  )
})

test_that("retained memory is O(n q) rather than O(n K)", {
  set.seed(1307)
  K <- 80L
  reference <- matrix(rnorm(K * 6L), nrow = K)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    reference,
    labels = rep(NA, K),
    control = list(unknown_prototypes = K)
  )
  n <- 5000L
  q <- 4L
  X <- matrix(rnorm(n * 6L), nrow = n)
  built <- rapid_uot_build_cost(
    X,
    prototypes,
    control = list(q = q, block_size = 256L),
    mode = "sparse"
  )

  expect_lte(built$diagnostics$stored_edges, n * q + K)
  expect_lte(built$diagnostics$peak_block_entries, 256L * K)
  expect_lt(as.numeric(object.size(built$cost)), n * K * 8)
  expect_identical(length(built$cost$row_ptr), n + 1L)
  expect_identical(length(built$cost$col_ptr), K + 1L)
  expect_false(built$diagnostics$dense_node_prototype_retained)
})

test_that("rectangular and stiff UOT problems remain finite", {
  set.seed(1308)
  reference <- matrix(rnorm(24L), nrow = 8L)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    reference,
    rep(c("a", "b"), each = 4L),
    control = list(prototypes_per_class = 2L)
  )
  X <- matrix(rnorm(33L) * 100, nrow = 11L)
  fit <- rapid_uot_solve(
    X,
    prototypes,
    control = list(
      q = 3L,
      epsilon = 1e-4,
      rho_source = 0.01,
      rho_target = 100,
      cost_cap = 1e4,
      max_iter = 5000L,
      tol = 1e-8
    )
  )

  expect_identical(dim(fit$coupling), c(11L, 4L))
  expect_true(all(is.finite(fit$coupling@x)))
  expect_true(all(is.finite(fit$row_mass)))
  expect_true(fit$diagnostics$converged)

  impossible <- matrix(1e200, nrow = 3L, ncol = 3L)
  expect_error(
    rapid_uot_build_cost(impossible, prototypes),
    "without a candidate"
  )
})

test_that("identical UOT problems reuse an exact warm result", {
  X <- matrix(seq_len(30L) / 10, nrow = 10L)
  prototypes <- manifoldalign:::.rapid_initialize_prototypes(
    X,
    rep(c("a", "b"), each = 5L),
    control = list(prototypes_per_class = 1L)
  )
  control <- list(q = 2L, epsilon = 0.2)
  cold <- rapid_uot_solve(X, prototypes, control = control)
  warm <- rapid_uot_solve(
    X, prototypes, control = control, warm_start = cold
  )

  expect_true(warm$diagnostics$warm_start_supplied)
  expect_true(warm$diagnostics$warm_start_reused)
  expect_identical(warm$solution, cold$solution)
  expect_equal(warm$coupling, cold$coupling)

  changed <- rapid_uot_solve(
    X + 0.01, prototypes, control = control, warm_start = cold
  )
  expect_true(changed$diagnostics$warm_start_supplied)
  expect_false(changed$diagnostics$warm_start_reused)
})
