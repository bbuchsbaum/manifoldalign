make_rapid_adapter_features <- function(latent, p, seed) {
  set.seed(seed)
  basis <- cbind(latent, latent[, 1L]^2, latent[, 2L]^2)
  basis %*% matrix(stats::rnorm(ncol(basis) * p), ncol(basis), p) +
    matrix(stats::rnorm(nrow(latent) * p, sd = 0.05), nrow(latent))
}

make_rapid_adapter_latent <- function(n, phase = 0) {
  index <- seq_len(n)
  cbind(
    x = index / max(n, 1L),
    y = sin(index / 5 + phase) + 0.15 * cos(index / 2)
  )
}

make_rapid_adapter_labels <- function(n, per_class = 3L) {
  truth <- rep(c("rural", "urban", "water"), length.out = n)
  labels <- rep(NA_character_, n)
  for (class in unique(truth)) {
    labels[head(which(truth == class), per_class)] <- class
  }
  list(truth = truth, labels = labels)
}

make_rapid_adapter_attributes <- function(n, phase = 0) {
  index <- seq_len(n)
  output <- data.frame(
    brightness = sin(index / 7 + phase),
    landform = factor(rep(c("flat", "ridge", "valley"), length.out = n)),
    protected = rep(c(TRUE, FALSE), length.out = n)
  )
  output$brightness[[min(5L, n)]] <- NA_real_
  output$landform[[min(9L, n)]] <- NA
  output
}

test_that("RAPID-MA adapter exposes native pair and multiset semantics", {
  latent_a <- make_rapid_adapter_latent(36L)
  latent_b <- make_rapid_adapter_latent(41L, 0.2)
  X_a <- make_rapid_adapter_features(latent_a, 8L, 6101L)
  X_b <- make_rapid_adapter_features(latent_b, 6L, 6102L)
  y_a <- make_rapid_adapter_labels(nrow(X_a))$labels
  y_b <- make_rapid_adapter_labels(nrow(X_b))$labels
  control <- list(max_iter = 1L, min_iter = 1L, seed = 11L)
  algo <- rapid_ma_aligner()

  expect_s3_class(algo, "rapid_ma_aligner")
  expect_identical(
    aligner_capabilities(algo),
    list(group = "perm", supports_multi = TRUE)
  )
  pair <- fit_pair(
    algo, X_a, X_b,
    labels = list(y_a, y_b),
    positions = list(latent_a, latent_b),
    ncomp = 5L,
    control = control,
    fine_match = FALSE
  )
  direct <- rapid_ma(
    list(i = X_a, j = X_b),
    labels = list(y_a, y_b),
    positions = list(latent_a, latent_b),
    ncomp = 5L,
    control = control
  )

  expect_s3_class(pair, "rapid_ma_pair_fit")
  expect_equal(pair$core$s, direct$s, tolerance = 1e-12)
  expect_equal(pair_loss(pair), direct$objective, tolerance = 1e-12)
  expect_identical(latent_dim(pair), 5L)

  third_latent <- make_rapid_adapter_latent(32L, -0.2)
  third <- make_rapid_adapter_features(third_latent, 7L, 6103L)
  native <- fit_many(
    algo,
    list(a = X_a, b = X_b, c = third),
    labels = list(
      make_rapid_adapter_labels(nrow(X_a))$labels,
      make_rapid_adapter_labels(nrow(X_b))$labels,
      make_rapid_adapter_labels(nrow(third))$labels
    ),
    positions = list(latent_a, latent_b, third_latent),
    ncomp = 4L,
    control = control
  )
  expect_s3_class(native, "rapid_ma")
  expect_identical(native$domain_names, c("a", "b", "c"))
  expect_identical(vapply(native$scores, nrow, integer(1)),
                   c(a = 36L, b = 41L, c = 32L))
})

test_that("OOS reapplication is exact and held-out interpolation is bounded", {
  latent_a <- make_rapid_adapter_latent(48L)
  latent_b <- make_rapid_adapter_latent(44L, 0.1)
  X_a <- make_rapid_adapter_features(latent_a, 11L, 6201L)
  X_b <- make_rapid_adapter_features(latent_b, 7L, 6202L)
  labels <- list(
    make_rapid_adapter_labels(48L)$labels,
    make_rapid_adapter_labels(44L)$labels
  )
  fit <- rapid_ma(
    list(north = X_a, south = X_b),
    labels = labels,
    positions = list(latent_a, latent_b),
    ncomp = 6L,
    control = list(max_iter = 2L, min_iter = 1L)
  )

  reapplied <- oos_predict(fit, X_a, side = "north", type = "all", k = 6L)
  expect_equal(reapplied$embedding, fit$scores$north, tolerance = 1e-12)
  expect_equal(
    reapplied$probabilities,
    fit$prediction_probabilities$north,
    tolerance = 1e-12
  )
  expect_identical(reapplied$class, fit$predictions$north)
  expect_true(all(reapplied$exact_reapplications))
  expect_identical(dim(reapplied$neighbours), c(48L, 6L))
  expect_equal(rowSums(reapplied$weights), rep(1, 48L), tolerance = 1e-12)

  perturbed <- X_a[1:5, , drop = FALSE] + 1e-3
  interpolated <- oos_predict(fit, perturbed, side = 1L, type = "all", k = 4L)
  expect_identical(dim(interpolated$embedding), c(5L, 6L))
  expect_identical(dim(interpolated$probabilities),
                   c(5L, ncol(fit$prototypes$class_prob)))
  expect_true(all(is.finite(interpolated$embedding)))
  expect_true(all(is.finite(interpolated$probabilities)))
  expect_equal(rowSums(interpolated$probabilities), rep(1, 5L), tolerance = 1e-12)
  expect_error(oos_predict(fit, X_a[, -1L], side = "north"), "expects 11")
})

test_that("OOS interpolation can retain position and attribute structure", {
  latent_a <- make_rapid_adapter_latent(48L)
  latent_b <- make_rapid_adapter_latent(44L, 0.1)
  X_a <- make_rapid_adapter_features(latent_a, 11L, 6251L)
  X_b <- make_rapid_adapter_features(latent_b, 7L, 6252L)
  attributes_a <- make_rapid_adapter_attributes(48L)
  attributes_b <- make_rapid_adapter_attributes(44L, 0.2)
  fit <- rapid_ma(
    list(north = X_a, south = X_b),
    labels = list(
      make_rapid_adapter_labels(48L)$labels,
      make_rapid_adapter_labels(44L)$labels
    ),
    positions = list(latent_a, latent_b),
    attributes = list(attributes_a, attributes_b),
    ncomp = 6L,
    control = list(max_iter = 2L, min_iter = 1L)
  )

  default <- oos_predict(fit, X_a[1:7, , drop = FALSE], side = "north",
                         type = "all", k = 6L)
  explicit_feature <- oos_predict(
    fit, X_a[1:7, , drop = FALSE], side = "north", type = "all", k = 6L,
    view_weights = c(feature = 1)
  )
  expect_equal(default$embedding, explicit_feature$embedding, tolerance = 0)
  expect_equal(default$probabilities, explicit_feature$probabilities,
               tolerance = 0)

  structured <- oos_predict(
    fit, X_a[1:7, , drop = FALSE], side = "north", type = "all", k = 6L,
    new_positions = latent_a[1:7, , drop = FALSE],
    new_attributes = attributes_a[1:7, , drop = FALSE],
    view_weights = c(feature = 1, position = 1, attribute = 1)
  )
  expect_equal(structured$embedding, fit$scores$north[1:7, , drop = FALSE],
               tolerance = 1e-12)
  expect_equal(
    structured$probabilities,
    fit$prediction_probabilities$north[1:7, , drop = FALSE],
    tolerance = 1e-12
  )
  expect_identical(structured$views_used,
                   c("feature", "position", "attribute"))
  expect_equal(unname(structured$view_weights), rep(1 / 3, 3),
               tolerance = 1e-12)
  expect_identical(vapply(structured$view_neighbours, nrow, integer(1)),
                   c(feature = 7L, position = 7L, attribute = 7L))

  repeated_x <- X_a[20L, , drop = FALSE] + 1e-3
  position_only <- oos_predict(
    fit, rbind(repeated_x, repeated_x), side = "north", type = "embedding",
    new_positions = latent_a[c(2L, 43L), , drop = FALSE],
    view_weights = c(position = 1), k = 4L
  )
  expect_false(isTRUE(all.equal(position_only[1L, ], position_only[2L, ])))

  reordered <- attributes_a[1:7, c("protected", "brightness", "landform")]
  expect_equal(
    oos_predict(
      fit, X_a[1:7, , drop = FALSE], side = "north", type = "probabilities",
      new_attributes = reordered, view_weights = c(attribute = 1)
    ),
    oos_predict(
      fit, X_a[1:7, , drop = FALSE], side = "north", type = "probabilities",
      new_attributes = attributes_a[1:7, , drop = FALSE],
      view_weights = c(attribute = 1)
    ),
    tolerance = 0
  )

  expect_error(
    oos_predict(fit, X_a[1:2, , drop = FALSE], side = "north",
                view_weights = c(position = 1)),
    "requires `new_positions`"
  )
  expect_error(
    oos_predict(
      fit, X_a[1:2, , drop = FALSE], side = "north",
      new_attributes = attributes_a[1:2, -1L, drop = FALSE],
      view_weights = c(attribute = 1)
    ),
    "fitted columns"
  )
  expect_error(
    oos_predict(fit, X_a[1:2, , drop = FALSE], side = "north",
                view_weights = c(feature = -1)),
    "nonnegative"
  )
})

test_that("bounded fine matching handles rectangles and fixes anchors", {
  source_latent <- make_rapid_adapter_latent(52L)
  target_latent <- make_rapid_adapter_latent(43L)
  source <- make_rapid_adapter_features(source_latent, 9L, 6301L)
  target <- make_rapid_adapter_features(target_latent, 6L, 6302L)
  source_labels <- make_rapid_adapter_labels(52L)$labels
  target_labels <- make_rapid_adapter_labels(43L)$labels
  anchors <- matrix(c(1L, 1L, 17L, 17L, 33L, 33L), ncol = 2L, byrow = TRUE)

  pair <- fit_pair(
    rapid_ma_aligner(), source, target,
    links = anchors,
    labels = list(source_labels, target_labels),
    positions = list(source_latent, target_latent),
    ncomp = 6L,
    control = list(max_iter = 2L, min_iter = 1L),
    match_control = list(
      candidate_cap = 12L,
      prototype_bucket_cap = 24L,
      structure_weight = 1,
      attribute_weight = 0
    )
  )
  matching <- pair$matching

  expect_s3_class(matching, "rapid_ma_matching")
  expect_s4_class(matching$assignment, "sparseMatrix")
  expect_identical(dim(matching$assignment), c(43L, 52L))
  expect_equal(matching$matches$target[anchors[, 1L]], anchors[, 2L])
  expect_true(all(matching$matches$anchored[anchors[, 1L]]))
  expect_true(all(matching$matches$confidence[anchors[, 1L]] == 1))
  expect_lte(matching$coverage, 43 / 52)
  expect_gt(matching$coverage, 0.7)
  expect_identical(
    length(matching$unmatched_source),
    52L - Matrix::nnzero(matching$assignment)
  )
  expect_lte(matching$candidate_edges, 52L * 12L)
  expect_false(matching$dense_pairwise_allocated)

  forward <- relative_transform(pair, from = "i", to = "j")
  reverse <- relative_transform(pair, from = "j", to = "i")
  expect_identical(dim(forward$op), c(43L, 52L))
  expect_identical(dim(reverse$op), c(52L, 43L))
  mapped <- apply_transform(forward, source)
  expect_identical(dim(mapped), c(43L, 9L))
})

test_that("fine matching accepts anchor-ID vectors and works across three domains", {
  latent <- make_rapid_adapter_latent(34L)
  domains <- list(
    a = make_rapid_adapter_features(latent, 7L, 6401L),
    b = make_rapid_adapter_features(latent, 8L, 6402L),
    c = make_rapid_adapter_features(latent, 5L, 6403L)
  )
  labels <- lapply(domains, function(x) make_rapid_adapter_labels(nrow(x))$labels)
  fit <- rapid_ma(
    domains,
    labels = labels,
    positions = list(latent, latent, latent),
    ncomp = 5L,
    control = list(max_iter = 1L, min_iter = 1L)
  )
  vec1 <- vec2 <- rep(NA_character_, 34L)
  vec1[c(2L, 18L, 30L)] <- vec2[c(2L, 18L, 30L)] <- c("u", "v", "w")
  matching <- rapid_ma_match(
    fit, from = "a", to = "c",
    anchors = list(vec1 = vec1, vec2 = vec2),
    candidate_cap = 16L,
    structure_weight = 1,
    attribute_weight = 0
  )

  expect_identical(matching$from, "a")
  expect_identical(matching$to, "c")
  expect_equal(matching$matches$target[c(2L, 18L, 30L)], c(2L, 18L, 30L))
  expect_equal(matching$anchors[, 1L], c(2L, 18L, 30L))
  expect_error(
    rapid_ma_match(fit, from = "a", to = "c", anchors = cbind(c(1L, 1L), c(2L, 3L))),
    "one-to-one"
  )
})
