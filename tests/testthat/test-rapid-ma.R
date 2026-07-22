rapid_ma_fit_list <- function(...) {
  manifoldalign:::rapid_ma.list(...)
}

rapid_control <- function(...) {
  manifoldalign::rapid_ma_control(...)
}

make_rapid_domain <- function(n, p, seed) {
  set.seed(seed)
  truth <- rep(c("forest", "water", "urban"), length.out = n)
  centers <- rbind(
    forest = c(-2, 0, 0),
    water = c(2, 0, 0),
    urban = c(0, 2.5, 0)
  )
  latent <- centers[truth, , drop = FALSE] + matrix(rnorm(n * 3L, sd = 0.35), n)
  loading <- matrix(rnorm(3L * p), 3L, p)
  X <- latent %*% loading + matrix(rnorm(n * p, sd = 0.15), n)
  position <- latent[, 1:2, drop = FALSE] + matrix(rnorm(n * 2L, sd = 0.1), n)
  attribute <- data.frame(
    elevation = latent[, 1L] + rnorm(n, sd = 0.1),
    moisture = latent[, 2L] + rnorm(n, sd = 0.1)
  )
  list(X = X, truth = truth, position = position, attribute = attribute)
}

sparse_labels <- function(truth, per_class = 5L) {
  labels <- rep(NA_character_, length(truth))
  for (class in unique(truth)) {
    index <- which(truth == class)
    labels[head(index, per_class)] <- class
  }
  labels
}

test_that("RAPID-MA fits rectangular structured domains with finite output", {
  d1 <- make_rapid_domain(72L, 14L, 1501L)
  d2 <- make_rapid_domain(57L, 9L, 1502L)
  labels <- list(sparse_labels(d1$truth, 6L), sparse_labels(d2$truth, 6L))
  fit <- rapid_ma_fit_list(
    list(source = d1$X, target = d2$X),
    labels = labels,
    positions = list(d1$position, d2$position),
    attributes = list(d1$attribute, d2$attribute),
    ncomp = 8L,
    control = list(
      max_iter = 3L,
      min_iter = 1L,
      supervised_weight = 1,
      uot = list(q = 8L, epsilon = 0.15, label_weight = 4)
    )
  )

  expect_s3_class(fit, "rapid_ma")
  expect_identical(vapply(fit$scores, nrow, integer(1)), c(source = 72L, target = 57L))
  expect_true(all(vapply(fit$scores, ncol, integer(1)) == 8L))
  expect_identical(dim(fit$s), c(129L, 8L))
  expect_true(all(is.finite(fit$s)))
  expect_true(all(vapply(fit$couplings, inherits, logical(1), "sparseMatrix")))
  expect_identical(dim(fit$couplings$source), c(72L, nrow(fit$prototypes$embedding)))
  expect_identical(dim(fit$couplings$target), c(57L, nrow(fit$prototypes$embedding)))
  expect_true(all(is.finite(fit$dropped_mass)))
  expect_true(fit$convergence$monotone_accepted)
  expect_gt(mean(fit$predictions$source == d1$truth), 0.8)
  expect_gt(mean(fit$predictions$target == d2$truth), 0.8)
})

test_that("accepted objective values are monotone and rollback state is explicit", {
  d1 <- make_rapid_domain(45L, 8L, 1503L)
  d2 <- make_rapid_domain(51L, 11L, 1504L)
  fit <- rapid_ma_fit_list(
    list(a = d1$X, b = d2$X),
    labels = list(sparse_labels(d1$truth), sparse_labels(d2$truth)),
    ncomp = 6L,
    control = list(max_iter = 3L, min_iter = 1L)
  )
  accepted <- fit$objective_history[fit$objective_history$accepted, ]

  expect_true(all(diff(accepted$objective) <= fit$control$objective_tolerance))
  expect_true(all(accepted$step[-1L] > 0 & accepted$step[-1L] <= 1))
  expect_true(fit$convergence$status %in% c("converged", "stalled", "max_iter"))
  expect_equal(tail(accepted$objective, 1L), fit$objective)
  expect_identical(
    fit$convergence$attempted_iterations,
    max(fit$objective_history$iteration)
  )
})

test_that("three-domain results are deterministic", {
  domains <- list(
    d1 = make_rapid_domain(34L, 7L, 1505L),
    d2 = make_rapid_domain(41L, 10L, 1506L),
    d3 = make_rapid_domain(29L, 6L, 1507L)
  )
  X <- lapply(domains, `[[`, "X")
  labels <- lapply(domains, function(x) sparse_labels(x$truth, 4L))
  control <- list(
    max_iter = 2L,
    min_iter = 1L,
    seed = 91L,
    uot = list(q = 8L, epsilon = 0.2)
  )

  first <- rapid_ma_fit_list(X, labels = labels, ncomp = 5L, control = control)
  second <- rapid_ma_fit_list(X, labels = labels, ncomp = 5L, control = control)

  expect_equal(second$s, first$s, tolerance = 1e-12)
  expect_equal(second$objective_history, first$objective_history, tolerance = 1e-12)
  expect_identical(second$predictions, first$predictions)
  for (name in names(first$couplings)) {
    expect_equal(second$couplings[[name]], first$couplings[[name]], tolerance = 1e-12)
  }
  expect_identical(first$domain_sizes, c(d1 = 34L, d2 = 41L, d3 = 29L))
})

test_that("fully unlabeled mode remains a structural embedding", {
  d1 <- make_rapid_domain(38L, 6L, 1508L)
  d2 <- make_rapid_domain(43L, 6L, 1509L)
  fit <- rapid_ma_fit_list(
    list(d1 = d1$X, d2 = d2$X),
    positions = list(d1$position, d2$position),
    ncomp = 5L,
    control = list(
      max_iter = 1L,
      min_iter = 1L,
      prototype = list(unknown_prototypes = 6L)
    )
  )

  expect_true(fit$prototypes$diagnostics$unknown_only)
  expect_true(all(fit$predictions$d1 == ".unknown"))
  expect_true(all(fit$predictions$d2 == ".unknown"))
  expect_equal(fit$objective_components[["supervised"]], 0)
  expect_true(all(is.finite(fit$s)))
  expect_true(length(fit$relation_barycenters) >= 1L)
})

test_that("hyperdesign label columns and custom relations use core semantics", {
  d1 <- make_rapid_domain(30L, 6L, 1510L)
  d2 <- make_rapid_domain(27L, 6L, 1511L)
  labels <- list(sparse_labels(d1$truth, 4L), sparse_labels(d2$truth, 4L))
  hd <- as_hyperdesign(
    list(d1$X, d2$X), labels = labels, label_name = "cover",
    domain_names = c("north", "south")
  )
  cycle <- function(n) Matrix::sparseMatrix(
    i = seq_len(n), j = c(2:n, 1L), x = 1, dims = c(n, n)
  )
  fit <- manifoldalign:::rapid_ma.hyperdesign(
    hd,
    y = cover,
    relations = list(
      north = list(topology = cycle(30L)),
      south = list(topology = cycle(27L))
    ),
    ncomp = 4L,
    control = list(max_iter = 1L, min_iter = 1L)
  )

  expect_identical(fit$domain_names, c("north", "south"))
  expect_identical(fit$labels, labels)
  expect_true(all(vapply(fit$relations, function(x) {
    "topology" %in% names(x$relations)
  }, logical(1))))
  expect_identical(dim(fit$s), c(57L, 4L))
})

test_that("RAPID-MA does not masquerade as a linear projector", {
  d1 <- make_rapid_domain(25L, 5L, 1512L)
  d2 <- make_rapid_domain(22L, 8L, 1513L)
  fit <- rapid_ma_fit_list(
    list(d1$X, d2$X),
    labels = list(sparse_labels(d1$truth, 3L), sparse_labels(d2$truth, 3L)),
    ncomp = 4L,
    control = list(max_iter = 1L, min_iter = 1L)
  )

  expect_false(inherits(fit, "multiblock_biprojector"))
  expect_null(fit$v)
  expect_null(fit$loadings)
  expect_true(fit$preprocessing$embedding_only)
  expect_false(fit$preprocessing$loadings_available)
  expect_identical(nrow(fit$domain_maps[[1L]]$coef), 4L)
})

test_that("domain failures identify the offending input", {
  X1 <- matrix(rnorm(60L), nrow = 20L)
  X2 <- matrix(rnorm(54L), nrow = 18L)
  bad_relation <- Matrix::Diagonal(17L)

  expect_error(
    rapid_ma_fit_list(
      list(good = X1, broken = X2),
      relations = list(NULL, list(topology = bad_relation)),
      control = list(max_iter = 1L)
    ),
    "Domain `broken` relation preparation"
  )
  X2[3L, 2L] <- Inf
  expect_error(
    rapid_ma_fit_list(list(good = X1, nonfinite = X2)),
    "Domain `nonfinite`"
  )
})

test_that("RAPID-MA control validates nested and outer settings", {
  control <- rapid_control(
    max_iter = 4L,
    map_type = "orthogonal",
    relation = list(feature_k = 9L, degree_cap = 12L),
    diffusion = list(steps = c(0L, 1L, 3L)),
    prototype = list(prototypes_per_class = 3L),
    uot = list(q = 5L),
    relmatch = list(relation_gate = "uniform")
  )

  expect_s3_class(control, "rapid_ma_control")
  expect_identical(control$relation$feature_k, 9L)
  expect_identical(control$diffusion$steps, c(0L, 1L, 3L))
  expect_identical(control$prototype$prototypes_per_class, 3L)
  expect_identical(control$uot$q, 5L)
  expect_identical(control$relmatch$relation_gate, "uniform")
  expect_error(rapid_control(backtrack_factor = 1), "must be < 1")
  expect_error(rapid_control(map_step = 1.1), "must be <= 1")
  expect_error(
    manifoldalign:::.rapid_resolve_control(list(not_a_setting = 1)),
    "Unknown RAPID-MA control"
  )
})
