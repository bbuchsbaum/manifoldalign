make_feature_domain <- function(n, p, seed) {
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
  colnames(X) <- paste0("map", seq_len(p))
  X
}

feature_fixture <- function(seed = 2101L, p = 6L) {
  X <- list(
    north = make_feature_domain(36L, p, seed),
    south = make_feature_domain(31L, p, seed + 1L)
  )
  fit <- rapid_ma(
    X,
    ncomp = 3L,
    control = rapid_ma_control(
      max_iter = 1L,
      min_iter = 1L,
      prototype = list(unknown_prototypes = 4L)
    )
  )
  list(X = X, fit = fit)
}

oracle_barycentre <- function(X, Z, weight = c("positive", "absolute", "shifted")) {
  weight <- match.arg(weight)
  W <- switch(
    weight,
    positive = pmax(X, 0),
    absolute = abs(X),
    shifted = sweep(X, 2L, apply(X, 2L, min), "-")
  )
  out <- matrix(NA_real_, ncol(X), ncol(Z))
  for (j in seq_len(ncol(X))) {
    total <- sum(W[, j])
    if (total > 0) {
      out[j, ] <- colSums(W[, j] * Z) / total
    }
  }
  out
}

test_that("feature scores match an explicit barycentre oracle", {
  fx <- feature_fixture()
  dual <- rapid_ma_feature_scores(fx$fit, fx$X)

  expect_s3_class(dual, "rapid_ma_feature_scores")
  expect_identical(dual$domain_names, c("north", "south"))
  expect_identical(dual$weight, "positive")
  for (name in names(fx$X)) {
    expected <- oracle_barycentre(fx$X[[name]], fx$fit$scores[[name]])
    expect_equal(unname(dual$scores[[name]]), expected, tolerance = 1e-12)
    expect_identical(rownames(dual$scores[[name]]), paste0("map", 1:6))
    expect_identical(colnames(dual$scores[[name]]),
                     colnames(fx$fit$scores[[name]]))
    expect_equal(
      unname(dual$mass[[name]]),
      unname(colSums(pmax(fx$X[[name]], 0))),
      tolerance = 1e-12
    )
    expect_identical(dual$zero_mass[[name]], integer(0))
  }

  stacked <- (dual$scores$north + dual$scores$south) / 2
  expect_equal(dual$pooled, stacked, tolerance = 1e-12)
  expected_dispersion <- sqrt((rowSums((dual$scores$north - stacked)^2) +
                                 rowSums((dual$scores$south - stacked)^2)) / 2)
  expect_equal(unname(dual$dispersion), unname(expected_dispersion),
               tolerance = 1e-12)
  expect_identical(unname(dual$n_domains), rep(2L, 6L))
  expect_output(print(dual), "RAPID-MA feature-side coordinates")
})

test_that("weight rules follow their definitions", {
  fx <- feature_fixture()
  for (weight in c("absolute", "shifted")) {
    dual <- rapid_ma_feature_scores(fx$fit, fx$X, weight = weight)
    for (name in names(fx$X)) {
      expected <- oracle_barycentre(fx$X[[name]], fx$fit$scores[[name]], weight)
      expect_equal(unname(dual$scores[[name]]), expected, tolerance = 1e-12)
    }
  }
  positive <- rapid_ma_feature_scores(fx$fit, fx$X, weight = "positive")
  absolute <- rapid_ma_feature_scores(fx$fit, fx$X, weight = "absolute")
  expect_false(isTRUE(all.equal(positive$scores$north, absolute$scores$north)))
})

test_that("indicator and scaled maps land where the barycentre says", {
  fx <- feature_fixture()
  Z <- fx$fit$scores$north
  n <- nrow(Z)
  indicator <- matrix(0, n, 3L)
  indicator[5L, 1L] <- 1
  indicator[c(2L, 9L), 2L] <- c(1, 3)
  indicator[, 3L] <- 7 * pmax(fx$X$north[, 1L], 0)
  colnames(indicator) <- c("single", "pair", "scaled")

  dual <- rapid_ma_feature_scores(fx$fit, list(north = indicator), pool = FALSE)
  expect_identical(names(dual$scores), "north")
  expect_null(dual$pooled)
  expect_equal(unname(dual$scores$north["single", ]), unname(Z[5L, ]),
               tolerance = 1e-12)
  expect_equal(
    unname(dual$scores$north["pair", ]),
    unname((Z[2L, ] + 3 * Z[9L, ]) / 4),
    tolerance = 1e-12
  )
  base <- rapid_ma_feature_scores(fx$fit, fx$X, pool = FALSE)
  expect_equal(unname(dual$scores$north["scaled", ]),
               unname(base$scores$north["map1", ]), tolerance = 1e-12)
})

test_that("sparse Matrix input matches dense input and stays sparse", {
  fx <- feature_fixture()
  sparse <- lapply(fx$X, function(X) {
    X[abs(X) < 0.5] <- 0
    Matrix::Matrix(X, sparse = TRUE)
  })
  dense <- lapply(sparse, as.matrix)
  for (weight in c("positive", "absolute", "shifted")) {
    a <- suppressWarnings(rapid_ma_feature_scores(fx$fit, sparse, weight = weight))
    b <- suppressWarnings(rapid_ma_feature_scores(fx$fit, dense, weight = weight))
    expect_equal(a$scores, b$scores, tolerance = 1e-12)
    expect_equal(a$pooled, b$pooled, tolerance = 1e-12)
    expect_identical(a$zero_mass, b$zero_mass)
    expect_true(all(is.finite(a$pooled[a$n_domains > 0L, ])))
  }
  W <- manifoldalign:::.rapid_feature_weights(sparse$north, "absolute")
  expect_s4_class(W, "dgCMatrix")
  expect_equal(as.matrix(W), abs(dense$north), tolerance = 1e-12,
               ignore_attr = TRUE)
})

test_that("zero-mass features are NA, reported, and excluded from pooling", {
  fx <- feature_fixture()
  X <- fx$X
  X$north[, 2L] <- 0
  X$north[, 4L] <- -abs(X$north[, 4L])
  X$south[, 4L] <- 0

  expect_warning(
    dual <- rapid_ma_feature_scores(fx$fit, X),
    "zero total weight.*north \\(2\\), south \\(1\\)"
  )
  expect_identical(dual$zero_mass$north, c(2L, 4L))
  expect_identical(dual$zero_mass$south, 4L)
  expect_true(all(is.na(dual$scores$north[c(2L, 4L), ])))
  expect_true(all(is.na(dual$scores$south[4L, ])))
  expect_equal(dual$pooled["map2", ], dual$scores$south["map2", ],
               tolerance = 1e-12)
  expect_true(all(is.na(dual$pooled["map4", ])))
  expect_identical(unname(dual$n_domains), c(2L, 1L, 2L, 0L, 2L, 2L))
  expect_true(is.na(dual$dispersion[["map2"]]))
  expect_true(is.na(dual$dispersion[["map4"]]))
  expect_false(is.na(dual$dispersion[["map1"]]))
})

test_that("hyperdesign and pair-fit inputs resolve domains by name", {
  fx <- feature_fixture()
  hd <- as_hyperdesign(unname(fx$X), domain_names = names(fx$X))
  from_hd <- rapid_ma_feature_scores(fx$fit, hd)
  from_list <- rapid_ma_feature_scores(fx$fit, fx$X)
  expect_equal(from_hd$scores, from_list$scores, tolerance = 1e-12)

  reversed <- rapid_ma_feature_scores(fx$fit, rev(fx$X))
  expect_identical(names(reversed$scores), c("north", "south"))
  expect_equal(reversed$scores, from_list$scores, tolerance = 1e-12)

  pair <- fit_pair(
    rapid_ma_aligner(),
    fx$X$north, fx$X$south,
    ncomp = 3L,
    control = rapid_ma_control(
      max_iter = 1L, min_iter = 1L,
      prototype = list(unknown_prototypes = 4L)
    ),
    fine_match = FALSE
  )
  dual <- rapid_ma_feature_scores(pair, list(fx$X$north, fx$X$south))
  expect_identical(dual$domain_names, c("i", "j"))
  expect_equal(unname(dual$scores$i),
               oracle_barycentre(fx$X$north, pair$core$scores$i),
               tolerance = 1e-12)
})

test_that("input contract violations are rejected explicitly", {
  fx <- feature_fixture()
  expect_error(rapid_ma_feature_scores(list(), fx$X), "must be a `rapid_ma`")
  expect_error(rapid_ma_feature_scores(fx$fit, fx$X$north), "list with one matrix")
  expect_error(rapid_ma_feature_scores(fx$fit, list(fx$X$north)),
               "one entry per fitted domain")
  expect_error(rapid_ma_feature_scores(fx$fit, list(east = fx$X$north)),
               "do not match fitted domains")
  expect_error(rapid_ma_feature_scores(fx$fit, list(north = NULL, south = NULL)),
               "at least one domain")
  expect_error(
    rapid_ma_feature_scores(fx$fit, list(north = fx$X$north[-1L, ])),
    "Domain `north`: `data` must have 36 rows"
  )
  expect_error(
    rapid_ma_feature_scores(fx$fit, list(north = as.data.frame(fx$X$north))),
    "must be a numeric matrix"
  )
  bad <- fx$X
  bad$south[2L, 3L] <- NaN
  expect_error(rapid_ma_feature_scores(fx$fit, bad), "finite")
  expect_error(rapid_ma_feature_scores(fx$fit, fx$X, pool = NA), "`pool`")

  uneven <- fx$X
  uneven$south <- uneven$south[, 1:4]
  expect_error(rapid_ma_feature_scores(fx$fit, uneven), "same number of features")
  ok <- rapid_ma_feature_scores(fx$fit, uneven, pool = FALSE)
  expect_identical(nrow(ok$scores$south), 4L)

  renamed <- fx$X
  colnames(renamed$south) <- paste0("other", 1:6)
  expect_error(rapid_ma_feature_scores(fx$fit, renamed), "identical feature names")
})
