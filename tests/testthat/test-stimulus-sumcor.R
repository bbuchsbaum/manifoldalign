ssc_fun <- function(name) get(name, envir = environment(stimulus_sumcor_align))

ssc_simulate <- function(n_events = 40L, n_rep = 2L, k = 3L,
                         dims = c(a = 20L, b = 24L, c = 18L),
                         noise = 0.4, seed = 1L, n_rest = 0L) {
  set.seed(seed)
  latent <- matrix(rnorm(n_events * k), n_events, k)
  event_ids <- sprintf("ev%03d", seq_len(n_events))
  subjects <- lapply(names(dims), function(nm) {
    p <- dims[[nm]]
    mix <- matrix(rnorm(k * p), k, p)
    rows <- rep(seq_len(n_events), times = n_rep)
    X <- latent[rows, , drop = FALSE] %*% mix +
      matrix(rnorm(length(rows) * p, sd = noise), length(rows), p)
    design <- data.frame(
      event_id = event_ids[rows],
      repeat_id = rep(seq_len(n_rep), each = n_events),
      run_id = paste0("run", rep(seq_len(n_rep), each = n_events)),
      obs = paste0(nm, "_", seq_along(rows)),
      stringsAsFactors = FALSE
    )
    if (n_rest > 0L) {
      X <- rbind(X, matrix(rnorm(n_rest * p), n_rest, p))
      design <- rbind(design, data.frame(
        event_id = NA_character_, repeat_id = NA_integer_, run_id = "run1",
        obs = paste0(nm, "_rest", seq_len(n_rest)), stringsAsFactors = FALSE
      ))
    }
    colnames(X) <- paste0(nm, "_v", seq_len(p))
    list(X = X, design = design)
  })
  names(subjects) <- names(dims)
  features <- latent %*% matrix(rnorm(k * 6L), k, 6L) +
    matrix(rnorm(n_events * 6L, sd = 0.2), n_events, 6L)
  rownames(features) <- event_ids
  list(
    X = lapply(subjects, `[[`, "X"),
    design = lapply(subjects, `[[`, "design"),
    latent = latent,
    event_ids = event_ids,
    features = features
  )
}

ssc_fit <- function(sim, ..., control = list(), ncomp = 3L) {
  base <- list(beta = 0, brain_rank = 12L, tol = 1e-13, grad_tol = 1e-10,
               max_iter = 500L)
  stimulus_sumcor_align(
    sim$X, design = sim$design, event_id = "event_id",
    repeat_id = "repeat_id", run_id = "run_id", observation_id = "obs",
    ncomp = ncomp, control = utils::modifyList(base, control), ...
  )
}

# Independent oracle for the identity basis: weighted centring, presentation
# covariance with shrinkage, overlap-normalised cross products.
ssc_identity_oracle <- function(X, design, rho, event_ids) {
  ev <- design$event_id
  keep <- !is.na(ev)
  X <- X[keep, , drop = FALSE]
  ev <- ev[keep]
  w <- rep(1 / length(event_ids), length(event_ids))
  names(w) <- event_ids
  counts <- table(ev)
  observed <- names(counts)
  m_s <- sum(w[observed])
  omega <- w[ev] / (as.numeric(counts[ev]) * m_s)
  mu <- colSums(X * omega)
  Xc <- sweep(X, 2L, mu)
  Chat <- crossprod(Xc * sqrt(omega))
  p <- ncol(X)
  C <- (1 - rho) * Chat + rho * sum(diag(Chat)) / p * diag(p)
  Ubar <- matrix(0, length(event_ids), p)
  for (e in observed) {
    Ubar[match(e, event_ids), ] <- colMeans(Xc[ev == e, , drop = FALSE])
  }
  list(C = C, Ubar = Ubar, w = w)
}

ssc_inv_sqrt <- function(C) {
  e <- eigen(C, symmetric = TRUE)
  e$vectors %*% (t(e$vectors) / sqrt(e$values))
}

test_that("two-subject fit reaches the exact CCA optimum and its constraints", {
  sim <- ssc_simulate(dims = c(a = 6L, b = 7L), seed = 11L)
  rho <- 0.1
  fit <- ssc_fit(sim, control = list(brain_basis = "identity",
                                     covariance_shrinkage = rho))
  o1 <- ssc_identity_oracle(sim$X$a, sim$design$a, rho, sim$event_ids)
  o2 <- ssc_identity_oracle(sim$X$b, sim$design$b, rho, sim$event_ids)
  expect_equal(fit$encoders$a$C, o1$C, tolerance = 1e-10, ignore_attr = TRUE)

  B12 <- crossprod(o1$Ubar * sqrt(o1$w), o2$Ubar * sqrt(o2$w))
  M <- ssc_inv_sqrt(o1$C) %*% B12 %*% ssc_inv_sqrt(o2$C)
  sv <- svd(M)
  expect_equal(fit$convergence$objective, sum(sv$d[1:3]), tolerance = 1e-8)
  expect_equal(fit$component_agreement, sv$d[1:3], tolerance = 1e-8)

  A_oracle <- ssc_inv_sqrt(o1$C) %*% sv$u[, 1:3]
  cross <- crossprod(fit$encoders$a$A, o1$C %*% A_oracle)
  expect_equal(abs(diag(cross)), rep(1, 3), tolerance = 1e-6)

  for (nm in c("a", "b")) {
    enc <- fit$encoders[[nm]]
    expect_equal(crossprod(enc$A, enc$C %*% enc$A), diag(3), tolerance = 1e-10)
  }
  expect_true(fit$convergence$converged)
  expect_true(fit$convergence$monotone)
})

test_that("the objective history is monotone and diagnostics are recorded", {
  sim <- ssc_simulate(seed = 12L)
  fit <- ssc_fit(sim)
  hist <- fit$convergence$history
  expect_true(all(diff(hist$objective) >= -1e-10 * max(1, abs(hist$objective))))
  expect_true(all(hist$max_constraint_error < 1e-10))
  expect_lt(utils::tail(hist$max_stationarity, 1L), 1e-10)
  expect_identical(unname(fit$supported_rank), rep(3L, 3L))
  expect_true(all(diff(fit$component_agreement) <= 1e-12))
  expect_output(print(fit), "Stimulus-weighted SUMCOR alignment")
})

test_that("training scores are the encoder applied to the training rows", {
  sim <- ssc_simulate(seed = 13L, n_rest = 4L)
  fit <- ssc_fit(sim)
  for (nm in names(sim$X)) {
    pred <- oos_predict(fit, sim$X[[nm]], side = nm)
    expect_identical(unname(pred), unname(fit$scores[[nm]]))
    expect_identical(rownames(fit$scores[[nm]]), sim$design[[nm]]$obs)
  }
  junk <- matrix(rnorm(3 * ncol(sim$X$a), sd = 50), 3, ncol(sim$X$a))
  colnames(junk) <- colnames(sim$X$a)
  alone <- oos_predict(fit, sim$X$a[1:5, ], side = "a")
  mixed <- oos_predict(fit, rbind(sim$X$a[1:5, ], junk), side = 1L)
  expect_equal(unname(mixed[1:5, ]), unname(alone), tolerance = 1e-12)
  expect_error(oos_predict(fit, sim$X$a[, -1], side = "a"), "columns")
  renamed <- sim$X$a
  colnames(renamed)[1] <- "other"
  expect_error(oos_predict(fit, renamed, side = "a"), "column names")
  expect_error(oos_predict(fit, sim$X$a, side = "zz"), "Unknown domain")
})

test_that("observations without an event are scored but do not change the fit", {
  plain <- ssc_simulate(seed = 14L)
  with_rest <- plain
  set.seed(140)
  for (nm in names(plain$X)) {
    rest <- matrix(rnorm(6L * ncol(plain$X[[nm]]), sd = 3), 6L)
    colnames(rest) <- colnames(plain$X[[nm]])
    with_rest$X[[nm]] <- rbind(plain$X[[nm]][1:10, ], rest, plain$X[[nm]][-(1:10), ])
    rest_design <- data.frame(event_id = NA_character_, repeat_id = NA_integer_,
                              run_id = "run1", obs = paste0(nm, "_rest", 1:6))
    with_rest$design[[nm]] <- rbind(plain$design[[nm]][1:10, ], rest_design,
                                    plain$design[[nm]][-(1:10), ])
  }
  fit_rest <- ssc_fit(with_rest)
  fit_plain <- ssc_fit(plain)
  for (nm in names(plain$X)) {
    expect_equal(fit_rest$encoders[[nm]]$weights, fit_plain$encoders[[nm]]$weights,
                 tolerance = 1e-10)
    expect_identical(nrow(fit_rest$scores[[nm]]), nrow(with_rest$X[[nm]]))
  }
  expect_identical(unname(fit_rest$provenance$observations_unlabelled),
                   rep(6L, 3L))
})

test_that("row order within a subject does not change the encoder", {
  sim <- ssc_simulate(seed = 15L)
  fit <- ssc_fit(sim)
  perm <- sample(nrow(sim$X$b))
  shuffled <- sim
  shuffled$X$b <- sim$X$b[perm, ]
  shuffled$design$b <- sim$design$b[perm, ]
  fit_shuffled <- ssc_fit(shuffled)
  for (nm in names(sim$X)) {
    expect_equal(fit_shuffled$encoders[[nm]]$weights, fit$encoders[[nm]]$weights,
                 tolerance = 1e-8)
  }
  expect_equal(fit_shuffled$scores$b[sim$design$b$obs, ], fit$scores$b,
               tolerance = 1e-8)
})

test_that("subject order does not change the per-subject encoders", {
  sim <- ssc_simulate(seed = 16L)
  fit <- ssc_fit(sim)
  reversed <- sim
  reversed$X <- rev(sim$X)
  reversed$design <- rev(sim$design)
  fit_rev <- ssc_fit(reversed)
  for (nm in names(sim$X)) {
    expect_equal(fit_rev$encoders[[nm]]$weights, fit$encoders[[nm]]$weights,
                 tolerance = 1e-6)
  }
})

test_that("logging does not change the computation", {
  sim <- ssc_simulate(seed = 17L)
  quiet <- ssc_fit(sim)
  loud <- suppressMessages(ssc_fit(sim, control = list(verbose = TRUE)))
  msgs <- suppressWarnings(capture_messages(
    ssc_fit(sim, control = list(verbose = TRUE, max_iter = 2L, tol = 0,
                                grad_tol = 0))
  ))
  expect_true(any(grepl("sweep 1", msgs)))
  for (nm in names(sim$X)) {
    expect_identical(loud$encoders[[nm]]$weights, quiet$encoders[[nm]]$weights)
  }
  expect_identical(loud$convergence$history, quiet$convergence$history)
})

test_that("an exhausted iteration budget is reported as non-convergence", {
  sim <- ssc_simulate(seed = 18L)
  expect_warning(
    fit <- ssc_fit(sim, control = list(max_iter = 1L, tol = 0, grad_tol = 0)),
    "did not converge"
  )
  expect_false(fit$convergence$converged)
  expect_identical(fit$convergence$status, "max_iter")
  expect_identical(fit$convergence$iterations, 1L)
})

test_that("the shared signal is recovered", {
  sim <- ssc_simulate(seed = 19L, noise = 0.3)
  fit <- ssc_fit(sim)
  cc <- cancor(fit$reference_scores, sim$latent)$cor
  expect_gt(min(cc), 0.95)
  for (nm in names(sim$X)) {
    cc_s <- cancor(fit$event_scores[[nm]], sim$latent)$cor
    expect_gt(min(cc_s), 0.9)
  }
})

test_that("beta = 0 is exactly the correspondence-only model", {
  sim <- ssc_simulate(seed = 20L)
  none <- ssc_fit(sim)
  ignored <- ssc_fit(sim, stimulus_features = list(acoustic = sim$features))
  expect_null(ignored$kernel$factor)
  for (nm in names(sim$X)) {
    expect_identical(ignored$encoders[[nm]]$weights, none$encoders[[nm]]$weights)
  }
  weighted <- ssc_fit(sim, stimulus_features = list(acoustic = sim$features),
                      control = list(beta = 1))
  expect_false(is.null(weighted$kernel$factor))
  expect_false(isTRUE(all.equal(weighted$encoders$a$weights,
                                none$encoders$a$weights)))
  expect_identical(weighted$kernel$families$family, "acoustic")
})

test_that("the stimulus kernel factor matches the defining formula", {
  sim <- ssc_simulate(seed = 21L)
  modality <- rep(c("music", "video"), each = 20L)
  events <- data.frame(event_id = sim$event_ids, modality = modality)
  F <- sim$features
  beta <- 0.7
  ridge <- 0.5
  fit <- ssc_fit(sim, events = events, stimulus_features = list(f = F),
                 control = list(beta = beta, feature_ridge = ridge))
  L <- fit$kernel$factor
  Fc <- sweep(F, 2L, colMeans(F))
  Fc <- sweep(Fc, 2L, apply(Fc, 2L, sd), "/")
  K <- matrix(0, nrow(F), nrow(F))
  kappa_max <- 0
  for (mod in c("music", "video")) {
    idx <- which(modality == mod)
    Fb <- Fc[idx, ]
    eta <- ridge * mean(svd(Fb)$d^2)
    Kb <- Fb %*% solve(crossprod(Fb) + eta * diag(ncol(Fb)), t(Fb))
    K[idx, idx] <- Kb
    kappa_max <- max(kappa_max, max(eigen(Kb, symmetric = TRUE)$values))
  }
  expect_equal(L %*% t(L), beta * K / kappa_max, tolerance = 1e-8)
  expect_true(all(abs((L %*% t(L))[modality == "music", modality == "video"]) < 1e-14))

  cross <- ssc_fit(sim, events = events, stimulus_features = list(f = F),
                   control = list(beta = beta, cross_modal = TRUE))
  expect_gt(max(abs((cross$kernel$factor %*% t(cross$kernel$factor))[
    modality == "music", modality == "video"])), 1e-3)
})

test_that("a kernel close to a multiple of the identity is flagged", {
  sim <- ssc_simulate(seed = 22L)
  F <- diag(length(sim$event_ids))
  rownames(F) <- sim$event_ids
  expect_warning(
    ssc_fit(sim, stimulus_features = list(identity = F),
            control = list(beta = 1)),
    "close to a multiple of the identity"
  )
})

test_that("missing events are masked by key and overlap is normalised", {
  sim <- ssc_simulate(seed = 23L)
  drop <- sim$design$c$event_id %in% sim$event_ids[31:40]
  sim$X$c <- sim$X$c[!drop, ]
  sim$design$c <- sim$design$c[!drop, ]
  fit <- ssc_fit(sim)
  expect_equal(fit$overlap["a", "c"], 0.75)
  expect_equal(fit$pair_weights["a", "c"], 1 / 0.75)
  expect_true(all(is.na(fit$event_scores$c[sim$event_ids[31:40], ])))
  expect_identical(fit$events$n_domains[31:40], rep(2L, 10L))
  expect_true(fit$convergence$converged)

  sparse <- ssc_simulate(seed = 23L)
  keep <- sparse$design$c$event_id %in% sim$event_ids[1:8]
  sparse$X$c <- sparse$X$c[keep, ]
  sparse$design$c <- sparse$design$c[keep, ]
  expect_error(ssc_fit(sparse, control = list(brain_rank = 5L)), "min_overlap")
})

test_that("event weights follow the declared rule", {
  sim <- ssc_simulate(seed = 24L)
  events <- data.frame(event_id = sim$event_ids,
                       modality = rep(c("music", "video"), c(10L, 30L)))
  fit <- ssc_fit(sim, events = events, control = list(event_weights = "modality"))
  w <- fit$events$weight
  expect_equal(sum(w[fit$events$modality == "music"]), 0.5)
  expect_equal(sum(w[fit$events$modality == "video"]), 0.5)
  custom <- stats::setNames(rep(1, 40L), sim$event_ids)
  custom[1:5] <- 3
  fit_custom <- ssc_fit(sim, control = list(event_weights = custom))
  expect_equal(fit_custom$events$weight, custom / sum(custom), ignore_attr = TRUE)
  expect_error(ssc_fit(sim, control = list(event_weights = custom[-1])),
               "lacks")
})

test_that("calibrating a training subject against the others reproduces it", {
  sim <- ssc_simulate(seed = 25L)
  fit <- ssc_fit(sim)
  cal <- calibrate_domain(fit, sim$X$a, design = sim$design$a,
                          domain_name = "a_again",
                          reference_domains = c("b", "c"))
  expect_equal(cal$encoders$a_again$weights, fit$encoders$a$weights,
               tolerance = 1e-6)
})

test_that("calibration adds a subject without touching the frozen reference", {
  full <- ssc_simulate(seed = 26L, dims = c(a = 20L, b = 24L, c = 18L, d = 22L))
  train <- full
  train$X <- full$X[c("a", "b", "c")]
  train$design <- full$design[c("a", "b", "c")]
  fit <- ssc_fit(train)
  cal <- calibrate_domain(fit, full$X$d, design = full$design$d,
                          domain_name = "d")
  for (nm in c("a", "b", "c")) {
    expect_identical(cal$encoders[[nm]], fit$encoders[[nm]])
    expect_identical(cal$scores[[nm]], fit$scores[[nm]])
  }
  expect_identical(cal$reference_scores, fit$reference_scores)
  expect_identical(cal$domain_roles[["d"]], "calibrated")
  expect_identical(cal$calibration$d$reference_domains, c("a", "b", "c"))
  pred <- oos_predict(cal, full$X$d, side = "d")
  expect_identical(unname(pred), unname(cal$scores$d))
  agreement <- diag(cor(cal$event_scores$d, cal$reference_scores))
  expect_gt(min(agreement), 0.9)
  expect_error(calibrate_domain(cal, full$X$d, design = full$design$d,
                                domain_name = "d"), "already exists")
  odd <- full$design$d
  odd$event_id[1:4] <- "novel"
  odd$repeat_id[1:4] <- 1:4
  expect_warning(calibrate_domain(fit, full$X$d, design = odd,
                                  domain_name = "d2"), "not in the fitted")
})

test_that("unsupported latent directions are completed and reported", {
  sim <- ssc_simulate(n_events = 3L, n_rep = 3L, seed = 27L,
                      dims = c(a = 8L, b = 9L))
  expect_warning(
    fit <- ssc_fit(sim, ncomp = 4L,
                   control = list(brain_rank = 6L, max_iter = 50L)),
    "supported"
  )
  expect_true(all(fit$supported_rank <= 3L))
  for (nm in c("a", "b")) {
    enc <- fit$encoders[[nm]]
    expect_equal(crossprod(enc$A, enc$C %*% enc$A), diag(4), tolerance = 1e-10)
  }
})

test_that("parcel PCA builds a block basis per parcel", {
  sim <- ssc_simulate(seed = 28L)
  parcels <- lapply(sim$X, function(X) rep(1:3, length.out = ncol(X)))
  fit <- ssc_fit(sim, parcels = parcels,
                 control = list(brain_basis = "parcel_pca", parcel_rank = 3L))
  P <- fit$encoders$a$basis$P
  expect_s4_class(P, "dgCMatrix")
  expect_identical(fit$encoders$a$basis$rank, 9L)
  for (j in seq_len(ncol(P))) {
    rows <- which(P[, j] != 0)
    expect_length(unique(parcels$a[rows]), 1L)
  }
  expect_error(ssc_fit(sim, control = list(brain_basis = "parcel_pca")),
               "requires `parcels`")
})

test_that("hyperdesign and list inputs give the same fit", {
  sim <- ssc_simulate(seed = 29L)
  hd <- structure(
    lapply(names(sim$X), function(nm) list(x = sim$X[[nm]],
                                           design = sim$design[[nm]])),
    names = names(sim$X),
    class = "hyperdesign"
  )
  fit_hd <- stimulus_sumcor_align(
    hd, event_id = "event_id", repeat_id = "repeat_id", observation_id = "obs",
    ncomp = 3L,
    control = list(beta = 0, brain_rank = 12L, tol = 1e-13, grad_tol = 1e-10,
                   max_iter = 500L)
  )
  fit_list <- ssc_fit(sim)
  for (nm in names(sim$X)) {
    expect_equal(fit_hd$encoders[[nm]]$weights, fit_list$encoders[[nm]]$weights,
                 tolerance = 1e-12)
  }
})

test_that("identifier and control contracts are enforced", {
  sim <- ssc_simulate(seed = 30L)
  expect_error(ssc_fit(sim, bogus = 1), "Unknown argument")
  expect_error(ssc_fit(sim, control = list(not_a_field = 1)), "Unknown")
  ctrl <- stimulus_sumcor_control()
  ctrl$beta <- -1
  expect_error(stimulus_sumcor_align(sim$X, design = sim$design, ncomp = 2L,
                                     control = ctrl), "beta")
  expect_error(ssc_fit(sim, ncomp = 30L), "exceeds the brain rank")
  dup <- sim
  dup$design$a$repeat_id[2] <- dup$design$a$repeat_id[1]
  dup$design$a$event_id[2] <- dup$design$a$event_id[1]
  expect_error(ssc_fit(dup), "must be unique")
  clash <- sim
  clash$design$a$modality <- "music"
  clash$design$b$modality <- "video"
  clash$design$c$modality <- "music"
  expect_error(ssc_fit(clash, modality = "modality"), "inconsistent")
  bad <- sim
  bad$X$a[1, 1] <- NA
  expect_error(ssc_fit(bad), "finite")
  expect_error(stimulus_sumcor_align(sim$X["a"], design = sim$design["a"]),
               "at least two")
  missing_col <- sim
  missing_col$design$b$event_id <- NULL
  expect_error(ssc_fit(missing_col), "missing")
})

test_that("the aligner declares an encoder contract and dispatches natively", {
  sim <- ssc_simulate(seed = 31L)
  algo <- stimulus_sumcor_aligner()
  caps <- aligner_capabilities(algo)
  expect_identical(caps$group, "encoder")
  expect_true(caps$supports_multi)
  expect_identical(caps$encoder, "linear")
  expect_true(caps$supports_predict)
  expect_true(caps$supports_calibration)
  expect_false(caps$requires_auxiliary_at_predict)
  expect_false(caps$invertible)
  expect_identical(aligner_capabilities(new_aligner("x", group = "O")),
                   list(group = "O", supports_multi = FALSE))
  ctrl <- list(beta = 0, brain_rank = 12L)
  via_many <- fit_many(algo, sim$X, design = sim$design, k = 2L, control = ctrl)
  via_align <- align_many(sim$X, algo, k = 2L, design = sim$design,
                          control = ctrl)
  direct <- stimulus_sumcor_align(sim$X, design = sim$design, ncomp = 2L,
                                  control = ctrl)
  expect_identical(via_many$encoders, direct$encoders)
  expect_identical(via_align$encoders, direct$encoders)
})
