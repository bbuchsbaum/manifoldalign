# Held-out-domain transfer evaluator and its ablation helpers.

# Shared latent per event, linear maps per domain, `reps` repeats per event.
# Domains named in `unrelated` get their own latent, independent of the shared
# one, so nothing about them transfers.
.cvt_sim <- function(n_events = 40, reps = 2, k = 2,
                     p = c(a = 8, b = 9, c = 7, d = 10), noise = 0.05,
                     unrelated = character(0), seed = 11) {
  set.seed(seed)
  latent <- matrix(rnorm(n_events * k), n_events, k)
  X <- list()
  design <- list()
  for (nm in names(p)) {
    lat <- if (nm %in% unrelated) {
      matrix(rnorm(n_events * k), n_events, k)
    } else latent
    rows <- rep(seq_len(n_events), times = reps)
    X[[nm]] <- lat[rows, , drop = FALSE] %*% matrix(rnorm(k * p[[nm]]), k,
                                                     p[[nm]]) +
      matrix(rnorm(length(rows) * p[[nm]], sd = noise), length(rows), p[[nm]])
    design[[nm]] <- data.frame(
      obs_id = paste0(nm, "_", seq_along(rows)),
      event_id = paste0("e", rows),
      repeat_id = rep(seq_len(reps), each = n_events),
      y = latent[rows, 1],
      cls = ifelse(latent[rows, 2] > 0, "pos", "neg"),
      stringsAsFactors = FALSE
    )
  }
  list(X = X, design = design, latent = latent)
}

.cvt_fit_sumcor <- function(data, design, ncomp = 2) {
  sumcor_align(data, design = design, event_id = "event_id",
               repeat_id = "repeat_id", observation_id = "obs_id",
               ncomp = ncomp,
               control = sumcor_control(beta = 0, basis_rank = 6))
}

# A sumcor method whose hooks record every row identifier they receive.
.cvt_recording_method <- function(log) {
  log$fit <- list()
  log$calibrate <- list()
  log$calibrate_columns <- list()
  log$calibrate_finite <- logical(0)
  list(
    fit = function(data, design) {
      log$fit[[length(log$fit) + 1L]] <- unlist(lapply(design, `[[`,
                                                        "obs_id"))
      log$fit_events <- unlist(lapply(design, `[[`, "event_id"))
      .cvt_fit_sumcor(data, design)
    },
    calibrate = function(fit, data, design, domain) {
      log$calibrate[[length(log$calibrate) + 1L]] <- design$obs_id
      log$calibrate_columns[[length(log$calibrate_columns) + 1L]] <-
        names(design)
      log$calibrate_finite <- c(log$calibrate_finite, all(is.finite(data)))
      calibrate_domain(fit, data, design = design, domain_name = domain)
    }
  )
}

# An exact oracle method: each domain's encoder is the least-squares map from
# its rows to the true latent of their events, so the shared space is known.
.cvt_oracle_method <- function(latent) {
  truth <- function(design) {
    latent[as.integer(sub("e", "", design$event_id)), , drop = FALSE]
  }
  lsq <- function(X, design) qr.solve(cbind(1, as.matrix(X)), truth(design))
  list(
    fit = function(data, design) {
      structure(list(W = Map(lsq, data, design)), class = "cvt_oracle_fit")
    },
    calibrate = function(fit, data, design, domain) {
      fit$W[[domain]] <- lsq(data, design)
      fit
    },
    encode = function(fit, x, domain) cbind(1, as.matrix(x)) %*% fit$W[[domain]]
  )
}

test_that("oracle: a known shared linear latent transfers almost perfectly", {
  sim <- .cvt_sim()
  res <- cv_alignment_transfer(sim$X,
                               list(oracle = .cvt_oracle_method(sim$latent)),
                               design = sim$design, event_id = "event_id",
                               target = "y")
  expect_s3_class(res, "cv_transfer_result")
  sc <- as.data.frame(res)
  expect_equal(nrow(sc), 4L)
  expect_equal(sc$heldout, c("a", "b", "c", "d"))
  expect_true(all(sc$status == "ok"))
  expect_true(all(sc$target_metric == "r2"))
  expect_true(all(sc$target_score > 0.995))
  expect_true(all(sc$recon_self_r2 > 0.995))
  expect_true(all(sc$recon_transfer_r2 > 0.995))
  expect_true(all(sc$recon_source_r2 > 0.995))
  expect_equal(sc$n_calibration + sc$n_evaluation, rep(80L, 4))
  expect_true(all(sc$n_transfer == sc$n_evaluation))
  expect_true(all(sc$n_target == sc$n_evaluation))

  cls <- cv_alignment_transfer(sim$X,
                               list(oracle = .cvt_oracle_method(sim$latent)),
                               design = sim$design, event_id = "event_id",
                               target = "cls")
  sc <- as.data.frame(cls)
  expect_true(all(sc$target_metric == "accuracy"))
  expect_true(all(sc$target_score >= 0.95))
  # The baseline is the majority source class rate on the evaluation rows.
  fold <- cls$folds[[1]]
  y_eval <- sim$design$a$cls[fold$evaluation]
  src <- unlist(lapply(sim$design[c("b", "c", "d")], `[[`, "cls"))
  majority <- names(which.max(table(src)))
  expect_equal(sc$target_baseline[1], mean(y_eval == majority))
})

test_that("sumcor_align transfers a shared linear latent out of the box", {
  sim <- .cvt_sim(n_events = 150)
  res <- cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                               design = sim$design, event_id = "event_id",
                               target = "y")
  sc <- as.data.frame(res)
  expect_true(all(sc$status == "ok"))
  expect_true(all(sc$target_score > 0.9))
  expect_true(all(sc$recon_self_r2 > 0.99))
  expect_true(all(sc$recon_transfer_r2 > 0.9))
  expect_true(all(sc$recon_source_r2 > 0.9))
})

test_that("negative control: an unrelated held-out domain scores at chance", {
  sim <- .cvt_sim(n_events = 80, unrelated = "d")
  folds <- list(list(heldout = "a", sources = c("b", "c")),
                list(heldout = "d", sources = c("a", "b", "c")))
  res <- cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                               design = sim$design, folds = folds,
                               event_id = "event_id", target = "cls")
  sc <- as.data.frame(res)
  expect_true(all(sc$status == "ok"))
  related <- sc[sc$heldout == "a", ]
  control <- sc[sc$heldout == "d", ]
  expect_gt(related$target_score, 0.9)
  # Chance for a balanced binary label is 0.5 (sd about 0.08 over 40 events).
  expect_lt(abs(control$target_score - 0.5), 0.2)
  expect_lt(control$recon_transfer_r2, 0.1)
  expect_lt(control$recon_source_r2, 0.1)
  # The held-out domain still reconstructs itself: compression is not transfer.
  expect_gt(control$recon_self_r2, 0.9)
})

test_that("calibration and readouts never see held-out evaluation rows", {
  sim <- .cvt_sim()
  log <- new.env()
  res <- cv_alignment_transfer(sim$X, list(rec = .cvt_recording_method(log)),
                               design = sim$design, event_id = "event_id",
                               target = "y")
  expect_true(all(res$scores$status == "ok"))
  for (f in seq_along(res$folds)) {
    fold <- res$folds[[f]]
    held_ids <- sim$design[[fold$heldout]]$obs_id
    expect_false(any(held_ids %in% log$fit[[f]]))
    expect_setequal(log$calibrate[[f]], held_ids[fold$calibration])
    expect_false(any(held_ids[fold$evaluation] %in% log$calibrate[[f]]))
    # The target column is withheld from calibration.
    expect_false("y" %in% log$calibrate_columns[[f]])
  }
})

test_that("poisoned evaluation rows leave calibration and readouts unchanged", {
  sim <- .cvt_sim()
  base <- cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                                design = sim$design, folds = "b",
                                event_id = "event_id", target = "y",
                                keep_models = TRUE)
  fold <- base$folds[[1]]
  explicit <- list(list(heldout = "b", calibration = fold$calibration))

  # Huge finite values: everything fitted on permitted partitions is
  # identical; only scores of the poisoned rows change.
  huge <- sim
  huge$X$b[fold$evaluation, ] <- 1e8
  poisoned <- cv_alignment_transfer(huge$X, list(sumcor = .cvt_fit_sumcor),
                                    design = huge$design, folds = explicit,
                                    event_id = "event_id", target = "y",
                                    keep_models = TRUE)
  expect_identical(poisoned$folds[[1]]$evaluation, fold$evaluation)
  m0 <- base$models[[1]]
  m1 <- poisoned$models[[1]]
  expect_identical(m1$calibrated$encoders$b, m0$calibrated$encoders$b)
  expect_identical(m1$readouts$heldout, m0$readouts$heldout)
  expect_identical(m1$readouts$target, m0$readouts$target)
  expect_identical(m1$readouts$sources, m0$readouts$sources)
  expect_false(isTRUE(all.equal(poisoned$scores$recon_self_r2,
                                base$scores$recon_self_r2)))

  # Missing values: fitting and calibration still succeed; the run fails only
  # where evaluation rows are encoded.
  na <- sim
  na$X$b[fold$evaluation, ] <- NA_real_
  log <- new.env()
  failed <- cv_alignment_transfer(na$X, list(rec = .cvt_recording_method(log)),
                                  design = na$design, folds = explicit,
                                  event_id = "event_id", target = "y")
  expect_equal(failed$scores$status, "failed")
  expect_equal(failed$scores$stage, "encode")
  expect_true(all(log$calibrate_finite))

  # Held-out labels: calibration-row labels are never used; evaluation-row
  # labels only enter the target score.
  relabel <- sim
  relabel$design$b$y[fold$calibration] <- 1e6
  same <- cv_alignment_transfer(relabel$X, list(sumcor = .cvt_fit_sumcor),
                                design = relabel$design, folds = explicit,
                                event_id = "event_id", target = "y")
  expect_identical(same$scores[, -1], base$scores[, -1])
  relabel <- sim
  relabel$design$b$y[fold$evaluation] <- -relabel$design$b$y[fold$evaluation]
  flipped <- cv_alignment_transfer(relabel$X, list(sumcor = .cvt_fit_sumcor),
                                   design = relabel$design, folds = explicit,
                                   event_id = "event_id", target = "y")
  expect_lt(flipped$scores$target_score, 0)
  expect_identical(flipped$scores$recon_transfer_r2,
                   base$scores$recon_transfer_r2)
})

test_that("event repeats and grouping columns never straddle the split", {
  sim <- .cvt_sim(n_events = 40)
  # Runs come in linked pairs: (1, 2), (3, 4), ... share their events.
  for (nm in names(sim$design)) {
    ev <- as.integer(sub("e", "", sim$design[[nm]]$event_id))
    pair <- (ev - 1L) %/% 5L
    sim$design[[nm]]$run <- 2L * pair + sim$design[[nm]]$repeat_id
  }
  res <- cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                               design = sim$design, event_id = "event_id",
                               group_by = "run", calibration_fraction = 0.4,
                               seed = 3)
  for (fold in res$folds) {
    d <- sim$design[[fold$heldout]]
    expect_length(intersect(d$event_id[fold$calibration],
                            d$event_id[fold$evaluation]), 0L)
    expect_length(intersect(d$run[fold$calibration], d$run[fold$evaluation]),
                  0L)
    expect_equal(fold$n_calibration_groups + fold$n_evaluation_groups, 8L)
  }
  # Without a target the target columns are NA, not a readout failure.
  expect_true(all(res$scores$status == "ok"))
  expect_true(all(is.na(res$scores$target_metric)))
  expect_true(all(is.na(res$scores$target_score)))
  expect_true(all(is.finite(res$scores$recon_transfer_r2)))

  # An explicit calibration set that splits a group is rejected.
  expect_error(
    cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                          design = sim$design, event_id = "event_id",
                          folds = list(list(heldout = "a", calibration = 1:40))),
    "splits"
  )
  # A held-out domain that is one connected group cannot be split.
  one <- sim
  one$design$a$block <- 1L
  expect_error(
    cv_alignment_transfer(one$X, list(sumcor = .cvt_fit_sumcor),
                          design = one$design, event_id = "event_id",
                          group_by = "block", folds = "a"),
    "single connected group"
  )
})

test_that("novel evaluation groups are removed from the sources", {
  sim <- .cvt_sim()
  log <- new.env()
  res <- cv_alignment_transfer(sim$X, list(rec = .cvt_recording_method(log)),
                               design = sim$design, folds = "c",
                               event_id = "event_id", target = "y",
                               evaluation_groups = "novel")
  fold <- res$folds[[1]]
  eval_events <- unique(sim$design$c$event_id[fold$evaluation])
  expect_length(intersect(log$fit_events, eval_events), 0L)
  expect_equal(res$scores$n_source_rows,
               3L * (80L - 2L * length(eval_events)))
  expect_equal(res$scores$status, "ok")
  expect_gt(res$scores$recon_transfer_r2, 0.9)
  expect_gt(res$scores$target_score, 0.9)
  expect_error(
    cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                          design = sim$design, evaluation_groups = "novel"),
    "requires `event_id`"
  )
})

test_that("methods without frozen calibration are reported as unsupported", {
  sim <- .cvt_sim(p = c(a = 6, b = 7, c = 5))
  fit_called <- FALSE
  methods <- list(
    grasp = grasp_aligner(),
    toy = function(data, design) {
      fit_called <<- TRUE
      structure(list(), class = "toy_transfer_fit")
    },
    aux = new_aligner("aux_transfer", group = "encoder", supports_multi = TRUE,
                      supports_calibration = TRUE,
                      requires_auxiliary_at_predict = TRUE),
    sumcor = .cvt_fit_sumcor
  )
  res <- cv_alignment_transfer(sim$X, methods, design = sim$design,
                               event_id = "event_id", target = "y")
  sc <- as.data.frame(res)
  expect_equal(nrow(sc), 3L * 4L)
  expect_true(all(sc$status[sc$method %in% c("grasp", "toy", "aux")] ==
                    "unsupported"))
  expect_true(all(sc$stage[sc$status == "unsupported"] == "capability"))
  expect_match(sc$message[sc$method == "grasp"][1], "supports_calibration")
  expect_match(sc$message[sc$method == "toy"][1], "calibrate_domain")
  expect_match(sc$message[sc$method == "aux"][1], "auxiliary")
  expect_true(fit_called)
  expect_true(all(is.na(sc$target_score[sc$status != "ok"])))
  expect_true(all(sc$status[sc$method == "sumcor"] == "ok"))

  # The sumcor aligner is supported through fit_many().
  res <- cv_alignment_transfer(
    sim$X,
    list(sa = list(fit = sumcor_aligner(),
                   args = list(ncomp = 2, event_id = "event_id",
                               repeat_id = "repeat_id",
                               control = sumcor_control(beta = 0,
                                                        basis_rank = 6)))),
    design = sim$design, event_id = "event_id", target = "y"
  )
  expect_true(all(res$scores$status == "ok"))
})

test_that("failures are kept as rows with their stage and message", {
  sim <- .cvt_sim(p = c(a = 6, b = 7, c = 5))
  methods <- list(
    boom = function(data, design) stop("fit exploded"),
    badcal = list(fit = .cvt_fit_sumcor,
                  calibrate = function(fit, data, design, domain) {
                    stop("calibration exploded")
                  }),
    badenc = list(fit = .cvt_fit_sumcor,
                  encode = function(fit, x, domain) matrix(0, 1, 2)),
    warn = function(data, design) {
      warning("watch out")
      .cvt_fit_sumcor(data, design)
    },
    sumcor = .cvt_fit_sumcor
  )
  res <- cv_alignment_transfer(sim$X, methods, design = sim$design,
                               event_id = "event_id", target = "y")
  sc <- as.data.frame(res)
  expect_equal(nrow(sc), 15L)
  by <- split(sc, sc$method)
  expect_true(all(by$boom$status == "failed"))
  expect_true(all(by$boom$stage == "fit"))
  expect_true(all(by$boom$message == "fit exploded"))
  expect_true(all(by$badcal$stage == "calibrate"))
  expect_true(all(by$badcal$message == "calibration exploded"))
  expect_true(all(by$badenc$stage == "encode"))
  expect_match(by$badenc$message[1], "returned 1 rows")
  expect_true(all(by$warn$status == "ok"))
  expect_true(all(by$warn$warnings == "watch out"))
  expect_true(all(by$sumcor$status == "ok"))
  expect_true(all(is.na(by$sumcor$warnings)))
  out <- capture.output(print(res))
  expect_true(any(grepl("failed at fit: fit exploded", out, fixed = TRUE)))
  expect_true(any(grepl("3 failed (of 3 folds)", out, fixed = TRUE)))
})

test_that("results are deterministic and the caller's RNG is restored", {
  sim <- .cvt_sim(p = c(a = 6, b = 7, c = 5))
  run <- function(seed) {
    cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                          design = sim$design, event_id = "event_id",
                          target = "y", seed = seed)
  }
  set.seed(99)
  before <- .Random.seed
  r1 <- run(5)
  expect_identical(.Random.seed, before)
  r2 <- run(5)
  expect_identical(r1$scores, r2$scores)
  expect_identical(r1$folds, r2$folds)
  r3 <- run(6)
  expect_false(identical(r1$folds, r3$folds))

  rm(".Random.seed", envir = globalenv())
  run(5)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(99)
})

test_that("sparse and hyperdesign inputs give the same scores", {
  sim <- .cvt_sim(p = c(a = 6, b = 7, c = 5))
  dense <- cv_alignment_transfer(sim$X, list(sumcor = .cvt_fit_sumcor),
                                 design = sim$design, event_id = "event_id",
                                 target = "y")
  sparse_X <- lapply(sim$X, function(x) Matrix::Matrix(x, sparse = TRUE))
  sparse <- cv_alignment_transfer(sparse_X, list(sumcor = .cvt_fit_sumcor),
                                  design = sim$design, event_id = "event_id",
                                  target = "y")
  expect_equal(sparse$scores, dense$scores, tolerance = 1e-8)

  skip_if_not_installed("multidesign")
  hd <- multidesign::hyperdesign(
    lapply(names(sim$X), function(nm) {
      multidesign::multidesign(sim$X[[nm]], sim$design[[nm]])
    }),
    block_names = names(sim$X)
  )
  from_hd <- cv_alignment_transfer(hd, list(sumcor = .cvt_fit_sumcor),
                                   event_id = "event_id", target = "y")
  expect_equal(from_hd$scores, dense$scores, tolerance = 1e-10)
  expect_error(
    cv_alignment_transfer(hd, list(sumcor = .cvt_fit_sumcor),
                          design = sim$design),
    "must be NULL for a hyperdesign"
  )
})

test_that("invalid specifications are rejected explicitly", {
  sim <- .cvt_sim(p = c(a = 6, b = 7, c = 5))
  go <- function(...) {
    cv_alignment_transfer(sim$X, design = sim$design, ...)
  }
  expect_error(go(methods = list(m = list(fit = .cvt_fit_sumcor, foo = 1))),
               "unknown specification element")
  expect_error(go(methods = list(.cvt_fit_sumcor, .cvt_fit_sumcor)),
               "unique, non-empty names")
  expect_error(go(methods = list(m = 1)), "function, an aligner")
  expect_error(go(methods = .cvt_fit_sumcor, folds = list(list(held = "a"))),
               "unknown or unnamed")
  expect_error(go(methods = .cvt_fit_sumcor, folds = "zz"), "`heldout`")
  expect_error(go(methods = .cvt_fit_sumcor,
                  folds = list(list(heldout = "a", sources = c("a", "b")))),
               "cannot also be a source")
  expect_error(go(methods = .cvt_fit_sumcor, target = "nope"),
               "missing from domain")
  expect_error(go(methods = .cvt_fit_sumcor, calibration_fraction = 1),
               "calibration_fraction")
  expect_error(go(methods = .cvt_fit_sumcor, target = "event_id",
                  event_id = "event_id"), "must not also be")
  expect_error(go(methods = .cvt_fit_sumcor, unknown_argument = 1),
               "unused argument")
})

test_that("explicit folds restrict the sources", {
  sim <- .cvt_sim()
  log <- new.env()
  res <- cv_alignment_transfer(
    sim$X, list(rec = .cvt_recording_method(log)), design = sim$design,
    folds = list(list(heldout = "d", sources = c("a", "b"))),
    event_id = "event_id", target = "y"
  )
  expect_equal(res$folds[[1]]$sources, c("a", "b"))
  expect_true(all(grepl("^[ab]_", log$fit[[1]])))
  expect_equal(res$scores$n_source_rows, 160L)
  expect_equal(res$scores$status, "ok")
})

# -- Ablation helpers --------------------------------------------------------------

test_that("event_order_kernel factors the Gaussian ordering kernel", {
  ids <- paste0("e", 1:12)
  ord <- c(3, 1, 7, 2, 9, 4, 5, 12, 6, 8, 11, 10)
  F <- event_order_kernel(ids, order = ord, bandwidth = 1.5)
  K <- exp(-outer(ord, ord, "-")^2 / (2 * 1.5^2))
  expect_equal(rownames(F), ids)
  expect_equal(unname(tcrossprod(F)), K, tolerance = 1e-8)

  # Keyed by event, not by row position.
  perm <- sample.int(12)
  Fp <- event_order_kernel(ids[perm], order = ord[perm], bandwidth = 1.5)
  expect_equal(tcrossprod(Fp)[ids, ids], tcrossprod(F), tolerance = 1e-8)

  blocks <- rep(c("r1", "r2"), each = 6)
  Fb <- event_order_kernel(ids, order = ord, bandwidth = 1.5, blocks = blocks)
  Kb <- tcrossprod(Fb)
  expect_equal(max(abs(Kb[1:6, 7:12])), 0)
  expect_equal(unname(Kb[1:6, 1:6]), K[1:6, 1:6], tolerance = 1e-8)

  low <- event_order_kernel(ids, order = ord, bandwidth = 1.5, rank = 2)
  expect_equal(ncol(low), 2L)
  expect_error(event_order_kernel(c("a", "a")), "unique")
  expect_error(event_order_kernel(ids, bandwidth = 0), "bandwidth")
  expect_error(event_order_kernel(ids, blocks = c(NA, blocks[-1])), "blocks")
})

test_that("an ordering kernel plugs into sumcor_align as event features", {
  sim <- .cvt_sim(p = c(a = 6, b = 7, c = 5))
  ids <- paste0("e", 1:40)
  F <- event_order_kernel(ids, bandwidth = 2, rank = 10)
  fit <- sumcor_align(sim$X, design = sim$design, event_id = "event_id",
                      repeat_id = "repeat_id", ncomp = 2,
                      event_features = list(order = F),
                      control = sumcor_control(beta = 0.5, basis_rank = 6,
                                               kernel_scale = FALSE))
  expect_s3_class(fit, "sumcor_align")
  expect_equal(fit$kernel$beta, 0.5)
})

test_that("circular-shift surrogates rotate the ordered sequence", {
  set.seed(4)
  feats <- matrix(cumsum(rnorm(60)), 30, 2,
                  dimnames = list(paste0("e", 1:30), c("f1", "f2")))
  s <- surrogate_event_features(feats, "circular_shift", min_shift = 5,
                                seed = 2)
  k <- attr(s, "surrogate")$shifts
  expect_true(k >= 5 && k <= 25)
  expect_equal(dimnames(s), dimnames(feats))
  expect_equal(unname(s[, , drop = FALSE]),
               unname(feats[((0:29 + k) %% 30) + 1, ]), ignore_attr = TRUE)
  # Circular autocorrelation is preserved exactly; the link to events is not.
  acf_circ <- function(x, lag) sum(x * x[((seq_along(x) - 1 + lag) %% 30) + 1])
  for (lag in 1:5) {
    expect_equal(acf_circ(s[, 1], lag), acf_circ(feats[, 1], lag))
  }
  expect_false(isTRUE(all.equal(unname(s[, 1]), unname(feats[, 1]))))

  # `order` controls the sequence; rows keep their identifiers.
  ord <- 30:1
  so <- surrogate_event_features(feats, "circular_shift", order = ord,
                                 seed = 2)
  k <- attr(so, "surrogate")$shifts
  rev_feats <- feats[30:1, ]
  expected <- rev_feats[((0:29 + k) %% 30) + 1, ]
  expect_equal(unname(so[30:1, ]), unname(expected), ignore_attr = TRUE)
  expect_error(surrogate_event_features(feats[1:3, ], min_shift = 2),
               "too short")
})

test_that("phase-randomised surrogates keep the periodogram and cross-spectrum", {
  set.seed(5)
  n <- 32
  base <- cumsum(rnorm(n))
  feats <- cbind(f1 = base, f2 = base + rnorm(n, sd = 0.3))
  rownames(feats) <- paste0("e", seq_len(n))
  for (m in c(n, n - 1)) {
    x <- feats[seq_len(m), ]
    s <- surrogate_event_features(x, "phase_randomize", seed = 1)
    fx <- stats::mvfft(unname(x))
    fs <- stats::mvfft(unname(s[, , drop = FALSE]))
    expect_equal(Mod(fs), Mod(fx), tolerance = 1e-8)
    expect_equal(fs[, 1] * Conj(fs[, 2]), fx[, 1] * Conj(fx[, 2]),
                 tolerance = 1e-8)
    expect_equal(unname(colMeans(s)), unname(colMeans(x)), tolerance = 1e-10)
    expect_lt(abs(stats::cor(s[, 1], x[, 1])), 0.9)
  }

  # Within blocks: each block's periodogram is preserved separately.
  blocks <- rep(1:2, each = n / 2)
  sb <- surrogate_event_features(feats, "phase_randomize", blocks = blocks,
                                 seed = 3)
  for (b in 1:2) {
    idx <- which(blocks == b)
    expect_equal(Mod(stats::fft(sb[idx, 1])), Mod(stats::fft(feats[idx, 1])),
                 tolerance = 1e-8)
  }
})

test_that("surrogates are seeded and restore the caller's RNG", {
  feats <- matrix(seq_len(40) / 7, 20, 2,
                  dimnames = list(paste0("e", 1:20), NULL))
  set.seed(42)
  before <- .Random.seed
  a <- surrogate_event_features(feats, "phase_randomize", seed = 8)
  expect_identical(.Random.seed, before)
  b <- surrogate_event_features(feats, "phase_randomize", seed = 8)
  expect_identical(a, b)
  c <- surrogate_event_features(feats, "phase_randomize", seed = 9)
  expect_false(identical(a, c))
  expect_error(surrogate_event_features(feats[1:2, ], "phase_randomize"),
               "at least three")
  expect_error(surrogate_event_features(feats, seed = "x"), "seed")
})
