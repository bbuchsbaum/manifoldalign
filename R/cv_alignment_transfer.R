# Held-out-domain transfer evaluation ---------------------------------------
#
# Estimand: how well an alignment fitted on source domains transfers to a
# domain it never saw. For each fold one domain is held out. Its rows are
# split, by connected groups, into a calibration part and an evaluation part.
#   1. The reference is fitted on the source domains only and then frozen.
#   2. The held-out domain is calibrated against the frozen reference from its
#      calibration rows only (methods without calibration are "unsupported").
#   3. Readouts are fitted on permitted partitions only: the target readout on
#      source rows, the held-out reconstruction readout D_h on calibration
#      rows, each source readout D_s on source rows outside evaluation events.
#   4. Scores are computed on the held-out evaluation rows.
# Every (fold, method) pair yields exactly one row with status "ok",
# "unsupported" or "failed"; no row is ever dropped.

#' Held-out-domain transfer evaluation for alignment methods
#'
#' Estimates how well an alignment fitted on source domains transfers to a
#' domain it never saw. Each fold holds out one whole domain. The held-out
#' domain's rows are split into a calibration part and an evaluation part such
#' that every connected group of rows stays on one side: all repeats of an
#' event (rows sharing an `event_id` value) and all rows sharing a value of any
#' `group_by` column. For every fold and method:
#' \enumerate{
#'   \item the method is fitted on the source domains only (the reference) and
#'     then frozen;
#'   \item the held-out domain is calibrated against the frozen reference using
#'     only its calibration rows. Methods that cannot add a domain to a frozen
#'     fit get an explicit `"unsupported"` row; they are never refitted with
#'     the held-out domain included;
#'   \item the held-out calibration and evaluation rows are encoded with the
#'     calibrated encoder, and source rows with the frozen reference;
#'   \item readouts are fitted on permitted partitions and scored on the
#'     held-out evaluation rows (see Scores).
#' }
#'
#' @section Methods:
#' `methods` is a named list. Each element is one of
#' \itemize{
#'   \item a function `function(data, design)` that fits the method on the
#'     source domains (a named list of row-subset matrices and the matching
#'     named list of design data frames) and returns a fitted object;
#'   \item an [aligner][new_aligner] object, fitted with
#'     `fit_many(aligner, data, design = design)`. Its
#'     [aligner_capabilities()] are checked first: an aligner that does not
#'     declare `supports_calibration = TRUE`, declares
#'     `supports_predict = FALSE`, or requires auxiliary inputs at prediction
#'     is reported as `"unsupported"` without being fitted;
#'   \item a list with element `fit` (a function or aligner as above) and
#'     optional elements `args` (a list of further arguments for the fit),
#'     `encode` (`function(fit, x, domain)` returning one latent row per row of
#'     `x`) and `calibrate` (`function(fit, data, design, domain)` returning
#'     the fit with `domain` added). Other elements are rejected.
#' }
#' By default rows are encoded with `oos_predict(fit, x, side = domain)` and
#' the held-out domain is calibrated with
#' `calibrate_domain(fit, data, design = design, domain_name = domain)`. A fit
#' whose class has no `calibrate_domain()` (or, without a custom `encode`, no
#' `oos_predict()`) method is reported as `"unsupported"`. [sumcor_align()]
#' fits work out of the box.
#'
#' @section Leakage controls:
#' The held-out domain's evaluation rows never reach fitting, calibration or
#' any readout fit; they are only encoded and scored. The `target` column is
#' removed from the calibration design, so calibration never sees held-out
#' labels. Source designs are passed unchanged (source labels are permitted
#' training data). With `evaluation_groups = "shared"` the source domains keep
#' their observations of the held-out evaluation events: this measures
#' transfer to a new domain for events the reference has seen. With
#' `evaluation_groups = "novel"` (requires `event_id`) source rows of the
#' held-out evaluation events are removed before fitting and before the target
#' readout is trained, which measures transfer to a new domain and new events.
#' In both modes, source rows of evaluation events are used only as prediction
#' inputs for `recon_transfer_r2` and as prediction targets for
#' `recon_source_r2`, and are always excluded from fitting the source readouts
#' `D_s`.
#'
#' @section Scores:
#' All readouts are ridge-regularised linear least-squares fits with an
#' intercept (ridge `readout_ridge` times the mean diagonal of the centred
#' latent Gram matrix).
#' \describe{
#'   \item{`target_score`}{the target readout is fitted on encoded source rows
#'     with a non-missing `target` and applied to the encoded held-out
#'     evaluation rows. A numeric target is scored by out-of-sample R-squared
#'     relative to the source training mean (`target_metric = "r2"`,
#'     `target_baseline = 0`). Any other target is a class label scored by
#'     accuracy (`target_metric = "accuracy"`) with a one-hot linear readout
#'     and arg-max prediction; `target_baseline` is the accuracy of always
#'     predicting the most frequent source class. A custom `target_readout`
#'     replaces the linear readout.}
#'   \item{`recon_self_r2`}{an explicit linear readout `D_h` from latent
#'     coordinates back to the held-out domain's native columns is fitted on
#'     the calibration rows; held-out evaluation rows are reconstructed from
#'     their own encodings. R-squared is pooled over columns, relative to the
#'     calibration column means.}
#'   \item{`recon_transfer_r2`}{(requires `event_id`) each evaluation row is
#'     predicted through `D_h` from the source domains' consensus encoding of
#'     its event: the mean over source domains of each domain's mean encoding
#'     of that event. Rows whose event no source domain observed are not
#'     scored; `n_transfer` counts the rows scored.}
#'   \item{`recon_source_r2`}{(requires `event_id`) for each source domain `s`
#'     an explicit readout `D_s` back to its native columns is fitted on its
#'     encoded rows outside the evaluation events; its rows of evaluation
#'     events are predicted from the held-out domain's mean evaluation encoding
#'     of the same event. R-squared is pooled over source domains and columns,
#'     relative to each `D_s` training mean.}
#' }
#' Every method is therefore scored in the same outcome spaces (the target and
#' each domain's native columns), whatever its latent dimension.
#'
#' @section Status:
#' Each fold and method gives one row. `status` is `"ok"`, `"unsupported"`
#' (the method cannot calibrate a new domain against a frozen fit or cannot
#' encode new rows; `message` gives the reason) or `"failed"` (an error was
#' raised; `stage` is one of `"fit"`, `"calibrate"`, `"encode"`,
#' `"readout"` and `message` keeps the error message). Unsupported rows have
#' `stage = "capability"`. Metrics are `NA` unless
#' the status is `"ok"`. Warnings raised by a method are recorded in the
#' `warnings` column instead of being re-signalled.
#'
#' @param data A named list of numeric matrices (rows are observations; base
#'   matrices, data frames or `Matrix` objects), one per domain, or a
#'   `hyperdesign`.
#' @param methods A named list of method specifications (see Methods). A
#'   single unnamed function or aligner is named `"method"`.
#' @param design For list `data`, a list of data frames (one per domain, rows
#'   matching `data`) holding identifier, grouping and target columns. Must be
#'   `NULL` for a `hyperdesign`, whose block designs are used.
#' @param folds `NULL` (default) for leave-one-domain-out over all domains; a
#'   character vector of domains to hold out, one fold each; or a list of fold
#'   specifications, each a list with `heldout` (one domain name) and optional
#'   `sources` (domain names used for fitting; default all others) and
#'   `calibration` (row indices or a logical vector over the held-out domain
#'   selecting its calibration rows; the rest are evaluation rows). An explicit
#'   `calibration` that splits a group is an error.
#' @param event_id Optional name of the design column holding correspondence
#'   keys shared across domains (events). All rows of one event stay on one
#'   side of the calibration split. Required for `recon_transfer_r2`,
#'   `recon_source_r2` and `evaluation_groups = "novel"`.
#' @param group_by Optional names of further design columns of the held-out
#'   domain whose values must not straddle the calibration split (for example
#'   a run or session column).
#' @param target Optional name of a design column holding the target variable
#'   for the target readout. It must be present in every domain.
#' @param calibration_fraction Fraction of the held-out domain's rows assigned
#'   (by whole groups, in random order) to calibration when a fold does not
#'   give `calibration` explicitly.
#' @param evaluation_groups `"shared"` (default) or `"novel"`; see Leakage
#'   controls.
#' @param target_readout Optional `function(z_train, y_train, z_test)`
#'   returning one prediction per row of `z_test`, replacing the default
#'   linear readout.
#' @param readout_ridge Nonnegative relative ridge for the linear readouts.
#' @param seed Integer seed. Group splits and each (fold, method) evaluation
#'   run under fixed seeds derived from it, and the caller's random number
#'   stream is restored on exit.
#' @param keep_models Logical; keep each fold's reference fit, calibrated fit
#'   and readouts in the result (can be large).
#'
#' @return An object of class `cv_transfer_result`: a list with `scores` (a
#'   data frame with one row per fold and method; columns `.fold`, `heldout`,
#'   `method`, `status`, `stage`, `message`, `warnings`, `target_metric`,
#'   `target_score`, `target_baseline`, `recon_self_r2`, `recon_transfer_r2`,
#'   `recon_source_r2`, `n_calibration`, `n_evaluation`, `n_target`,
#'   `n_transfer`, `n_source_rows` and `.seed`), `folds` (the held-out domain,
#'   sources, and calibration and evaluation row indices of each fold),
#'   `settings`, `models` (when `keep_models = TRUE`) and `call`. Use
#'   `as.data.frame()` for the scores.
#'
#' @seealso [cv_alignment_rows()] for row-holdout evaluation within domains,
#'   [calibrate_domain()], [sumcor_align()].
#' @examples
#' set.seed(1)
#' n_events <- 30
#' latent <- matrix(rnorm(n_events * 2), n_events, 2)
#' make_domain <- function(p) {
#'   rows <- rep(seq_len(n_events), times = 2)
#'   X <- latent[rows, ] %*% matrix(rnorm(2 * p), 2, p) +
#'     matrix(rnorm(length(rows) * p, sd = 0.2), length(rows), p)
#'   list(X = X, design = data.frame(
#'     event_id = paste0("e", rows), repeat_id = rep(1:2, each = n_events),
#'     y = latent[rows, 1]
#'   ))
#' }
#' doms <- list(a = make_domain(8), b = make_domain(9), c = make_domain(7))
#' fit_sumcor <- function(data, design) {
#'   sumcor_align(data, design = design, event_id = "event_id",
#'                repeat_id = "repeat_id", ncomp = 2,
#'                control = sumcor_control(beta = 0, basis_rank = 5))
#' }
#' res <- cv_alignment_transfer(
#'   lapply(doms, `[[`, "X"), methods = list(sumcor = fit_sumcor),
#'   design = lapply(doms, `[[`, "design"),
#'   event_id = "event_id", target = "y"
#' )
#' res
#' as.data.frame(res)[, c("heldout", "status", "target_score",
#'                        "recon_transfer_r2")]
#' @export
cv_alignment_transfer <- function(data,
                                  methods,
                                  design = NULL,
                                  folds = NULL,
                                  event_id = NULL,
                                  group_by = NULL,
                                  target = NULL,
                                  calibration_fraction = 0.5,
                                  evaluation_groups = c("shared", "novel"),
                                  target_readout = NULL,
                                  readout_ridge = 1e-6,
                                  seed = 1L,
                                  keep_models = FALSE) {
  call <- match.call()
  evaluation_groups <- match.arg(evaluation_groups)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      abs(seed - round(seed)) > 1e-8) {
    stop("`seed` must be a single integer.", call. = FALSE)
  }
  seed <- as.integer(round(seed))
  if (!is.numeric(calibration_fraction) || length(calibration_fraction) != 1L ||
      !is.finite(calibration_fraction) || calibration_fraction <= 0 ||
      calibration_fraction >= 1) {
    stop("`calibration_fraction` must be a single number in (0, 1).",
         call. = FALSE)
  }
  if (!is.numeric(readout_ridge) || length(readout_ridge) != 1L ||
      !is.finite(readout_ridge) || readout_ridge < 0) {
    stop("`readout_ridge` must be a single nonnegative number.", call. = FALSE)
  }
  if (!is.logical(keep_models) || length(keep_models) != 1L ||
      is.na(keep_models)) {
    stop("`keep_models` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(target_readout) && !is.function(target_readout)) {
    stop("`target_readout` must be NULL or a function.", call. = FALSE)
  }
  event_id <- .cvt_column_names(event_id, "event_id", single = TRUE)
  group_by <- .cvt_column_names(group_by, "group_by", single = FALSE)
  target <- .cvt_column_names(target, "target", single = TRUE)
  if (identical(evaluation_groups, "novel") && is.null(event_id)) {
    stop("`evaluation_groups = \"novel\"` requires `event_id`.", call. = FALSE)
  }
  if (!is.null(target) && target %in% c(event_id, group_by)) {
    stop("`target` must not also be an `event_id` or `group_by` column.",
         call. = FALSE)
  }

  ctx <- .cvt_resolve_data(data, design)
  methods <- .cvt_resolve_methods(methods)
  .cvt_check_columns(ctx, event_id, "event_id", ctx$domain_names)
  .cvt_check_columns(ctx, target, "target", ctx$domain_names)
  # Assign with `[` so NULL settings stay present as named elements: `$`
  # then matches them exactly instead of partially matching another field.
  ctx[c("event_id", "target", "target_type", "evaluation_groups",
        "readout_ridge", "target_readout", "events")] <- list(
    event_id,
    target,
    .cvt_target_type(ctx, target),
    evaluation_groups,
    readout_ridge,
    target_readout,
    if (is.null(event_id)) NULL else
      lapply(ctx$design, function(d) as.character(d[[event_id]]))
  )

  fold_list <- .with_local_seed(seed, .cvt_resolve_folds(
    folds, ctx, event_id, group_by, calibration_fraction
  ))

  rows <- list()
  models <- list()
  for (f in seq_along(fold_list)) {
    fold <- fold_list[[f]]
    for (m in seq_along(methods)) {
      run_seed <- as.integer((seed + 7919 * f + 104729 * m) %%
                               .Machine$integer.max)
      res <- .with_local_seed(run_seed, .cvt_evaluate(
        methods[[m]], fold, ctx, keep_models
      ))
      row <- res$row
      row$.fold <- f
      row$heldout <- fold$heldout
      row$method <- names(methods)[[m]]
      row$.seed <- run_seed
      rows[[length(rows) + 1L]] <- row
      if (keep_models) {
        models[[length(models) + 1L]] <- c(
          list(.fold = f, heldout = fold$heldout, method = names(methods)[[m]]),
          res$models
        )
      }
    }
  }

  scores <- .cvt_rows_to_frame(rows)
  structure(
    list(
      scores = scores,
      folds = fold_list,
      settings = list(
        methods = names(methods),
        event_id = event_id,
        group_by = group_by,
        target = target,
        target_type = ctx$target_type,
        calibration_fraction = calibration_fraction,
        evaluation_groups = evaluation_groups,
        readout_ridge = readout_ridge,
        custom_target_readout = !is.null(target_readout),
        seed = seed
      ),
      models = if (keep_models) models else NULL,
      call = call
    ),
    class = "cv_transfer_result"
  )
}

# -- Input resolution ----------------------------------------------------------

.cvt_column_names <- function(x, name, single) {
  if (is.null(x)) return(NULL)
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) ||
      anyDuplicated(x) || (single && length(x) != 1L)) {
    stop("`", name, "` must be ",
         if (single) "one non-empty column name" else
           "a character vector of unique column names", ".", call. = FALSE)
  }
  x
}

.cvt_as_data_matrix <- function(x, domain) {
  if (is.data.frame(x)) x <- as.matrix(x)
  ok <- (is.matrix(x) && is.numeric(x)) || methods::is(x, "dMatrix")
  if (!ok) {
    stop("Domain `", domain, "` must be a numeric matrix, data frame or ",
         "Matrix.", call. = FALSE)
  }
  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop("Domain `", domain, "` needs at least two rows and one column.",
         call. = FALSE)
  }
  x
}

.cvt_resolve_data <- function(data, design) {
  if (inherits(data, "hyperdesign")) {
    if (!is.null(design)) {
      stop("`design` must be NULL for a hyperdesign; block designs are used.",
           call. = FALSE)
    }
    resolved <- resolve_hyperdesign(data, keep_sparse = TRUE)
    X <- lapply(resolved$domains, `[[`, "x")
    design <- lapply(resolved$domains, `[[`, "design")
    names(X) <- names(design) <- resolved$domain_names
  } else if (is.list(data) && !is.data.frame(data)) {
    nms <- names(data)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
      stop("List `data` must have unique, non-empty domain names.",
           call. = FALSE)
    }
    X <- data
    if (is.null(design)) {
      design <- vector("list", length(X))
      names(design) <- nms
    } else {
      if (!is.list(design) || is.data.frame(design) ||
          length(design) != length(X)) {
        stop("`design` must be a list with one data frame per domain.",
             call. = FALSE)
      }
      if (!is.null(names(design))) {
        if (!setequal(names(design), nms)) {
          stop("`design` names must match the domain names of `data`.",
               call. = FALSE)
        }
        design <- design[nms]
      } else {
        names(design) <- nms
      }
    }
  } else {
    stop("`data` must be a named list of matrices or a hyperdesign.",
         call. = FALSE)
  }
  domain_names <- names(X)
  if (length(domain_names) < 2L) {
    stop("At least two domains are needed.", call. = FALSE)
  }
  if (anyDuplicated(domain_names)) {
    stop("Domain names must be unique.", call. = FALSE)
  }
  for (d in domain_names) {
    X[[d]] <- .cvt_as_data_matrix(X[[d]], d)
    n <- nrow(X[[d]])
    if (is.null(design[[d]])) {
      design[[d]] <- data.frame(row.names = seq_len(n))
    } else {
      design[[d]] <- as.data.frame(design[[d]], stringsAsFactors = FALSE)
    }
    if (nrow(design[[d]]) != n) {
      stop("Domain `", d, "`: design has ", nrow(design[[d]]), " rows but ",
           "the data have ", n, ".", call. = FALSE)
    }
  }
  list(X = X, design = design, domain_names = domain_names)
}

.cvt_check_columns <- function(ctx, columns, name, domains) {
  for (col in columns) {
    for (d in domains) {
      if (!col %in% names(ctx$design[[d]])) {
        stop("Design column `", col, "` (", name, ") is missing from domain `",
             d, "`.", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

.cvt_target_type <- function(ctx, target) {
  if (is.null(target)) return(NA_character_)
  values <- lapply(ctx$design, `[[`, target)
  numeric <- vapply(values, function(v) is.numeric(v) && !is.factor(v),
                    logical(1))
  if (all(numeric)) return("numeric")
  if (any(numeric)) {
    stop("`target` column `", target, "` is numeric in some domains and ",
         "categorical in others.", call. = FALSE)
  }
  "categorical"
}

.cvt_method_fields <- c("fit", "args", "encode", "calibrate")

.cvt_resolve_methods <- function(methods) {
  if (is.function(methods) || inherits(methods, "aligner")) {
    methods <- list(method = methods)
  }
  if (!is.list(methods) || !length(methods)) {
    stop("`methods` must be a non-empty named list of method specifications.",
         call. = FALSE)
  }
  nms <- names(methods)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    stop("`methods` must have unique, non-empty names.", call. = FALSE)
  }
  out <- lapply(nms, function(nm) {
    spec <- methods[[nm]]
    if (is.function(spec) || inherits(spec, "aligner")) {
      spec <- list(fit = spec)
    }
    if (!is.list(spec)) {
      stop("Method `", nm, "` must be a function, an aligner, or a list.",
           call. = FALSE)
    }
    fields <- names(spec)
    if (length(spec) && (is.null(fields) || any(!nzchar(fields)))) {
      stop("Method `", nm, "`: specification elements must be named.",
           call. = FALSE)
    }
    unknown <- setdiff(fields, .cvt_method_fields)
    if (length(unknown)) {
      stop("Method `", nm, "`: unknown specification element(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    if (!(is.function(spec$fit) || inherits(spec$fit, "aligner"))) {
      stop("Method `", nm, "`: `fit` must be a function or an aligner.",
           call. = FALSE)
    }
    if (is.null(spec$args)) spec$args <- list()
    if (!is.list(spec$args)) {
      stop("Method `", nm, "`: `args` must be a list.", call. = FALSE)
    }
    for (hook in c("encode", "calibrate")) {
      if (!is.null(spec[[hook]]) && !is.function(spec[[hook]])) {
        stop("Method `", nm, "`: `", hook, "` must be a function.",
             call. = FALSE)
      }
    }
    spec
  })
  names(out) <- nms
  out
}

# -- Folds and group splits -----------------------------------------------------

# Connected groups of rows: rows sharing a non-missing value of any key column
# belong to the same group. Labels are propagated to the minimum row index
# until stable.
.cvt_row_groups <- function(design, keys) {
  n <- nrow(design)
  comp <- seq_len(n)
  key_values <- lapply(keys, function(k) {
    v <- design[[k]]
    if (is.factor(v)) v <- as.character(v)
    v
  })
  repeat {
    before <- comp
    for (v in key_values) {
      ok <- !is.na(v)
      if (any(ok)) {
        comp[ok] <- stats::ave(comp[ok], as.character(v[ok]), FUN = min)
      }
    }
    if (identical(before, comp)) break
  }
  match(comp, unique(comp))
}

.cvt_split_groups <- function(groups, fraction, domain) {
  ids <- unique(groups)
  if (length(ids) < 2L) {
    stop("Held-out domain `", domain, "` forms a single connected group; it ",
         "cannot be split into calibration and evaluation rows.",
         call. = FALSE)
  }
  sizes <- tabulate(groups, nbins = length(ids))
  ordered <- ids[sample.int(length(ids))]
  target_rows <- fraction * length(groups)
  cal <- integer(0)
  total <- 0
  for (g in ordered[-length(ordered)]) {
    if (length(cal) && total >= target_rows) break
    cal <- c(cal, g)
    total <- total + sizes[[g]]
  }
  sort(which(groups %in% cal))
}

.cvt_fold_fields <- c("heldout", "sources", "calibration")

.cvt_resolve_folds <- function(folds, ctx, event_id, group_by, fraction) {
  domains <- ctx$domain_names
  if (is.null(folds)) {
    folds <- lapply(domains, function(d) list(heldout = d))
  } else if (is.character(folds)) {
    folds <- lapply(folds, function(d) list(heldout = d))
  } else if (!is.list(folds) || !length(folds)) {
    stop("`folds` must be NULL, a character vector of domain names, or a ",
         "list of fold specifications.", call. = FALSE)
  }
  lapply(seq_along(folds), function(i) {
    spec <- folds[[i]]
    if (!is.list(spec)) {
      stop("Fold ", i, " must be a list with element `heldout`.",
           call. = FALSE)
    }
    unknown <- setdiff(names(spec), .cvt_fold_fields)
    if (length(unknown) || (length(spec) && is.null(names(spec)))) {
      stop("Fold ", i, ": unknown or unnamed element(s); allowed: ",
           paste(.cvt_fold_fields, collapse = ", "), ".", call. = FALSE)
    }
    held <- spec$heldout
    if (!is.character(held) || length(held) != 1L || !held %in% domains) {
      stop("Fold ", i, ": `heldout` must name one domain of `data`.",
           call. = FALSE)
    }
    sources <- if (is.null(spec$sources)) setdiff(domains, held) else
      spec$sources
    if (!is.character(sources) || !length(sources) || anyNA(sources) ||
        anyDuplicated(sources) || length(setdiff(sources, domains))) {
      stop("Fold ", i, ": `sources` must name distinct domains of `data`.",
           call. = FALSE)
    }
    if (held %in% sources) {
      stop("Fold ", i, ": the held-out domain `", held, "` cannot also be a ",
           "source.", call. = FALSE)
    }
    .cvt_check_columns(ctx, group_by, "group_by", held)
    design <- ctx$design[[held]]
    n <- nrow(design)
    groups <- .cvt_row_groups(design, c(event_id, group_by))
    if (is.null(spec$calibration)) {
      cal <- .cvt_split_groups(groups, fraction, held)
    } else {
      cal <- spec$calibration
      if (is.logical(cal)) {
        if (length(cal) != n || anyNA(cal)) {
          stop("Fold ", i, ": a logical `calibration` needs one non-missing ",
               "value per row of `", held, "`.", call. = FALSE)
        }
        cal <- which(cal)
      } else if (!is.numeric(cal) || anyNA(cal) ||
                 any(abs(cal - round(cal)) > 1e-8) || any(cal < 1) ||
                 any(cal > n) || anyDuplicated(cal)) {
        stop("Fold ", i, ": `calibration` must hold distinct row indices of ",
             "`", held, "`.", call. = FALSE)
      }
      cal <- sort(as.integer(round(cal)))
      if (!length(cal) || length(cal) == n) {
        stop("Fold ", i, ": `calibration` must leave both calibration and ",
             "evaluation rows.", call. = FALSE)
      }
      split <- intersect(groups[cal], groups[-cal])
      if (length(split)) {
        stop("Fold ", i, ": `calibration` splits ", length(split), " group(s) ",
             "of rows linked by `event_id`/`group_by` across calibration and ",
             "evaluation.", call. = FALSE)
      }
    }
    eval_rows <- setdiff(seq_len(n), cal)
    list(
      heldout = held,
      sources = sources,
      calibration = cal,
      evaluation = eval_rows,
      n_calibration_groups = length(unique(groups[cal])),
      n_evaluation_groups = length(unique(groups[eval_rows]))
    )
  })
}

# -- Per-fold evaluation ----------------------------------------------------------

.cvt_metric_names <- c("target_score", "target_baseline", "recon_self_r2",
                       "recon_transfer_r2", "recon_source_r2")
.cvt_count_names <- c("n_calibration", "n_evaluation", "n_target",
                      "n_transfer", "n_source_rows")

.cvt_empty_row <- function(fold) {
  row <- list(status = NA_character_, stage = NA_character_,
              message = NA_character_, warnings = NA_character_,
              target_metric = NA_character_)
  for (nm in .cvt_metric_names) row[[nm]] <- NA_real_
  for (nm in .cvt_count_names) row[[nm]] <- NA_integer_
  row$n_calibration <- length(fold$calibration)
  row$n_evaluation <- length(fold$evaluation)
  row
}

.cvt_has_method <- function(generic, fit) {
  ns <- asNamespace("manifoldalign")
  for (cls in class(fit)) {
    if (!is.null(utils::getS3method(generic, cls, optional = TRUE,
                                    envir = ns))) {
      return(TRUE)
    }
  }
  FALSE
}

.cvt_default_encode <- function(fit, x, domain) {
  oos_predict(fit, x, side = domain)
}

.cvt_default_calibrate <- function(fit, data, design, domain) {
  calibrate_domain(fit, data, design = design, domain_name = domain)
}

.cvt_evaluate <- function(spec, fold, ctx, keep_models) {
  row <- .cvt_empty_row(fold)
  warn <- new.env(parent = emptyenv())
  warn$messages <- character(0)
  models <- NULL
  finish <- function(status, stage = NA_character_, message = NA_character_) {
    row$status <- status
    row$stage <- stage
    row$message <- message
    if (length(warn$messages)) {
      row$warnings <- paste(unique(warn$messages), collapse = " | ")
    }
    list(row = row, models = models)
  }
  attempt <- function(expr) {
    withCallingHandlers(
      tryCatch(list(ok = TRUE, value = expr),
               error = function(e) list(ok = FALSE,
                                        message = conditionMessage(e))),
      warning = function(w) {
        warn$messages <- c(warn$messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }

  held <- fold$heldout
  sources <- fold$sources
  cal <- fold$calibration
  ev <- fold$evaluation
  X_h <- ctx$X[[held]]
  eval_events <- character(0)
  if (!is.null(ctx$event_id)) {
    e <- ctx$events[[held]][ev]
    eval_events <- unique(e[!is.na(e)])
  }
  # Source rows available for fitting the reference and the target readout.
  train_rows <- lapply(sources, function(s) {
    n <- nrow(ctx$X[[s]])
    if (identical(ctx$evaluation_groups, "novel")) {
      which(!(ctx$events[[s]] %in% eval_events))
    } else {
      seq_len(n)
    }
  })
  names(train_rows) <- sources
  row$n_source_rows <- sum(lengths(train_rows))
  if (any(lengths(train_rows) == 0L)) {
    return(finish("failed", "fit",
                  "a source domain has no rows outside the evaluation events"))
  }

  if (inherits(spec$fit, "aligner")) {
    caps <- aligner_capabilities(spec$fit)
    if (is.null(spec$calibrate) && !isTRUE(caps$supports_calibration)) {
      return(finish("unsupported", "capability",
                    "aligner does not declare `supports_calibration = TRUE`"))
    }
    if (is.null(spec$encode) && isFALSE(caps$supports_predict)) {
      return(finish("unsupported", "capability",
                    "aligner declares `supports_predict = FALSE`"))
    }
    if (isTRUE(caps$requires_auxiliary_at_predict)) {
      return(finish("unsupported", "capability",
                    "aligner requires auxiliary inputs at prediction"))
    }
  }

  src_data <- lapply(sources, function(s) {
    ctx$X[[s]][train_rows[[s]], , drop = FALSE]
  })
  src_design <- lapply(sources, function(s) {
    ctx$design[[s]][train_rows[[s]], , drop = FALSE]
  })
  names(src_data) <- names(src_design) <- sources

  # 1. Reference on source domains only.
  fitted <- attempt({
    if (inherits(spec$fit, "aligner")) {
      do.call(fit_many, c(list(spec$fit, src_data, design = src_design),
                          spec$args))
    } else {
      do.call(spec$fit, c(list(src_data, src_design), spec$args))
    }
  })
  if (!fitted$ok) return(finish("failed", "fit", fitted$message))
  fit_ref <- fitted$value

  calibrate <- spec$calibrate
  if (is.null(calibrate)) {
    if (!.cvt_has_method("calibrate_domain", fit_ref)) {
      return(finish("unsupported", "capability", paste0(
        "no calibrate_domain() method for class ",
        paste(class(fit_ref), collapse = "/"),
        "; the held-out domain cannot be added to a frozen fit")))
    }
    calibrate <- .cvt_default_calibrate
  }
  encode <- spec$encode
  if (is.null(encode)) {
    if (!.cvt_has_method("oos_predict", fit_ref)) {
      return(finish("unsupported", "capability", paste0(
        "no oos_predict() method for class ",
        paste(class(fit_ref), collapse = "/"))))
    }
    encode <- .cvt_default_encode
  }

  # 2. Calibrate the held-out domain from its calibration rows only. The
  # target column is withheld from calibration.
  cal_design <- ctx$design[[held]][cal, , drop = FALSE]
  if (!is.null(ctx$target)) {
    cal_design <- cal_design[, setdiff(names(cal_design), ctx$target),
                             drop = FALSE]
  }
  calibrated <- attempt(calibrate(fit_ref, X_h[cal, , drop = FALSE],
                                  cal_design, held))
  if (!calibrated$ok) return(finish("failed", "calibrate", calibrated$message))
  fit_cal <- calibrated$value

  # 3. Encode: held-out rows with the calibrated fit, source rows with the
  # frozen reference.
  encoded <- attempt({
    Z_cal <- .cvt_encode(encode, fit_cal, X_h[cal, , drop = FALSE], held)
    Z_eval <- .cvt_encode(encode, fit_cal, X_h[ev, , drop = FALSE], held)
    Z_src <- lapply(sources, function(s) {
      .cvt_encode(encode, fit_ref, ctx$X[[s]], s)
    })
    names(Z_src) <- sources
    k <- unique(c(ncol(Z_cal), ncol(Z_eval), vapply(Z_src, ncol, integer(1))))
    if (length(k) != 1L) {
      stop("Encodings of different domains have different latent dimensions.",
           call. = FALSE)
    }
    list(cal = Z_cal, eval = Z_eval, src = Z_src)
  })
  if (!encoded$ok) return(finish("failed", "encode", encoded$message))
  Z <- encoded$value

  # 4. Readouts on permitted partitions, scored on evaluation rows.
  scored <- attempt(.cvt_score(ctx, fold, Z, train_rows, eval_events))
  if (!scored$ok) return(finish("failed", "readout", scored$message))
  for (nm in names(scored$value$row)) row[[nm]] <- scored$value$row[[nm]]
  if (keep_models) {
    models <- list(reference = fit_ref, calibrated = fit_cal,
                   readouts = scored$value$readouts)
  }
  finish("ok")
}

.cvt_encode <- function(encode, fit, x, domain) {
  Z <- encode(fit, x, domain)
  if (inherits(Z, "Matrix") || is.data.frame(Z)) Z <- as.matrix(Z)
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("Encoding of domain `", domain, "` is not a numeric matrix.",
         call. = FALSE)
  }
  if (nrow(Z) != nrow(x)) {
    stop("Encoding of domain `", domain, "` returned ", nrow(Z), " rows for ",
         nrow(x), " input rows.", call. = FALSE)
  }
  if (any(!is.finite(Z))) {
    stop("Encoding of domain `", domain, "` contains non-finite values.",
         call. = FALSE)
  }
  Z
}

# -- Readouts -------------------------------------------------------------------

# Ridge least squares with an intercept: Y ~ Z. Y may be a (sparse) Matrix;
# centring Z alone suffices for the cross-product because Zc has zero column
# sums.
.cvt_linear_fit <- function(Z, Y, ridge) {
  if (nrow(Z) < 1L) stop("A readout has no training rows.", call. = FALSE)
  z_mean <- colMeans(Z)
  Zc <- sweep(Z, 2L, z_mean, "-")
  G <- crossprod(Zc)
  lambda <- ridge * mean(diag(G))
  if (!is.finite(lambda) || lambda <= 0) lambda <- ridge
  y_mean <- if (inherits(Y, "Matrix")) Matrix::colMeans(Y) else colMeans(Y)
  coef <- solve(G + diag(lambda, ncol(Z)), as.matrix(crossprod(Zc, Y)))
  list(z_mean = z_mean, y_mean = as.numeric(y_mean), coef = coef)
}

.cvt_linear_predict <- function(model, Z) {
  P <- sweep(Z, 2L, model$z_mean, "-") %*% model$coef
  sweep(P, 2L, model$y_mean, "+")
}

# Pooled sums for R-squared: residuals against predictions and against a
# baseline fitted on the permitted partition.
.cvt_r2_sums <- function(Y, P, baseline) {
  Y <- as.matrix(Y)
  c(sse = sum((Y - P)^2), sst = sum(sweep(Y, 2L, baseline, "-")^2))
}

.cvt_r2 <- function(sums) {
  if (is.null(sums) || !is.finite(sums[["sst"]]) || sums[["sst"]] <= 0) {
    return(NA_real_)
  }
  unname(1 - sums[["sse"]] / sums[["sst"]])
}

.cvt_event_means <- function(Z, events) {
  keep <- !is.na(events)
  if (!any(keep)) return(NULL)
  ev <- events[keep]
  rs <- rowsum(Z[keep, , drop = FALSE], ev, reorder = TRUE)
  rs / as.vector(table(ev)[rownames(rs)])
}

.cvt_score <- function(ctx, fold, Z, train_rows, eval_events) {
  held <- fold$heldout
  sources <- fold$sources
  cal <- fold$calibration
  ev <- fold$evaluation
  ridge <- ctx$readout_ridge
  X_h <- ctx$X[[held]]
  out <- list()
  readouts <- list()

  # (a) Target readout: fitted on encoded source training rows only.
  if (!is.null(ctx$target)) {
    z_train <- do.call(rbind, lapply(sources, function(s) {
      Z$src[[s]][train_rows[[s]], , drop = FALSE]
    }))
    y_train <- unlist(lapply(sources, function(s) {
      v <- ctx$design[[s]][[ctx$target]][train_rows[[s]]]
      if (is.factor(v)) as.character(v) else v
    }), use.names = FALSE)
    y_eval <- ctx$design[[held]][[ctx$target]][ev]
    if (is.factor(y_eval)) y_eval <- as.character(y_eval)
    keep_train <- !is.na(y_train)
    keep_eval <- !is.na(y_eval)
    z_train <- z_train[keep_train, , drop = FALSE]
    y_train <- y_train[keep_train]
    z_test <- Z$eval[keep_eval, , drop = FALSE]
    y_test <- y_eval[keep_eval]
    out$n_target <- length(y_test)
    numeric_target <- identical(ctx$target_type, "numeric")
    out$target_metric <- if (numeric_target) "r2" else "accuracy"
    if (!length(y_train)) stop("No source rows have a target value.",
                               call. = FALSE)
    if (length(y_test)) {
      if (!is.null(ctx$target_readout)) {
        pred <- ctx$target_readout(z_train, y_train, z_test)
        if (length(pred) != nrow(z_test)) {
          stop("`target_readout` returned ", length(pred), " predictions for ",
               nrow(z_test), " rows.", call. = FALSE)
        }
        readouts$target <- "custom"
      } else if (numeric_target) {
        model <- .cvt_linear_fit(z_train, matrix(y_train, ncol = 1L), ridge)
        pred <- as.numeric(.cvt_linear_predict(model, z_test))
        readouts$target <- model
      } else {
        classes <- sort(unique(as.character(y_train)), method = "radix")
        onehot <- outer(as.character(y_train), classes, "==") * 1
        model <- .cvt_linear_fit(z_train, onehot, ridge)
        scores <- .cvt_linear_predict(model, z_test)
        pred <- classes[max.col(scores, ties.method = "first")]
        model$classes <- classes
        readouts$target <- model
      }
      if (numeric_target) {
        pred <- as.numeric(pred)
        sse <- sum((y_test - pred)^2)
        sst <- sum((y_test - mean(y_train))^2)
        out$target_score <- if (sst > 0) 1 - sse / sst else NA_real_
        out$target_baseline <- 0
      } else {
        out$target_score <- mean(as.character(pred) == as.character(y_test))
        tab <- table(as.character(y_train))
        majority <- names(tab)[which.max(tab)]
        out$target_baseline <- mean(as.character(y_test) == majority)
      }
    }
  }

  # (b) Held-out reconstruction readout D_h, fitted on calibration rows only.
  D_h <- .cvt_linear_fit(Z$cal, X_h[cal, , drop = FALSE], ridge)
  readouts$heldout <- D_h
  X_eval <- as.matrix(X_h[ev, , drop = FALSE])
  out$recon_self_r2 <- .cvt_r2(.cvt_r2_sums(
    X_eval, .cvt_linear_predict(D_h, Z$eval), D_h$y_mean))

  if (!is.null(ctx$event_id)) {
    eval_ev <- ctx$events[[held]][ev]

    # Sources -> held-out: consensus source encoding of each evaluation event.
    per_source <- lapply(sources, function(s) {
      .cvt_event_means(Z$src[[s]], ctx$events[[s]])
    })
    all_ev <- unique(unlist(lapply(per_source, rownames)))
    sums <- matrix(0, length(all_ev), ncol(Z$eval),
                   dimnames = list(all_ev, NULL))
    counts <- stats::setNames(numeric(length(all_ev)), all_ev)
    for (M in per_source) {
      if (is.null(M)) next
      sums[rownames(M), ] <- sums[rownames(M), , drop = FALSE] + M
      counts[rownames(M)] <- counts[rownames(M)] + 1
    }
    has <- !is.na(eval_ev) & eval_ev %in% all_ev
    out$n_transfer <- sum(has)
    if (any(has)) {
      key <- eval_ev[has]
      consensus <- sums[key, , drop = FALSE] / counts[key]
      out$recon_transfer_r2 <- .cvt_r2(.cvt_r2_sums(
        X_eval[has, , drop = FALSE], .cvt_linear_predict(D_h, consensus),
        D_h$y_mean))
    }

    # Held-out -> sources: D_s fitted on each source's rows outside the
    # evaluation events, applied to the held-out mean evaluation encoding.
    heldout_means <- .cvt_event_means(Z$eval, eval_ev)
    total <- NULL
    readouts$sources <- list()
    for (s in sources) {
      ev_s <- ctx$events[[s]]
      fit_rows <- which(!(ev_s %in% eval_events))
      test_rows <- which(!is.na(ev_s) & ev_s %in% rownames(heldout_means))
      if (!length(fit_rows) || !length(test_rows)) next
      D_s <- .cvt_linear_fit(Z$src[[s]][fit_rows, , drop = FALSE],
                             ctx$X[[s]][fit_rows, , drop = FALSE], ridge)
      readouts$sources[[s]] <- D_s
      P <- .cvt_linear_predict(D_s, heldout_means[ev_s[test_rows], ,
                                                  drop = FALSE])
      part <- .cvt_r2_sums(ctx$X[[s]][test_rows, , drop = FALSE], P,
                           D_s$y_mean)
      total <- if (is.null(total)) part else total + part
    }
    out$recon_source_r2 <- .cvt_r2(total)
  }
  list(row = out, readouts = readouts)
}

.cvt_rows_to_frame <- function(rows) {
  cols <- c(".fold", "heldout", "method", "status", "stage", "message",
            "warnings", "target_metric", .cvt_metric_names, .cvt_count_names,
            ".seed")
  df <- lapply(cols, function(col) {
    vals <- lapply(rows, `[[`, col)
    unlist(lapply(vals, function(v) if (is.null(v)) NA else v),
           use.names = FALSE)
  })
  names(df) <- cols
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (col in c(.cvt_count_names, ".fold", ".seed")) {
    df[[col]] <- as.integer(df[[col]])
  }
  for (col in .cvt_metric_names) df[[col]] <- as.numeric(df[[col]])
  for (col in c("heldout", "method", "status", "stage", "message", "warnings",
                "target_metric")) {
    df[[col]] <- as.character(df[[col]])
  }
  df
}

# -- Result methods ---------------------------------------------------------------

#' @method as.data.frame cv_transfer_result
#' @export
as.data.frame.cv_transfer_result <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
  df <- x$scores
  if (!is.null(row.names)) rownames(df) <- row.names
  df
}

#' @method print cv_transfer_result
#' @export
print.cv_transfer_result <- function(x, ...) {
  scores <- x$scores
  cat("\n=== Held-out Domain Transfer Result ===\n")
  cat("\nFolds: ", length(x$folds), " (held out: ",
      paste(vapply(x$folds, `[[`, character(1), "heldout"), collapse = ", "),
      ")\n", sep = "")
  cat("Evaluation groups: ", x$settings$evaluation_groups, "\n", sep = "")
  metrics <- c("target_score", "recon_self_r2", "recon_transfer_r2",
               "recon_source_r2")
  for (m in unique(scores$method)) {
    s <- scores[scores$method == m, , drop = FALSE]
    n <- nrow(s)
    ok <- s$status == "ok"
    cat("\nMethod ", m, ": ", sum(ok), " ok, ",
        sum(s$status == "unsupported"), " unsupported, ",
        sum(s$status == "failed"), " failed (of ", n, " folds)\n", sep = "")
    if (any(ok)) {
      for (metric in metrics) {
        vals <- s[[metric]][ok]
        if (all(is.na(vals))) next
        label <- if (identical(metric, "target_score")) {
          paste0(metric, " (", s$target_metric[ok][1L], ")")
        } else metric
        cat("  ", label, ": mean = ", format(round(mean(vals, na.rm = TRUE), 4)),
            " over ", sum(!is.na(vals)), " of ", n, " folds\n", sep = "")
      }
    }
    bad <- s[!ok, , drop = FALSE]
    for (i in seq_len(nrow(bad))) {
      cat("  fold ", bad$.fold[i], " (", bad$heldout[i], "): ", bad$status[i],
          " at ", bad$stage[i], ": ", bad$message[i], "\n", sep = "")
    }
    n_warn <- sum(!is.na(s$warnings))
    if (n_warn) cat("  ", n_warn, " fold(s) raised warnings; see `warnings`.\n",
                    sep = "")
  }
  cat("\n=======================================\n")
  invisible(x)
}

# -- Ablation helpers -------------------------------------------------------------

#' Event-similarity kernel from event ordering
#'
#' Builds a low-rank factor of a Gaussian similarity kernel over event
#' positions, `k(i, j) = exp(-(o_i - o_j)^2 / (2 h^2))`, optionally block
#' diagonal so events in different blocks (for example runs) have zero
#' similarity. The factor `F` satisfies `F F' = K` up to the eigenvalues
#' dropped by `rank` and `tol`. It is an ablation control: an event metric
#' that knows only *when* events occurred, not *what* they were. It can be
#' passed as one family of `event_features` to [sumcor_align()]; use
#' `sumcor_control(kernel_scale = FALSE)` there so the kernel's spectrum is
#' not flattened by column scaling.
#'
#' @param event_ids Unique event identifiers.
#' @param order Numeric event positions, one per event. Defaults to
#'   `seq_along(event_ids)`.
#' @param bandwidth Positive kernel bandwidth `h`, in units of `order`.
#' @param blocks Optional block labels, one per event; the kernel is zero
#'   between blocks. `NA` is not allowed.
#' @param rank Optional maximum number of components kept per block.
#' @param tol Relative eigenvalue threshold below which components are dropped.
#' @return A numeric matrix with one row per event (row names `event_ids`, in
#'   input order) and columns `order_1`, `order_2`, ...; attribute
#'   `eigenvalues` holds the kept eigenvalues.
#' @examples
#' F <- event_order_kernel(paste0("e", 1:20), bandwidth = 2,
#'                         blocks = rep(1:2, each = 10))
#' K <- tcrossprod(F)
#' round(K[1:3, 1:3], 3)
#' K[1, 11] # different blocks
#' @export
event_order_kernel <- function(event_ids, order = NULL, bandwidth = 1,
                               blocks = NULL, rank = NULL, tol = 1e-10) {
  event_ids <- .cvt_event_ids(event_ids)
  n <- length(event_ids)
  if (is.null(order)) order <- seq_len(n)
  if (!is.numeric(order) || length(order) != n || any(!is.finite(order))) {
    stop("`order` must be finite numbers, one per event.", call. = FALSE)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a single positive number.", call. = FALSE)
  }
  if (!is.null(rank) && (!is.numeric(rank) || length(rank) != 1L ||
                         !is.finite(rank) || rank < 1 ||
                         abs(rank - round(rank)) > 1e-8)) {
    stop("`rank` must be NULL or a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0 ||
      tol >= 1) {
    stop("`tol` must be a single number in [0, 1).", call. = FALSE)
  }
  block_idx <- .cvt_blocks(blocks, n)
  parts <- list()
  values <- numeric(0)
  for (b in block_idx) {
    o <- order[b]
    K <- exp(-outer(o, o, "-")^2 / (2 * bandwidth^2))
    eg <- eigen(K, symmetric = TRUE)
    keep <- which(eg$values > tol * max(eg$values))
    if (!is.null(rank)) keep <- utils::head(keep, as.integer(round(rank)))
    V <- eg$vectors[, keep, drop = FALSE]
    # Deterministic signs: largest-magnitude entry of each vector positive.
    pivot <- apply(abs(V), 2L, which.max)
    flip <- sign(V[cbind(pivot, seq_along(keep))])
    flip[flip == 0] <- 1
    V <- sweep(V, 2L, flip * sqrt(eg$values[keep]), "*")
    block_F <- matrix(0, n, length(keep))
    block_F[b, ] <- V
    parts[[length(parts) + 1L]] <- block_F
    values <- c(values, eg$values[keep])
  }
  out <- do.call(cbind, parts)
  dimnames(out) <- list(event_ids, paste0("order_", seq_len(ncol(out))))
  attr(out, "eigenvalues") <- values
  attr(out, "bandwidth") <- bandwidth
  out
}

#' Surrogate event features that keep autocorrelation but break event identity
#'
#' Creates an ablation control for event-feature metrics. Within each block,
#' the feature sequence (rows taken in event order) is replaced by a surrogate
#' with the same autocorrelation, and the result is written back to the
#' original rows, so each event identifier now carries another position's
#' features. `"circular_shift"` rotates the whole sequence by one random
#' offset (identical for all columns), which preserves every column's circular
#' autocorrelation and all cross-column relations exactly.
#' `"phase_randomize"` multiplies the discrete Fourier transform by random
#' phases shared by all columns, preserving each column's mean and periodogram
#' (hence its circular autocorrelation) and all cross-spectra exactly.
#'
#' @param features Numeric matrix with one row per event (row names are event
#'   identifiers and are kept).
#' @param method `"circular_shift"` or `"phase_randomize"`.
#' @param order Optional numeric event positions, one per row; rows are
#'   processed in this order. Defaults to the row order.
#' @param blocks Optional block labels, one per row; surrogates are formed
#'   within each block separately.
#' @param min_shift Smallest circular offset allowed (and smallest distance
#'   from a full rotation); applies to `"circular_shift"`.
#' @param seed Optional integer seed. When given, the caller's random number
#'   stream is restored afterwards.
#' @return A matrix with the dimensions and dimnames of `features`; attribute
#'   `surrogate` records the method and, for circular shifts, each block's
#'   offset.
#' @examples
#' feats <- matrix(cumsum(rnorm(40)), 20, 2,
#'                 dimnames = list(paste0("e", 1:20), c("f1", "f2")))
#' shifted <- surrogate_event_features(feats, "circular_shift", seed = 1)
#' attr(shifted, "surrogate")$shifts
#' phased <- surrogate_event_features(feats, "phase_randomize", seed = 1)
#' all.equal(Mod(stats::fft(phased[, 1])), Mod(stats::fft(feats[, 1])))
#' @export
surrogate_event_features <- function(features,
                                     method = c("circular_shift",
                                                "phase_randomize"),
                                     order = NULL, blocks = NULL,
                                     min_shift = 1L, seed = NULL) {
  method <- match.arg(method)
  if (is.data.frame(features)) features <- as.matrix(features)
  if (!is.matrix(features) || !is.numeric(features) || nrow(features) < 2L ||
      any(!is.finite(features))) {
    stop("`features` must be a finite numeric matrix with at least two rows.",
         call. = FALSE)
  }
  n <- nrow(features)
  if (is.null(order)) order <- seq_len(n)
  if (!is.numeric(order) || length(order) != n || any(!is.finite(order))) {
    stop("`order` must be finite numbers, one per row.", call. = FALSE)
  }
  if (!is.numeric(min_shift) || length(min_shift) != 1L ||
      !is.finite(min_shift) || min_shift < 1 ||
      abs(min_shift - round(min_shift)) > 1e-8) {
    stop("`min_shift` must be a single positive integer.", call. = FALSE)
  }
  min_shift <- as.integer(round(min_shift))
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                         !is.finite(seed))) {
    stop("`seed` must be NULL or a single integer.", call. = FALSE)
  }
  block_idx <- .cvt_blocks(blocks, n)
  run <- function() {
    out <- features
    shifts <- integer(0)
    for (b in block_idx) {
      rows <- b[base::order(order[b], method = "radix")]
      m <- length(rows)
      X <- features[rows, , drop = FALSE]
      if (identical(method, "circular_shift")) {
        lo <- min_shift
        hi <- m - min_shift
        if (hi < lo) {
          stop("A block of ", m, " events is too short for `min_shift = ",
               min_shift, "`.", call. = FALSE)
        }
        k <- if (hi == lo) lo else lo - 1L + sample.int(hi - lo + 1L, 1L)
        shifts <- c(shifts, k)
        out[rows, ] <- X[((seq_len(m) - 1L + k) %% m) + 1L, , drop = FALSE]
      } else {
        if (m < 3L) {
          stop("Phase randomisation needs at least three events per block.",
               call. = FALSE)
        }
        half <- (m - 1L) %/% 2L
        phi <- stats::runif(half, 0, 2 * pi)
        rot <- complex(m, real = 1)
        rot[1L + seq_len(half)] <- exp(1i * phi)
        rot[m + 1L - seq_len(half)] <- exp(-1i * phi)
        FX <- stats::mvfft(X) * rot
        out[rows, ] <- Re(stats::mvfft(FX, inverse = TRUE)) / m
      }
    }
    attr(out, "surrogate") <- list(
      method = method,
      shifts = if (identical(method, "circular_shift")) shifts else NULL
    )
    out
  }
  if (is.null(seed)) run() else .with_local_seed(as.integer(seed), run())
}

.cvt_event_ids <- function(event_ids) {
  if (is.factor(event_ids)) event_ids <- as.character(event_ids)
  if (!is.character(event_ids) || length(event_ids) < 1L ||
      anyNA(event_ids) || any(!nzchar(event_ids)) || anyDuplicated(event_ids)) {
    stop("`event_ids` must be unique, non-empty identifiers.", call. = FALSE)
  }
  event_ids
}

.cvt_blocks <- function(blocks, n) {
  if (is.null(blocks)) return(list(seq_len(n)))
  if (length(blocks) != n || anyNA(blocks)) {
    stop("`blocks` must have one non-missing label per event.", call. = FALSE)
  }
  labels <- as.character(blocks)
  lapply(unique(labels), function(l) which(labels == l))
}
