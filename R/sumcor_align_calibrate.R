# Event-keyed SUMCOR: prediction, frozen calibration, adapter -----------------

.sumcor_domain_index <- function(fit, side) {
  if (length(side) != 1L || is.na(side)) {
    stop("`side` must identify exactly one domain.", call. = FALSE)
  }
  if (is.character(side)) {
    idx <- match(side, fit$domain_names)
    if (is.na(idx)) stop("Unknown domain `", side, "`.", call. = FALSE)
    return(idx)
  }
  if (!is.numeric(side) || abs(side - round(side)) > 1e-8 ||
      side < 1 || side > length(fit$domain_names)) {
    stop("`side` must be a domain name or index.", call. = FALSE)
  }
  as.integer(side)
}

#' Encode new observations with a SUMCOR encoder
#'
#' Applies the fitted encoder `f_s(x) = phi_s(x) A_s` of one domain to new
#' rows. Training scores were produced by the same function, so re-encoding a
#' domain's training observations reproduces `fit$scores` exactly. Rows are
#' encoded independently: the result for a row does not depend on the other
#' rows in `newX`.
#'
#' @param fit_or_transform A `sumcor_align` fit.
#' @param newX New observations (rows) with the same columns as the domain's
#'   training data. Column names, when both are present, must match.
#' @param side Domain name or index.
#' @param ... Unused; unknown arguments are rejected.
#' @return A numeric matrix with one row per row of `newX` and `ncomp` columns.
#' @export
oos_predict.sumcor_align <- function(fit_or_transform, newX, side, ...) {
  .sumcor_reject_dots(...)
  fit <- fit_or_transform
  m <- .sumcor_domain_index(fit, side)
  enc <- fit$encoders[[m]]
  X <- .sumcor_as_numeric_matrix(newX, "newX")
  if (ncol(X) != enc$n_columns) {
    stop("`newX` has ", ncol(X), " columns; domain `", fit$domain_names[[m]],
         "` was fitted with ", enc$n_columns, ".", call. = FALSE)
  }
  if (!is.null(colnames(X)) && !is.null(enc$column_names) &&
      !identical(colnames(X), enc$column_names)) {
    stop("`newX` column names do not match those of domain `",
         fit$domain_names[[m]], "`.", call. = FALSE)
  }
  Z <- .sumcor_encode(enc, X)
  rownames(Z) <- rownames(X)
  Z
}

#' Calibrate a new domain against a frozen alignment reference
#'
#' Adds a domain to a fitted alignment without changing any existing encoder
#' or reference score. The new domain's preprocessing and encoder are fitted
#' from its own calibration data only.
#'
#' @param fit A fitted alignment object.
#' @param data The new domain's calibration data.
#' @param ... Method-specific arguments.
#' @return The fit with the new domain added.
#' @export
calibrate_domain <- function(fit, data, ...) UseMethod("calibrate_domain")

#' @export
calibrate_domain.default <- function(fit, data, ...) {
  stop("calibrate_domain() is not implemented for objects of class ",
       paste(class(fit), collapse = ", "), ".", call. = FALSE)
}

#' Calibrate a new domain for SUMCOR alignment
#'
#' Fits the new domain's basis and repeat-level covariance from its own
#' calibration observations, then solves its encoder in one polar step against
#' the frozen event scores of the reference domains: with the reference fixed
#' this block subproblem is solved exactly. Existing encoders, scores and
#' `reference_scores` are left unchanged. Events the fit never saw are
#' excluded and counted. The identifier columns of the original fit are used.
#'
#' Calibrating a training domain's own data against the other training
#' domains reproduces that domain's fitted encoder whenever the fit converged,
#' which makes this a direct stationarity check.
#'
#' @param fit A `sumcor_align` fit.
#' @param data The new domain's observation-by-column matrix.
#' @param design A data frame with the identifier columns used by the fit.
#' @param domain_name Name for the new domain; must not already exist.
#' @param reference_domains Fitted domains that form the frozen reference.
#'   Defaults to every domain with role `"reference"`.
#' @param feature_groups Column groups for the new domain when the fit used
#'   `basis = "grouped_pca"`.
#' @param ... Unused; unknown arguments are rejected.
#' @return The `sumcor_align` fit with the new domain appended and its role
#'   recorded as `"calibrated"`. Calibration provenance is stored in
#'   `fit$calibration[[domain_name]]`.
#' @export
calibrate_domain.sumcor_align <- function(
  fit,
  data,
  design,
  domain_name,
  reference_domains = NULL,
  feature_groups = NULL,
  ...
) {
  .sumcor_reject_dots(...)
  domain_name <- .sumcor_name(domain_name, "domain_name", allow_null = FALSE)
  if (domain_name %in% fit$domain_names) {
    stop("Domain `", domain_name, "` already exists in the fit.", call. = FALSE)
  }
  if (is.null(reference_domains)) {
    reference_domains <- fit$domain_names[fit$domain_roles == "reference"]
  }
  reference_domains <- as.character(reference_domains)
  if (!length(reference_domains) ||
      length(setdiff(reference_domains, fit$domain_names))) {
    stop("`reference_domains` must name fitted domains.", call. = FALSE)
  }
  control <- fit$control
  dom <- .sumcor_resolve_domain(data, design, fit$ids, domain_name)

  universe <- fit$events$event_id
  obs <- dom$observations
  known <- obs$used & obs$event_id %in% universe
  unknown_events <- unique(obs$event_id[obs$used & !known])
  if (length(unknown_events)) {
    warning(length(unknown_events), " event(s) of domain `", domain_name,
            "` are not in the fitted event set and were excluded.",
            call. = FALSE)
  }
  obs$used <- known
  dom$observations <- obs
  if (sum(known) < 2L) {
    stop("Domain `", domain_name, "` has fewer than two observations of ",
         "fitted events.", call. = FALSE)
  }

  w <- fit$events$weight
  prep <- .sumcor_prepare_domain(dom, fit$events, w, control, feature_groups,
                                 domain_name)
  if (fit$ncomp > prep$p) {
    stop("Domain `", domain_name, "` has basis rank ", prep$p, ", below ",
         "`ncomp = ", fit$ncomp, "`.", call. = FALSE)
  }

  overlap <- stats::setNames(numeric(length(reference_domains)),
                             reference_domains)
  acc <- NULL
  for (ref in reference_domains) {
    Z <- fit$event_scores[[ref]]
    mask_t <- !is.na(Z[, 1L])
    m_st <- sum(w[prep$mask & mask_t])
    overlap[[ref]] <- m_st
    if (m_st <= 0 || m_st < control$min_overlap) next
    Z[!mask_t, ] <- 0
    term <- (Z * sqrt(w)) / m_st
    acc <- if (is.null(acc)) term else acc + term
  }
  if (is.null(acc)) {
    stop("Domain `", domain_name, "` shares less than `min_overlap` event ",
         "weight with every reference domain.", call. = FALSE)
  }
  G <- crossprod(prep$Ut, .sumcor_apply_H(fit$kernel$factor, unname(acc)))
  pol <- .sumcor_polar(G, NULL, control$rank_tol)
  A <- prep$C_inv_half %*% pol$Q
  if (pol$rank < fit$ncomp) {
    warning("Only ", pol$rank, " of ", fit$ncomp, " latent directions are ",
            "supported for calibrated domain `", domain_name, "`.",
            call. = FALSE)
  }

  enc <- .sumcor_encoder(prep, A, colnames(dom$X), control)
  scores <- .sumcor_encode(enc, dom$X)
  rownames(scores) <- obs$observation_id

  fit$domain_names <- c(fit$domain_names, domain_name)
  fit$domain_roles[[domain_name]] <- "calibrated"
  fit$scores[[domain_name]] <- scores
  fit$event_scores[[domain_name]] <- .sumcor_event_scores(prep, A, universe)
  fit$encoders[[domain_name]] <- enc
  fit$observations[[domain_name]] <- obs
  fit$supported_rank[[domain_name]] <- pol$rank
  fit$calibration[[domain_name]] <- list(
    reference_domains = reference_domains,
    overlap = overlap,
    events_used = sum(prep$mask),
    observations_used = sum(known),
    unknown_events = unknown_events,
    supported_rank = pol$rank,
    singular_values = pol$singular_values,
    objective_contribution = sum(diag(crossprod(G, pol$Q)))
  )
  fit
}

#' Construct a SUMCOR aligner descriptor
#'
#' The aligner is a native multi-domain encoder: [align_many()] and
#' [fit_many()] call [sumcor_align()] directly and never synchronise pairwise
#' transforms. Its capability record declares a linear encoder with
#' prediction and new-domain calibration, no auxiliary inputs at prediction,
#' and no inverse.
#'
#' @return An object of class `c("sumcor_aligner", "aligner")`.
#' @export
sumcor_aligner <- function() {
  new_aligner(
    "sumcor",
    group = "encoder",
    supports_multi = TRUE,
    encoder = "linear",
    supports_predict = TRUE,
    supports_calibration = TRUE,
    requires_auxiliary_at_predict = FALSE,
    invertible = FALSE
  )
}

#' @export
fit_many.sumcor_aligner <- function(algo, domains, k = NULL, ncomp = NULL,
                                    ...) {
  if (is.null(ncomp)) ncomp <- if (is.null(k)) 10L else k
  sumcor_align(domains, ncomp = ncomp, ...)
}
