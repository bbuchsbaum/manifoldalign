#' Generate a Synthetic Two-Domain Spiral Fixture
#'
#' Creates a deterministic two-domain spiral fixture inspired by the qualitative
#' examples used in manifold-alignment papers. It is not a digitized or exact
#' reproduction of a published KEMA benchmark.
#'
#' @param n_per_domain Number of samples per domain.
#' @param noise_level Gaussian noise standard deviation.
#' @param seed Random seed for reproducibility.
#' @return A list with domain data, labels, and combined strata.
#' @examples
#' data <- generate_spiral_validation_data(n_per_domain = 50, seed = 123)
#' str(data)
#' @export
generate_spiral_validation_data <- function(n_per_domain = 100,
                                            noise_level = 0.1,
                                            seed = 42) {
  set.seed(seed)

  t <- seq(0, 4 * pi, length.out = n_per_domain)
  x1 <- cbind(
    t * cos(t) + stats::rnorm(n_per_domain, 0, noise_level),
    t * sin(t) + stats::rnorm(n_per_domain, 0, noise_level)
  )
  x2 <- cbind(
    0.8 * t * cos(t + pi / 4) + stats::rnorm(n_per_domain, 0, noise_level),
    0.8 * t * sin(t + pi / 4) + stats::rnorm(n_per_domain, 0, noise_level)
  )
  labels1 <- ifelse(t < 2 * pi, "early", "late")
  labels2 <- ifelse(t < 2 * pi, "early", "late")

  strata <- list(
    list(x = x1, design = data.frame(labels = labels1)),
    list(x = x2, design = data.frame(labels = labels2))
  )

  list(
    strata = strata,
    labels = c(labels1, labels2),
    domain1 = strata[[1L]],
    domain2 = strata[[2L]]
  )
}

#' Withdrawn KEMA Paper-Eigenvalue Validator
#'
#' This compatibility shim is retained so historical callers receive an
#' explicit failure. The former routine did not reproduce a documented paper
#' fixture and solved a different generalized problem whose right-hand side was
#' only `lambda * I`; its claimed comparison with paper eigenvalues was therefore
#' not valid accuracy evidence.
#'
#' @param expected_eigenvals Deprecated and ignored.
#' @param tolerance Deprecated and ignored.
#' @param n_per_domain Deprecated and ignored.
#' @return This function emits a deprecation warning and then errors.
#' @seealso [run_kema_validation_suite()]
#' @examples
#' \dontrun{
#' validate_kema_eigenvalues()
#' }
#' @export
validate_kema_eigenvalues <- function(expected_eigenvals = c(0.82, 0.41),
                                      tolerance = 0.1,
                                      n_per_domain = 100) {
  .Deprecated(
    msg = paste0(
      "validate_kema_eigenvalues() is deprecated and withdrawn because it ",
      "did not reproduce a documented paper fixture or the implemented KEMA ",
      "generalized eigenproblem. Use run_kema_validation_suite()."
    )
  )
  stop(
    "Historical KEMA paper-eigenvalue validation is invalid. Use the ",
    "fidelity-gated residual and backend-agreement checks in ",
    "run_kema_validation_suite().",
    call. = FALSE
  )
}

#' Withdrawn KEMA Out-of-Sample Reconstruction Validator
#'
#' This compatibility shim is retained temporarily so existing callers receive
#' an explicit failure instead of a fabricated validation score. The former
#' implementation did not fit KEMA or evaluate held-out predictions; it returned
#' `expected_error` plus random noise and has therefore been withdrawn.
#'
#' @param expected_error Deprecated and ignored.
#' @param tolerance Deprecated and ignored.
#' @param test_fraction Deprecated and ignored.
#' @return This function emits a deprecation warning and then errors.
#' @details
#' No direct reconstruction replacement is currently available because KEMA's
#' out-of-sample projection has not passed an independent reconstruction oracle.
#' Use [cv_alignment_rows()] for leakage-safe held-out alignment scoring.
#' Historical results from this function are invalid.
#' @seealso [cv_alignment_rows()]
#' @examples
#' \dontrun{
#' validate_out_of_sample_reconstruction()
#' }
#' @export
validate_out_of_sample_reconstruction <- function(expected_error = 0.14,
                                                  tolerance = 0.05,
                                                  test_fraction = 0.2) {
  .Deprecated(
    msg = paste0(
      "validate_out_of_sample_reconstruction() is deprecated and withdrawn ",
      "because it did not fit KEMA or evaluate held-out predictions. ",
      "Use cv_alignment_rows() for leakage-safe held-out alignment scoring."
    )
  )
  stop(
    "KEMA out-of-sample reconstruction validation is withdrawn until a ",
    "verified OOS predictor and independent reconstruction oracle are available.",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.kema_validation_hyperdesign <- function(data) {
  multidesign::hyperdesign(list(
    domain1 = multidesign::multidesign(
      data$domain1$x,
      data.frame(labels = factor(data$domain1$design$labels))
    ),
    domain2 = multidesign::multidesign(
      data$domain2$x,
      data.frame(labels = factor(data$domain2$design$labels))
    )
  ))
}

#' Run Enforced KEMA Numerical Validation
#'
#' Fits deterministic KEMA fixtures and records the same numerical contracts
#' required of every returned fit: generalized-eigen backward residual,
#' B-orthogonality, numerical-zero exclusion, and agreement between matrix and
#' operator backends up to subspace indeterminacy.
#'
#' @param verbose Whether to print a concise report.
#' @return A list containing `fidelity_validation`, `backend_agreement`,
#'   `withdrawn_validations`, and `overall_success`.
#' @details
#' Withdrawn paper-eigenvalue, regression-solver, and fabricated reconstruction
#' checks are reported for auditability but do not contribute to
#' `overall_success`.
#' @examples
#' \donttest{
#' results <- run_kema_validation_suite(verbose = FALSE)
#' results$overall_success
#' }
#' @export
run_kema_validation_suite <- function(verbose = TRUE) {
  if (verbose) {
    cat("Running enforced KEMA numerical validation\n")
  }

  fidelity_validation <- tryCatch({
    data <- generate_spiral_validation_data(
      n_per_domain = 30,
      noise_level = 0.05,
      seed = 42
    )
    hd <- .kema_validation_hyperdesign(data)
    fit <- kema(
      hd,
      y = labels,
      ncomp = 2,
      knn = 5,
      kernel = kernlab::rbfdot(sigma = 0.5),
      lambda = 1e-3,
      backend = "full_exact"
    )

    success <- isTRUE(fit$fidelity$passed) &&
      fit$fidelity$max_rel_residual <= 1e-6 &&
      fit$fidelity$max_b_orth_offdiag <= 1e-6 &&
      all(fit$eigenvalues$values > fit$fidelity$eigenvalue_zero_tol)
    list(
      success = success,
      backend = fit$backend,
      eigensolver = fit$fidelity$eigensolver,
      max_rel_residual = fit$fidelity$max_rel_residual,
      max_b_orth_offdiag = fit$fidelity$max_b_orth_offdiag,
      eigenvalue_zero_tol = fit$fidelity$eigenvalue_zero_tol,
      eigenvalues = fit$eigenvalues$values
    )
  }, error = function(e) {
    list(success = FALSE, message = conditionMessage(e))
  })

  backend_agreement <- tryCatch({
    data <- generate_spiral_validation_data(
      n_per_domain = 20,
      noise_level = 0.05,
      seed = 91
    )
    hd <- .kema_validation_hyperdesign(data)
    fit_backend <- function(backend) {
      kema(
        hd,
        y = labels,
        ncomp = 2,
        knn = 4,
        kernel = kernlab::rbfdot(sigma = 0.5),
        lambda = 1e-3,
        backend = backend
      )
    }
    full <- fit_backend("full_exact")
    operator <- fit_backend("operator_exact")

    Q_full <- qr.Q(qr(as.matrix(full$s)))
    Q_operator <- qr.Q(qr(as.matrix(operator$s)))
    correlations <- svd(
      crossprod(Q_full, Q_operator),
      nu = 0,
      nv = 0
    )$d
    eigenvalue_error <- max(abs(
      full$eigenvalues$values - operator$eigenvalues$values
    ))
    success <- isTRUE(full$fidelity$passed) &&
      isTRUE(operator$fidelity$passed) &&
      min(correlations) >= 1 - 1e-6 &&
      eigenvalue_error <= 1e-6
    list(
      success = success,
      canonical_correlations = correlations,
      max_eigenvalue_error = eigenvalue_error,
      full_max_rel_residual = full$fidelity$max_rel_residual,
      operator_max_rel_residual = operator$fidelity$max_rel_residual
    )
  }, error = function(e) {
    list(success = FALSE, message = conditionMessage(e))
  })

  withdrawn_validations <- list(
    paper_eigenvalue_ratios = list(
      status = "withdrawn",
      reason = paste0(
        "The former validator did not reproduce a documented paper fixture ",
        "and solved the wrong generalized right-hand side."
      )
    ),
    regression_solver_consistency = list(
      status = "withdrawn",
      reason = paste0(
        "solver='regression' never selected a distinct implementation, so ",
        "agreement with solver='exact' was tautological."
      )
    ),
    out_of_sample_reconstruction = list(
      status = "withdrawn",
      reason = paste0(
        "The former routine returned expected_error plus random noise without ",
        "fitting KEMA or evaluating held-out predictions."
      ),
      replacement = "cv_alignment_rows"
    )
  )

  overall_success <- isTRUE(fidelity_validation$success) &&
    isTRUE(backend_agreement$success)
  if (verbose) {
    cat(
      "  fidelity gate: ",
      if (isTRUE(fidelity_validation$success)) "PASS" else "FAIL",
      "\n",
      sep = ""
    )
    cat(
      "  backend agreement: ",
      if (isTRUE(backend_agreement$success)) "PASS" else "FAIL",
      "\n",
      sep = ""
    )
  }

  list(
    fidelity_validation = fidelity_validation,
    backend_agreement = backend_agreement,
    withdrawn_validations = withdrawn_validations,
    overall_success = overall_success
  )
}
