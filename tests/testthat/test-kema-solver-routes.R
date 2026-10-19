# KEMA quotient solver routes: preconditioned LOBPCG, the certified dense
# fallback, and the matrix-free operator_exact backend. Every route is checked
# against an independent dense oracle of the original KEMA pencil (Eq. 6):
#   K (L + mu L_s) K alpha = theta (K L_d K + lambda I) alpha.

# Independent dense oracle: forms the original n x n pencil with base R and
# solves it through the Cholesky factor of the SPD right-hand side. It shares
# no code with the quotient construction or eigencore.
kema_dense_pencil_oracle <- function(X_list, y, kernel, knn, sigma, mu,
                                     lambda, ncomp) {
  W <- as.matrix(Matrix::bdiag(lapply(X_list, function(x) {
    g <- adjoin::graph_weights(
      x,
      weight_mode = "normalized",
      neighbor_mode = "knn",
      k = knn,
      type = "normal",
      sigma = sigma
    )
    A <- as.matrix(adjoin::adjacency(g))
    A <- (A + t(A)) / 2
    diag(A) <- 0
    A
  })))
  labelled <- !is.na(y)
  yy <- as.character(y)
  Ws <- outer(yy, yy, "==") & outer(labelled, labelled, "&")
  Wd <- outer(yy, yy, "!=") & outer(labelled, labelled, "&")
  Ws[is.na(Ws)] <- FALSE
  Wd[is.na(Wd)] <- FALSE
  Ws <- Ws * 1
  Wd <- Wd * 1
  diag(Ws) <- 0
  diag(Wd) <- 0
  laplacian <- function(A) diag(rowSums(A)) - A
  M <- laplacian(W) + mu * laplacian(Ws)
  Ld <- laplacian(Wd)
  K <- as.matrix(Matrix::bdiag(lapply(X_list, function(x) {
    as.matrix(kernlab::kernelMatrix(kernel, x))
  })))
  K <- (K + t(K)) / 2
  A <- K %*% M %*% K
  B <- K %*% Ld %*% K + lambda * diag(nrow(K))
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  R <- chol(B)
  C <- backsolve(R, t(backsolve(R, A, transpose = TRUE)), transpose = TRUE)
  C <- (C + t(C)) / 2
  e <- eigen(C, symmetric = TRUE)
  ord <- order(e$values)
  values <- e$values[ord]
  zero_tol <- sqrt(.Machine$double.eps) * max(abs(values))
  take <- which(values > zero_tol)[seq_len(ncomp)]
  alpha <- backsolve(R, e$vectors[, ord[take], drop = FALSE])
  list(values = values[take], scores = K %*% alpha, all_values = values)
}

# Cosines of the principal angles between two column spaces.
kema_subspace_cosines <- function(S1, S2) {
  Q1 <- qr.Q(qr(as.matrix(S1)))
  Q2 <- qr.Q(qr(as.matrix(S2)))
  svd(crossprod(Q1, Q2), nu = 0, nv = 0)$d
}

kema_route_problem <- function(n = c(40L, 50L), p = 3L, seed = 20261019L) {
  set.seed(seed)
  X_list <- list(
    source = matrix(rnorm(n[1] * p), n[1], p),
    target = matrix(rnorm(n[2] * p) + 0.3, n[2], p)
  )
  y_list <- lapply(n, function(m) {
    y <- rep(c("a", "b", "c"), length.out = m)
    y[seq(2L, m, by = 3L)] <- NA
    y
  })
  hd <- manifoldalign:::as_hyperdesign(
    X_list,
    labels = y_list,
    label_name = "lbl",
    domain_names = names(X_list)
  )
  list(hd = hd, X_list = X_list, y = unlist(y_list, use.names = FALSE))
}

test_that("every KEMA quotient route matches the dense original-pencil oracle", {
  skip_if_not_installed("eigencore", minimum_version = "2.0.0")
  prob <- kema_route_problem()
  kernel <- kernlab::rbfdot(sigma = 2)
  oracle <- kema_dense_pencil_oracle(
    prob$X_list, prob$y, kernel,
    knn = 4, sigma = 0.73, mu = 0.5, lambda = 1e-3, ncomp = 3
  )
  # The oracle's spectrum has the graph-component null direction and then a
  # gap; the requested pairs are well separated from it.
  expect_gt(oracle$values[1], 100 * sqrt(.Machine$double.eps) *
              max(abs(oracle$all_values)))

  fit_route <- function(backend, control = NULL) {
    kema(prob$hd, y = lbl, preproc = NULL, ncomp = 3, knn = 4,
         sigma = 0.73, u = 0.5, kernel = kernel, lambda = 1e-3,
         backend = backend, backend_control = control)
  }
  routes <- list(
    dense_full = list(backend = "full_exact", control = NULL,
                      solver = "eigencore_dense_full", dense = TRUE),
    dense_lobpcg = list(
      backend = "full_exact",
      control = list(dense_exact_max_dim = 3L, dense_fallback_max_dim = 3L),
      solver = "eigencore_lobpcg_deflated", dense = TRUE
    ),
    operator = list(
      backend = "operator_exact",
      control = list(dense_exact_max_dim = 3L),
      solver = "eigencore_lobpcg_operator", dense = FALSE
    )
  )
  for (nm in names(routes)) {
    route <- routes[[nm]]
    fit <- fit_route(route$backend, route$control)
    expect_identical(fit$fidelity$eigensolver, route$solver, label = nm)
    expect_identical(fit$fidelity$dense_quotient_formed, route$dense,
                     label = nm)
    expect_true(fit$fidelity$passed, label = nm)
    expect_true(fit$fidelity$eigencore_certificate_passed, label = nm)
    expect_lte(fit$fidelity$max_rel_residual, 1e-6)
    expect_lte(fit$fidelity$max_b_orth_offdiag, 1e-6)
    expect_equal(fit$eigenvalues$values, oracle$values, tolerance = 1e-7,
                 label = nm)
    expect_gt(min(kema_subspace_cosines(fit$s, oracle$scores)), 1 - 1e-7,
              label = nm)
  }
})

test_that("operator_exact never forms the dense quotient matrices", {
  skip_if_not_installed("eigencore", minimum_version = "2.0.0")
  prob <- kema_route_problem()
  kernel <- kernlab::rbfdot(sigma = 2)
  fit_backend <- function(backend) {
    kema(prob$hd, y = lbl, preproc = NULL, ncomp = 3, knn = 4,
         kernel = kernel, lambda = 1e-3, backend = backend,
         backend_control = list(dense_exact_max_dim = 3L))
  }
  full <- fit_backend("full_exact")
  testthat::local_mocked_bindings(
    .kema_dense_quotient_pencil = function(...) {
      stop("dense KEMA quotient formed")
    },
    .kema_quotient_H = function(...) {
      stop("sparse block-diagonal quotient basis formed")
    },
    .package = "manifoldalign"
  )
  operator <- fit_backend("operator_exact")
  expect_error(fit_backend("full_exact"), "dense KEMA quotient formed")

  expect_identical(operator$backend, "operator_exact")
  expect_false(operator$fidelity$dense_quotient_formed)
  expect_identical(
    operator$fidelity$preconditioner,
    "block_jacobi_shifted_cholesky"
  )
  expect_equal(operator$eigenvalues$values, full$eigenvalues$values,
               tolerance = 1e-8)
  expect_gt(min(kema_subspace_cosines(operator$s, full$s)), 1 - 1e-8)
})

test_that("matrix-free quotient products equal the dense quotient pencil", {
  prob <- kema_route_problem(n = c(17L, 23L))
  K_list <- lapply(prob$X_list, function(x) {
    K <- as.matrix(kernlab::kernelMatrix(kernlab::rbfdot(sigma = 0.4), x))
    methods::as((K + t(K)) / 2, "dgCMatrix")
  })
  graphs <- manifoldalign:::build_original_class_graphs(prob$y)
  Ld <- Matrix::Diagonal(x = Matrix::rowSums(graphs$Wd)) - graphs$Wd
  W <- Matrix::bdiag(lapply(c(17L, 23L), function(m) {
    A <- Matrix::sparseMatrix(i = 1:(m - 1), j = 2:m, x = 1, dims = c(m, m))
    A + Matrix::t(A)
  }))
  M <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W +
    0.5 * (Matrix::Diagonal(x = Matrix::rowSums(graphs$Ws)) - graphs$Ws)
  quotient <- manifoldalign:::build_kema_kernel_quotient(
    K_list, full_form = TRUE,
    component_basis = manifoldalign:::kema_graph_component_basis(W)
  )
  lambda <- 1e-2
  dense <- manifoldalign:::.kema_dense_quotient_pencil(quotient, M, Ld, lambda)
  operator <- manifoldalign:::.kema_operator_quotient_pencil(
    quotient, M, Ld, lambda
  )
  set.seed(3)
  X <- matrix(rnorm(quotient$rank * 4), quotient$rank, 4)
  expect_equal(operator$A_apply(X), dense$A %*% X, tolerance = 1e-12)
  expect_equal(operator$B_apply(X), dense$B %*% X, tolerance = 1e-12)

  # The low-rank metric square root reproduces B and inverts exactly.
  root <- operator$constraints_to_standard
  inv_root <- operator$to_original
  expect_equal(root(root(X)), dense$B %*% X, tolerance = 1e-10)
  expect_equal(inv_root(root(X)), X, tolerance = 1e-10)
  # Whitened operator C = B^{-1/2} A B^{-1/2} is the standard form.
  C_dense <- inv_root(t(inv_root(dense$A)))
  C_applied <- operator$C$apply(X)
  expect_equal(as.matrix(C_applied), C_dense %*% X, tolerance = 1e-9)

  # Chunked per-domain Gram blocks equal the dense congruence, whatever the
  # chunk width (the transient S U product is at most n_i x chunk).
  U <- quotient$U_blocks[[2]]
  d <- quotient$d_blocks[[2]]
  rows <- 18:40
  S <- M[rows, rows]
  reference <- diag(d) %*% crossprod(U, as.matrix(S %*% U)) %*% diag(d)
  for (chunk in c(1L, 5L, 256L)) {
    expect_equal(
      manifoldalign:::.kema_quotient_block_gram(U, d, S, diag_add = 0.25,
                                                chunk = chunk),
      reference + 0.25 * diag(length(d)),
      tolerance = 1e-12
    )
    expect_equal(
      manifoldalign:::.kema_quotient_block_diag(U, d, S, chunk = chunk),
      diag(reference),
      tolerance = 1e-12
    )
  }
})

test_that("LOBPCG uses guard vectors within the deflated quotient", {
  expect_identical(manifoldalign:::.kema_lobpcg_block_size(8, 549), 16L)
  expect_identical(manifoldalign:::.kema_lobpcg_block_size(2, 100), 6L)
  # Never wider than a third of the deflated space, never narrower than ncomp.
  expect_identical(manifoldalign:::.kema_lobpcg_block_size(3, 10), 3L)
  expect_identical(manifoldalign:::.kema_lobpcg_block_size(4, 9), 4L)
})

test_that("RAPID scaling KEMA fixture at 550 nodes is certified (regression)", {
  # Regenerated scaling evidence: every backend failed this exact benchmark
  # call ("eigencore failed to certify the KEMA quotient solve", LOBPCG max
  # backward error 1.065e-05) because unpreconditioned LOBPCG exhausted 500
  # iterations on a quotient whose spectrum spans [0, 1e4] while the
  # targeted eigenvalues are below 0.5. The pencil is well posed.
  skip_on_cran()
  skip_if_not_installed("eigencore", minimum_version = "2.0.0")
  fixture <- rapid_ma_benchmark_fixture(
    n_common = 250L, extra_target = 50L, labels_per_class = 4L, seed = 6001L
  )
  data <- lapply(fixture$data, as.matrix)
  hd <- manifoldalign:::as_hyperdesign(
    data, labels = fixture$labels, label_name = "benchmark_label",
    domain_names = names(data)
  )
  oracle <- kema_dense_pencil_oracle(
    lapply(data, function(x) scale(x, center = TRUE, scale = FALSE)),
    unlist(fixture$labels, use.names = FALSE),
    kernlab::rbfdot(sigma = 0.73),
    knn = 5, sigma = 0.73, mu = 0.5, lambda = 0.01, ncomp = 8
  )
  for (backend in c("auto", "operator_exact")) {
    fit <- kema(hd, y = benchmark_label, ncomp = 8L, knn = 5L,
                solver = "exact", sample_frac = 1, lambda = 0.01,
                backend = backend)
    expect_true(fit$fidelity$passed, label = backend)
    expect_true(fit$fidelity$eigencore_certificate_passed, label = backend)
    expect_equal(fit$fidelity$quotient_dimension, 550L)
    expect_false(fit$fidelity$eigensolver == "eigencore_dense_full",
                 label = backend)
    expect_equal(fit$eigenvalues$values, oracle$values, tolerance = 1e-7,
                 label = backend)
    expect_gt(min(kema_subspace_cosines(fit$s, oracle$scores)), 1 - 1e-6,
              label = backend)
  }
})

test_that("an uncertified partial solve falls back to the dense solve", {
  skip_if_not_installed("eigencore", minimum_version = "2.0.0")
  prob <- kema_route_problem()
  kernel <- kernlab::rbfdot(sigma = 2)
  # One LOBPCG iteration cannot certify; the dense generalized solve then
  # supplies the result, which passes the same certificate and fidelity gate.
  fit <- kema(prob$hd, y = lbl, preproc = NULL, ncomp = 3, knn = 4,
              kernel = kernel, lambda = 1e-3, backend = "full_exact",
              backend_control = list(dense_exact_max_dim = 3L,
                                     eigencore_maxit = 1L))
  expect_identical(fit$fidelity$eigensolver, "eigencore_dense_full")
  expect_length(fit$fidelity$partial_attempts, 1L)
  expect_false(fit$fidelity$partial_attempts[[1]]$certified)
  expect_true(fit$fidelity$eigencore_certificate_passed)
  expect_true(fit$fidelity$passed)

  # Without an admissible dense fallback the uncertified solve is an error,
  # never a returned fit.
  expect_error(
    kema(prob$hd, y = lbl, preproc = NULL, ncomp = 3, knn = 4,
         kernel = kernel, lambda = 1e-3, backend = "full_exact",
         backend_control = list(dense_exact_max_dim = 3L,
                                dense_fallback_max_dim = 3L,
                                eigencore_maxit = 1L)),
    "failed to certify the KEMA quotient solve"
  )
})

test_that("the dense-memory guard is per backend and skips before allocation", {
  skip_if_not_installed("eigencore", minimum_version = "2.0.0")
  est <- manifoldalign:::.kema_dense_memory_estimate
  # The RAPID scaling task at 5,000 common nodes (11,000 samples): the dense
  # quotient backends exceed the default 8 GiB guard, operator_exact does not.
  expect_gt(est(c(5000, 6000), c(5000, 6000), 24, "full_exact", 8), 8 * 1024^3)
  expect_lt(est(c(5000, 6000), c(5000, 6000), 24, "operator_exact", 8),
            8 * 1024^3)
  # The operator estimate grows with sum(r_i^2), the dense one with r_total^2.
  expect_lt(est(c(1e4, 1e4), c(1e4, 1e4), 0, "operator_exact", 2),
            est(c(1e4, 1e4), c(1e4, 1e4), 0, "full_exact", 2))

  prob <- kema_route_problem(n = c(120L, 150L))
  kernel <- kernlab::rbfdot(sigma = 2)
  n_lab <- sum(!is.na(prob$y))
  dense_est <- est(c(120, 150), c(120, 150), n_lab, "full_exact", 2)
  operator_est <- est(c(120, 150), c(120, 150), n_lab, "operator_exact", 2)
  expect_lt(operator_est, dense_est)
  cap <- (operator_est + dense_est) / 2

  fit <- kema(prob$hd, y = lbl, preproc = NULL, ncomp = 2, knn = 4,
              kernel = kernel, lambda = 1e-3,
              backend_control = list(max_dense_bytes = cap))
  expect_identical(fit$backend, "operator_exact")
  expect_identical(fit$backend_candidates, "operator_exact")
  skipped <- Filter(function(a) identical(a$skipped, "max_dense_bytes"),
                    fit$fidelity_history)
  expect_setequal(vapply(skipped, `[[`, character(1), "backend"),
                  c("full_exact", "reduced_exact"))
  expect_true(all(!vapply(skipped, `[[`, logical(1), "solved")))
  expect_true(fit$fidelity$passed)
  expect_equal(unname(fit$fidelity$dense_memory_estimates[["operator_exact"]]),
               operator_est)

  expect_error(
    kema(prob$hd, y = lbl, preproc = NULL, ncomp = 2, knn = 4,
         kernel = kernel, lambda = 1e-3, backend = "full_exact",
         backend_control = list(max_dense_bytes = cap)),
    "dense O\\(n\\^2\\).*cheapest candidate backend 'full_exact'.*max_dense_bytes"
  )
})
