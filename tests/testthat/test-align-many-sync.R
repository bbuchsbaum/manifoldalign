# Orientation, synchronisation and diagnostics contract of align_many().
#
# Contract (see ?new_align_transform): an edge transform g_ij maps domain i
# onto domain j (X_i %*% g_ij = X_j for "O"/"GL", g_ij %*% X_i = X_j for
# "perm"); per-domain transforms T_i map each domain onto one consensus, so
# apply_transform(T_i, X_i) must agree across domains. The tests below check
# the transformed observations directly, not only the reported residuals.

# -- Test adapters ------------------------------------------------------------

# Exact linear adapter: g_ij solves X_i g = X_j in least squares.
.exact_linear_aligner <- function(group) {
  new_aligner("test_exact_linear", group = group)
}

.register_exact_linear <- function() {
  ns <- asNamespace("manifoldalign")
  registerS3method("fit_pair", "test_exact_linear_aligner",
    function(algo, X_i, X_j, links = NULL, ...) {
      structure(list(op = qr.solve(as.matrix(X_i), as.matrix(X_j)),
                     group = algo$group, links = links),
                class = "test_exact_linear_fit")
    }, envir = ns)
  registerS3method("relative_transform", "test_exact_linear_fit",
    function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
      from <- match.arg(from)
      op <- if (from == "i") fit$op else solve(fit$op)
      new_align_transform(fit$group, op, from = from, to = to, k = ncol(fit$op))
    }, envir = ns)
}

# Exact row-matching adapter: rows are matched by their first column.
.exact_perm_aligner <- function() new_aligner("test_exact_perm", group = "perm")

.register_exact_perm <- function() {
  ns <- asNamespace("manifoldalign")
  registerS3method("fit_pair", "test_exact_perm_aligner",
    function(algo, X_i, X_j, links = NULL, ...) {
      op <- outer(X_j[, 1], X_i[, 1], function(a, b) as.numeric(abs(a - b) < 1e-12))
      structure(list(op = op), class = "test_exact_perm_fit")
    }, envir = ns)
  registerS3method("relative_transform", "test_exact_perm_fit",
    function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
      from <- match.arg(from)
      op <- if (from == "i") fit$op else t(fit$op)
      new_align_transform("perm", op, from = from, to = to)
    }, envir = ns)
}

.rot3 <- function(a, b, c) {
  Rz <- matrix(c(cos(a), sin(a), 0, -sin(a), cos(a), 0, 0, 0, 1), 3)
  Rx <- matrix(c(1, 0, 0, 0, cos(b), sin(b), 0, -sin(b), cos(b)), 3)
  Ry <- matrix(c(cos(c), 0, -sin(c), 0, 1, 0, sin(c), 0, cos(c)), 3)
  Rz %*% Rx %*% Ry
}

.max_embedding_gap <- function(Z) {
  max(vapply(Z[-1], function(z) max(abs(as.matrix(z) - as.matrix(Z[[1]]))), numeric(1)))
}

# -- Orthogonal group ----------------------------------------------------------

test_that("O(3) sync aligns non-commuting rotations and diagnostics certify it", {
  .register_exact_linear()
  set.seed(11)
  X <- matrix(rnorm(60), 20, 3)
  R2 <- .rot3(0.3, 0.7, -0.4)
  R3 <- .rot3(-1.0, 0.2, 0.9)
  expect_gt(max(abs(R2 %*% R3 - R3 %*% R2)), 0.1)  # genuinely non-commuting
  domains <- list(X, X %*% R2, X %*% R3)

  res <- align_many(domains, .exact_linear_aligner("O"), graph = "complete")
  expect_lt(.max_embedding_gap(res$embeddings), 1e-10)
  expect_lt(res$diagnostics$mean_edge_residual, 1e-10)
  expect_lt(res$diagnostics$mean_cycle_residual, 1e-10)

  # The relation T_i = g_ij T_j holds edge by edge.
  for (idx in seq_len(nrow(res$graph))) {
    i <- res$graph$i[idx]; j <- res$graph$j[idx]
    g <- relative_transform(res$pair_fits[[idx]], "i", "j")
    expect_equal(g$op %*% res$transforms[[j]]$op, res$transforms[[i]]$op,
                 tolerance = 1e-10)
  }
})

test_that("O(3) local refinement keeps an exact solution exact", {
  .register_exact_linear()
  set.seed(12)
  X <- matrix(rnorm(60), 20, 3)
  domains <- list(X, X %*% .rot3(0.3, 0.7, -0.4), X %*% .rot3(-1, 0.2, 0.9),
                  X %*% .rot3(0.5, -0.6, 0.1))
  res <- align_many(domains, .exact_linear_aligner("O"), graph = "complete",
                    sync = list(refine = list(gpa = TRUE)))
  expect_lt(.max_embedding_gap(res$embeddings), 1e-10)
  expect_lt(res$diagnostics$mean_edge_residual, 1e-10)
})

test_that("O(3) refinement repairs a perturbed consensus", {
  set.seed(13)
  N <- 4
  R <- list(diag(3), .rot3(0.3, 0.7, -0.4), .rot3(-1, 0.2, 0.9), .rot3(0.5, -0.6, 0.1))
  E <- data.frame(i = c(1, 1, 1, 2, 2, 3), j = c(2, 3, 4, 3, 4, 4))
  G <- lapply(seq_len(nrow(E)), function(e) {
    new_align_transform("O", R[[E$i[e]]] %*% t(R[[E$j[e]]]), from = E$i[e], to = E$j[e], k = 3)
  })
  W <- rep(1, nrow(E))
  # Perturb the exact consensus of domain 3 and refine.
  T0 <- lapply(seq_len(N), function(i) new_align_transform("O", R[[i]], from = i, to = "consensus", k = 3))
  T0[[3]] <- new_align_transform("O", .rot3(0.05, 0, 0) %*% R[[3]], from = 3, to = "consensus", k = 3)
  before <- manifoldalign:::compute_diagnostics_("O", T0, G, E, W)$mean_edge_residual
  T1 <- manifoldalign:::refine_orthogonal_local_(T0, G, W, E, N = N)
  after <- manifoldalign:::compute_diagnostics_("O", T1, G, E, W)$mean_edge_residual
  expect_lt(after, before / 2)
  # Repeated sweeps converge back to an exact (gauge-rotated) solution.
  for (s in 1:30) T1 <- manifoldalign:::refine_orthogonal_local_(T1, G, W, E, N = N)
  expect_lt(manifoldalign:::compute_diagnostics_("O", T1, G, E, W)$mean_edge_residual, 1e-6)
})

# -- General linear group ------------------------------------------------------

test_that("GL sync aligns diagonal and sheared maps; residual certifies it", {
  .register_exact_linear()
  set.seed(14)
  X <- matrix(rnorm(40), 20, 2)
  domains <- list(X, X %*% diag(c(2, 0.5)), X %*% matrix(c(1, 0.3, -0.2, 1.5), 2))

  exact <- align_many(domains, .exact_linear_aligner("GL"), graph = "complete",
                      sync = list(gl = list(ridge = 0)))
  expect_lt(.max_embedding_gap(exact$embeddings), 1e-10)
  expect_lt(exact$diagnostics$mean_edge_residual, 1e-10)
  expect_lt(exact$diagnostics$mean_cycle_residual, 1e-10)
  # Gauge: domain 1 keeps its coordinates.
  expect_equal(exact$embeddings[[1]], X, tolerance = 1e-10)

  ridged <- align_many(domains, .exact_linear_aligner("GL"), graph = "complete")
  expect_lt(.max_embedding_gap(ridged$embeddings), 1e-4)
  expect_lt(ridged$diagnostics$mean_edge_residual, 1e-4)
})

test_that("GL diagnostics flag a consensus that misaligns the data", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  D <- diag(c(2, 0.5))
  E <- data.frame(i = 1, j = 2)
  G <- list(new_align_transform("GL", D, from = 1, to = 2, k = 2))  # X D = X_2
  good <- list(new_align_transform("GL", diag(2), from = 1, to = "consensus"),
               new_align_transform("GL", solve(D), from = 2, to = "consensus"))
  bad <- list(good[[1]], new_align_transform("GL", D, from = 2, to = "consensus"))
  d_good <- manifoldalign:::compute_diagnostics_("GL", good, G, E, 1)
  d_bad <- manifoldalign:::compute_diagnostics_("GL", bad, G, E, 1)
  expect_lt(d_good$mean_edge_residual, 1e-12)
  expect_gt(d_bad$mean_edge_residual, 0.5)
  # Embeddings agree exactly for the good consensus only.
  expect_equal(apply_transform(good[[2]], X %*% D), apply_transform(good[[1]], X))
  expect_gt(max(abs(apply_transform(bad[[2]], X %*% D) - X)), 1)
})

# -- Graph shapes and cycle diagnostics ----------------------------------------

test_that("star and mst graphs with four domains run and align", {
  .register_exact_linear()
  set.seed(15)
  X <- matrix(rnorm(90), 30, 3)
  domains <- list(X, X %*% .rot3(0.3, 0.7, -0.4), X %*% .rot3(-1, 0.2, 0.9),
                  X %*% .rot3(0.5, -0.6, 0.1))
  for (graph in c("star", "mst")) {
    res <- align_many(domains, .exact_linear_aligner("O"), graph = graph)
    expect_equal(nrow(res$graph), 3L)
    expect_lt(.max_embedding_gap(res$embeddings), 1e-10)
    expect_lt(res$diagnostics$mean_edge_residual, 1e-10)
    expect_null(res$diagnostics$cycle_residuals)
  }
})

test_that("cycle diagnostics handle absent edges and reversed edge orientation", {
  R <- list(diag(3), .rot3(0.3, 0.7, -0.4), .rot3(-1, 0.2, 0.9), .rot3(0.5, -0.6, 0.1))
  # Triangle 1-2-3 with the 1-3 edge stored as 3 -> 1, plus a pendant 4.
  E <- data.frame(i = c(1, 2, 3, 1), j = c(2, 3, 1, 4))
  G <- lapply(seq_len(nrow(E)), function(e) {
    new_align_transform("O", R[[E$i[e]]] %*% t(R[[E$j[e]]]), from = E$i[e], to = E$j[e], k = 3)
  })
  T_list <- lapply(1:4, function(i) new_align_transform("O", R[[i]], from = i, to = "consensus", k = 3))
  d <- manifoldalign:::compute_diagnostics_("O", T_list, G, E, rep(1, 4))
  expect_equal(nrow(d$cycle_residuals), 1L)
  expect_equal(unlist(d$cycle_residuals[1, c("i", "j", "k")]), c(i = 1, j = 2, k = 3))
  expect_lt(d$mean_cycle_residual, 1e-12)
  expect_lt(d$mean_edge_residual, 1e-12)

  # A broken edge is detected by the cycle residual.
  G[[2]] <- new_align_transform("O", .rot3(0.4, 0, 0) %*% G[[2]]$op, from = 2, to = 3, k = 3)
  d2 <- manifoldalign:::compute_diagnostics_("O", T_list, G, E, rep(1, 4))
  expect_gt(d2$mean_cycle_residual, 0.05)
})

# -- Permutation group ---------------------------------------------------------

.perm_fixture <- function() {
  set.seed(16)
  P <- cbind(seq(0.5, 3, by = 0.5), matrix(rnorm(12), 6, 2))  # unique keys
  list(
    domains = list(P, P[c(4, 1, 6, 2, 5, 3), ], P[c(5, 2, 6, 3), ], P[c(2, 6, 1, 4, 3, 5), ]),
    P = P
  )
}

.expect_perm_alignment <- function(res, P, tol = 1e-8, present_rows = NULL) {
  for (i in seq_along(res$transforms)) {
    Ti <- as.matrix(res$transforms[[i]]$op)
    expect_equal(nrow(Ti), nrow(P))
    present <- rowSums(Ti) > 0.5
    if (!is.null(present_rows[[i]])) expect_identical(which(present), present_rows[[i]])
    Zi <- as.matrix(res$embeddings[[i]])
    expect_equal(Zi[present, , drop = FALSE], P[present, , drop = FALSE], tolerance = tol)
    expect_true(all(abs(rowSums(Ti)[!present]) < tol))
  }
  expect_lt(res$diagnostics$mean_edge_residual, max(tol, 1e-6))
}

test_that("permutation star, global and spectral sync recover hard maps with unequal sizes", {
  .register_exact_perm()
  fx <- .perm_fixture()
  for (mode in c("star", "global", "spectral")) {
    res <- align_many(fx$domains, .exact_perm_aligner(), graph = "complete",
                      sync = list(perm = list(mode = mode)))
    .expect_perm_alignment(res, fx$P)
  }
  skip_if_not_installed("clue")
  res <- align_many(fx$domains, .exact_perm_aligner(), graph = "complete",
                    sync = list(perm = list(mode = "spectral", rounding = "hungarian")))
  .expect_perm_alignment(res, fx$P)
})

test_that("permutation sync works on identity maps and without hub edges", {
  .register_exact_perm()
  fx <- .perm_fixture()
  same <- rep(list(fx$P), 5)
  # Default for N > 4 is the spectral backend.
  res <- align_many(same, .exact_perm_aligner(), graph = "complete")
  .expect_perm_alignment(res, fx$P)

  # Chain 1-2-3-4: star consensus composes along the path to the hub. Rows
  # 1 and 4 of P are absent from domain 3, so no edge links them to domain 4
  # and they must stay unassigned there rather than be matched spuriously.
  present_rows <- list(1:6, 1:6, c(2L, 3L, 5L, 6L), c(2L, 3L, 5L, 6L))
  G <- list()
  E <- data.frame(i = c(1, 2, 3), j = c(2, 3, 4))
  fits <- lapply(seq_len(nrow(E)), function(e) {
    fit_pair(.exact_perm_aligner(), fx$domains[[E$i[e]]], fx$domains[[E$j[e]]])
  })
  G <- lapply(fits, relative_transform, from = "i", to = "j")
  for (mode in c("star", "global", "spectral")) {
    T_list <- switch(mode,
      star = manifoldalign:::permutation_sync_star(G, rep(1, 3), E, N = 4),
      global = manifoldalign:::permutation_sync(G, rep(1, 3), E, N = 4),
      spectral = manifoldalign:::permutation_sync_spectral(G, rep(1, 3), E, N = 4))
    Z <- lapply(1:4, function(i) apply_transform(T_list[[i]], fx$domains[[i]]))
    res <- list(transforms = T_list, embeddings = Z,
                diagnostics = manifoldalign:::compute_diagnostics_("perm", T_list, G, E, rep(1, 3)))
    .expect_perm_alignment(res, fx$P, tol = 1e-6, present_rows = present_rows)
  }
})

test_that("permutation sync rejects inconsistent domain sizes", {
  E <- data.frame(i = c(1, 2), j = c(2, 3))
  G <- list(diag(3), matrix(1, 3, 4))  # implies n_2 = 3 and n_2 = 4
  expect_error(manifoldalign:::permutation_sync_star(G, c(1, 1), E, N = 3),
               "implies 4 rows for domain 2")
})

# -- align_many argument handling ---------------------------------------------

test_that("align_many passes per-edge links and rejects unknown edges", {
  .register_exact_linear()
  set.seed(17)
  X <- matrix(rnorm(30), 10, 3)
  domains <- list(X, X %*% .rot3(0.3, 0.7, -0.4), X %*% .rot3(-1, 0.2, 0.9))
  L12 <- data.frame(a = 1:3, b = 1:3)
  res <- align_many(domains, .exact_linear_aligner("O"), graph = "complete",
                    links = list("1:2" = L12))
  keys <- paste(res$graph$i, res$graph$j, sep = ":")
  expect_identical(res$pair_fits[[match("1:2", keys)]]$links, L12)
  expect_null(res$pair_fits[[match("2:3", keys)]]$links)

  expect_error(align_many(domains, .exact_linear_aligner("O"), graph = "complete",
                          links = list("2:1" = L12)),
               "opposite orientation")
  expect_error(align_many(domains, .exact_linear_aligner("O"), graph = "star",
                          links = list("2:3" = L12)),
               "not edges of the dataset graph")
  expect_error(align_many(domains, .exact_linear_aligner("O"), links = L12),
               "edge keys")
  two <- align_many(domains[1:2], .exact_linear_aligner("O"), links = L12)
  expect_identical(two$pair_fits[[1]]$links, L12)
})

test_that("align_many warns with the cause when a domain cannot be embedded", {
  ns <- asNamespace("manifoldalign")
  registerS3method("fit_pair", "test_basis_aligner",
    function(algo, X_i, X_j, links = NULL, ...) structure(list(), class = "test_basis_fit"),
    envir = ns)
  registerS3method("relative_transform", "test_basis_fit",
    function(fit, from = "i", to = "j", ...) new_align_transform("O", diag(2), from = from, to = to, k = 2),
    envir = ns)
  domains <- replicate(3, matrix(rnorm(15), 5, 3), simplify = FALSE)
  w <- testthat::capture_warnings(
    res <- align_many(domains, new_aligner("test_basis", group = "O"), graph = "complete")
  )
  expect_length(w, 3L)
  expect_match(w[1], "could not embed domain 1 .*dimension mismatch")
  expect_true(all(vapply(res$embeddings, is.null, logical(1))))
})

test_that("align_many rejects groups and consensus modes it cannot honour", {
  domains <- replicate(3, matrix(rnorm(15), 5, 3), simplify = FALSE)
  expect_error(align_many(domains, new_aligner("test_kernel", group = "kernel")),
               "group 'kernel'")
  expect_error(align_many(domains, new_aligner("test_enc", group = "encoder")),
               "group 'encoder'")
  expect_error(align_many(domains, new_aligner("test_o", group = "O"), consensus = "tree"),
               "not implemented")
  .register_exact_linear()
  expect_error(align_many(domains, .exact_linear_aligner("GL"),
                          sync = list(gl = list(iters = 5))),
               "only applies")
})

test_that("align_many graph scoring does not disturb the caller's RNG", {
  .register_exact_linear()
  set.seed(18)
  X <- matrix(rnorm(90), 30, 3)
  domains <- list(X, X %*% .rot3(0.3, 0.7, -0.4), X %*% .rot3(-1, 0.2, 0.9),
                  X %*% .rot3(0.5, -0.6, 0.1))
  set.seed(99); expected <- runif(3)
  set.seed(99)
  r1 <- align_many(domains, .exact_linear_aligner("O"), graph = "mst")
  expect_identical(runif(3), expected)
  r2 <- align_many(domains, .exact_linear_aligner("O"), graph = "mst")
  expect_identical(r1$graph, r2$graph)
})

test_that("native fit_many methods map k to ncomp explicitly", {
  calls <- list()
  record <- function(name) {
    force(name)
    function(data, ...) {
      calls[[name]] <<- list(...)
      structure(list(), class = "mock_native_fit")
    }
  }
  local_mocked_bindings(
    kema = record("kema"),
    cone_align_multiple = record("cone"),
    global_geo_align = record("global_geo"),
    spectral_mnn_align = record("spectral_mnn"),
    ssma_align = record("ssma"),
    generalized_procrustes = record("gpa"),
    .package = "manifoldalign"
  )
  hd <- structure(list(a = list(x = diag(2)), b = list(x = diag(2))), class = "hyperdesign")
  algos <- list(kema = kema_aligner(), cone = cone_multi_aligner(),
                global_geo = global_geo_aligner(), spectral_mnn = spectral_mnn_aligner(),
                ssma = ssma_aligner())
  for (nm in names(algos)) {
    align_many(hd, algos[[nm]], k = 3, y = "label")
    expect_identical(calls[[nm]]$ncomp, 3)
    expect_null(calls[[nm]]$k)
    expect_null(calls[[nm]]$knn)
    expect_identical(calls[[nm]]$y, "label")
    expect_error(fit_many(algos[[nm]], hd, k = 3, ncomp = 4), "disagree")
  }
  align_many(hd, gpa_aligner(), y = "label")
  expect_null(calls$gpa$k)
  expect_error(align_many(hd, gpa_aligner(), k = 2), "no latent dimension")
})

test_that("aligner constructors declare the full capability contract", {
  fields <- c("group", "supports_multi", "encoder", "supports_predict",
              "supports_calibration", "requires_auxiliary_at_predict", "invertible")
  algos <- list(kema_aligner(), gpa_aligner(), cone_multi_aligner(),
                global_geo_aligner(), spectral_mnn_aligner(), ssma_aligner(),
                grasp_aligner())
  for (a in algos) {
    caps <- aligner_capabilities(a)
    expect_identical(names(caps), fields)
    expect_true(caps$encoder %in% c("linear", "nonlinear", "none"))
    for (f in fields[-(1:3)]) expect_true(is.logical(caps[[f]]) && length(caps[[f]]) == 1L)
  }
  expect_identical(aligner_capabilities(grasp_aligner())$encoder, "none")
  expect_false(aligner_capabilities(grasp_aligner())$supports_predict)
})

# -- GRASP adapter orientation -----------------------------------------------

test_that("GRASP relative transform maps domain i coordinates onto domain j", {
  # Hand-built fit: M is a non-symmetric signed permutation with C_j M = C_i.
  M <- matrix(0, 3, 3); M[cbind(c(2, 3, 1), 1:3)] <- c(1, -1, 1)
  expect_false(isTRUE(all.equal(M, t(M))))
  C_j <- matrix(rnorm(12), 4, 3)
  fit <- structure(list(rotation = M, coefficients = list(i = C_j %*% M, j = C_j),
                        k = 3, loss = 0), class = "grasp_pair_fit")
  g_ij <- relative_transform(fit, "i", "j")
  g_ji <- relative_transform(fit, "j", "i")
  expect_equal(apply_transform(g_ij, fit$coefficients$i), fit$coefficients$j)
  expect_equal(apply_transform(g_ji, fit$coefficients$j), fit$coefficients$i)
})

test_that("GRASP pair fit stores coefficients consistent with its transform", {
  skip_on_cran()
  set.seed(19)
  X1 <- matrix(rnorm(80), 40, 2)
  perm <- sample(40)
  fit <- fit_pair(grasp_aligner(), X1, X1[perm, ], ncomp = 4, q_descriptors = 20)
  g <- relative_transform(fit, "i", "j")
  C <- fit$coefficients
  resid_ij <- sum((apply_transform(g, C$i) - C$j)^2)
  expect_equal(fit$loss, sum((C$i - C$j %*% fit$rotation)^2), tolerance = 1e-10)
  expect_equal(resid_ij, fit$loss, tolerance = 1e-8)
})

# -- Sparse inputs and eigen-solver fallbacks ---------------------------------

test_that("apply_transform preserves Matrix inputs and matches dense results", {
  set.seed(20)
  Xs <- Matrix::rsparsematrix(6, 3, 0.5)
  Xd <- as.matrix(Xs)
  tr <- new_align_transform("GL", matrix(rnorm(9), 3, 3), from = 1, to = 2)
  out <- apply_transform(tr, Xs)
  expect_s4_class(out, "Matrix")
  expect_equal(as.matrix(out), Xd %*% tr$op)
  expect_true(is.matrix(apply_transform(tr, Xd)))

  P <- Matrix::sparseMatrix(i = c(2, 1, 3), j = 1:3, x = 1, dims = c(3, 3))
  trp <- new_align_transform("perm", P, from = 1, to = 2)
  outp <- apply_transform(trp, Matrix::rsparsematrix(3, 4, 0.5))
  expect_s4_class(outp, "sparseMatrix")
})

test_that("resolve_hyperdesign and as_hyperdesign keep sparse blocks on request", {
  Xs <- Matrix::rsparsematrix(5, 3, 0.4)
  hd <- as_hyperdesign(list(Xs, as.matrix(Xs)))
  expect_s4_class(hd$domain1$x, "dgCMatrix")
  expect_equal(nrow(hd$domain1$design), 5L)
  expect_equal(nrow(hd$domain2$design), 5L)
  dense <- manifoldalign:::resolve_hyperdesign(hd)
  expect_true(is.matrix(dense$domains[[1]]$x))
  sparse <- manifoldalign:::resolve_hyperdesign(hd, keep_sparse = TRUE)
  expect_s4_class(sparse$domains[[1]]$x, "dgCMatrix")
  expect_equal(manifoldalign:::feature_block_indices(sparse$domains),
               list(domain1 = 1:3, domain2 = 4:6))
})

test_that("graph scoring accepts sparse domains", {
  set.seed(21)
  doms <- replicate(4, Matrix::rsparsematrix(20, 5, 0.5), simplify = FALSE)
  S <- manifoldalign:::compute_dataset_pair_scores_(doms, NULL)
  expect_equal(dim(S), c(4L, 4L))
  expect_true(all(is.finite(S)))
})

test_that("iterative and dense top eigenvectors span the same subspace", {
  skip_if_not_installed("RSpectra")
  set.seed(22)
  A <- matrix(rnorm(400 * 400), 400)
  M <- Matrix::forceSymmetric(Matrix::Matrix(crossprod(A) / 400))
  V_dense <- manifoldalign:::.sync_top_eigenvectors(M, 3)
  V_iter <- manifoldalign:::.sync_top_eigenvectors(M, 3, dense_max = 0L)
  expect_equal(V_dense %*% t(V_dense), V_iter %*% t(V_iter), tolerance = 1e-6)
})

test_that("every exported aligner declares the full capability contract", {
  fields <- c("group", "supports_multi", "encoder", "supports_predict",
              "supports_calibration", "requires_auxiliary_at_predict",
              "invertible")
  constructors <- grep("_aligner$", getNamespaceExports("manifoldalign"),
                       value = TRUE)
  constructors <- setdiff(constructors, "new_aligner")
  for (ctor in constructors) {
    caps <- aligner_capabilities(getExportedValue("manifoldalign", ctor)())
    expect_identical(names(caps), fields, info = ctor)
    expect_true(caps$encoder %in% c("linear", "nonlinear", "none"), info = ctor)
  }
})
