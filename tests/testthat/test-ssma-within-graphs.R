.within_graph_fixture <- function() {
  set.seed(411)
  x <- matrix(rnorm(64), 16, 4)
  list(a = x, b = x + matrix(rnorm(64, sd = .08), 16, 4))
}
.within_graph_control <- function(...) {
  ssma_align_control(knn = 4, sigma = .8, compression = 'none',
                     auto_correspondences = FALSE, normalization = 'domain_scores',
                     domain_tol = 1e-7, domain_maxiter = 100000,
                     seed = 412, verbose = FALSE, ...)
}
.within_graph_constraints <- function(x) {
  list(linear_correspondence(names(x), lapply(x, function(z) z[1:8, , drop = FALSE]), strength = .3))
}

test_that('supplied SSMA graphs match an independent edge quadratic oracle', {
  x <- .within_graph_fixture()
  A <- matrix(0, 16, 16)
  for (i in 1:15) A[i, i + 1] <- A[i + 1, i] <- i / 16
  g <- .ssma_resolve_within_graphs(list(A, Matrix::Matrix(A, sparse = TRUE)),
                                lapply(x, function(z) list(x = z)), names(x))
  set.seed(413); B <- matrix(rnorm(8), 4, 2)
  for (i in 1:2) {
    z <- x[[i]] %*% B
    direct <- sum(vapply(1:15, function(j) (j / 16) * sum((z[j, ] - z[j + 1, ])^2), numeric(1)))
    L <- diag(rowSums(A)) - A
    expect_equal(sum(diag(crossprod(z, L %*% z))), direct, tolerance = 1e-12)
    expect_equal(as.matrix(g[[i]]$A), A)
    expect_equal(g[[i]]$degree, rowSums(A))
    expect_s4_class(g[[i]]$A, 'sparseMatrix')
  }
})

test_that('native and supplied SSMA sample graphs give the same domain-score fit', {
  x <- .within_graph_fixture(); ctrl <- .within_graph_control()
  G <- lapply(x, function(z) .ssma_build_within_graph(z, ctrl$knn, ctrl$sigma)$A)
  native <- ssma_align(x, preproc = NULL, ncomp = 2, constraints = .within_graph_constraints(x), control = ctrl)
  supplied <- ssma_align(x, preproc = NULL, ncomp = 2, constraints = .within_graph_constraints(x), within_graphs = G, control = ctrl)
  expect_true(native$domain_diagnostics$converged)
  expect_true(supplied$domain_diagnostics$converged)
  expect_equal(supplied$domain_diagnostics$objective, native$domain_diagnostics$objective, tolerance = 1e-10)
  expect_equal(tcrossprod(supplied$s), tcrossprod(native$s), tolerance = 1e-8)
  expect_identical(native$within_graphs_source, 'feature_knn')
  expect_identical(supplied$within_graphs_source, 'provided')
  expect_length(supplied$within_graphs_diagnostics, 2)
})

test_that('supplied SSMA graphs bypass kNN and preserve row-permutation equivariance', {
  x <- .within_graph_fixture(); ctrl <- .within_graph_control()
  G <- lapply(x, function(z) .ssma_build_within_graph(z, ctrl$knn, ctrl$sigma)$A)
  local_mocked_bindings(.ssma_build_within_graph = function(...) stop('unexpected graph construction'))
  f <- ssma_align(x, preproc = NULL, ncomp = 2, constraints = .within_graph_constraints(x), within_graphs = G, control = ctrl)
  perm <- list(rev(seq_len(16)), c(9:16, 1:8))
  xp <- Map(function(z, p) z[p, , drop = FALSE], x, perm)
  gp <- Map(function(g, p) g[p, p, drop = FALSE], G, perm)
  fp <- ssma_align(xp, preproc = NULL, ncomp = 2, constraints = .within_graph_constraints(x), within_graphs = gp, control = ctrl)
  expect_equal(tcrossprod(fp$v), tcrossprod(f$v), tolerance = 1e-8)
  expect_error(ssma_align(x, preproc = NULL, ncomp = 2, control = ctrl), 'unexpected graph construction')
})

test_that('supplied SSMA graphs resolve names and reject ambiguous sample identity', {
  x <- list(a = matrix(1:12, 4, 3), b = matrix(1:15, 5, 3))
  d <- lapply(x, function(z) list(x = z)); g <- list(a = matrix(0, 4, 4), b = matrix(0, 5, 5))
  out <- .ssma_resolve_within_graphs(g[2:1], d, names(x))
  expect_equal(vapply(out, function(z) nrow(z$A), integer(1)), c(4L, 5L))
  expect_equal(out[[1]]$degree, rep(0, 4))
  bad <- g; names(bad) <- c('a', 'other')
  expect_error(.ssma_resolve_within_graphs(bad, d, names(x)), 'domain names')
  names(bad) <- c('a', 'a')
  expect_error(.ssma_resolve_within_graphs(bad, d, names(x)), 'domain names')
  dimnames(g$a) <- list(letters[1:4], letters[1:4])
  expect_error(.ssma_resolve_within_graphs(g, d, names(x)), 'row/column names')
  rownames(d$a$x) <- letters[1:4]
  expect_silent(.ssma_resolve_within_graphs(g, d, names(x)))
  g$a <- g$a[4:1, 4:1]
  expect_error(.ssma_resolve_within_graphs(g, d, names(x)), 'row/column names')
})

test_that('supplied SSMA graphs reject malformed weights and serial conflicts', {
  x <- .within_graph_fixture(); d <- lapply(x, function(z) list(x = z))
  g <- list(matrix(0, 16, 16), matrix(0, 16, 16))
  resolve <- function(a) .ssma_resolve_within_graphs(list(a, g[[2]]), d, names(x))
  expect_error(.ssma_resolve_within_graphs(g[1], d, names(x)), 'one adjacency')
  expect_error(resolve(matrix(0, 15, 15)), 'size mismatch')
  expect_error(resolve(matrix('x', 16, 16)), 'numeric')
  for (v in c(-1, NA_real_, Inf)) {
    a <- g[[1]]; a[1, 2] <- a[2, 1] <- v
    expect_error(resolve(a), 'finite and nonnegative')
  }
  a <- g[[1]]; a[1, 2] <- 1
  expect_error(resolve(a), 'symmetric')
  a[2, 1] <- 1; a[1, 1] <- 1
  expect_error(resolve(a), 'zero diagonal')
  ctrl <- .within_graph_control()
  expect_error(ssma_align(x, within_graphs = g, serial_index = list(1:16, 1:16), control = ctrl), 'cannot be combined')
  ctrl$serial_control <- ssma_serial_control(enabled = TRUE)
  expect_error(ssma_align(x, within_graphs = g, control = ctrl), 'cannot be combined')
})
