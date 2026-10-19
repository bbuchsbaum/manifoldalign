# Fixed seeds and independent quadratic references: no neural data required.
.lc_fixture <- function() {
  set.seed(5821)
  x <- list(a = matrix(rnorm(24 * 5), 24, 5), b = matrix(rnorm(27 * 7), 27, 7))
  h <- list(matrix(rnorm(4 * 5), 4, 5), matrix(rnorm(4 * 7), 4, 7))
  list(x = x, h = h, link = linear_correspondence(c("a", "b"), h,
       weights = c(1, 0, 2, 3), strength = 2, name = "shared"),
       rows = data.frame(domain_i = 1L, index_i = 1:12, domain_j = 2L,
                         index_j = 2:13, weight = seq_len(12) / 12))
}
.lc_fit <- function(f, constraints = list(f$link), solver = "reduced", compression = "none", ...) {
  ssma_align(f$x, correspondences = f$rows, constraints = constraints,
             preproc = NULL, ncomp = 3,
             control = ssma_align_control(knn = 8, sigma = 2, solver = solver,
               compression = compression, operator_tol = 1e-11,
               operator_maxiter = 5000, ...))
}
.lc_projector <- function(x) { q <- qr.Q(qr(x)); tcrossprod(q) }
.lc_weights <- function(fit) as.matrix(fit$v)

# This oracle does not use any SSMA constraint assembly or operator code.
.lc_oracle <- function(maps, weights, strength = 1) {
  h <- if (length(maps) == 2) cbind(as.matrix(maps[[1]]), -as.matrix(maps[[2]])) else as.matrix(maps[[1]])
  out <- matrix(0, ncol(h), ncol(h))
  for (i in seq_len(nrow(h))) out <- out + strength * weights[i] * outer(h[i, ], h[i, ])
  out
}

test_that("linear constraints validate coordinates, identities and confidence weights", {
  f <- .lc_fixture()
  expect_s3_class(f$link, "linear_alignment_constraint")
  expect_equal(sum(f$link$effective_weights), 1)
  expect_identical(ssma_align_control(compression = "none")$compression, "none")
  expect_error(ssma_align_control(compression = "unknown"), "arg")
  expect_error(ssma_align_control(auto_correspondences = NA), "auto_correspondences")
  for (bad in list(-1, NA_real_, Inf, c(1, 2), "x", 1i, matrix(1, 2, 2))) {
    expect_error(linear_correspondence(c("a", "b"), f$h, weights = bad), "weights")
  }
  for (bad in list(-1, NA_real_, Inf, c(1, 2), 1i, matrix(1, 1, 1))) {
    expect_error(linear_penalty("a", f$h[[1]], strength = bad), "strength")
  }
  expect_error(linear_correspondence(c("a", "a"), f$h), "distinct")
  expect_error(linear_correspondence(c("a", "b"), list(f$h[[1]], f$h[[2]][1:3, ])), "same number")
  expect_error(linear_penalty("a", matrix(NA_real_, 2, 3)), "finite")
  expect_error(linear_penalty("a", matrix(1, 0, 3)), "positive dimensions")
  expect_error(linear_penalty("a", matrix(1i, 2, 3)), "numeric")
  expect_error(linear_penalty("a", diag(2), weights = 0), "positive total")
  expect_error(linear_penalty("a", diag(2), space = "raw"), "processed")
  expect_equal(linear_penalty("a", diag(2), weights = 0, reduction = "sum")$effective_weights, c(0, 0))
  huge <- linear_penalty("a", diag(2), weights = c(1e308, 1e308))
  expect_equal(huge$effective_weights, c(.5, .5))
  expect_error(linear_penalty("a", diag(2), weights = 1e308, strength = 2, reduction = "sum"), "overflow")
  h <- f$h
  rownames(h[[1]]) <- letters[1:4]
  expect_error(linear_correspondence(c("a", "b"), h), "row names")
  rownames(h[[2]]) <- rev(letters[1:4])
  expect_error(linear_correspondence(c("a", "b"), h), "row names")
  rownames(h[[2]]) <- letters[1:4]
  expect_s3_class(linear_correspondence(c("a", "b"), h), "linear_alignment_constraint")
  colnames(h[[1]]) <- rep("same", 5)
  expect_error(linear_correspondence(c("a", "b"), h), "column names")
})

test_that("linear constraint quadratic and operator match an independent outer-product oracle", {
  f <- .lc_fixture()
  domains <- lapply(f$x, function(x) list(x = x))
  cts <- .resolve_linear_constraints(list(f$link), domains, names(domains))
  basis <- lapply(f$x, function(x) list(P = NULL, T = x))
  terms <- .compile_linear_constraints(cts, basis)
  ranks <- c(5, 7)
  oracle <- .lc_oracle(f$h, c(1, 0, 2, 3) / 6, 2)
  actual <- .add_linear_constraint_matrix(matrix(0, 12, 12), terms, ranks)
  expect_equal(actual, oracle, tolerance = 1e-12)
  set.seed(65)
  b <- matrix(rnorm(12 * 3), 12, 3)
  got <- do.call(rbind, .linear_constraint_apply(terms, list(b[1:5, ], b[6:12, ])))
  expect_equal(got, oracle %*% b, tolerance = 1e-12)
  residual <- f$h[[1]] %*% b[1:5, ] - f$h[[2]] %*% b[6:12, ]
  energy <- 2 * sum((c(1, 0, 2, 3) / 6) * rowSums(residual^2))
  expect_equal(sum(b * got), energy, tolerance = 1e-12)
  expect_gte(min(eigen(oracle, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
  d <- .linear_constraint_diagnostics(terms, list(b[1:5, ], b[6:12, ]))[[1]]
  expect_equal(d$objective, energy, tolerance = 1e-12)
  expect_equal(d$discarded_energy_fraction, c(a = 0, b = 0))
})

test_that("sparse maps, unary penalties and split weights preserve the same quadratic", {
  f <- .lc_fixture()
  sparse <- lapply(f$h, function(h) Matrix::Matrix(h, sparse = TRUE))
  link <- linear_correspondence(c("a", "b"), sparse, weights = c(1, 0, 2, 3), strength = 2)
  expect_s4_class(link$maps[[1]], "sparseMatrix")
  penalty <- linear_penalty("a", Matrix::Diagonal(5), strength = .4, name = "ridge")
  domains <- lapply(f$x, function(x) list(x = x))
  basis <- lapply(f$x, function(x) list(P = NULL, T = x))
  compile <- function(x) .compile_linear_constraints(.resolve_linear_constraints(x, domains, names(domains)), basis)
  terms <- compile(list(link, penalty))
  oracle <- .lc_oracle(f$h, c(1, 0, 2, 3) / 6, 2)
  oracle[1:5, 1:5] <- oracle[1:5, 1:5] + diag(.4 / 5, 5)
  expect_equal(.add_linear_constraint_matrix(matrix(0, 12, 12), terms, c(5, 7)), oracle, tolerance = 1e-12)
  expanded <- lapply(f$h, function(h) h[c(1, 1, 2, 3, 4), , drop = FALSE])
  split <- linear_correspondence(c("a", "b"), expanded, weights = c(.5, .5, 0, 2, 3), strength = 2)
  expect_equal(.add_linear_constraint_matrix(matrix(0, 12, 12), compile(list(split)), c(5, 7)),
               .lc_oracle(f$h, c(1, 0, 2, 3) / 6, 2), tolerance = 1e-12)
})

test_that("compressed constraints equal projection of the full objective and quantify discarded directions", {
  f <- .lc_fixture()
  set.seed(91)
  q <- lapply(c(5, 7), function(p) qr.Q(qr(matrix(rnorm(p * 3), p, 3))))
  basis <- lapply(seq_along(q), function(i) list(P = q[[i]], T = f$x[[i]] %*% q[[i]]))
  domains <- lapply(f$x, function(x) list(x = x))
  terms <- .compile_linear_constraints(.resolve_linear_constraints(list(f$link), domains, names(domains)), basis)
  Q <- as.matrix(Matrix::bdiag(q))
  oracle <- crossprod(Q, .lc_oracle(f$h, c(1, 0, 2, 3) / 6, 2) %*% Q)
  expect_equal(.add_linear_constraint_matrix(matrix(0, 6, 6), terms, c(3, 3)), oracle, tolerance = 1e-12)
  for (i in 1:2) {
    residual <- f$h[[i]] - f$h[[i]] %*% tcrossprod(q[[i]])
    expect_equal(unname(terms[[1]]$discarded_energy_fraction[i]), sum(residual^2) / sum(f$h[[i]]^2), tolerance = 1e-12)
  }
})

test_that("SSMA base eigen fallback returns the smallest spectrum before truncation", {
  fit <- .ssma_solve_reduced_geigen(diag(1:8), diag(8), 2, 1e-8, 1e-6)
  expect_equal(fit$eigenvalues, (1:2) / (1 + 1e-6), tolerance = 1e-12)
  expect_equal(.lc_projector(fit$alpha), diag(c(1, 1, rep(0, 6))), tolerance = 1e-12)
})

test_that("SSMA constraints preserve baseline with absent or zero-strength penalties", {
  f <- .lc_fixture()
  zero <- f$link; zero$strength <- 0
  a <- .lc_fit(f, NULL)
  b <- .lc_fit(f, list())
  c <- .lc_fit(f, list(zero))
  expect_equal(a$s, b$s, tolerance = 1e-12)
  expect_equal(.lc_projector(a$s), .lc_projector(c$s), tolerance = 1e-10)
  expect_equal(a$generalized_eigenvalues, c$generalized_eigenvalues, tolerance = 1e-10)
  expect_equal(a$domain_names, c("a", "b"))
  expect_equal(a$reduced_ranks, c(5L, 7L))
  expect_true(all(vapply(a$compression_bases, is.null, logical(1))))
  expect_equal(c$constraint_diagnostics[[1]]$objective, 0)
})

test_that("SSMA supplied maps enforce names, dimensions and mutable object validation", {
  f <- .lc_fixture()
  bad <- f$link; bad$domains[1] <- "missing"
  expect_error(.lc_fit(f, list(bad)), "Unknown constraint domain")
  bad <- f$link; bad$maps[[1]] <- bad$maps[[1]][, 1:4]
  expect_error(.lc_fit(f, list(bad)), "processed features")
  bad <- f$link; bad$weights[1] <- -1
  expect_error(.lc_fit(f, list(bad)), "weights")
  expect_error(.lc_fit(f, list(f$link, f$link)), "names must be unique")
  expect_error(.lc_fit(f, f$link), "list of")
  colnames(f$x$a) <- paste0("a", 1:5)
  colnames(f$link$maps[[1]]) <- rev(colnames(f$x$a))
  expect_error(.lc_fit(f), "feature names in order")
  colnames(f$link$maps[[1]]) <- colnames(f$x$a)
  expect_s3_class(.lc_fit(f), "ssma_align")
  names(f$x) <- c("a", "a")
  expect_error(.lc_fit(f), "unique")
})

test_that("SSMA constraints match observation-link energy without changing normalization", {
  f <- .lc_fixture()
  rows <- c(2, 4, 6, 8)
  map <- linear_correspondence(c("a", "b"), list(f$x$a[rows, ], f$x$b[rows + 1, ]),
                              weights = c(.2, .7, 1.4, 2), reduction = "sum")
  extra <- data.frame(domain_i = 1L, index_i = rows, domain_j = 2L,
                      index_j = rows + 1, weight = c(.2, .7, 1.4, 2))
  a <- .lc_fit(f, list(map))
  f$rows <- rbind(f$rows, extra)
  b <- .lc_fit(f, list())
  expect_equal(a$generalized_eigenvalues, b$generalized_eigenvalues, tolerance = 1e-9)
  expect_equal(.lc_projector(a$s), .lc_projector(b$s), tolerance = 1e-8)
})

test_that("SSMA constrained dense and operator solvers agree on unequal domains", {
  f <- .lc_fixture()
  terms <- list(f$link, linear_penalty("a", diff(diag(5)), strength = .2))
  a <- .lc_fit(f, terms)
  b <- .lc_fit(f, terms, solver = "operator")
  expect_equal(a$generalized_eigenvalues, b$generalized_eigenvalues, tolerance = 1e-7)
  expect_equal(.lc_projector(a$s), .lc_projector(b$s), tolerance = 1e-6)
  expect_equal(a$constraint_diagnostics[[1]]$objective, b$constraint_diagnostics[[1]]$objective, tolerance = 1e-7)
})

test_that("SSMA feature constraints work without automatically discovering sample links", {
  f <- .lc_fixture()
  fit <- ssma_align(f$x, constraints = list(f$link), preproc = NULL, ncomp = 3,
                    control = ssma_align_control(compression = "none", auto_correspondences = FALSE))
  expect_equal(nrow(fit$correspondences), 0L)
  expect_equal(fit$domain_names, c("a", "b"))
  expect_true(all(is.finite(fit$s)))
  expect_gt(fit$constraint_diagnostics[[1]]$objective, 0)
})

test_that("SSMA constraints are equivariant to feature permutations and orthogonal bases", {
  f <- .lc_fixture(); a <- .lc_fit(f)
  g <- f
  perms <- list(c(5, 2, 1, 3, 4), c(7, 1, 3, 2, 4, 6, 5))
  for (i in 1:2) {
    g$x[[i]] <- f$x[[i]][, perms[[i]]]
    g$link$maps[[i]] <- f$link$maps[[i]][, perms[[i]]]
  }
  b <- .lc_fit(g)
  expect_equal(.lc_projector(a$s), .lc_projector(b$s), tolerance = 1e-8)
  set.seed(251)
  for (i in 1:2) {
    p <- ncol(f$x[[i]]); q <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
    g$x[[i]] <- f$x[[i]] %*% q
    g$link$maps[[i]] <- f$link$maps[[i]] %*% q
  }
  c <- .lc_fit(g)
  expect_equal(.lc_projector(a$s), .lc_projector(c$s), tolerance = 1e-8)
  expect_equal(a$generalized_eigenvalues, c$generalized_eigenvalues, tolerance = 1e-9)
})

test_that("SSMA constraints preserve domain identity when domains and rows are reordered", {
  f <- .lc_fixture(); a <- .lc_fit(f)
  g <- f; g$x <- rev(g$x)
  g$rows$domain_i <- 2L; g$rows$domain_j <- 1L
  b <- .lc_fit(g)
  order_back <- c(28:51, 1:27)
  expect_equal(.lc_projector(a$s), .lc_projector(b$s[order_back, ]), tolerance = 1e-8)
  set.seed(56); pa <- sample(24); pb <- sample(27)
  g <- f; g$x$a <- f$x$a[pa, ]; g$x$b <- f$x$b[pb, ]
  g$rows$index_i <- match(f$rows$index_i, pa)
  g$rows$index_j <- match(f$rows$index_j, pb)
  b <- .lc_fit(g)
  expect_equal(.lc_projector(a$s), .lc_projector(b$s[c(order(pa), 24 + order(pb)), ]), tolerance = 1e-8)
})

test_that("SSMA reduced feature maps compose and diagnostics report actual fitted residuals", {
  f <- .lc_fixture()
  fit <- .lc_fit(f, compression = "svd", rank_per_domain = 2)
  expect_equal(fit$reduced_ranks, c(2L, 2L))
  expect_equal(dim(fit$compression_bases[[1]]), c(5L, 2L))
  w <- .lc_weights(fit)
  expect_equal(as.matrix(fit$s), rbind(f$x$a %*% w[1:5, ], f$x$b %*% w[6:12, ]), tolerance = 1e-10)
  residual <- f$h[[1]] %*% w[1:5, ] - f$h[[2]] %*% w[6:12, ]
  expected <- 2 * sum((c(1, 0, 2, 3) / 6) * rowSums(residual^2))
  expect_equal(fit$constraint_diagnostics[[1]]$objective, expected, tolerance = 1e-10)
  expect_gt(min(fit$constraint_diagnostics[[1]]$discarded_energy_fraction), 0)
})

test_that("stronger contradictory constraints reduce their optimized total residual", {
  f <- .lc_fixture()
  # Opposing maps cannot both agree except by suppressing their response.
  other <- linear_correspondence(c("a", "b"), list(f$h[[1]], -f$h[[2]]), name = "opposing")
  weak <- .lc_fit(f, list(f$link, other))
  f$link$strength <- 20; other$strength <- 10
  strong <- .lc_fit(f, list(f$link, other))
  energy <- function(fit) sum(vapply(fit$constraint_diagnostics, function(x) x$objective, numeric(1)))
  expect_lt(energy(strong) / 10, energy(weak))
  expect_true(all(is.finite(strong$s)))
})

test_that("SSMA constrained eigenproblem matches a tiny independently assembled oracle", {
  set.seed(923)
  x <- list(matrix(rnorm(9 * 4), 9, 4), matrix(rnorm(10 * 3), 10, 3))
  domains <- lapply(x, function(a) list(x = a))
  graphs <- lapply(c(9, 10), function(n) {
    a <- diag(1, n)[, c(2:n, 1)]
    a <- a + t(a)
    list(A = Matrix::Matrix(a, sparse = TRUE), degree = rowSums(a))
  })
  # Explicit edge-incidence construction, independent of SSMA's Laplacian code.
  L <- matrix(0, 19, 19)
  for (i in 1:2) {
    off <- c(0, 9)[i]; n <- nrow(x[[i]])
    for (j in seq_len(n)) {
      e <- numeric(19); e[off + j] <- 1; e[off + (j %% n) + 1] <- -1
      L <- L + tcrossprod(e)
    }
  }
  e <- numeric(19); e[1] <- 1; e[10] <- -1
  L <- L + .7 * tcrossprod(e)
  h <- list(matrix(rnorm(3 * 4), 3, 4), matrix(rnorm(3 * 3), 3, 3))
  constraint <- linear_correspondence(c("a", "b"), h, weights = c(1, 2, 4), strength = .3)
  resolved <- .resolve_linear_constraints(list(constraint), domains, c("a", "b"))
  actual <- .ssma_build_reduced_system(domains, graphs, Matrix::Matrix(L, sparse = TRUE),
                                       128, 1e-6, compression = "none", constraints = resolved)
  X <- as.matrix(Matrix::bdiag(x))
  A <- crossprod(X, L %*% X) + .lc_oracle(h, c(1, 2, 4) / 7, .3)
  B <- crossprod(X, X * 2)
  expect_equal(actual$A_r, A, tolerance = 1e-12)
  expect_equal(actual$B_r, B, tolerance = 1e-12)
  fit <- .ssma_solve_reduced_geigen(actual$A_r, actual$B_r, 2, 1e-8, 1e-6)
  # Symmetric inverse square root from a separate eigen decomposition.
  eb <- eigen(B + diag(1e-6, 7), symmetric = TRUE)
  invhalf <- eb$vectors %*% diag(1 / sqrt(eb$values)) %*% t(eb$vectors)
  oracle <- eigen(invhalf %*% A %*% invhalf, symmetric = TRUE)
  ix <- order(oracle$values)[1:2]
  ref <- invhalf %*% oracle$vectors[, ix]
  expect_equal(fit$eigenvalues, oracle$values[ix], tolerance = 1e-10)
  expect_equal(.lc_projector(fit$alpha), .lc_projector(ref), tolerance = 1e-9)
  expect_equal(A %*% fit$alpha, (B + diag(1e-6, 7)) %*% fit$alpha %*% diag(fit$eigenvalues), tolerance = 1e-9)
})

test_that("true feature links complement sparse sample anchors on held-out synthetic samples", {
  set.seed(781)
  n <- 70; p <- 6
  x <- matrix(rnorm(n * p), n, p) %*% diag(seq(2, .5, length.out = p))
  rotation <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  y <- x %*% rotation + matrix(rnorm(n * p, sd = .12), n, p)
  rows <- data.frame(domain_i = 1L, index_i = 1:3, domain_j = 2L, index_j = 1:3)
  true <- linear_correspondence(c("a", "b"), list(diag(p), rotation), strength = 1e5)
  wrong <- linear_correspondence(c("a", "b"), list(diag(p), -rotation), strength = 1e5)
  fit <- function(ct) ssma_align(list(a = x, b = y), correspondences = rows,
    constraints = ct, preproc = NULL, ncomp = 3,
    control = ssma_align_control(compression = "none", sigma = 2, knn = 12))
  error <- function(z) {
    a <- z$s[4:n, ]; b <- z$s[n + (4:n), ]
    sum((a - b)^2) / (sum(a^2) + sum(b^2))
  }
  baseline <- fit(list()); aligned <- fit(list(true)); negative <- fit(list(wrong))
  expect_lt(error(aligned), .5 * error(baseline))
  expect_lt(error(aligned), .2 * error(negative))
})

test_that("linear constraints remain sparse with very wide externally reduced interfaces", {
  p <- 100000L
  h <- Matrix::sparseMatrix(i = 1:3, j = p - 2:0, x = 1, dims = c(3L, p))
  link <- linear_correspondence(c("a", "b"), list(h, h), weights = c(1, 2, 3))
  link$ids <- 1:2
  terms <- .compile_linear_constraints(list(link), list(list(P = NULL), list(P = NULL)))
  expect_s4_class(terms[[1]]$maps[[1]], "sparseMatrix")
  expect_lt(as.numeric(object.size(terms)), 2e6)
  a <- matrix(0, p, 1); b <- a
  a[p - 2:0, 1] <- c(1, 3, 5)
  b[p - 2:0, 1] <- c(3, 2, 1)
  out <- .linear_constraint_apply(terms, list(a, b))
  expected <- (c(1, 2, 3) / 6) * c(-2, 1, 4)
  expect_equal(out[[1]][p - 2:0, 1], expected, tolerance = 1e-12)
  expect_equal(out[[2]][p - 2:0, 1], -expected, tolerance = 1e-12)
  expect_equal(sum(abs(out[[1]][1:(p - 3), 1])), 0)
})

test_that("constraints compose additively across three domains and honor processed coordinates", {
  f <- .lc_fixture()
  set.seed(751); f$x$c <- matrix(rnorm(21 * 3), 21, 3)
  h3 <- matrix(rnorm(4 * 3), 4, 3)
  third <- linear_correspondence(c("b", "c"), list(f$h[[2]], h3), name = "second-pair")
  penalty <- linear_penalty("c", diff(diag(3)), name = "third-smooth")
  terms <- list(f$link, third, penalty)
  domains <- lapply(f$x, function(x) list(x = x))
  basis <- lapply(f$x, function(x) list(P = NULL))
  compiled <- .compile_linear_constraints(.resolve_linear_constraints(terms, domains, names(domains)), basis)
  A <- matrix(0, 15, 15)
  A[1:12, 1:12] <- .lc_oracle(f$h, c(1, 0, 2, 3) / 6, 2)
  A[6:15, 6:15] <- A[6:15, 6:15] + .lc_oracle(list(f$h[[2]], h3), rep(.25, 4))
  A[13:15, 13:15] <- A[13:15, 13:15] + .lc_oracle(list(diff(diag(3))), c(.5, .5))
  expect_equal(.add_linear_constraint_matrix(matrix(0, 15, 15), compiled, c(5, 7, 3)), A, tolerance = 1e-12)
  blocks <- list(matrix(rnorm(5 * 2), 5), matrix(rnorm(7 * 2), 7), matrix(rnorm(3 * 2), 3))
  expect_equal(do.call(rbind, .linear_constraint_apply(compiled, blocks)), A %*% do.call(rbind, blocks), tolerance = 1e-12)
  # Same processed-space maps under implicit versus explicit preprocessing.
  process <- function(x) sweep(x, 2, colMeans(x)) * 2
  g <- .lc_fixture(); h <- g; h$x <- lapply(g$x, process)
  a <- ssma_align(g$x, preproc = process, correspondences = g$rows,
                  constraints = list(g$link), ncomp = 3,
                  control = ssma_align_control(compression = "none", knn = 8, sigma = 2))
  b <- .lc_fit(h)
  expect_equal(.lc_projector(a$s), .lc_projector(b$s), tolerance = 1e-9)
})
