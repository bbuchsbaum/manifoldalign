# Regression tests for the PARROT proximal solver, the shared Sinkhorn kernel
# and the compiled-kernel dispatch (backlog items B3.17 / B3.19 / call_cpp_or).

library(testthat)
library(Matrix)
library(manifoldalign)

ns_fun <- function(name) get(name, envir = asNamespace("manifoldalign"))

cpp_symbol_loaded <- function(symbol) {
  isTRUE(ns_fun(".cpp_symbol_available")(symbol))
}

parrot_hd <- function(X1, X2, anchors1, anchors2) {
  hd <- list(
    domain1 = list(x = X1, design = data.frame(node_id = seq_len(nrow(X1)), anchors = anchors1)),
    domain2 = list(x = X2, design = data.frame(node_id = seq_len(nrow(X2)), anchors = anchors2))
  )
  class(hd) <- "hyperdesign"
  hd
}

ring_pair <- function(n = 16, seed = 7, noise = 0.03, offset = 0) {
  set.seed(seed)
  theta <- 2 * pi * seq_len(n) / n
  X1 <- cbind(cos(theta), sin(theta), sin(2 * theta)) +
    matrix(rnorm(n * 3, sd = 0.02), n, 3)
  X2 <- X1 + matrix(rnorm(n * 3, sd = noise), n, 3)
  list(X1 = X1 + offset, X2 = X2 + offset)
}

# Independent log-domain Sinkhorn oracle (long run, scalar loops).
sinkhorn_oracle <- function(C, a, b, eps, iters = 20000) {
  f <- numeric(nrow(C))
  g <- numeric(ncol(C))
  lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
  for (it in seq_len(iters)) {
    for (j in seq_len(ncol(C))) g[j] <- eps * log(b[j]) - eps * lse((f - C[, j]) / eps)
    for (i in seq_len(nrow(C))) f[i] <- eps * log(a[i]) - eps * lse((g - C[i, ]) / eps)
  }
  exp(outer(f, g, "+") / eps - C / eps)
}

# -------------------------------------------------------------------------
# Shared Sinkhorn kernel
# -------------------------------------------------------------------------

test_that("sinkhorn_unified honours stabilized = TRUE for large epsilon", {
  a <- rep(1 / 4, 4)
  b <- rep(1 / 5, 5)
  # exp(-1000) underflows: the plain-kernel path would return a zero plan.
  P <- sinkhorn_unified(matrix(1000, 4, 5), a, b, epsilon = 1, stabilized = TRUE)
  expect_equal(P, outer(a, b), tolerance = 1e-12)

  set.seed(1)
  C <- matrix(runif(12, 900, 1000), 3, 4)
  a <- c(0.2, 0.3, 0.5)
  b <- c(0.1, 0.4, 0.25, 0.25)
  P <- sinkhorn_unified(C, a, b, epsilon = 2, max_iter = 5000, tol = 1e-12,
                        stabilized = TRUE)
  expect_true(all(is.finite(P)))
  expect_equal(sum(P), 1, tolerance = 1e-10)
  expect_equal(P, sinkhorn_oracle(C, a, b, 2, iters = 3000), tolerance = 1e-8)
})

test_that("sinkhorn_unified_potentials reports convergence and consistent potentials", {
  set.seed(2)
  C <- matrix(runif(30, 0, 7), 5, 6)
  a <- rep(1 / 5, 5)
  b <- rep(1 / 6, 6)
  for (stab in c(TRUE, FALSE)) {
    res <- sinkhorn_unified_potentials(C, a, b, epsilon = 0.5, max_iter = 5000,
                                       tol = 1e-12, stabilized = stab)
    expect_true(res$converged)
    expect_lt(res$marginal_error, 1e-12)
    # No post-hoc rebalancing: the plan is exactly the one implied by the potentials.
    expect_equal(res$pi, exp(-C / 0.5 + outer(as.numeric(res$log_u), as.numeric(res$log_v), "+")),
                 tolerance = 1e-12)
    expect_equal(res$pi, sinkhorn_oracle(C, a, b, 0.5, iters = 3000), tolerance = 1e-8)
  }

  short <- sinkhorn_unified_potentials(C, a, b, epsilon = 0.01, max_iter = 2,
                                       tol = 1e-12, stabilized = TRUE)
  expect_false(short$converged)
  expect_equal(short$iterations, 2L)
  expect_gt(short$marginal_error, 1e-12)
})

test_that("PARROT Sinkhorn: R and C++ solve the same unscaled problem", {
  skip_if_not(cpp_symbol_loaded("_manifoldalign_sinkhorn_unified_potentials"))
  solve <- ns_fun("solve_sinkhorn_stabilized")
  last_backend <- ns_fun(".last_cpp_backend")

  set.seed(3)
  for (eps in c(0.05, 0.3, 2)) {
    C <- matrix(runif(7 * 9, 0, 5), 7, 9)   # cost range != 1
    r <- solve(C, eps, max_iter = 5000, tol = 1e-13, use_cpp = FALSE)
    cc <- solve(C, eps, max_iter = 5000, tol = 1e-13, use_cpp = TRUE)
    expect_identical(last_backend("solve_sinkhorn_stabilized"), "cpp")
    expect_true(r$converged)
    expect_true(cc$converged)
    expect_equal(r$plan, cc$plan, tolerance = 1e-9)
    expect_equal(r$iterations, cc$iterations)
    # Raw-cost convention: matches the independent oracle on the unscaled cost.
    expect_equal(r$plan, sinkhorn_oracle(C, rep(1 / 7, 7), rep(1 / 9, 9), eps, iters = 2000),
                 tolerance = 1e-8)
  }

  # Metamorphic: scaling cost and epsilon together leaves the plan unchanged,
  # scaling only the cost does not (no hidden rescaling to [0, 1]).
  C <- matrix(runif(36, 0, 1), 6, 6)
  p1 <- solve(C, 0.1, 5000, 1e-13)$plan
  p3 <- solve(3 * C, 0.3, 5000, 1e-13)$plan
  p3b <- solve(3 * C, 0.1, 5000, 1e-13)$plan
  expect_equal(p1, p3, tolerance = 1e-9)
  expect_gt(max(abs(p1 - p3b)), 1e-3)
})

test_that("PARROT Sinkhorn respects max_iter and reports non-convergence", {
  solve_r <- ns_fun("solve_sinkhorn_stabilized_r")
  set.seed(4)
  C <- matrix(runif(64, 0, 5), 8, 8)
  res <- solve_r(C, 0.01, max_iter = 3, tol = 1e-10)
  expect_identical(res$iterations, 3L)
  expect_false(res$converged)
  expect_gt(res$marginal_error, 1e-10)
  # Returned as computed: rows exact after the last u-update, columns not.
  expect_lt(max(abs(rowSums(res$plan) - 1 / 8)), 1e-12)
})

# -------------------------------------------------------------------------
# Objective and proximal loop
# -------------------------------------------------------------------------

make_problem <- function(lambda_e = 0.3, lambda_n = 0.2, lambda_p = 0.4) {
  rp <- ring_pair(n = 10, seed = 11)
  nets <- ns_fun("extract_parrot_networks")(list(list(x = rp$X1), list(x = rp$X2)))
  vec <- c(1:3, rep(NA, 7))
  ai <- list(vec1 = vec, vec2 = vec, idx1 = 1:3, idx2 = 1:3, n1 = 10, n2 = 10)
  rwr <- ns_fun("compute_parrot_rwr")(nets, ai, 0.15, 100, 1e-8)
  C <- ns_fun("compute_parrot_cost")(nets, rwr, ai, 0.2, 0.15, 0.1)
  list(prob = ns_fun(".parrot_problem")(nets, C, ai, lambda_e, lambda_n, lambda_p),
       nets = nets, ai = ai, rwr = rwr, C = C)
}

test_that("PARROT objective gradient matches finite differences", {
  pr <- make_problem()
  prob <- pr$prob
  obj <- ns_fun(".parrot_objective")
  grad_f <- ns_fun(".parrot_smooth_gradient")

  set.seed(5)
  S <- matrix(runif(100, 0.5, 1.5), 10, 10)
  S <- S / sum(S)
  G <- grad_f(prob, S) + (prob$lambda_n + prob$lambda_p) * (log(S) + 1)
  for (k in 1:3) {
    D <- matrix(rnorm(100), 10, 10) * 1e-3
    h <- 1e-4
    fd <- (obj(prob, S + h * D) - obj(prob, S - h * D)) / (2 * h)
    expect_equal(fd, sum(G * D), tolerance = 1e-6)
  }
})

test_that("PARROT proximal loop iterates and decreases the objective", {
  pr <- make_problem()
  tr <- ns_fun("solve_parrot_proximal")(
    pr$nets, pr$rwr, pr$ai, 0.3, 0.2, 0.4, 0.05, 0.2, 0.15, 0.1, "sinkhorn",
    max_iter = 5000, max_outer = 15, tol = 1e-8
  )
  expect_gt(tr$outer_iterations, 1L)
  expect_length(tr$objective_history, tr$outer_iterations + 1L)
  expect_true(all(diff(tr$objective_history) <= 1e-12 * max(abs(tr$objective_history))))
  expect_lt(tail(tr$objective_history, 1), tr$objective_history[1])
  # Inner solves are reported honestly; this sharp problem gets close but the
  # residual is part of the result rather than hidden by rebalancing.
  expect_lt(tr$sinkhorn$marginal_error, 1e-4)
  expect_equal(max(abs(c(rowSums(tr$transport_plan) - 0.1, colSums(tr$transport_plan) - 0.1))),
               tr$sinkhorn$marginal_error, tolerance = 1e-6)

  # It is not the old single entropic OT on the RWR cost.
  single <- ns_fun("solve_sinkhorn_stabilized")(pr$C, 0.05 + 0.2 + 0.4, 2000, 1e-12)$plan
  expect_gt(max(abs(tr$transport_plan - single)), 1e-4)
  # And the returned plan beats that single solve on the actual objective.
  expect_lt(tr$final_objective, ns_fun(".parrot_objective")(pr$prob, single))
})

test_that("parrot() returns objective history and the anchor prior moves mass", {
  rp <- ring_pair(n = 16, seed = 21)
  # Declared anchor pairs (1,9), (2,10), (3,11), (4,12) contradict the
  # features (which match i to i), so only the anchor prior can pull mass there.
  anchors1 <- c(1:4, rep(NA, 12))
  anchors2 <- c(rep(NA, 8), 1:4, rep(NA, 4))
  hd <- parrot_hd(rp$X1, rp$X2, anchors1, anchors2)

  fit0 <- parrot(hd, anchors = anchors, lambda_e = 0, lambda_n = 0.05, lambda_p = 0,
                 max_iter = 5000, preproc = NULL)
  fit1 <- parrot(hd, anchors = anchors, lambda_e = 0, lambda_n = 0.05, lambda_p = 2,
                 max_iter = 5000, preproc = NULL)

  for (fit in list(fit0, fit1)) {
    expect_true(is.numeric(fit$objective_history))
    expect_gt(length(fit$objective_history), 2L)
    expect_true(all(diff(fit$objective_history) <= 1e-12 * max(abs(fit$objective_history))))
    expect_gt(fit$convergence$outer_iterations, 1L)
  }

  anchor_mass <- function(S) sum(as.matrix(S)[cbind(1:4, 9:12)])
  # A strong anchor prior concentrates mass on the anchor pairs (the old code
  # only used lambda_p as extra entropy, which spreads mass instead).
  expect_gt(anchor_mass(fit1$transport_plan), anchor_mass(fit0$transport_plan) + 0.05)
  expect_gt(anchor_mass(fit1$transport_plan), 0.9 * 4 / 16)
})

test_that("parrot() applies the default center() preprocessing", {
  rp <- ring_pair(n = 14, seed = 31, offset = 5)
  anchors <- c(1:3, rep(NA, 11))
  hd_raw <- parrot_hd(rp$X1, rp$X2, anchors, anchors)
  hd_centered <- parrot_hd(scale(rp$X1, scale = FALSE), scale(rp$X2, scale = FALSE),
                           anchors, anchors)

  fit_default <- parrot(hd_raw, anchors = anchors, max_iter = 5000)
  fit_null <- parrot(hd_raw, anchors = anchors, max_iter = 5000, preproc = NULL)
  fit_manual <- parrot(hd_centered, anchors = anchors, max_iter = 5000, preproc = NULL)

  # center() is fitted per domain: identical to centering by hand ...
  expect_equal(as.matrix(fit_default$transport_plan), as.matrix(fit_manual$transport_plan),
               tolerance = 1e-8)
  # ... and different from no preprocessing on offset data.
  expect_gt(max(abs(fit_default$transport_plan - fit_null$transport_plan)), 1e-4)
  expect_s3_class(fit_default$preproc, "pre_processor")
})

test_that("parrot() rejects solver = 'exact' and reports Sinkhorn non-convergence", {
  rp <- ring_pair(n = 12, seed = 41)
  anchors <- c(1:3, rep(NA, 9))
  hd <- parrot_hd(rp$X1, rp$X2, anchors, anchors)

  expect_error(parrot(hd, anchors = anchors, solver = "exact"), "not implemented")
  expect_error(parrot(hd, anchors = anchors, solver = "bogus"), "should be one of")

  expect_warning(
    fit <- parrot(hd, anchors = anchors, max_iter = 2, tol = 1e-12),
    "did not reach tol"
  )
  expect_false(fit$convergence$sinkhorn_converged)
  expect_gt(fit$convergence$sinkhorn_marginal_error, 1e-12)
})

test_that("R and C++ PARROT pipelines agree tightly", {
  skip_if_not(cpp_symbol_loaded("_manifoldalign_compute_parrot_cost_cpp"))
  rp <- ring_pair(n = 15, seed = 51)
  anchors <- c(1:3, rep(NA, 12))
  hd <- parrot_hd(rp$X1, rp$X2, anchors, anchors)
  fr <- parrot(hd, anchors = anchors, max_iter = 5000, tol = 1e-8, use_cpp = FALSE)
  fc <- parrot(hd, anchors = anchors, max_iter = 5000, tol = 1e-8, use_cpp = TRUE)
  last_backend <- ns_fun(".last_cpp_backend")
  expect_identical(last_backend("solve_sinkhorn_stabilized"), "cpp")
  expect_identical(last_backend("compute_parrot_cost"), "cpp")
  expect_equal(as.matrix(fr$transport_plan), as.matrix(fc$transport_plan), tolerance = 1e-6)
  expect_equal(fr$objective_history, fc$objective_history, tolerance = 1e-8)
})

# -------------------------------------------------------------------------
# Edge gradient: sparse reformulation vs the original O(|E1||E2|) loop
# -------------------------------------------------------------------------

edge_gradient_loop_oracle <- function(S, A1, A2, X1, X2) {
  grad <- matrix(0, nrow(S), ncol(S))
  e1 <- which(as.matrix(A1) != 0, arr.ind = TRUE)
  e1 <- e1[e1[, 1] < e1[, 2], , drop = FALSE]
  e2 <- which(as.matrix(A2) != 0, arr.ind = TRUE)
  e2 <- e2[e2[, 1] < e2[, 2], , drop = FALSE]
  for (p in seq_len(nrow(e1))) {
    i <- e1[p, 1]; j <- e1[p, 2]
    d1 <- sum((X1[i, ] - X1[j, ])^2)
    for (q in seq_len(nrow(e2))) {
      a <- e2[q, 1]; b <- e2[q, 2]
      d2 <- sum((X2[a, ] - X2[b, ])^2)
      diff <- 2 * (d1 - d2)
      grad[i, a] <- grad[i, a] + S[j, b] * diff
      grad[j, b] <- grad[j, b] + S[i, a] * diff
    }
  }
  grad
}

test_that("edge-consistency gradient: sparse R and C++ match the loop oracle", {
  set.seed(61)
  n1 <- 9
  n2 <- 7
  X1 <- matrix(rnorm(n1 * 3), n1, 3)
  X2 <- matrix(rnorm(n2 * 3), n2, 3)
  A1 <- Matrix::rsparsematrix(n1, n1, nnz = 20, rand.x = function(k) runif(k, 0.1, 2))
  A1 <- A1 + Matrix::t(A1)
  A2 <- Matrix::rsparsematrix(n2, n2, nnz = 12)          # asymmetric, signed
  S <- matrix(runif(n1 * n2), n1, n2)
  nets <- list(list(adjacency = A1, features = X1), list(adjacency = A2, features = X2))

  oracle <- edge_gradient_loop_oracle(S, A1, A2, X1, X2)
  expect_equal(ns_fun("compute_edge_gradient")(S, nets, use_cpp = FALSE), oracle,
               tolerance = 1e-12)
  skip_if_not(cpp_symbol_loaded("_manifoldalign_compute_edge_gradient_cpp"))
  expect_equal(ns_fun("compute_edge_gradient")(S, nets, use_cpp = TRUE), oracle,
               tolerance = 1e-12)
  expect_identical(ns_fun(".last_cpp_backend")("compute_edge_gradient"), "cpp")

  # Single-edge graphs (the old R loop dropped to a vector here).
  B1 <- Matrix::sparseMatrix(i = 1, j = 2, x = 1, dims = c(n1, n1))
  B2 <- Matrix::sparseMatrix(i = 3, j = 5, x = 1, dims = c(n2, n2))
  nets1 <- list(list(adjacency = B1, features = X1), list(adjacency = B2, features = X2))
  expect_equal(ns_fun("compute_edge_gradient")(S, nets1, use_cpp = TRUE),
               edge_gradient_loop_oracle(S, B1, B2, X1, X2), tolerance = 1e-12)
  expect_equal(ns_fun("compute_edge_gradient")(S, nets1, use_cpp = FALSE),
               edge_gradient_loop_oracle(S, B1, B2, X1, X2), tolerance = 1e-12)
})

# -------------------------------------------------------------------------
# call_cpp_or(): fall back only when the symbol is unavailable
# -------------------------------------------------------------------------

test_that("call_cpp_or re-signals compiled errors instead of rerunning R", {
  skip_if_not(cpp_symbol_loaded("_manifoldalign_compute_squared_distances_cpp"))
  sqd <- ns_fun("compute_squared_distances_cpp_fallback")
  last_backend <- ns_fun(".last_cpp_backend")
  # Mismatched feature dimensions: the compiled kernel errors; that error must surface.
  expect_error(sqd(matrix(1, 3, 2), matrix(1, 4, 3)))
  expect_identical(last_backend("compute_squared_distances"), "cpp")

  ran_fallback <- FALSE
  call_cpp_or <- ns_fun("call_cpp_or")
  expect_error(
    call_cpp_or("_manifoldalign_compute_squared_distances_cpp",
                call_cpp = function() stop("boom from compiled code"),
                fallback = function() { ran_fallback <<- TRUE; 0 }),
    "boom from compiled code"
  )
  expect_false(ran_fallback)

  # A genuinely missing symbol falls back and records the R backend.
  out <- call_cpp_or("_manifoldalign_no_such_symbol",
                     call_cpp = function() stop("should not run"),
                     fallback = function() "r-result",
                     kernel = "missing_kernel")
  expect_identical(out, "r-result")
  expect_identical(last_backend("missing_kernel"), "r")
})

test_that("R fallbacks are numerically equivalent to the compiled kernels", {
  skip_if_not(cpp_symbol_loaded("_manifoldalign_compute_squared_distances_cpp"))
  last_backend <- ns_fun(".last_cpp_backend")
  run_both <- function(kernel, f) {
    old <- options(manifoldalign.force_r_fallback = FALSE)
    on.exit(options(old))
    cpp <- f()
    expect_identical(last_backend(kernel), "cpp")
    options(manifoldalign.force_r_fallback = TRUE)
    r <- f()
    expect_identical(last_backend(kernel), "r")
    list(cpp = cpp, r = r)
  }

  set.seed(71)
  X <- matrix(rnorm(30), 10, 3)
  Y <- matrix(rnorm(21), 7, 3)
  res <- run_both("compute_squared_distances",
                  function() ns_fun("compute_squared_distances_cpp_fallback")(X, Y))
  expect_equal(res$r, res$cpp, tolerance = 1e-12)
  # Large offsets make the expanded-square formula round below zero; both clamp.
  Xo <- X + 1e6
  res <- run_both("compute_squared_distances",
                  function() ns_fun("compute_squared_distances_cpp_fallback")(Xo, Xo))
  expect_true(all(res$cpp >= 0))
  expect_true(all(res$r >= 0))

  C <- matrix(runif(48, 0, 4), 6, 8)
  res <- run_both("solve_sinkhorn_stabilized",
                  function() ns_fun("solve_sinkhorn_stabilized_cpp_fallback")(C, 0.2, 3000, 1e-12))
  expect_equal(res$r$plan, res$cpp$plan, tolerance = 1e-9)
  expect_equal(res$r$iterations, res$cpp$iterations)

  pr <- make_problem()
  W1t <- as.matrix(Matrix::t(pr$nets[[1]]$transition))
  E <- diag(10)[, 1:3]
  res <- run_both("compute_rwr_vectorized",
                  function() ns_fun("compute_rwr_vectorized_cpp_fallback")(W1t, E, 0.15, 200, 1e-10))
  expect_equal(res$r, res$cpp, tolerance = 1e-12)

  W1 <- as.matrix(pr$nets[[1]]$transition)
  W2t <- t(as.matrix(pr$nets[[2]]$transition))
  res <- run_both("solve_sylvester_rwr",
                  function() ns_fun("solve_sylvester_rwr_cpp_fallback")(W1, W2t, pr$C, 0.15, 0.1, 1e-10, 200))
  expect_equal(res$r, res$cpp, tolerance = 1e-10)

  res <- run_both("compute_parrot_cost",
                  function() ns_fun("compute_parrot_cost_cpp_fallback")(pr$nets, pr$rwr, pr$ai, 0.2, 0.15, 0.1))
  expect_equal(res$r, res$cpp, tolerance = 1e-10)
})

# -------------------------------------------------------------------------
# Validation data: anchors must point at the same node in both networks
# -------------------------------------------------------------------------

test_that("generate_parrot_validation_data maps anchors through the inverse permutation", {
  skip_if_not_installed("multidesign")
  set.seed(91)
  # permute_fraction = 1 leaves no unpermuted node, so anchors are drawn from
  # all nodes (the fallback branch) and include permuted ones.
  vd <- generate_parrot_validation_data(n_nodes = 24, n_anchors = 8, noise_level = 0,
                                        structure = "random", permute_fraction = 1)
  hd <- vd$data
  X1 <- hd[[1]]$x
  X2 <- hd[[2]]$x
  a1 <- hd[[1]]$design$anchors
  a2 <- hd[[2]]$design$anchors
  perm <- vd$ground_truth$permutation
  expect_true(any(perm[vd$ground_truth$anchor_indices$net1] !=
                    vd$ground_truth$anchor_indices$net1))
  for (val in a1[!is.na(a1)]) {
    i <- which(a1 == val)
    j <- which(a2 == val)
    expect_length(j, 1L)
    # Noise-free: the anchored rows hold the same node.
    expect_equal(X2[j, ], X1[i, ])
    expect_identical(j, vd$ground_truth$inverse_perm[i])
  }
  expect_identical(vd$ground_truth$anchor_indices$net2,
                   vd$ground_truth$inverse_perm[vd$ground_truth$anchor_indices$net1])
})
