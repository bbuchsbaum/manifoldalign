# Exactness and feasibility of the FPGW linear oracles (network simplex) and
# of the Frank-Wolfe solver for each FPGW variant.

ns_ot <- function(...) manifoldalign:::network_simplex_ot_cpp(...)

# Dense lpSolve reference for small problems:
#   min <C, P>  s.t.  rows (dir_r) a, cols (dir_c) b, [sum(P) == mass]
lp_reference <- function(C, a, b, dir_r = "=", dir_c = "=", mass = NULL) {
  n <- nrow(C); m <- ncol(C)
  A <- rbind(
    t(vapply(seq_len(n), function(i) as.numeric(row(C) == i), numeric(n * m))),
    t(vapply(seq_len(m), function(j) as.numeric(col(C) == j), numeric(n * m)))
  )
  dirs <- c(rep(dir_r, n), rep(dir_c, m))
  rhs <- c(a, b)
  if (!is.null(mass)) {
    A <- rbind(A, rep(1, n * m)); dirs <- c(dirs, "="); rhs <- c(rhs, mass)
  }
  res <- lpSolve::lp("min", as.vector(C), A, dirs, rhs)
  stopifnot(res$status == 0)
  res$objval
}

test_that("network simplex matches an LP oracle and is exact on the marginals", {
  skip_if_not_installed("lpSolve")
  set.seed(11)
  for (rep in 1:60) {
    n <- sample(1:8, 1); m <- sample(1:8, 1)
    C <- matrix(rnorm(n * m), n, m) * 10^runif(1, -3, 3)  # includes negative costs
    a <- rexp(n); a <- a / sum(a)
    b <- rexp(m); b <- b / sum(b)
    P <- ns_ot(C, a, b)
    ref <- lp_reference(C, a, b)
    expect_lte(abs(sum(C * P) - ref), 1e-9 * max(1, abs(ref)))
    expect_true(all(P >= 0))
    expect_lt(max(abs(rowSums(P) - a)), 1e-14)
    expect_lt(max(abs(colSums(P) - b)), 1e-14)
    # A basic (vertex) solution has at most n + m - 1 non-zeros; the old
    # multiplicative rebalancing smeared mass off the vertex.
    expect_lte(sum(P > 0), n + m - 1)
  }
})

test_that("network simplex keeps small masses (no fixed 1e6 quantisation)", {
  a <- c(4e-7, 1 - 4e-7)
  b <- c(0.5, 0.5)
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  P <- ns_ot(C, a, b)
  # The old adapter rounded 4e-7 * 1e6 to 0 and then rebalanced.
  # Exact up to one quantisation unit (2^-50 of the total mass, ~8.9e-16)
  expect_lt(max(abs(rowSums(P) - a)), 2e-15)
  expect_lt(max(abs(colSums(P) - b)), 2e-15)
  expect_lt(abs(P[1, 1] - 4e-7), 2e-15)
  expect_equal(P[1, 2], 0)

  # Masses spanning many orders of magnitude, plus an exact zero
  a <- c(1e-12, 3e-9, 0, 0.25, 0.75 - 1e-12 - 3e-9)
  b <- c(0.1, 0.2, 0.3, 0.4)
  set.seed(5)
  P <- ns_ot(matrix(runif(20), 5, 4), a, b)
  expect_true(all(abs(rowSums(P) - a) <= 2e-15))
  expect_lt(max(abs(colSums(P) - b)), 2e-15)
  expect_true(all(P[3, ] == 0))
})

test_that("network simplex solves odd sizes and degenerate inputs exactly", {
  expect_equal(ns_ot(matrix(3, 1, 1), 1, 1), matrix(1, 1, 1))
  P <- ns_ot(matrix(c(1, 2, 3), 1, 3), 1, c(0.2, 0.3, 0.5))
  expect_lt(max(abs(as.vector(P) - c(0.2, 0.3, 0.5))), 2e-15)
  # constant cost: any feasible plan is optimal; marginals still exact
  P <- ns_ot(matrix(7, 3, 5), rep(1/3, 3), rep(1/5, 5))
  expect_lt(max(abs(rowSums(P) - 1/3)), 2e-15)
  expect_lt(max(abs(colSums(P) - 1/5)), 2e-15)
  expect_error(ns_ot(matrix(c(1, NA, 1, 1), 2, 2), c(.5, .5), c(.5, .5)), "finite")
  expect_error(ns_ot(matrix(1, 2, 2), c(.5, .5), c(.7, .7)), "equal mass")
})

test_that("mass-constrained oracle is exact, respects bounds and mass", {
  skip_if_not_installed("lpSolve")
  set.seed(21)
  for (rep in 1:30) {
    n <- sample(2:7, 1); m <- sample(2:7, 1)
    C <- matrix(rnorm(n * m), n, m)
    a <- rexp(n); a <- a / sum(a) * runif(1, 0.5, 1.5)
    b <- rexp(m); b <- b / sum(b) * runif(1, 0.5, 1.5)
    rho <- runif(1, 0.05, 1) * min(sum(a), sum(b))
    if (rep == 1) rho <- min(sum(a), sum(b))  # full budget with unequal totals
    G <- manifoldalign:::partial_ot_mass(C, a, b, rho)
    ref <- lp_reference(C, a, b, "<=", "<=", mass = rho)
    expect_lte(abs(sum(C * G) - ref), 1e-9 * max(1, abs(ref)))
    expect_lt(abs(sum(G) - rho), 1e-12)
    expect_lte(max(rowSums(G) - a), 1e-14)
    expect_lte(max(colSums(G) - b), 1e-14)
    expect_true(all(G >= 0))
    # the lpSolve reference implementation agrees
    G_lp <- manifoldalign:::partial_ot_mass_lp_r(C, a, b, rho)
    expect_equal(sum(C * G_lp), ref, tolerance = 1e-8)
  }
  expect_error(manifoldalign:::partial_ot_mass(matrix(1, 2, 2), c(.5, .5), c(.5, .5), 1.5),
               "exceeds")
})

test_that("penalized-variant oracle solves the <= polytope LP exactly", {
  skip_if_not_installed("lpSolve")
  set.seed(31)
  for (rep in 1:30) {
    n <- sample(2:7, 1); m <- sample(2:7, 1)
    C <- matrix(rnorm(n * m), n, m)  # mixed signs: only negative entries pay off
    a <- rep(1 / n, n); b <- rep(1 / m, m)
    G <- manifoldalign:::partial_ot_tv(C, a, b)
    ref <- lp_reference(C, a, b, "<=", "<=")
    expect_lte(abs(sum(C * G) - ref), 1e-9 * max(1, abs(ref)))
    expect_true(all(G >= 0))
    expect_lte(max(rowSums(G) - a), 1e-14)
    expect_lte(max(colSums(G) - b), 1e-14)
  }
  # Non-negative costs: transporting nothing is optimal
  G <- manifoldalign:::partial_ot_tv(matrix(runif(9) + 0.1, 3, 3), rep(1/3, 3), rep(1/3, 3))
  expect_equal(sum(G), 0)
})

fw_problem <- function(seed, n1 = 6, n2 = 5) {
  set.seed(seed)
  X <- matrix(rnorm(n1 * 2), n1, 2)
  Y <- matrix(rnorm(n2 * 2), n2, 2)
  Cx <- as.matrix(dist(X)); Cx <- Cx / max(Cx)
  Cy <- as.matrix(dist(Y)); Cy <- Cy / max(Cy)
  list(C = manifoldalign:::compute_feature_cost(X, Y), Cx = Cx, Cy = Cy,
       p = rep(1 / n1, n1), q = rep(1 / n2, n2))
}

test_that("Frank-Wolfe never increases the objective and stays feasible (all variants)", {
  for (seed in 1:5) {
    pr <- fw_problem(seed)
    P_bal <- outer(pr$p, pr$q)
    for (cfg in list(list(lambda = 0, rho = NULL, P0 = P_bal),
                     list(lambda = 0, rho = 0.6, P0 = 0.6 * P_bal),
                     list(lambda = 0.2, rho = NULL, P0 = P_bal),
                     list(lambda = 0.02, rho = NULL, P0 = 0.5 * P_bal))) {
      sol <- manifoldalign:::fw_fpgw(pr$C, pr$Cx, pr$Cy, pr$p, pr$q,
                                     omega1 = 0.3, lambda = cfg$lambda, rho = cfg$rho,
                                     P0 = cfg$P0, max_iter = 300, tol = 1e-10)
      tr <- attr(sol, "obj_trace")
      expect_true(all(diff(tr) <= 1e-12 * max(1, abs(tr))),
                  label = paste("monotone trace", sol$variant, seed))
      P <- sol$P
      expect_true(all(P >= 0))
      expect_lte(max(rowSums(P) - pr$p), 1e-12)
      expect_lte(max(colSums(P) - pr$q), 1e-12)
      if (!is.null(cfg$rho)) expect_lt(abs(sum(P) - cfg$rho), 1e-12)
      if (sol$variant == "balanced") {
        expect_lt(max(abs(rowSums(P) - pr$p)), 1e-12)
        expect_lt(max(abs(colSums(P) - pr$q)), 1e-12)
      }
      if (sol$converged) expect_lte(sol$gap, 1e-10)
    }
  }
})

test_that("penalized variant: gradient and step use the correct mass-penalty terms", {
  pr <- fw_problem(8)
  lambda <- 0.3
  omega1 <- 0.3
  obj <- function(P) {
    M <- outer(as.vector(pr$Cx^2 %*% rowSums(P)), as.vector(pr$Cy^2 %*% colSums(P)), "+") -
      2 * pr$Cx %*% P %*% t(pr$Cy)
    omega1 * sum(pr$C * P) + (1 - omega1) * sum(M * P) + lambda * (2 - 2 * sum(P)^2)
  }
  P0 <- 0.4 * outer(pr$p, pr$q)
  sol <- manifoldalign:::fw_fpgw(pr$C, pr$Cx, pr$Cy, pr$p, pr$q, omega1 = omega1,
                                 lambda = lambda, rho = NULL, P0 = P0, max_iter = 0)
  # Gradient at P0 matches finite differences of the documented objective
  G <- sol$grad
  h <- 1e-6
  for (k in c(1, 7, 13, 30)) {
    E <- matrix(0, nrow(P0), ncol(P0)); E[k] <- h
    expect_equal(G[k], (obj(P0 + E) - obj(P0 - E)) / (2 * h), tolerance = 1e-6)
  }
  # Reported objective equals the documented objective
  expect_equal(sol$objective, obj(P0), tolerance = 1e-12)
  # Full solve: objective at the returned plan is below the start
  sol2 <- manifoldalign:::fw_fpgw(pr$C, pr$Cx, pr$Cy, pr$p, pr$q, omega1 = omega1,
                                  lambda = lambda, rho = NULL, P0 = P0, max_iter = 200,
                                  tol = 1e-10)
  expect_equal(sol2$objective, obj(sol2$P), tolerance = 1e-12)
  expect_lt(sol2$objective, obj(P0))
})

test_that("mass-constrained FPGW never reports convergence at an infeasible mass", {
  skip_if_not_installed("multidesign")
  # Identical, degenerate domains: all costs are zero, so the FW gap is zero
  # immediately. The old solver started from a full-mass plan and returned it
  # with converged = TRUE although rho = 0.5.
  X <- matrix(1, 6, 2)
  design <- data.frame(id = 1:6)
  hd <- multidesign::hyperdesign(list(
    a = multidesign::multidesign(X, design),
    b = multidesign::multidesign(X, design)
  ))
  for (eps in c(0, 0.01)) {
    res <- fpgw(hd, rho = 0.5, epsilon = eps)
    P <- res$transport_plans[[1]]
    expect_lt(abs(sum(P) - 0.5), 1e-12)
    expect_lte(max(rowSums(P) - 1/6), 1e-14)
    expect_true(res$converged[1, 2])
  }
})

test_that("fw_fpgw refuses an infeasible starting plan", {
  pr <- fw_problem(3)
  P_full <- outer(pr$p, pr$q)
  expect_error(
    manifoldalign:::fw_fpgw(pr$C, pr$Cx, pr$Cy, pr$p, pr$q, omega1 = 0.5,
                            lambda = 0, rho = 0.5, P0 = P_full),
    "infeasible"
  )
  expect_error(
    manifoldalign:::fw_fpgw(pr$C, pr$Cx, pr$Cy, pr$p, pr$q, omega1 = 0.5,
                            lambda = 0, rho = NULL, P0 = 0.5 * P_full),
    "infeasible"
  )
})

test_that("fpgw validates its arguments", {
  skip_if_not_installed("multidesign")
  set.seed(1)
  design <- data.frame(id = 1:5)
  hd <- multidesign::hyperdesign(list(
    a = multidesign::multidesign(matrix(rnorm(10), 5, 2), design),
    b = multidesign::multidesign(matrix(rnorm(10), 5, 2), design)
  ))
  expect_error(fpgw(hd, omega1 = 1.5), "omega1")
  expect_error(fpgw(hd, omega1 = -0.1), "omega1")
  expect_error(fpgw(hd, rho = 1.5), "rho")
  expect_error(fpgw(hd, rho = -0.2), "rho")
  expect_error(fpgw(hd, epsilon = -1), "epsilon")
  expect_error(fpgw(hd, max_iter = 0), "max_iter")
  expect_error(fpgw(hd, tol = 0), "tol")
  expect_error(fpgw(hd, lambda = -1), "lambda")
  expect_error(fpgw(hd, metric = "chebyshev"), "metric")
  expect_error(fpgw(hd, omgea1 = 0.5), "unused argument")
  expect_error(fpgw(hd, rho = 0.5, lambda = 0.1), "not both")
})

test_that("LP fallbacks guard their size and fail loudly", {
  skip_if_not_installed("lpSolve")
  old <- options(manifoldalign.lp_max_variables = 10)
  on.exit(options(old), add = TRUE)
  expect_error(manifoldalign:::classical_ot_lp_r(matrix(0, 4, 4), rep(.25, 4), rep(.25, 4)),
               "lp_max_variables")
  expect_error(manifoldalign:::partial_ot_mass_lp_r(matrix(0, 4, 4), rep(.25, 4), rep(.25, 4), .5),
               "lp_max_variables")
  options(old)
  # Infeasible LP (unequal masses) is an error, not a warning with garbage
  expect_error(manifoldalign:::classical_ot_lp_r(matrix(0, 2, 2), c(.5, .5), c(1, 1)),
               "lpSolve failed")
  P <- manifoldalign:::classical_ot_lp_r(matrix(c(0, 1, 1, 0), 2), c(.5, .5), c(.5, .5))
  expect_equal(P, diag(0.5, 2))
})

test_that("dense pairwise OT memory estimate and guard", {
  # 30 domains x 2,000 rows: 13.92 GB of plans alone, plus 0.96 GB of
  # within-domain distances and the working matrices of one pair
  est <- manifoldalign:::.ot_dense_bytes(rep(2000, 30))
  expect_equal(est, 8 * (435 * 4e6 + 30 * 4e6 + 6 * 4e6))
  expect_error(manifoldalign:::.ot_check_dense_memory(rep(2000, 30), "fpgw"),
               "ot_max_dense_bytes")
  expect_silent(manifoldalign:::.ot_check_dense_memory(rep(50, 3), "fpgw"))
  skip_if_not_installed("multidesign")
  design <- data.frame(id = 1:5)
  hd <- multidesign::hyperdesign(list(
    a = multidesign::multidesign(matrix(rnorm(10), 5, 2), design),
    b = multidesign::multidesign(matrix(rnorm(10), 5, 2), design)
  ))
  old <- options(manifoldalign.ot_max_dense_bytes = 100)
  on.exit(options(old), add = TRUE)
  expect_error(fpgw(hd), "ot_max_dense_bytes")
})
