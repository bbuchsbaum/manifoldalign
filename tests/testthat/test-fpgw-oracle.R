library(testthat)
library(manifoldalign)

fpgw_oracle_tolerance <- function(...) {
  values <- abs(as.numeric(unlist(list(...), recursive = TRUE)))
  2e-8 * max(1, values[is.finite(values)], na.rm = TRUE)
}

fpgw_oracle_units <- function(x, quantum, name) {
  raw <- x / quantum
  units <- as.integer(round(raw))
  if (any(abs(raw - units) > 1e-10)) {
    stop(name, " must be an integer multiple of quantum.", call. = FALSE)
  }
  units
}

# Exhaustive integer-flow enumeration is independent of both the package's
# LEMON transport backend and lpSolve. With integer capacities, total
# unimodularity guarantees that a linear transport optimum is among these
# integer vertices.
fpgw_enumerate_tiny_transport <- function(p, q, quantum, rho = NULL) {
  row_cap <- fpgw_oracle_units(p, quantum, "p")
  col_cap <- fpgw_oracle_units(q, quantum, "q")
  balanced <- is.null(rho)
  total <- if (balanced) {
    if (sum(row_cap) != sum(col_cap)) {
      stop("Balanced marginal totals differ.", call. = FALSE)
    }
    sum(row_cap)
  } else {
    fpgw_oracle_units(rho, quantum, "rho")
  }
  if (total > min(sum(row_cap), sum(col_cap))) {
    stop("Requested mass exceeds capacity.", call. = FALSE)
  }

  n <- length(row_cap)
  m <- length(col_cap)
  values <- integer(n * m)
  plans <- list()

  visit <- function(index, remaining, row_used, col_used) {
    if (index > length(values)) {
      feasible <- remaining == 0L
      if (balanced) {
        feasible <- feasible &&
          identical(row_used, row_cap) &&
          identical(col_used, col_cap)
      }
      if (feasible) {
        plans[[length(plans) + 1L]] <<- matrix(values, n, m) * quantum
      }
      return(invisible(NULL))
    }

    if (
      remaining < 0L ||
        remaining > sum(row_cap - row_used) ||
        remaining > sum(col_cap - col_used)
    ) {
      return(invisible(NULL))
    }

    i <- ((index - 1L) %% n) + 1L
    j <- ((index - 1L) %/% n) + 1L
    upper <- min(
      remaining,
      row_cap[[i]] - row_used[[i]],
      col_cap[[j]] - col_used[[j]]
    )
    if (upper < 0L) {
      return(invisible(NULL))
    }
    for (amount in 0:upper) {
      values[[index]] <<- amount
      next_rows <- row_used
      next_cols <- col_used
      next_rows[[i]] <- next_rows[[i]] + amount
      next_cols[[j]] <- next_cols[[j]] + amount
      visit(index + 1L, remaining - amount, next_rows, next_cols)
    }
    values[[index]] <<- 0L
    invisible(NULL)
  }

  visit(1L, total, integer(n), integer(m))
  if (!length(plans)) {
    stop("No feasible tiny transport plans were enumerated.", call. = FALSE)
  }
  plans
}

fpgw_exact_linear_oracle <- function(cost, p, q, quantum, rho = NULL) {
  plans <- fpgw_enumerate_tiny_transport(p, q, quantum, rho)
  objective <- vapply(plans, function(plan) sum(cost * plan), numeric(1))
  best <- min(objective)
  tolerance <- fpgw_oracle_tolerance(cost, best)
  list(
    objective = best,
    plans = plans[objective <= best + tolerance],
    n_feasible = length(plans),
    tolerance = tolerance
  )
}

fpgw_four_index_objective <- function(P, feature_cost, Cx, Cy, omega1) {
  feature <- sum(feature_cost * P)
  structure <- 0
  for (i in seq_len(nrow(P))) {
    for (j in seq_len(ncol(P))) {
      for (k in seq_len(nrow(P))) {
        for (ell in seq_len(ncol(P))) {
          structure <- structure +
            (Cx[i, k] - Cy[j, ell])^2 * P[i, j] * P[k, ell]
        }
      }
    }
  }
  c(
    feature = feature,
    structure = structure,
    total = omega1 * feature + (1 - omega1) * structure
  )
}

test_that("compiled FPGW transport oracles match exhaustive weighted tiny problems", {
  cost <- matrix(c(
    0.2, 1.1,
    1.3, 0.1,
    0.9, 0.4
  ), 2, 3)
  p <- c(0.6, 0.4)
  q <- c(0.5, 0.3, 0.2)
  quantum <- 0.1

  expect_true(manifoldalign:::has_lemon_network_simplex())

  balanced_oracle <- fpgw_exact_linear_oracle(cost, p, q, quantum)
  balanced_cpp <- manifoldalign:::network_simplex_ot_cpp(cost, p, q)
  balanced_dispatch <- manifoldalign:::classical_ot_lp(cost, p, q)
  marginal_tol <- fpgw_oracle_tolerance(p, q)

  expect_gt(balanced_oracle$n_feasible, 1L)
  expect_equal(
    sum(cost * balanced_cpp),
    balanced_oracle$objective,
    tolerance = balanced_oracle$tolerance
  )
  expect_equal(balanced_dispatch, balanced_cpp, tolerance = marginal_tol)
  expect_equal(rowSums(balanced_cpp), p, tolerance = marginal_tol)
  expect_equal(colSums(balanced_cpp), q, tolerance = marginal_tol)

  rho <- 0.6
  partial_oracle <- fpgw_exact_linear_oracle(cost, p, q, quantum, rho)
  partial_cpp <- manifoldalign:::partial_ot_mass_rcpp(cost, p, q, rho)
  partial_dispatch <- manifoldalign:::partial_ot_mass(cost, p, q, rho)

  expect_gt(partial_oracle$n_feasible, 1L)
  expect_equal(
    sum(cost * partial_cpp),
    partial_oracle$objective,
    tolerance = partial_oracle$tolerance
  )
  expect_equal(partial_dispatch, partial_cpp, tolerance = marginal_tol)
  expect_equal(sum(partial_cpp), rho, tolerance = marginal_tol)
  expect_lte(max(rowSums(partial_cpp) - p), marginal_tol)
  expect_lte(max(colSums(partial_cpp) - q), marginal_tol)
  expect_gte(min(partial_cpp), -marginal_tol)
})

test_that("FPGW objective matches a four-index oracle at fusion-weight endpoints", {
  source_coordinate <- c(0, 1, 3)
  target_permutation <- c(3, 1, 2)
  target_coordinate <- source_coordinate[target_permutation]
  Cx <- as.matrix(stats::dist(source_coordinate))
  Cy <- as.matrix(stats::dist(target_coordinate))
  feature_cost <- outer(
    source_coordinate,
    target_coordinate,
    function(x, y) (x - y)^2
  )
  p <- q <- rep(1 / 3, 3)
  P0 <- outer(p, q)
  truth <- matrix(0, 3, 3)
  for (i in seq_along(source_coordinate)) {
    truth[i, which(target_coordinate == source_coordinate[[i]])] <- 1 / 3
  }

  initial_components <- fpgw_four_index_objective(
    P0, feature_cost, Cx, Cy, omega1 = 0.5
  )
  expect_gt(initial_components[["feature"]], 0)
  expect_gt(initial_components[["structure"]], 0)

  for (omega1 in c(0, 0.4, 1)) {
    solution <- manifoldalign:::fw_fpgw(
      feature_cost, Cx, Cy, p, q,
      omega1 = omega1, lambda = 0, rho = NULL,
      P0 = P0, max_iter = 20, tol = 1e-11
    )
    initial <- fpgw_four_index_objective(
      P0, feature_cost, Cx, Cy, omega1
    )
    final <- fpgw_four_index_objective(
      solution$P, feature_cost, Cx, Cy, omega1
    )
    trace <- attr(solution, "obj_trace")
    tolerance <- fpgw_oracle_tolerance(
      feature_cost, Cx, Cy, initial, final
    )

    expect_true(solution$converged, info = paste("omega1", omega1))
    expect_equal(solution$P, truth, tolerance = tolerance)
    expect_equal(trace[[1L]], initial[["total"]], tolerance = tolerance)
    expect_equal(tail(trace, 1L), final[["total"]], tolerance = tolerance)
    expect_equal(solution$distance^2, final[["total"]], tolerance = tolerance)
    expect_equal(final[["total"]], 0, tolerance = tolerance)
    expect_equal(rowSums(solution$P), p, tolerance = tolerance)
    expect_equal(colSums(solution$P), q, tolerance = tolerance)
  }
})

test_that("weighted partial FPGW is exact and target-permutation equivariant", {
  feature_cost <- matrix(c(
    0.2, 1.1,
    1.3, 0.1,
    0.9, 0.4
  ), 2, 3)
  p <- c(0.6, 0.4)
  q <- c(0.5, 0.3, 0.2)
  rho <- 0.6
  quantum <- 0.1
  P0 <- matrix(c(
    0.2, 0.2,
    0.1, 0.0,
    0.1, 0.0
  ), 2, 3)
  exact <- fpgw_exact_linear_oracle(
    feature_cost, p, q, quantum, rho
  )
  zero_x <- matrix(0, 2, 2)
  zero_y <- matrix(0, 3, 3)
  solution <- manifoldalign:::fw_fpgw(
    feature_cost, zero_x, zero_y, p, q,
    omega1 = 1, lambda = 0, rho = rho,
    P0 = P0, max_iter = 20, tol = 1e-11
  )
  tolerance <- fpgw_oracle_tolerance(feature_cost, p, q, rho)

  expect_true(solution$converged)
  expect_equal(
    sum(feature_cost * solution$P),
    exact$objective,
    tolerance = exact$tolerance
  )
  expect_equal(sum(solution$P), rho, tolerance = tolerance)
  expect_lte(max(rowSums(solution$P) - p), tolerance)
  expect_lte(max(colSums(solution$P) - q), tolerance)
  expect_gte(min(solution$P), -tolerance)
  expect_true(all(diff(attr(solution, "obj_trace")) <= tolerance))

  target_order <- c(3, 1, 2)
  permuted <- manifoldalign:::fw_fpgw(
    feature_cost[, target_order, drop = FALSE],
    zero_x,
    zero_y[target_order, target_order, drop = FALSE],
    p,
    q[target_order],
    omega1 = 1, lambda = 0, rho = rho,
    P0 = P0[, target_order, drop = FALSE],
    max_iter = 20, tol = 1e-11
  )

  expect_true(permuted$converged)
  expect_equal(
    permuted$P,
    solution$P[, target_order, drop = FALSE],
    tolerance = tolerance
  )
  expect_equal(permuted$distance, solution$distance, tolerance = tolerance)
})

test_that("withdrawn penalized FPGW cannot silently return balanced transport", {
  cost <- matrix(c(0, 1, 1, 0), 2, 2)
  geometry <- matrix(c(0, 1, 1, 0), 2, 2)
  marginal <- c(0.5, 0.5)
  balanced_start <- outer(marginal, marginal)

  expect_error(
    manifoldalign:::fw_fpgw(
      cost, geometry, geometry, marginal, marginal,
      omega1 = 0.5, lambda = 0.1, rho = NULL,
      P0 = balanced_start, max_iter = 2
    ),
    "lambda > 0.*withdrawn.*rho"
  )
  expect_error(
    manifoldalign:::partial_ot_tv(cost, marginal, marginal, lambda = 0.1),
    "withdrawn.*rho"
  )
})

test_that("partial FPGW requires a feasible exact-mass starting coupling", {
  cost <- matrix(c(0, 1, 1, 0), 2, 2)
  geometry <- matrix(c(0, 1, 1, 0), 2, 2)
  marginal <- c(0.5, 0.5)
  rho <- 0.6
  balanced_start <- outer(marginal, marginal)

  expect_error(
    manifoldalign:::fw_fpgw(
      cost, geometry, geometry, marginal, marginal,
      omega1 = 0.5, lambda = 0, rho = rho,
      P0 = balanced_start, max_iter = 2
    ),
    "P0.*feasible.*rho"
  )

  feasible_start <- rho * balanced_start
  solution <- manifoldalign:::fw_fpgw(
    cost, geometry, geometry, marginal, marginal,
    omega1 = 0.5, lambda = 0, rho = rho,
    P0 = feasible_start, max_iter = 20, tol = 1e-11
  )
  tolerance <- fpgw_oracle_tolerance(cost, marginal, rho)

  expect_equal(sum(solution$P), rho, tolerance = tolerance)
  expect_lte(max(rowSums(solution$P) - marginal), tolerance)
  expect_lte(max(colSums(solution$P) - marginal), tolerance)
  expect_true(all(diff(attr(solution, "obj_trace")) <= tolerance))
})
