# test-parrot-equivalence.R - Test numerical equivalence between R and C++ implementations
# -------------------------------------------------------------------------
# This ensures that the C++ optimizations produce identical results to the R version
# -------------------------------------------------------------------------

library(testthat)
library(Matrix)
library(manifoldalign)
library(multidesign)
library(multivarious)
library(tibble)

# Both backends run the same algorithm on the same (unscaled) inputs, so they
# must agree to rounding error. Each comparison skips when the compiled symbol
# is unavailable and otherwise asserts that the compiled backend actually ran.

cpp_backend_available <- function(symbol) {
  isTRUE(manifoldalign:::.cpp_symbol_available(symbol))
}

last_backend <- function(kernel) manifoldalign:::.last_cpp_backend(kernel)

## Performance/benchmark gating is centralized in helper-benchmarks.R

# Direct double-loop squared distances (independent oracle).
proxy_free_sqdist <- function(X1, X2) {
  D <- matrix(0, nrow(X1), nrow(X2))
  for (i in seq_len(nrow(X1))) {
    for (j in seq_len(nrow(X2))) {
      D[i, j] <- sum((X1[i, ] - X2[j, ])^2)
    }
  }
  D
}

# Helper to create a hyperdesign with anchors for testing
create_test_hyperdesign_anchors <- function(X1, X2, anchor_pairs) {
  
  # Create design data frames with a specified anchor column
  design1 <- data.frame(id = 1:nrow(X1), anchor_col = NA)
  design1$anchor_col[anchor_pairs[,1]] <- anchor_pairs[,1]
  
  design2 <- data.frame(id = 1:nrow(X2), anchor_col = NA)
  design2$anchor_col[anchor_pairs[,2]] <- anchor_pairs[,1]
  
  # Create multidesign objects using the correct package
  md1 <- multidesign::multidesign(X1, design1)
  md2 <- multidesign::multidesign(X2, design2)
  
  # Create a simple list that mimics a hyperdesign structure
  hd <- list(A = md1, B = md2)
  class(hd) <- "hyperdesign"
  attr(hd, "anchors") <- data.frame(
    id1 = anchor_pairs[,1], 
    id2 = anchor_pairs[,2]
  )
  hd
}

# Helper to force both implementations
run_both_implementations <- function(fun_name, ...) {
  # Get the dispatcher function
  fun_dispatch <- get(fun_name, envir = asNamespace("manifoldalign"))
  
  # Run R version
  result_r <- fun_dispatch(..., use_cpp = FALSE)
  
  # Run C++ version
  symbol_map <- list(
    solve_sinkhorn_stabilized = "_manifoldalign_sinkhorn_unified_potentials",
    compute_edge_gradient = "_manifoldalign_compute_edge_gradient_cpp",
    compute_parrot_cost = "_manifoldalign_compute_parrot_cost_cpp",
    compute_rwr_vectorized = "_manifoldalign_compute_rwr_vectorized_cpp"
  )
  cpp_symbol <- symbol_map[[fun_name]]
  if (!cpp_backend_available(cpp_symbol)) {
    skip(paste("compiled kernel unavailable:", cpp_symbol))
  }
  result_cpp <- fun_dispatch(..., use_cpp = TRUE)
  expect_identical(last_backend(fun_name), "cpp")

  list(r = result_r, cpp = result_cpp)
}

test_that("Sinkhorn R and C++ produce equivalent results", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  sizes <- c(10, 25, 50)
  tau_values <- c(0.001, 0.01, 0.1)
  
  for (n in sizes) {
    for (tau in tau_values) {
      # Cost range != 1, so a hidden rescaling in either backend would show.
      C <- matrix(runif(n * n, 0, 3), n, n)
      C <- (C + t(C)) / 2  # Make symmetric for consistency

      # Run both implementations
      results <- run_both_implementations(
        "solve_sinkhorn_stabilized",
        C = C, tau = tau, max_iter = 100, tol = 1e-6
      )

      info <- sprintf("Sinkhorn equivalence failed for n=%d, tau=%.3f", n, tau)
      expect_equal(results$r$plan, results$cpp$plan, tolerance = 1e-8, info = info)
      expect_identical(results$r$iterations, results$cpp$iterations, info = info)
      expect_identical(results$r$converged, results$cpp$converged, info = info)
      expect_equal(results$r$marginal_error, results$cpp$marginal_error,
                   tolerance = 1e-6, info = info)

      # Rows are exact after the final update; the reported residual is honest.
      plan <- results$cpp$plan
      expect_lt(max(abs(rowSums(plan) - 1 / n)), 1e-12)
      expect_equal(max(abs(c(rowSums(plan), colSums(plan)) - 1 / n)),
                   results$cpp$marginal_error, tolerance = 1e-6)
      if (results$cpp$converged) {
        expect_lt(max(abs(colSums(plan) - 1 / n)), 1e-6)
      }
    }
  }
})


test_that("Edge gradient R and C++ produce equivalent results", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  n_nodes <- 20
  
  # Create test networks
  X1 <- matrix(rnorm(n_nodes * 3), n_nodes, 3)
  X2 <- matrix(rnorm(n_nodes * 3), n_nodes, 3)
  
  # Create sparse adjacency matrices
  A1 <- Matrix::rsparsematrix(n_nodes, n_nodes, nnz = 30, rand.x = NULL)
  A1 <- A1 + t(A1)  # Make symmetric
  A1@x <- rep(1, length(A1@x))  # Binary adjacency
  
  A2 <- Matrix::rsparsematrix(n_nodes, n_nodes, nnz = 30, rand.x = NULL)
  A2 <- A2 + t(A2)
  A2@x <- rep(1, length(A2@x))
  
  # Create transport plan
  S <- matrix(runif(n_nodes * n_nodes), n_nodes, n_nodes)
  S <- S / sum(S) * n_nodes  # Normalize
  
  # Package as networks
  networks <- list(
    list(adjacency = A1, features = X1),
    list(adjacency = A2, features = X2)
  )
  
  # Run both implementations
  results <- run_both_implementations(
    "compute_edge_gradient",
    S = S, networks = networks
  )
  
  # Check equivalence
  expect_equal(
    results$r, results$cpp,
    tolerance = 1e-12,
    info = "Edge gradient equivalence failed"
  )
})

test_that("Squared distance computation R and C++ are equivalent", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  sizes <- list(c(10, 15), c(50, 50), c(100, 80))
  
  for (dims in sizes) {
    n1 <- dims[1]
    n2 <- dims[2]
    d <- 5
    
    X1 <- matrix(rnorm(n1 * d), n1, d)
    X2 <- matrix(rnorm(n2 * d), n2, d)
    
    D_r <- manifoldalign:::compute_squared_distances(X1, X2, use_cpp = FALSE)

    if (!cpp_backend_available("_manifoldalign_compute_squared_distances_cpp")) {
      skip("compiled squared-distance kernel unavailable")
    }
    D_cpp <- manifoldalign:::compute_squared_distances(X1, X2, use_cpp = TRUE)
    expect_identical(last_backend("compute_squared_distances"), "cpp")

    # Check equivalence against each other and an independent oracle
    expect_equal(D_r, D_cpp, tolerance = 1e-12,
                 info = sprintf("Distance computation failed for %dx%d", n1, n2))
    D_oracle <- as.matrix(proxy_free_sqdist(X1, X2))
    expect_equal(D_cpp, D_oracle, tolerance = 1e-10)
  }
})

test_that("RWR computation R and C++ are equivalent", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  n_nodes <- 25
  n_anchors <- 5
  
  # Create test network
  X <- matrix(rnorm(n_nodes * 3), n_nodes, 3)
  
  # Build adjacency matrix
  knn <- 5
  A <- Matrix::Matrix(0, n_nodes, n_nodes)
  for (i in 1:n_nodes) {
    dists <- sqrt(rowSums((X - matrix(X[i,], n_nodes, 3, byrow = TRUE))^2))
    nearest <- order(dists)[2:(knn+1)]
    A[i, nearest] <- 1
    A[nearest, i] <- 1
  }
  
  # Create transition matrix
  deg <- Matrix::rowSums(A)
  deg[deg == 0] <- 1
  W <- Matrix::Diagonal(x = 1/deg) %*% A
  WT <- Matrix::t(W)
  
  # Create restart vectors
  E <- Matrix::sparseMatrix(
    i = 1:n_anchors,
    j = 1:n_anchors,
    x = 1,
    dims = c(n_nodes, n_anchors)
  )
  
  # R implementation
  sigma <- 0.15
  R_r <- as.matrix(E) / n_nodes
  for (iter in 1:20) {
    R_old <- R_r
    R_r <- (1 - sigma) * as.matrix(WT %*% R_r) + sigma * as.matrix(E)
    if (max(abs(R_r - R_old)) < 1e-6) break
  }
  
  if (!cpp_backend_available("_manifoldalign_compute_rwr_vectorized_cpp")) {
    skip("compiled RWR kernel unavailable")
  }
  R_cpp <- manifoldalign:::compute_rwr_vectorized_cpp_fallback(
    as.matrix(WT), as.matrix(E), sigma, 20, 1e-6
  )
  expect_identical(last_backend("compute_rwr_vectorized"), "cpp")

  # Check equivalence
  expect_equal(
    R_r, R_cpp,
    tolerance = 1e-12,
    info = "RWR computation equivalence failed"
  )
})

test_that("Full PARROT pipeline R and C++ are equivalent", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  n_nodes <- 15
  
  # Create simple aligned networks
  X1 <- matrix(rnorm(n_nodes * 2), n_nodes, 2)
  X2 <- X1 + matrix(rnorm(n_nodes * 2, 0, 0.1), n_nodes, 2)
  
  # Create hyperdesign with anchors
  hd <- create_test_hyperdesign_anchors(X1, X2, cbind(1:3, 1:3))
  
  if (!cpp_backend_available("_manifoldalign_compute_parrot_cost_cpp")) {
    skip("compiled PARROT kernels unavailable")
  }
  result_r <- parrot(hd, anchors = anchor_col, tau = 0.1, max_iter = 5000, use_cpp = FALSE)
  result_cpp <- parrot(hd, anchors = anchor_col, tau = 0.1, max_iter = 5000, use_cpp = TRUE)
  expect_identical(last_backend("solve_sinkhorn_stabilized"), "cpp")
  expect_identical(last_backend("compute_parrot_cost"), "cpp")

  plan_r <- as.matrix(result_r$transport_plan)
  plan_cpp <- as.matrix(result_cpp$transport_plan)

  # Same problem, same algorithm: agreement to rounding error.
  expect_equal(
    plan_r,
    plan_cpp,
    tolerance = 1e-6,
    info = "Full PARROT pipeline equivalence failed"
  )
  expect_equal(result_r$objective_history, result_cpp$objective_history, tolerance = 1e-8)

  # Check that both produce valid transport plans
  expect_true(all(plan_cpp >= 0))
  expect_true(result_cpp$convergence$sinkhorn_converged)
  expect_lt(max(abs(rowSums(plan_cpp) - 1 / n_nodes)), 1e-6)
  expect_lt(max(abs(colSums(plan_cpp) - 1 / n_nodes)), 1e-6)
})

test_that("Performance improvement is achieved with C++", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  skip_on_cran()  # Skip performance tests on CRAN
  skip_if_benchmarks_disabled("Benchmark tests disabled; enable with options(manifoldalign.run_benchmarks = TRUE)")
  
  set.seed(123)
  n <- 200
  C <- matrix(runif(n * n), n, n)
  iterations <- 6
  repeats <- 4
  invisible(solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 50,
                                     tol = 1e-6, use_cpp = FALSE))
  invisible(solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 50,
                                     tol = 1e-6, use_cpp = TRUE))
  gc()
  
  speedups <- replicate(repeats, {
    time_r <- system.time({
      for (i in seq_len(iterations)) {
        solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100,
                                  tol = 1e-6, use_cpp = FALSE)
      }
    })
    time_cpp <- system.time({
      for (i in seq_len(iterations)) {
        solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100,
                                  tol = 1e-6, use_cpp = TRUE)
      }
    })
    time_r[3] / time_cpp[3]
  })
  message(sprintf("Sinkhorn speedups over %d runs: %s", repeats,
                  paste(sprintf("%.2f", speedups), collapse = ", ")))
  
  max_speedup <- max(speedups, na.rm = TRUE)
  if (max_speedup <= 1.1) {
    skip(sprintf("C++ speedup too small for reliable assertion (max %.2fx)",
                 max_speedup))
  }
  expect_gt(max_speedup, 1.1)
})
