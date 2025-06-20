# benchmark_parrot.R - Performance benchmarking for PARROT optimizations
# -------------------------------------------------------------------------
# Compare R vs C++ implementations across different problem sizes
# -------------------------------------------------------------------------

library(manifoldalign)
library(microbenchmark)
library(ggplot2)

# Source the test helper if needed
source("../testthat/test-parrot.R")

#' Benchmark PARROT performance across implementations
#' 
#' @param sizes Vector of network sizes to test
#' @param n_reps Number of repetitions for timing
#' @return Data frame with benchmark results
benchmark_parrot_scalability <- function(sizes = c(10, 25, 50, 100, 200), 
                                       n_reps = 10) {
  results <- list()
  
  for (n in sizes) {
    message(sprintf("Benchmarking size n=%d...", n))
    
    # Create test data
    set.seed(123)
    test_case <- create_aligned_networks(n, n_features = 5, noise_level = 0.1)
    
    # Benchmark both implementations
    mb <- microbenchmark(
      R = {
        set_parrot_use_rcpp(FALSE)
        parrot(test_case$hyperdesign, anchors = anchors, 
               tau = 0.05, max_iter = 30, tol = 1e-4)
      },
      Cpp = {
        set_parrot_use_rcpp(TRUE)
        parrot(test_case$hyperdesign, anchors = anchors,
               tau = 0.05, max_iter = 30, tol = 1e-4)
      },
      times = n_reps
    )
    
    # Store results
    results[[as.character(n)]] <- data.frame(
      size = n,
      implementation = mb$expr,
      time = mb$time / 1e9  # Convert to seconds
    )
  }
  
  # Combine results
  do.call(rbind, results)
}

#' Benchmark individual components
#' 
#' @param n Network size
#' @return List of component timings
benchmark_parrot_components <- function(n = 100) {
  message(sprintf("Component benchmarking for n=%d", n))
  
  # Create test data
  set.seed(123)
  C <- matrix(runif(n * n), n, n)
  X1 <- matrix(rnorm(n * 10), n, 10)
  X2 <- matrix(rnorm(n * 10), n, 10)
  S <- matrix(runif(n * n), n, n)
  S <- S / sum(S) * n
  
  # Create sparse adjacencies
  A1 <- Matrix::rsparsematrix(n, n, nnz = 5*n)
  A1 <- A1 + t(A1)
  A2 <- Matrix::rsparsematrix(n, n, nnz = 5*n)
  A2 <- A2 + t(A2)
  
  networks <- list(
    list(adjacency = A1, features = X1),
    list(adjacency = A2, features = X2)
  )
  
  # Benchmark Sinkhorn
  sinkhorn_bench <- microbenchmark(
    R = solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, 
                                 tol = 1e-6, use_rcpp = FALSE),
    Cpp = solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100,
                                   tol = 1e-6, use_rcpp = TRUE),
    times = 20
  )
  
  # Benchmark edge gradient
  edge_bench <- microbenchmark(
    R = compute_edge_gradient(S, networks, use_rcpp = FALSE),
    Cpp = compute_edge_gradient(S, networks, use_rcpp = TRUE),
    times = 20
  )
  
  # Benchmark distance computation
  dist_bench <- microbenchmark(
    R = {
      X1_sq <- rowSums(X1^2)
      X2_sq <- rowSums(X2^2)
      outer(X1_sq, X2_sq, "+") - 2 * X1 %*% t(X2)
    },
    Cpp = manifoldalign:::compute_squared_distances_cpp(X1, X2),
    times = 20
  )
  
  list(
    sinkhorn = summary(sinkhorn_bench),
    edge_gradient = summary(edge_bench),
    distance = summary(dist_bench)
  )
}

#' Generate performance report
#' 
#' @param results Benchmark results from benchmark_parrot_scalability
#' @param output_file Path to save plot
generate_parrot_performance_report <- function(results, 
                                             output_file = "parrot_performance.pdf") {
  # Calculate speedups
  speedup_data <- results %>%
    group_by(size) %>%
    summarise(
      time_r = mean(time[implementation == "R"]),
      time_cpp = mean(time[implementation == "Cpp"]),
      speedup = time_r / time_cpp
    )
  
  # Create plots
  p1 <- ggplot(results, aes(x = size, y = time, color = implementation)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    scale_y_log10() +
    labs(title = "PARROT Performance: R vs C++",
         x = "Network Size (nodes)",
         y = "Time (seconds, log scale)",
         color = "Implementation") +
    theme_minimal()
  
  p2 <- ggplot(speedup_data, aes(x = size, y = speedup)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    labs(title = "C++ Speedup over R",
         x = "Network Size (nodes)",
         y = "Speedup Factor") +
    theme_minimal()
  
  # Save plots
  pdf(output_file, width = 10, height = 6)
  print(p1)
  print(p2)
  dev.off()
  
  # Print summary
  cat("\nPARROT Performance Summary\n")
  cat("==========================\n")
  print(speedup_data)
  cat("\nAverage speedup:", mean(speedup_data$speedup), "x\n")
  
  invisible(speedup_data)
}

# Run benchmarks if executed directly
if (sys.nframe() == 0) {
  message("Running PARROT performance benchmarks...")
  
  # Overall scalability
  results <- benchmark_parrot_scalability()
  speedup_summary <- generate_parrot_performance_report(results)
  
  # Component analysis
  comp_results <- benchmark_parrot_components(n = 100)
  
  message("\nComponent-wise speedups:")
  for (comp in names(comp_results)) {
    r_time <- comp_results[[comp]]$mean[comp_results[[comp]]$expr == "R"]
    cpp_time <- comp_results[[comp]]$mean[comp_results[[comp]]$expr == "Cpp"]
    speedup <- r_time / cpp_time
    message(sprintf("  %s: %.1fx speedup", comp, speedup))
  }
  
  message("\nBenchmarking complete!")
}