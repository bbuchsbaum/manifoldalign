#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION") &&
      grepl("^Package:\\s*manifoldalign", paste(readLines("DESCRIPTION", n = 5L), collapse = "\\n"))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(manifoldalign)
})

pairwise_cosine <- function(S1, S2) {
  denom <- sqrt(rowSums(S1^2)) * sqrt(rowSums(S2^2)) + 1e-9
  rowSums(S1 * S2) / denom
}

mean_pair_alignment <- function(scores, n_per_domain) {
  m <- length(n_per_domain)
  starts <- cumsum(c(1L, head(n_per_domain, -1L)))
  blocks <- lapply(seq_len(m), function(i) {
    s <- starts[i]
    e <- s + n_per_domain[i] - 1L
    scores[s:e, , drop = FALSE]
  })
  pairs <- utils::combn(seq_len(m), 2, simplify = FALSE)
  vals <- vapply(pairs, function(p) {
    mean(pairwise_cosine(blocks[[p[1]]], blocks[[p[2]]]))
  }, numeric(1))
  mean(vals)
}

latent_recovery <- function(scores_block, latent) {
  C <- suppressWarnings(stats::cor(scores_block, latent))
  if (anyNA(C)) return(NA_real_)
  sum(svd(C)$d) / ncol(latent)
}

generate_hyperdesign <- function(n = 1000L,
                                 m = 3L,
                                 latent_dim = 4L,
                                 noise_sd = 0.08,
                                 seed = 77L) {
  set.seed(seed)
  dims <- seq(10L, by = 2L, length.out = m)
  Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
  domains <- lapply(seq_along(dims), function(i) {
    Wi <- matrix(rnorm(latent_dim * dims[i]), latent_dim, dims[i])
    Xi <- Z %*% Wi + matrix(rnorm(n * dims[i], sd = noise_sd), n, dims[i])
    Xi <- scale(Xi)
    list(x = Xi, design = data.frame(sample_id = seq_len(n)))
  })
  names(domains) <- paste0("domain", seq_along(domains))
  class(domains) <- c("hyperdesign", "list")
  list(hd = domains, latent = Z)
}

load_baseline_method <- function(path) {
  if (!file.exists(path)) {
    stop("Baseline file not found: ", path, call. = FALSE)
  }
  old_env <- new.env(parent = asNamespace("manifoldalign"))
  sys.source(path, envir = old_env)
  old_env$coupled_diagonalization.hyperdesign
}

run_one_fit <- function(method_fn, hd, latent, params, seed = 1L) {
  set.seed(seed)
  timing <- system.time({
    fit <- do.call(method_fn, c(list(data = hd), params))
  })
  n_per_domain <- vapply(hd, function(d) nrow(d$x), integer(1))
  first_block <- fit$s[seq_len(n_per_domain[1]), , drop = FALSE]
  list(
    elapsed_sec = as.numeric(timing[["elapsed"]]),
    final_cost = as.numeric(fit$final_cost),
    pair_cosine = mean_pair_alignment(fit$s, n_per_domain),
    latent_recovery = latent_recovery(first_block, latent),
    fit = fit
  )
}

run_sweep_method <- function(mode,
                             hd,
                             latent,
                             mu_grid,
                             base_params,
                             baseline_method,
                             current_method,
                             seed_base) {
  stopifnot(mode %in% c("baseline_uncached", "current_uncached", "current_cached"))
  fn <- if (mode == "baseline_uncached") baseline_method else current_method

  total_elapsed <- 0
  costs <- numeric(length(mu_grid))
  pair_vals <- numeric(length(mu_grid))
  latent_vals <- numeric(length(mu_grid))
  cache <- NULL

  for (idx in seq_along(mu_grid)) {
    mu <- mu_grid[[idx]]
    params <- base_params
    params$mu_coupling <- mu

    if (mode == "current_cached" && idx > 1) {
      params$spectral_cache <- cache
    }

    out <- run_one_fit(
      method_fn = fn,
      hd = hd,
      latent = latent,
      params = params,
      seed = seed_base + idx
    )

    if (mode == "current_cached" && idx == 1) {
      cache <- out$fit$diagnostics$spectral_cache
    }

    total_elapsed <- total_elapsed + out$elapsed_sec
    costs[idx] <- out$final_cost
    pair_vals[idx] <- out$pair_cosine
    latent_vals[idx] <- out$latent_recovery
  }

  data.frame(
    method = mode,
    total_elapsed_sec = total_elapsed,
    mean_final_cost = mean(costs),
    mean_pair_cosine = mean(pair_vals),
    mean_latent_recovery = mean(latent_vals),
    stringsAsFactors = FALSE
  )
}

benchmark_cache_sweep <- function(baseline_path,
                                  n_runs = 3L,
                                  out_raw = "dev-scripts/benchmarks/coupled_diagonalization_cache_sweep_results.csv",
                                  out_summary = "dev-scripts/benchmarks/coupled_diagonalization_cache_sweep_summary.csv") {
  baseline_method <- load_baseline_method(baseline_path)
  current_method <- manifoldalign:::coupled_diagonalization.hyperdesign

  mu_grid <- c(0.5, 1.0, 2.0, 4.0, 8.0)
  base_params <- list(
    ncomp = 5,
    ncomp_per_domain = 14,
    knn = 10,
    max_iter = 20,
    step_size = 0.3,
    tol = 1e-5,
    preproc = NULL,
    verbose = FALSE
  )

  raw <- do.call(rbind, lapply(seq_len(n_runs), function(run_id) {
    data_obj <- generate_hyperdesign(seed = 77L + run_id)
    hd <- data_obj$hd
    latent <- data_obj$latent

    run_rows <- do.call(rbind, lapply(
      c("baseline_uncached", "current_uncached", "current_cached"),
      function(mode) {
        out <- run_sweep_method(
          mode = mode,
          hd = hd,
          latent = latent,
          mu_grid = mu_grid,
          base_params = base_params,
          baseline_method = baseline_method,
          current_method = current_method,
          seed_base = 1000 + 100 * run_id
        )
        out$run <- run_id
        out
      }
    ))
    run_rows
  }))

  methods <- unique(raw$method)
  summary <- do.call(rbind, lapply(methods, function(mth) {
    sub <- raw[raw$method == mth, , drop = FALSE]
    data.frame(
      method = mth,
      runs = nrow(sub),
      mean_total_elapsed_sec = mean(sub$total_elapsed_sec),
      mean_final_cost = mean(sub$mean_final_cost),
      mean_pair_cosine = mean(sub$mean_pair_cosine),
      mean_latent_recovery = mean(sub$mean_latent_recovery),
      stringsAsFactors = FALSE
    )
  }))

  b <- summary[summary$method == "baseline_uncached", "mean_total_elapsed_sec"]
  u <- summary[summary$method == "current_uncached", "mean_total_elapsed_sec"]
  c <- summary[summary$method == "current_cached", "mean_total_elapsed_sec"]

  speed <- data.frame(
    speedup_baseline_over_current_uncached = b / u,
    speedup_baseline_over_current_cached = b / c,
    speedup_current_uncached_over_cached = u / c,
    stringsAsFactors = FALSE
  )

  utils::write.csv(raw, out_raw, row.names = FALSE)
  utils::write.csv(summary, out_summary, row.names = FALSE)

  cat("Coupled Diagonalization Cache Sweep Benchmark\n")
  cat("==========================================\n")
  print(summary, row.names = FALSE)
  cat("\nSpeedups:\n")
  print(speed, row.names = FALSE)
  cat(sprintf("\nWrote raw results: %s\n", out_raw))
  cat(sprintf("Wrote summary: %s\n", out_summary))

  invisible(list(raw = raw, summary = summary, speedups = speed))
}

baseline_file <- Sys.getenv("CD_BASELINE_FILE", unset = "")
if (!nzchar(baseline_file)) {
  stop(
    "Set CD_BASELINE_FILE to a saved pre-refactor coupled_diagonalization file ",
    "(e.g., /tmp/coupled_diagonalization_before_*.R).",
    call. = FALSE
  )
}

benchmark_cache_sweep(baseline_file)
