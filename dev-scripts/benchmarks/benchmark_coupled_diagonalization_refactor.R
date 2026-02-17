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

generate_hyperdesign <- function(n,
                                 m,
                                 latent_dim = 4L,
                                 noise_sd = 0.08,
                                 seed = 42L) {
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

run_trial <- function(method_fn, hd, latent, params, seed = 1L, return_fit = FALSE) {
  set.seed(seed)
  timing <- system.time({
    fit <- do.call(method_fn, c(list(data = hd), params))
  })
  n_per_domain <- vapply(hd, function(d) nrow(d$x), integer(1))
  first_block <- fit$s[seq_len(n_per_domain[1]), , drop = FALSE]
  metrics <- data.frame(
    elapsed_sec = as.numeric(timing[["elapsed"]]),
    final_cost = as.numeric(fit$final_cost),
    converged = as.logical(fit$converged),
    iterations = as.integer(fit$iterations),
    mean_pair_cosine = mean_pair_alignment(fit$s, n_per_domain),
    latent_recovery = latent_recovery(first_block, latent),
    stringsAsFactors = FALSE
  )
  if (return_fit) {
    return(list(metrics = metrics, fit = fit))
  }
  metrics
}

summarize_by_method <- function(df) {
  methods <- unique(df$method)
  do.call(rbind, lapply(methods, function(mth) {
    sub <- df[df$method == mth, , drop = FALSE]
    data.frame(
      method = mth,
      runs = nrow(sub),
      mean_elapsed_sec = mean(sub$elapsed_sec),
      mean_final_cost = mean(sub$final_cost),
      mean_iterations = mean(sub$iterations),
      converged_rate = mean(sub$converged),
      mean_pair_cosine = mean(sub$mean_pair_cosine),
      mean_latent_recovery = mean(sub$latent_recovery, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

run_scenario <- function(scenario, baseline_method, current_method, n_runs, use_cache = TRUE) {
  bench_data <- generate_hyperdesign(n = scenario$n, m = scenario$m, seed = scenario$seed)
  hd <- bench_data$hd
  latent <- bench_data$latent

  params <- list(
    ncomp = scenario$ncomp,
    ncomp_per_domain = scenario$ncomp_per_domain,
    mu_coupling = scenario$mu_coupling,
    knn = scenario$knn,
    max_iter = scenario$max_iter,
    step_size = scenario$step_size,
    tol = scenario$tol,
    preproc = NULL,
    verbose = FALSE
  )

  baseline_runs <- do.call(rbind, lapply(seq_len(n_runs), function(i) {
    out <- run_trial(baseline_method, hd, latent, params, seed = scenario$seed + i)
    out$scenario <- scenario$name
    out$method <- "baseline_pre_refactor"
    out$run <- i
    out
  }))

  current_runs <- do.call(rbind, lapply(seq_len(n_runs), function(i) {
    out <- run_trial(current_method, hd, latent, params, seed = scenario$seed + i)
    out$scenario <- scenario$name
    out$method <- "current_refactor_uncached"
    out$run <- i
    out
  }))
  
  cached_runs <- NULL
  if (isTRUE(use_cache)) {
    cache <- NULL
    cached_runs <- do.call(rbind, lapply(seq_len(n_runs), function(i) {
      params_i <- params
      if (!is.null(cache)) {
        params_i$spectral_cache <- cache
      }
      out_i <- run_trial(
        method_fn = current_method,
        hd = hd,
        latent = latent,
        params = params_i,
        seed = scenario$seed + i,
        return_fit = TRUE
      )
      if (is.null(cache)) {
        cache <- out_i$fit$diagnostics$spectral_cache
      }
      out <- out_i$metrics
      out$scenario <- scenario$name
      out$method <- "current_refactor_cached"
      out$run <- i
      out
    }))
  }

  raw <- rbind(
    baseline_runs,
    current_runs,
    if (!is.null(cached_runs)) cached_runs
  )
  summary <- summarize_by_method(raw)
  summary$scenario <- scenario$name

  b <- summary[summary$method == "baseline_pre_refactor", , drop = FALSE]
  compare_methods <- setdiff(summary$method, "baseline_pre_refactor")
  deltas <- do.call(rbind, lapply(compare_methods, function(mth) {
    c <- summary[summary$method == mth, , drop = FALSE]
    data.frame(
      scenario = scenario$name,
      compared_method = mth,
      speedup_baseline_over_method = b$mean_elapsed_sec / c$mean_elapsed_sec,
      final_cost_delta_method_minus_baseline = c$mean_final_cost - b$mean_final_cost,
      pair_cosine_delta_method_minus_baseline = c$mean_pair_cosine - b$mean_pair_cosine,
      latent_recovery_delta_method_minus_baseline = c$mean_latent_recovery - b$mean_latent_recovery,
      stringsAsFactors = FALSE
    )
  }))

  list(raw = raw, summary = summary, deltas = deltas)
}

benchmark_refactor <- function(baseline_path,
                               n_runs = 3L,
                               use_cache_for_current = TRUE,
                               out_csv = "dev-scripts/benchmarks/coupled_diagonalization_refactor_results.csv",
                               out_summary = "dev-scripts/benchmarks/coupled_diagonalization_refactor_summary.csv",
                               out_deltas = "dev-scripts/benchmarks/coupled_diagonalization_refactor_deltas.csv") {
  baseline_method <- load_baseline_method(baseline_path)
  current_method <- manifoldalign:::coupled_diagonalization.hyperdesign

  scenarios <- list(
    list(
      name = "default_quality",
      n = 220L,
      m = 3L,
      ncomp = 4L,
      ncomp_per_domain = 10L,
      mu_coupling = 2.0,
      knn = 10L,
      max_iter = 70L,
      step_size = 0.3,
      tol = 1e-5,
      seed = 42L
    ),
    list(
      name = "high_q_speed",
      n = 1000L,
      m = 3L,
      ncomp = 5L,
      ncomp_per_domain = 14L,
      mu_coupling = 2.0,
      knn = 10L,
      max_iter = 20L,
      step_size = 0.3,
      tol = 1e-5,
      seed = 77L
    )
  )

  pieces <- lapply(scenarios, run_scenario,
                   baseline_method = baseline_method,
                   current_method = current_method,
                   n_runs = n_runs,
                   use_cache = use_cache_for_current)

  raw <- do.call(rbind, lapply(pieces, `[[`, "raw"))
  summary <- do.call(rbind, lapply(pieces, `[[`, "summary"))
  deltas <- do.call(rbind, lapply(pieces, `[[`, "deltas"))

  utils::write.csv(raw, out_csv, row.names = FALSE)
  utils::write.csv(summary, out_summary, row.names = FALSE)
  utils::write.csv(deltas, out_deltas, row.names = FALSE)

  cat("Coupled Diagonalization Refactor Benchmark\n")
  cat("========================================\n")
  cat(sprintf("Current cache mode enabled: %s\n\n", if (isTRUE(use_cache_for_current)) "yes" else "no"))
  print(summary[order(summary$scenario, summary$method), ], row.names = FALSE)
  cat("\nScenario deltas (method - baseline where noted):\n")
  print(deltas, row.names = FALSE)
  cat(sprintf("\nWrote raw results: %s\n", out_csv))
  cat(sprintf("Wrote summary: %s\n", out_summary))
  cat(sprintf("Wrote deltas: %s\n", out_deltas))

  invisible(list(raw = raw, summary = summary, deltas = deltas))
}

baseline_file <- Sys.getenv("CD_BASELINE_FILE", unset = "")
if (!nzchar(baseline_file)) {
  stop(
    "Set CD_BASELINE_FILE to a saved pre-refactor coupled_diagonalization file ",
    "(e.g., /tmp/coupled_diagonalization_before_*.R).",
    call. = FALSE
  )
}

benchmark_refactor(baseline_file, use_cache_for_current = TRUE)
