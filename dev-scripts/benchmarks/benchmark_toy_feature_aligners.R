#!/usr/bin/env Rscript

# benchmark_toy_feature_aligners.R
# ------------------------------------------------------------
# Compare feature-domain alignment methods on the shared synthetic
# benchmark registry used by the toy benchmark tests.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION") &&
      grepl("^Package:\\s*manifoldalign", paste(readLines("DESCRIPTION", n = 5L), collapse = "\n"))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(manifoldalign)
})

`%||%` <- function(a, b) if (is.null(a) || identical(a, "")) b else a

write_csv_safe <- function(df, path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  utils::write.csv(df, file = path, row.names = FALSE)
  invisible(path)
}

identity_correspondences <- function(n) {
  data.frame(
    domain_i = 1L,
    index_i = seq_len(n),
    domain_j = 2L,
    index_j = seq_len(n),
    stringsAsFactors = FALSE
  )
}

make_pair_hyperdesign <- function(scenario) {
  manifoldalign:::synthetic_alignment_to_hyperdesign(
    scenario,
    views = c("view1", "view2"),
    label = scenario$true_labels,
    label_col = "label",
    id_col = "sample_id"
  )
}

extract_first_block_scores <- function(fit, scenario) {
  if (is.null(fit$s)) {
    stop("Fit does not expose a score matrix `s`.", call. = FALSE)
  }
  n <- nrow(scenario$view1)
  as.matrix(fit$s[seq_len(n), , drop = FALSE])
}

score_fit <- function(fit, scenario, method) {
  emb <- extract_first_block_scores(fit, scenario)
  latent <- scenario$latent
  latent_mat <- if (is.matrix(latent)) latent else matrix(as.numeric(latent), ncol = 1L)

  corr_mat <- suppressWarnings(stats::cor(latent_mat, emb, use = "pairwise.complete.obs"))
  max_abs_latent_cor <- if (length(corr_mat)) max(abs(corr_mat), na.rm = TRUE) else NA_real_
  if (!is.finite(max_abs_latent_cor)) {
    max_abs_latent_cor <- NA_real_
  }

  labels <- as.character(scenario$true_labels)
  centroid_mat <- vapply(unique(labels), function(lbl) {
    colMeans(emb[labels == lbl, , drop = FALSE])
  }, numeric(ncol(emb)))
  if (is.null(dim(centroid_mat))) {
    centroid_mat <- matrix(centroid_mat, ncol = 1L)
    colnames(centroid_mat) <- unique(labels)
  }
  dists <- sapply(seq_len(ncol(centroid_mat)), function(j) {
    rowSums((emb - matrix(centroid_mat[, j], nrow(emb), ncol(emb), byrow = TRUE))^2)
  })
  pred <- colnames(centroid_mat)[max.col(-dists, ties.method = "first")]
  centroid_acc <- mean(pred == labels)

  c(
    max_abs_latent_cor = max_abs_latent_cor,
    centroid_acc = centroid_acc,
    n_rows = nrow(emb),
    n_cols = ncol(emb)
  )
}

build_methods <- function() {
  methods <- list()

  if (requireNamespace("multivarious", quietly = TRUE)) {
    methods$spectral_mnn <- function(scenario) {
      corr <- identity_correspondences(nrow(scenario$view1))
      spectral_mnn_align(
        list(view1 = scenario$view1, view2 = scenario$view2),
        correspondences = corr,
        preproc = multivarious::center(),
        ncomp = 2,
        control = spectral_mnn_align_control(
          knn = min(10L, nrow(scenario$view1) - 1L),
          q_descriptors = min(10L, nrow(scenario$view1) - 1L),
          max_anchors = nrow(scenario$view1),
          refine_rounds = 0L,
          verbose = FALSE
        )
      )
    }

    methods$ssma <- function(scenario) {
      corr <- identity_correspondences(nrow(scenario$view1))
      ssma_align(
        list(view1 = scenario$view1, view2 = scenario$view2),
        correspondences = corr,
        preproc = multivarious::center(),
        ncomp = 2,
        control = ssma_align_control(
          knn = min(10L, nrow(scenario$view1) - 1L),
          rank_per_domain = min(16L, nrow(scenario$view1)),
          verbose = FALSE
        )
      )
    }
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE)) {
    methods$lowrank <- function(scenario) {
      hd <- make_pair_hyperdesign(scenario)
      lv <- levels(scenario$true_labels)
      S <- diag(length(lv))
      dimnames(S) <- list(lv, lv)
      lowrank_align(
        hd,
        label,
        simfun = createSimFun(S),
        preproc = multivarious::center(),
        ncomp = 2,
        lambda = 0.01,
        solver = "operator"
      )
    }
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE) &&
      requireNamespace("kernlab", quietly = TRUE)) {
    methods$kema <- function(scenario) {
      hd <- make_pair_hyperdesign(scenario)
      n <- nrow(scenario$view1)
      kema(
        hd,
        label,
        kernel = kernlab::rbfdot(sigma = 1.0),
        preproc = multivarious::center(),
        ncomp = 2,
        solver = "exact",
        sigma = 1.0,
        sample_frac = if (n > 40L) 0.4 else 1.0
      )
    }
  }

  methods
}

main <- function() {
  here <- dirname(normalizePath(sys.frame(1)$ofile %||% "dev-scripts/benchmarks/benchmark_toy_feature_aligners.R", mustWork = FALSE))
  out_results <- file.path(here, "toy_feature_aligners_results.csv")
  out_summary <- file.path(here, "toy_feature_aligners_summary.csv")

  profile <- match.arg(Sys.getenv("MANIFOLDALIGN_TOY_PROFILE", "fast"), c("full", "fast"))
  seeds <- as.integer(strsplit(Sys.getenv("MANIFOLDALIGN_TOY_SEEDS", "1,2,3"), ",")[[1]])
  seeds <- seeds[is.finite(seeds)]
  scenarios <- manifoldalign:::synthetic_alignment_scenarios(profile = profile)$scenario
  methods <- build_methods()

  if (!length(methods)) {
    stop("No benchmark methods available. Install multidesign/multivarious/kernlab as needed.", call. = FALSE)
  }

  cat("Toy feature aligner benchmark\n")
  cat("----------------------------\n")
  cat("profile: ", profile, "\n", sep = "")
  cat("seeds: ", paste(seeds, collapse = ","), "\n", sep = "")
  cat("scenarios: ", paste(scenarios, collapse = ", "), "\n", sep = "")
  cat("methods: ", paste(names(methods), collapse = ", "), "\n\n", sep = "")

  results <- manifoldalign:::run_synthetic_alignment_benchmark(
    methods = methods,
    scenarios = scenarios,
    profile = profile,
    seeds = seeds,
    score_fn = score_fit,
    verbose = TRUE
  )

  summary_df <- manifoldalign:::summarize_synthetic_alignment_benchmark(
    results[is.na(results$error), , drop = FALSE],
    metric_cols = c("runtime_sec", "max_abs_latent_cor", "centroid_acc")
  )

  write_csv_safe(results, out_results)
  write_csv_safe(summary_df, out_summary)

  cat("Wrote:\n")
  cat("  - ", out_results, "\n", sep = "")
  cat("  - ", out_summary, "\n\n", sep = "")

  print(summary_df)

  invisible(list(results = results, summary = summary_df))
}

if (sys.nframe() == 0) {
  main()
}
