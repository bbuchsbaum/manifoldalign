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
  if (!is.null(fit$s)) {
    n <- nrow(scenario$view1)
    return(as.matrix(fit$s[seq_len(n), , drop = FALSE]))
  }

  if (!is.null(fit$A_est)) {
    emb <- t(as.matrix(fit$A_est))
    if (nrow(emb) != nrow(scenario$view1)) {
      stop("Consensus embedding row count does not match scenario rows.", call. = FALSE)
    }
    return(emb)
  }

  stop("Fit does not expose `s` or `A_est` for benchmarking.", call. = FALSE)
}

score_toy_feature_fit <- function(fit, scenario, method) {
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

benchmark_method_spec <- function(fit, meta = list()) {
  list(fit = fit, meta = meta)
}

method_meta <- function(method_family,
                        supervision_regime,
                        side_information,
                        dimensionality_constraint = "unrestricted",
                        variant_family = method_family,
                        variant = "default",
                        backend = NA_character_,
                        tuned = FALSE,
                        kernel = NA_character_,
                        landmarking = NA_character_,
                        scale_selection = NA_character_) {
  list(
    method_family = method_family,
    supervision_regime = supervision_regime,
    side_information = side_information,
    dimensionality_constraint = dimensionality_constraint,
    variant_family = variant_family,
    variant = variant,
    backend = backend,
    tuned = isTRUE(tuned),
    kernel = kernel,
    landmarking = landmarking,
    scale_selection = scale_selection
  )
}

make_multiscale_method <- function(name,
                                   backend = c("spectral", "randomized_dwt"),
                                   scale_selection = c("energy", "supervised", "hybrid"),
                                   tuned = FALSE,
                                   supervision_weight = 0.7,
                                   max_basis_per_scale = 12L) {
  backend <- match.arg(backend)
  scale_selection <- match.arg(scale_selection)

  benchmark_method_spec(
    fit = function(scenario) {
      corr <- identity_correspondences(nrow(scenario$view1))
      base_rank <- min(16L, nrow(scenario$view1))
      control_args <- list(
        enabled = TRUE,
        backend = backend,
        rank_per_domain = base_rank,
        knn = min(10L, nrow(scenario$view1) - 1L),
        eigen_k_max = 32L,
        scale_selection = scale_selection,
        supervision_weight = supervision_weight,
        max_basis_per_scale = max_basis_per_scale,
        tune = tuned,
        verbose = FALSE
      )

      if (identical(backend, "randomized_dwt")) {
        control_args$randomized_sketch_size <- min(16L, nrow(scenario$view1))
        control_args$randomized_oversample <- 4L
        control_args$randomized_power_iters <- 1L
      }

      if (isTRUE(tuned)) {
        control_args$candidate_cross_edge_weight <- c(0.5, 1, 2)
        control_args$candidate_rank_per_domain <- sort(unique(pmax(4L, c(base_rank %/% 2L, base_rank))))
        control_args$candidate_max_levels <- c(3L, 4L, 5L)
      }

      multiscale_manifold_align(
        list(view1 = scenario$view1, view2 = scenario$view2),
        correspondences = corr,
        preproc = multivarious::center(),
        ncomp = 2,
        control = do.call(multiscale_manifold_align_control, control_args),
        verbose = FALSE
      )
    },
    meta = method_meta(
      method_family = "multiscale",
      supervision_regime = "exact_correspondence",
      side_information = "row_id_correspondence",
      dimensionality_constraint = "unrestricted",
      variant_family = "multiscale",
      variant = name,
      backend = backend,
      tuned = tuned,
      scale_selection = scale_selection
    )
  )
}

build_toy_feature_methods <- function(include_multiscale_ablations = TRUE) {
  methods <- list()

  if (requireNamespace("multivarious", quietly = TRUE)) {
    methods$spectral_mnn_exact <- benchmark_method_spec(
      fit = function(scenario) {
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
      },
      meta = method_meta(
        method_family = "spectral_mnn",
        supervision_regime = "exact_correspondence",
        side_information = "row_id_correspondence",
        dimensionality_constraint = "unrestricted",
        variant = "exact"
      )
    )

    methods$ssma_exact <- benchmark_method_spec(
      fit = function(scenario) {
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
      },
      meta = method_meta(
        method_family = "ssma",
        supervision_regime = "exact_correspondence",
        side_information = "row_id_correspondence",
        dimensionality_constraint = "unrestricted",
        variant = "exact"
      )
    )

    if (isTRUE(include_multiscale_ablations)) {
      methods$multiscale_spectral_energy <- make_multiscale_method(
        name = "spectral_energy",
        backend = "spectral",
        scale_selection = "energy",
        tuned = FALSE,
        supervision_weight = 0
      )
      methods$multiscale_spectral_supervised <- make_multiscale_method(
        name = "spectral_supervised",
        backend = "spectral",
        scale_selection = "supervised",
        tuned = FALSE,
        supervision_weight = 1
      )
      methods$multiscale_spectral_hybrid <- make_multiscale_method(
        name = "spectral_hybrid",
        backend = "spectral",
        scale_selection = "hybrid",
        tuned = FALSE,
        supervision_weight = 0.7
      )
      methods$multiscale_spectral_hybrid_tuned <- make_multiscale_method(
        name = "spectral_hybrid_tuned",
        backend = "spectral",
        scale_selection = "hybrid",
        tuned = TRUE,
        supervision_weight = 0.7
      )
      methods$multiscale_randomized_hybrid_tuned <- make_multiscale_method(
        name = "randomized_hybrid_tuned",
        backend = "randomized_dwt",
        scale_selection = "hybrid",
        tuned = TRUE,
        supervision_weight = 0.7
      )
    }
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE)) {
    methods$lowrank_label <- benchmark_method_spec(
      fit = function(scenario) {
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
      },
      meta = method_meta(
        method_family = "lowrank",
        supervision_regime = "label_supervision",
        side_information = "labels",
        dimensionality_constraint = "unrestricted",
        variant = "label"
      )
    )
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE) &&
      requireNamespace("genpca", quietly = TRUE) &&
      requireNamespace("PRIMME", quietly = TRUE)) {
    methods$gpca_label <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_pair_hyperdesign(scenario)
        gpca_align(
          hd,
          label,
          preproc = multivarious::center(),
          ncomp = 2,
          u = 0.5,
          lambda = 0.1,
          control = gpca_align_control(
            knn = min(10L, nrow(scenario$view1) - 1L),
            verbose = FALSE
          )
        )
      },
      meta = method_meta(
        method_family = "gpca",
        supervision_regime = "label_supervision",
        side_information = "labels",
        dimensionality_constraint = "unrestricted",
        variant = "label"
      )
    )
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE) &&
      requireNamespace("kernlab", quietly = TRUE)) {
    methods$kema_linear_label <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_pair_hyperdesign(scenario)
        kema(
          hd,
          label,
          kernel = kernlab::vanilladot(),
          preproc = multivarious::center(),
          ncomp = 2,
          solver = "exact",
          sample_frac = 1
        )
      },
      meta = method_meta(
        method_family = "kema",
        supervision_regime = "label_supervision",
        side_information = "labels",
        dimensionality_constraint = "unrestricted",
        variant = "linear_label",
        kernel = "linear",
        landmarking = "full"
      )
    )

    methods$kema_rbf_full_label <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_pair_hyperdesign(scenario)
        sigma_hat <- mean(vapply(
          list(scenario$view1, scenario$view2),
          function(X) choose_sigma(X),
          numeric(1)
        ))
        kema(
          hd,
          label,
          kernel = kernlab::rbfdot(sigma = sigma_hat),
          preproc = multivarious::center(),
          ncomp = 2,
          solver = "exact",
          sigma = sigma_hat,
          sample_frac = 1
        )
      },
      meta = method_meta(
        method_family = "kema",
        supervision_regime = "label_supervision",
        side_information = "labels",
        dimensionality_constraint = "unrestricted",
        variant = "rbf_full_label",
        kernel = "rbf",
        landmarking = "full"
      )
    )

    methods$kema_rbf_reduced_label <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_pair_hyperdesign(scenario)
        sigma_hat <- mean(vapply(
          list(scenario$view1, scenario$view2),
          function(X) choose_sigma(X),
          numeric(1)
        ))
        n <- nrow(scenario$view1)
        kema(
          hd,
          label,
          kernel = kernlab::rbfdot(sigma = sigma_hat),
          preproc = multivarious::center(),
          ncomp = 2,
          solver = "exact",
          sigma = sigma_hat,
          sample_frac = if (n > 40L) 0.5 else 0.75
        )
      },
      meta = method_meta(
        method_family = "kema",
        supervision_regime = "label_supervision",
        side_information = "labels",
        dimensionality_constraint = "unrestricted",
        variant = "rbf_reduced_label",
        kernel = "rbf",
        landmarking = "reduced"
      )
    )
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE)) {
    methods$genproc_rowid_exact <- benchmark_method_spec(
      fit = function(scenario) {
        if (ncol(scenario$view1) != ncol(scenario$view2)) {
          stop(
            "generalized_procrustes requires equal feature dimensions across domains",
            call. = FALSE
          )
        }
        hd <- make_pair_hyperdesign(scenario)
        generalized_procrustes(
          hd,
          sample_id,
          preproc = multivarious::center(),
          verbose = FALSE
        )
      },
      meta = method_meta(
        method_family = "generalized_procrustes",
        supervision_regime = "exact_correspondence",
        side_information = "row_id_correspondence",
        dimensionality_constraint = "equal_feature_dim",
        variant = "rowid_exact"
      )
    )
  }

  methods
}

summarize_toy_feature_regime_ranks <- function(summary_df) {
  metric_cols <- c(
    runtime_sec = "min",
    max_abs_latent_cor = "max",
    centroid_acc = "max"
  )
  needed <- c(
    "scenario", "profile", "supervision_regime", "dimensionality_constraint",
    "method", "method_family"
  )
  if (!all(needed %in% names(summary_df))) {
    stop("Summary data frame is missing required regime columns.", call. = FALSE)
  }

  strata <- unique(summary_df[c("scenario", "profile", "supervision_regime", "dimensionality_constraint")])
  rows <- vector("list", nrow(strata))

  for (i in seq_len(nrow(strata))) {
    key <- strata[i, , drop = FALSE]
    idx <- rep(TRUE, nrow(summary_df))
    for (nm in names(key)) {
      idx <- idx & summary_df[[nm]] == key[[nm]]
    }
    block <- summary_df[idx, , drop = FALSE]
    rank_block <- block

    rank_cols <- character(0)
    for (metric in names(metric_cols)) {
      metric_col <- paste0(metric, "_mean")
      if (!metric_col %in% names(block)) {
        next
      }
      vals <- block[[metric_col]]
      dir <- metric_cols[[metric]]
      rank_vals <- rep(NA_real_, length(vals))
      ok <- is.finite(vals)
      if (any(ok)) {
        rank_vals[ok] <- if (identical(dir, "min")) {
          rank(vals[ok], ties.method = "average")
        } else {
          rank(-vals[ok], ties.method = "average")
        }
      }
      rank_name <- paste0(metric, "_rank")
      rank_block[[rank_name]] <- rank_vals
      rank_cols <- c(rank_cols, rank_name)
    }

    if (length(rank_cols)) {
      rank_block$overall_rank <- rowMeans(rank_block[rank_cols], na.rm = TRUE)
    } else {
      rank_block$overall_rank <- NA_real_
    }

    rows[[i]] <- rank_block[order(rank_block$overall_rank, rank_block$method), , drop = FALSE]
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
