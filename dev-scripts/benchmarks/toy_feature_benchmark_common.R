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

scenario_view_names <- function(scenario) {
  view_names <- scenario$view_names %||% grep("^view[0-9]+$", names(scenario), value = TRUE)
  if (!length(view_names)) {
    stop("Scenario does not expose any view matrices.", call. = FALSE)
  }
  view_names
}

scenario_primary_view <- function(scenario) {
  scenario_view_names(scenario)[[1L]]
}

scenario_domains <- function(scenario, views = scenario_view_names(scenario)) {
  out <- lapply(views, function(view_name) as.matrix(scenario[[view_name]]))
  names(out) <- views
  out
}

scenario_row_ids <- function(scenario, views = scenario_view_names(scenario)) {
  manifoldalign:::.scenario_row_ids(scenario, views = views)
}

scenario_shared_row_ids <- function(scenario, views = scenario_view_names(scenario)) {
  row_ids <- scenario_row_ids(scenario, views = views)
  shared <- Reduce(intersect, row_ids)
  if (!length(shared)) {
    stop("Scenario does not contain any row IDs shared across all selected views.", call. = FALSE)
  }
  shared
}

scenario_correspondences <- function(scenario, views = scenario_view_names(scenario)) {
  row_ids <- scenario_row_ids(scenario, views = views)
  pairs <- utils::combn(seq_along(views), 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    ids_i <- row_ids[[pair[1L]]]
    ids_j <- row_ids[[pair[2L]]]
    shared <- intersect(ids_i, ids_j)
    if (!length(shared)) {
      return(NULL)
    }
    idx_i <- match(shared, ids_i)
    idx_j <- match(shared, ids_j)
    keep <- !is.na(idx_i) & !is.na(idx_j)
    data.frame(
      domain_i = pair[1L],
      index_i = idx_i[keep],
      domain_j = pair[2L],
      index_j = idx_j[keep],
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      domain_i = integer(0),
      index_i = integer(0),
      domain_j = integer(0),
      index_j = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

scenario_anchor_correspondence <- function(scenario, views = scenario_view_names(scenario)) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix is required to build anchor correspondence matrices.", call. = FALSE)
  }

  row_ids <- scenario_row_ids(scenario, views = views)
  shared <- scenario_shared_row_ids(scenario, views = views)
  mats <- lapply(row_ids, function(ids) {
    idx <- match(shared, ids)
    keep <- !is.na(idx)
    Matrix::sparseMatrix(
      i = idx[keep],
      j = seq_len(sum(keep)),
      x = 1,
      dims = c(length(ids), length(shared))
    )
  })
  names(mats) <- views
  mats
}

make_scenario_hyperdesign <- function(scenario) {
  manifoldalign:::synthetic_alignment_to_hyperdesign(
    scenario,
    views = scenario_view_names(scenario),
    label_col = "label",
    id_col = "sample_id"
  )
}

make_scenario_hyperdesign_pair <- function(scenario,
                                           views = scenario_view_names(scenario)[1:2]) {
  manifoldalign:::synthetic_alignment_to_hyperdesign(
    scenario,
    views = views,
    label_col = "label",
    id_col = "sample_id"
  )
}

make_scenario_hyperdesign_anchored <- function(scenario,
                                               views = scenario_view_names(scenario)[1:2],
                                               anchor_frac = 0.5,
                                               seed = 42L) {
  row_ids <- manifoldalign:::.scenario_row_ids(scenario, views = views)
  shared <- Reduce(intersect, row_ids)

  set.seed(seed)
  n_anchors <- max(2L, round(length(shared) * anchor_frac))
  anchor_ids <- sort(sample(shared, n_anchors))

  domains <- lapply(views, function(vn) {
    X <- as.matrix(scenario[[vn]])
    ids <- row_ids[[vn]]
    label_obj <- if (!is.null(scenario$observed_labels_by_view)) {
      scenario$observed_labels_by_view[[vn]]
    } else {
      scenario$true_labels
    }
    design_df <- data.frame(
      sample_id = ids,
      label = as.character(label_obj),
      anchor_id = ifelse(ids %in% anchor_ids, ids, NA_integer_),
      stringsAsFactors = FALSE
    )
    multidesign::multidesign(X, design_df)
  })
  names(domains) <- views
  multidesign::hyperdesign(domains)
}

extract_first_block_scores <- function(fit, scenario) {
  primary_view <- scenario_primary_view(scenario)
  if (!is.null(fit$s)) {
    n <- nrow(scenario[[primary_view]])
    return(as.matrix(fit$s[seq_len(n), , drop = FALSE]))
  }

  if (!is.null(fit$A_est)) {
    emb <- t(as.matrix(fit$A_est))
    if (nrow(emb) != nrow(scenario[[primary_view]])) {
      stop("Consensus embedding row count does not match scenario rows.", call. = FALSE)
    }
    return(emb)
  }

  stop("Fit does not expose `s` or `A_est` for benchmarking.", call. = FALSE)
}

score_toy_feature_fit <- function(fit, scenario, method) {
  emb <- extract_first_block_scores(fit, scenario)
  primary_view <- scenario_primary_view(scenario)
  latent <- if (!is.null(scenario$latent_by_view)) scenario$latent_by_view[[primary_view]] else scenario$latent
  latent_mat <- if (is.matrix(latent)) latent else matrix(as.numeric(latent), ncol = 1L)

  corr_mat <- suppressWarnings(stats::cor(latent_mat, emb, use = "pairwise.complete.obs"))
  max_abs_latent_cor <- if (length(corr_mat)) max(abs(corr_mat), na.rm = TRUE) else NA_real_
  if (!is.finite(max_abs_latent_cor)) {
    max_abs_latent_cor <- NA_real_
  }

  labels_obj <- if (!is.null(scenario$true_labels_by_view)) scenario$true_labels_by_view[[primary_view]] else scenario$true_labels
  labels <- as.character(labels_obj)
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
      domains <- scenario_domains(scenario)
      corr <- scenario_correspondences(scenario, views = names(domains))
      base_rank <- min(16L, min(vapply(domains, nrow, integer(1))))
      control_args <- list(
        enabled = TRUE,
        backend = backend,
        rank_per_domain = base_rank,
        knn = min(10L, min(vapply(domains, nrow, integer(1))) - 1L),
        eigen_k_max = 32L,
        scale_selection = scale_selection,
        supervision_weight = supervision_weight,
        max_basis_per_scale = max_basis_per_scale,
        tune = tuned,
        verbose = FALSE
      )

      if (identical(backend, "randomized_dwt")) {
        control_args$randomized_sketch_size <- min(16L, min(vapply(domains, nrow, integer(1))))
        control_args$randomized_oversample <- 4L
        control_args$randomized_power_iters <- 1L
      }

      if (isTRUE(tuned)) {
        control_args$candidate_cross_edge_weight <- c(0.5, 1, 2)
        control_args$candidate_rank_per_domain <- sort(unique(pmax(4L, c(base_rank %/% 2L, base_rank))))
        control_args$candidate_max_levels <- c(3L, 4L, 5L)
      }

      multiscale_manifold_align(
        domains,
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
        domains <- scenario_domains(scenario)
        corr <- scenario_correspondences(scenario, views = names(domains))
        base_n <- min(vapply(domains, nrow, integer(1)))
        spectral_mnn_align(
          domains,
          correspondences = corr,
          preproc = multivarious::center(),
          ncomp = 2,
          control = spectral_mnn_align_control(
            knn = min(10L, base_n - 1L),
            q_descriptors = min(10L, base_n - 1L),
            max_anchors = base_n,
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
        domains <- scenario_domains(scenario)
        corr <- scenario_correspondences(scenario, views = names(domains))
        ssma_align(
          domains,
          correspondences = corr,
          preproc = multivarious::center(),
          ncomp = 2,
          control = ssma_align_control(
            knn = min(10L, min(vapply(domains, nrow, integer(1))) - 1L),
            rank_per_domain = min(16L, min(vapply(domains, nrow, integer(1)))),
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
    methods$global_geo_rowid_exact <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_scenario_hyperdesign(scenario)
        domains <- scenario_domains(scenario)
        global_geo_align(
          hd,
          correspondences = scenario_correspondences(scenario, views = names(domains)),
          preproc = NULL,
          ncomp = 2,
          control = global_geo_align_control(
            knn = min(10L, min(vapply(domains, nrow, integer(1))) - 1L),
            n_landmarks_total = max(12L, min(48L, sum(vapply(domains, nrow, integer(1))) %/% 2L)),
            k_embed = 8L,
            scale_method = "none",
            ridge_eps = 1e-4,
            verbose = FALSE,
            seed = 1L
          )
        )
      },
      meta = method_meta(
        method_family = "global_geo",
        supervision_regime = "exact_correspondence",
        side_information = "row_id_correspondence",
        dimensionality_constraint = "unrestricted",
        variant = "rowid_exact"
      )
    )

    methods$lowrank_label <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_scenario_hyperdesign(scenario)
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
      requireNamespace("Matrix", quietly = TRUE)) {
    methods$coupled_diag_rowid_exact <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_scenario_hyperdesign(scenario)
        domains <- scenario_domains(scenario)
        n_rows <- min(vapply(domains, nrow, integer(1)))
        coupled_diagonalization(
          hd,
          correspondence = scenario_anchor_correspondence(scenario, views = names(domains)),
          preproc = multivarious::center(),
          ncomp = 2,
          ncomp_per_domain = min(8L, max(4L, n_rows - 1L)),
          mu_coupling = 0.5,
          knn = min(10L, n_rows - 1L),
          max_iter = 80L,
          tol = 1e-5,
          verbose = FALSE
        )
      },
      meta = method_meta(
        method_family = "coupled_diagonalization",
        supervision_regime = "exact_correspondence",
        side_information = "row_id_correspondence",
        dimensionality_constraint = "unrestricted",
        variant = "rowid_exact"
      )
    )
  }

  if (requireNamespace("multivarious", quietly = TRUE)) {
    methods$mma_consensus_unsupervised <- benchmark_method_spec(
      fit = function(scenario) {
        domains <- scenario_domains(scenario)
        mma_align_multiple(
          domains,
          preproc = multivarious::center(),
          ncomp = 2,
          knn = min(10L, min(vapply(domains, nrow, integer(1))) - 1L),
          match_to = "consensus",
          consensus_centers = "min",
          em_max_iter = 25L,
          em_tol = 1e-4,
          final_assignment = "none"
        )
      },
      meta = method_meta(
        method_family = "mma",
        supervision_regime = "unsupervised_multiset",
        side_information = "none",
        dimensionality_constraint = "unrestricted",
        variant = "consensus_unsupervised"
      )
    )
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE) &&
      requireNamespace("genpca", quietly = TRUE) &&
      requireNamespace("PRIMME", quietly = TRUE)) {
    methods$gpca_label <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_scenario_hyperdesign(scenario)
        gpca_align(
          hd,
          label,
          preproc = multivarious::center(),
          ncomp = 2,
          u = 0.5,
          lambda = 0.1,
          control = gpca_align_control(
            knn = min(10L, min(vapply(scenario_domains(scenario), nrow, integer(1))) - 1L),
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
        hd <- make_scenario_hyperdesign(scenario)
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
        hd <- make_scenario_hyperdesign(scenario)
        domains <- scenario_domains(scenario)
        sigma_hat <- mean(vapply(
          domains,
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
        hd <- make_scenario_hyperdesign(scenario)
        domains <- scenario_domains(scenario)
        sigma_hat <- mean(vapply(
          domains,
          function(X) choose_sigma(X),
          numeric(1)
        ))
        n <- min(vapply(domains, nrow, integer(1)))
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

  methods$multiset_uot_side_geometry <- benchmark_method_spec(
    fit = function(scenario) {
      if (is.null(scenario$latent_by_view)) {
        stop("multiset_uot benchmark requires per-view latent geometry.", call. = FALSE)
      }

      views <- scenario_view_names(scenario)
      latent_views <- lapply(views, function(view) {
        latent <- scenario$latent_by_view[[view]]
        as.matrix(latent)
      })
      names(latent_views) <- views

      template_view <- views[[1L]]
      template_Y <- latent_views[[template_view]]
      datasets <- lapply(views, function(view) {
        X <- latent_views[[view]]
        list(
          X = X + matrix(stats::rnorm(length(X), sd = 0.02), nrow = nrow(X), ncol = ncol(X)),
          alpha = rep(1, nrow(X))
        )
      })

      fit_uot <- multiset_uot_align(
        datasets = datasets,
        template = list(Y = template_Y, beta = rep(1, nrow(template_Y))),
        epsilon = 0.5,
        rho1 = 1,
        rho2 = 1,
        lambda_anat = 1,
        lambda_feat = 0,
        neighbor_mode = "knn",
        k_neighbors = min(10L, nrow(template_Y) - 1L),
        max_outer = 3,
        max_inner = 200,
        tol_inner = 1e-6,
        tol_outer = 1e-4,
        parallel = FALSE,
        verbose = FALSE
      )

      mapped <- multiset_uot_map(
        fit_uot,
        datasets = datasets,
        signals = lapply(datasets, function(d) t(d$X))
      )
      scores <- do.call(rbind, lapply(mapped, function(sig) t(as.matrix(sig))))

      list(
        s = scores,
        fit = fit_uot,
        mapped = mapped,
        template = fit_uot$template
      )
    },
    meta = method_meta(
      method_family = "multiset_uot",
      supervision_regime = "transport_side_geometry",
      side_information = "shared_geometry",
      dimensionality_constraint = "shared_geometry_dim",
      variant = "latent_geometry_template"
    )
  )

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE)) {
    methods$genproc_rowid_exact <- benchmark_method_spec(
      fit = function(scenario) {
        dims <- vapply(scenario_domains(scenario), ncol, integer(1))
        if (length(unique(dims)) != 1L) {
          stop(
            "generalized_procrustes requires equal feature dimensions across domains",
            call. = FALSE
          )
        }
        hd <- make_scenario_hyperdesign(scenario)
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

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE)) {
    methods$cone_align_unsupervised <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_scenario_hyperdesign_pair(scenario)
        domains <- scenario_domains(scenario, views = scenario_view_names(scenario)[1:2])
        base_n <- min(vapply(domains, nrow, integer(1)))
        cone_align(
          hd,
          preproc = multivarious::center(),
          ncomp = min(10L, base_n - 1L),
          sigma = 0.73,
          lambda = 0.1,
          use_laplacian = TRUE,
          solver = "linear",
          knn = min(10L, base_n - 1L),
          max_iter = 30L,
          tol = 0.01
        )
      },
      meta = method_meta(
        method_family = "cone_align",
        supervision_regime = "unsupervised_pairwise",
        side_information = "none",
        dimensionality_constraint = "unrestricted",
        variant = "unsupervised"
      )
    )

    methods$grasp_unsupervised <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_scenario_hyperdesign_pair(scenario)
        domains <- scenario_domains(scenario, views = scenario_view_names(scenario)[1:2])
        base_n <- min(vapply(domains, nrow, integer(1)))
        grasp(
          hd,
          preproc = multivarious::center(),
          ncomp = min(20L, base_n - 1L),
          q_descriptors = min(50L, base_n - 1L),
          sigma = 0.73,
          lambda = 0.1,
          use_laplacian = TRUE,
          solver = "linear"
        )
      },
      meta = method_meta(
        method_family = "grasp",
        supervision_regime = "unsupervised_pairwise",
        side_information = "none",
        dimensionality_constraint = "unrestricted",
        variant = "unsupervised"
      )
    )
  }

  if (requireNamespace("multidesign", quietly = TRUE) &&
      requireNamespace("multivarious", quietly = TRUE)) {
    methods$parrot_partial_anchor <- benchmark_method_spec(
      fit = function(scenario) {
        hd <- make_scenario_hyperdesign_anchored(scenario, anchor_frac = 0.5)
        domains <- scenario_domains(scenario, views = scenario_view_names(scenario)[1:2])
        base_n <- min(vapply(domains, nrow, integer(1)))
        parrot(
          hd,
          anchors = anchor_id,
          preproc = multivarious::center(),
          ncomp = min(10L, base_n - 1L),
          sigma = 0.15,
          lambda = 0.1,
          tau = 0.05,
          solver = "sinkhorn",
          max_iter = 100L,
          tol = 1e-6,
          use_cpp = FALSE,
          verbose = FALSE
        )
      },
      meta = method_meta(
        method_family = "parrot",
        supervision_regime = "partial_anchor",
        side_information = "anchor_subset",
        dimensionality_constraint = "unrestricted",
        variant = "partial_anchor"
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
