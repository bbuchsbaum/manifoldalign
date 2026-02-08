#' Control settings for `global_geo_align.hyperdesign()`
#'
#' @param knn Number of within-dataset nearest neighbors used to build each
#'   domain graph.
#' @param metric Distance metric for within-dataset graph edges (`"euclidean"`
#'   or `"cosine"`).
#' @param corr_edge_weight Weight assigned to correspondence edges.
#' @param n_landmarks_total Total number of landmarks across all datasets.
#' @param landmarks_per_dataset Optional fixed landmark count per dataset.
#' @param landmark_corr_frac Fraction of landmarks reserved for correspondence
#'   nodes when available.
#' @param k_embed Intermediate geometry rank used in landmark MDS. If `NULL`,
#'   defaults to `max(3 * ncomp, ncomp + 5)` inside the solver.
#' @param scale_method Dataset scale harmonization strategy.
#' @param scale_graph_reg Ridge penalty for global log-scale fitting.
#' @param scale_anchor_max Maximum matched anchors per dataset pair used during
#'   scale estimation.
#' @param max_pairs_per_label Maximum correspondences sampled per shared label
#'   when `y` is used to generate correspondences.
#' @param ridge_eps Ridge regularization base value for linear map estimation.
#' @param ridge_relative Whether to scale ridge by per-dataset feature energy.
#' @param ginv_method Method for solving `(X^T X + eps I)^{-1} C`.
#' @param eig_tol Minimum eigenvalue threshold retained in spectral steps.
#' @param embed_block_size Column block size used for landmark out-of-sample
#'   embedding.
#' @param seed Random seed for landmark/correspondence sampling.
#' @param verbose Logical toggle for progress messages.
#'
#' @return A list of class `global_geo_align_control`.
#' @export
global_geo_align_control <- function(
  knn = 15L,
  metric = c("euclidean", "cosine"),
  corr_edge_weight = 0,
  n_landmarks_total = 1000L,
  landmarks_per_dataset = NULL,
  landmark_corr_frac = 0.3,
  k_embed = NULL,
  scale_method = c("none", "pairwise_eta", "robust_median"),
  scale_graph_reg = 1e-4,
  scale_anchor_max = 128L,
  max_pairs_per_label = 100L,
  ridge_eps = 1e-4,
  ridge_relative = TRUE,
  ginv_method = c("auto", "woodbury", "direct"),
  eig_tol = 1e-8,
  embed_block_size = 2000L,
  seed = 1L,
  verbose = TRUE
) {
  metric <- match.arg(metric)
  scale_method <- match.arg(scale_method)
  ginv_method <- match.arg(ginv_method)

  if (!is.numeric(knn) || length(knn) != 1L || is.na(knn) || knn < 1) {
    stop("`knn` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(corr_edge_weight) || length(corr_edge_weight) != 1L ||
      is.na(corr_edge_weight) || corr_edge_weight < 0) {
    stop("`corr_edge_weight` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(n_landmarks_total) || length(n_landmarks_total) != 1L ||
      is.na(n_landmarks_total) || n_landmarks_total < 2) {
    stop("`n_landmarks_total` must be >= 2.", call. = FALSE)
  }
  if (!is.null(landmarks_per_dataset) &&
      (!is.numeric(landmarks_per_dataset) || length(landmarks_per_dataset) != 1L ||
       is.na(landmarks_per_dataset) || landmarks_per_dataset < 1)) {
    stop("`landmarks_per_dataset` must be NULL or a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(landmark_corr_frac) || length(landmark_corr_frac) != 1L ||
      is.na(landmark_corr_frac) || landmark_corr_frac < 0 || landmark_corr_frac > 1) {
    stop("`landmark_corr_frac` must be in [0, 1].", call. = FALSE)
  }
  if (!is.null(k_embed) &&
      (!is.numeric(k_embed) || length(k_embed) != 1L || is.na(k_embed) || k_embed < 2)) {
    stop("`k_embed` must be NULL or an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(scale_graph_reg) || length(scale_graph_reg) != 1L ||
      is.na(scale_graph_reg) || scale_graph_reg < 0) {
    stop("`scale_graph_reg` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(scale_anchor_max) || length(scale_anchor_max) != 1L ||
      is.na(scale_anchor_max) || scale_anchor_max < 2) {
    stop("`scale_anchor_max` must be >= 2.", call. = FALSE)
  }
  if (!is.numeric(max_pairs_per_label) || length(max_pairs_per_label) != 1L ||
      is.na(max_pairs_per_label) || max_pairs_per_label < 1) {
    stop("`max_pairs_per_label` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(ridge_eps) || length(ridge_eps) != 1L || is.na(ridge_eps) || ridge_eps <= 0) {
    stop("`ridge_eps` must be a positive scalar.", call. = FALSE)
  }
  if (!is.logical(ridge_relative) || length(ridge_relative) != 1L || is.na(ridge_relative)) {
    stop("`ridge_relative` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(eig_tol) || length(eig_tol) != 1L || is.na(eig_tol) || eig_tol < 0) {
    stop("`eig_tol` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(embed_block_size) || length(embed_block_size) != 1L ||
      is.na(embed_block_size) || embed_block_size < 1) {
    stop("`embed_block_size` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      knn = as.integer(knn),
      metric = metric,
      corr_edge_weight = as.numeric(corr_edge_weight),
      n_landmarks_total = as.integer(n_landmarks_total),
      landmarks_per_dataset = if (is.null(landmarks_per_dataset)) NULL else as.integer(landmarks_per_dataset),
      landmark_corr_frac = as.numeric(landmark_corr_frac),
      k_embed = if (is.null(k_embed)) NULL else as.integer(k_embed),
      scale_method = scale_method,
      scale_graph_reg = as.numeric(scale_graph_reg),
      scale_anchor_max = as.integer(scale_anchor_max),
      max_pairs_per_label = as.integer(max_pairs_per_label),
      ridge_eps = as.numeric(ridge_eps),
      ridge_relative = ridge_relative,
      ginv_method = ginv_method,
      eig_tol = as.numeric(eig_tol),
      embed_block_size = as.integer(embed_block_size),
      seed = as.integer(seed),
      verbose = verbose
    ),
    class = "global_geo_align_control"
  )
}

resolve_global_geo_align_control <- function(control) {
  defaults <- global_geo_align_control()

  if (missing(control) || is.null(control)) {
    return(defaults)
  }
  if (inherits(control, "global_geo_align_control")) {
    merged <- modifyList(defaults, control)
    return(do.call(global_geo_align_control, merged))
  }
  if (!is.list(control)) {
    stop("`control` must be NULL, a named list, or global_geo_align_control().", call. = FALSE)
  }
  if (is.null(names(control)) || any(names(control) == "")) {
    stop("`control` must be a named list.", call. = FALSE)
  }

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown control argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  merged <- modifyList(defaults, control)
  do.call(global_geo_align_control, merged)
}

#' Global-Geometry Multi-Dataset Alignment
#'
#' Aligns multiple domains by preserving global graph geometry induced by
#' within-domain kNN neighborhoods and cross-domain correspondence edges.
#'
#' @param data A hyperdesign object or list of matrices.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `multiblock_biprojector` object with subclass `global_geo_align`.
#' @export
global_geo_align <- function(data, ...) {
  UseMethod("global_geo_align")
}

#' @rdname global_geo_align
#' @method global_geo_align hyperdesign
#' @export
global_geo_align.hyperdesign <- function(
  data,
  y = NULL,
  correspondences = NULL,
  feature_correspondence = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = global_geo_align_control(),
  ...
) {
  ctrl <- resolve_global_geo_align_control(control)

  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1) {
    stop("`ncomp` must be a positive scalar.", call. = FALSE)
  }

  resolved <- resolve_hyperdesign(data)
  domains <- resolved$domains
  domain_names <- resolved$domain_names
  n_domains <- resolved$n_domains

  if (n_domains < 2L) {
    stop("global_geo_align requires at least two domains.", call. = FALSE)
  }

  if (isTRUE(ctrl$verbose)) {
    message("global_geo_align: preprocessing ", n_domains, " domains.")
  }

  pdata <- domains
  proclist <- vector("list", n_domains)
  names(proclist) <- domain_names

  is_stateful_preproc <- inherits(preproc, "prepper") || inherits(preproc, "pre_processor")
  broadcast <- is.function(preproc) || is_stateful_preproc
  if (broadcast) {
    if (is_stateful_preproc) {
      preproc <- replicate(
        n_domains,
        unserialize(serialize(preproc, connection = NULL)),
        simplify = FALSE
      )
    } else {
      preproc <- rep(list(preproc), n_domains)
    }
  } else if (is.list(preproc) && length(preproc) != n_domains) {
    stop("Length of `preproc` must match number of domains.", call. = FALSE)
  }

  for (i in seq_len(n_domains)) {
    Xi <- as.matrix(domains[[i]]$x)
    if (!all(is.finite(Xi))) {
      stop("Domain ", i, " contains non-finite values.", call. = FALSE)
    }
    if (!is.null(preproc) && !is.null(preproc[[i]])) {
      pre_i <- preproc[[i]]
      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        pre_i <- unserialize(serialize(pre_i, connection = NULL))
      }
      templ <- multivarious::prep(pre_i)
      Xi_tf <- multivarious::init_transform(templ, Xi)
      pdata[[i]]$x <- as.matrix(Xi_tf)
      proc_attr <- attr(Xi_tf, "preproc")
      proclist[[i]] <- if (is.null(proc_attr)) templ else proc_attr
    } else {
      pdata[[i]]$x <- Xi
      proclist[[i]] <- NULL
    }
  }
  names(pdata) <- domain_names

  sab_info <- NULL
  pdata_geom <- pdata
  if (!is.null(feature_correspondence)) {
    sab_info <- .apply_spatial_anchor_basis(
      pdata = pdata,
      domain_names = domain_names,
      spec = feature_correspondence,
      seed = ctrl$seed
    )
    if (!is.null(sab_info)) {
      pdata_geom <- sab_info$pdata_tilde
    }
    if (!is.null(sab_info) && isTRUE(ctrl$verbose)) {
      message(
        "global_geo_align: applied spatial anchor basis (K=",
        ncol(pdata_geom[[1]]$x),
        ")."
      )
    }
  }

  proc <- NULL
  if (!all(vapply(proclist, is.null, logical(1)))) {
    sample_idx <- compute_block_indices(pdata)
    sample_idx_list <- split(sample_idx, row(sample_idx))
    proc <- multivarious::concat_pre_processors(proclist, sample_idx_list)
  }

  domain_sizes <- vapply(pdata_geom, function(d) nrow(d$x), integer(1))
  if (any(domain_sizes < 2L)) {
    stop("Each domain must contain at least two samples.", call. = FALSE)
  }

  y_expr <- if (missing(y)) NULL else rlang::enexpr(y)

  corr_df <- .resolve_correspondences(
    pdata = pdata_geom,
    domain_names = domain_names,
    domain_sizes = domain_sizes,
    y_expr = y_expr,
    correspondences = correspondences,
    max_pairs_per_label = ctrl$max_pairs_per_label,
    seed = ctrl$seed
  )
  if (!nrow(corr_df)) {
    stop("No correspondences available. Supply `correspondences` or `y`.", call. = FALSE)
  }

  scales <- .estimate_domain_scales(
    pdata = pdata_geom,
    corr_df = corr_df,
    method = ctrl$scale_method,
    reg = ctrl$scale_graph_reg,
    anchor_max = ctrl$scale_anchor_max,
    seed = ctrl$seed
  )

  offsets <- c(0L, cumsum(domain_sizes))
  total_n <- sum(domain_sizes)

  edge_df <- .build_global_graph_edges(
    pdata = pdata_geom,
    domain_sizes = domain_sizes,
    offsets = offsets,
    corr_df = corr_df,
    scales = scales,
    knn = ctrl$knn,
    metric = ctrl$metric,
    corr_edge_weight = ctrl$corr_edge_weight
  )

  if (!nrow(edge_df)) {
    stop("No graph edges were constructed. Check `knn` and correspondences.", call. = FALSE)
  }

  graph <- igraph::graph_from_data_frame(
    edge_df,
    directed = FALSE,
    vertices = data.frame(name = as.character(seq_len(total_n)))
  )
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
                            edge.attr.comb = list(weight = "min", "ignore"))

  landmark_ids <- .select_landmarks(
    n_landmarks_total = ctrl$n_landmarks_total,
    landmarks_per_dataset = ctrl$landmarks_per_dataset,
    domain_sizes = domain_sizes,
    offsets = offsets,
    corr_df = corr_df,
    corr_frac = ctrl$landmark_corr_frac,
    seed = ctrl$seed
  )

  if (isTRUE(ctrl$verbose)) {
    message("global_geo_align: running Dijkstra from ", length(landmark_ids), " landmarks.")
  }

  D_L_all <- igraph::distances(
    graph,
    v = as.character(landmark_ids),
    to = igraph::V(graph),
    weights = igraph::E(graph)$weight,
    algorithm = "dijkstra"
  )

  if (any(!is.finite(D_L_all))) {
    stop(
      "Global graph is disconnected for one or more landmarks. ",
      "Add correspondences or increase graph connectivity.",
      call. = FALSE
    )
  }

  D_LL <- D_L_all[, landmark_ids, drop = FALSE]

  k_embed_default <- max(3L * as.integer(ncomp), as.integer(ncomp) + 5L)
  k_embed <- if (is.null(ctrl$k_embed)) k_embed_default else as.integer(ctrl$k_embed)
  mds <- .landmark_mds(D_LL, k_embed = k_embed, eig_tol = ctrl$eig_tol)

  if (mds$k < ncomp) {
    warning("Only ", mds$k, " positive landmark components found; reducing ncomp.", call. = FALSE)
    ncomp <- mds$k
  }

  Y <- .embed_from_landmarks(
    D_L_all = D_L_all,
    mds = mds,
    block_size = ctrl$embed_block_size
  )

  map_fit <- .solve_linear_maps_low_rank(
    pdata = pdata_geom,
    Y = Y,
    offsets = offsets,
    ncomp = as.integer(ncomp),
    ridge_eps = ctrl$ridge_eps,
    ridge_relative = ctrl$ridge_relative,
    ginv_method = ctrl$ginv_method,
    eig_tol = ctrl$eig_tol,
    feature_penalty = if (!is.null(sab_info) && !is.null(sab_info$smoothness)) {
      replicate(
        n_domains,
        sab_info$smoothness$lambda * sab_info$smoothness$laplacian,
        simplify = FALSE
      )
    } else {
      NULL
    }
  )

  if (!is.null(sab_info)) {
    anchor_blocks <- feature_block_indices(pdata_geom)
    A_anchor_list <- lapply(anchor_blocks, function(idx) {
      map_fit$loadings[idx, , drop = FALSE]
    })
    A_voxel_list <- lapply(seq_along(A_anchor_list), function(i) {
      as.matrix(sab_info$B_list[[i]] %*% A_anchor_list[[i]])
    })
    score_list <- lapply(seq_along(A_voxel_list), function(i) {
      as.matrix(pdata[[i]]$x) %*% A_voxel_list[[i]]
    })
    map_fit$anchor_loadings <- do.call(rbind, A_anchor_list)
    map_fit$loadings <- do.call(rbind, A_voxel_list)
    map_fit$scores <- do.call(rbind, score_list)
  }

  feature_blocks <- feature_block_indices(pdata)
  extras <- list(
    correspondences = corr_df,
    domain_scales = stats::setNames(scales, domain_names),
    landmarks = landmark_ids,
    geometry_rank = mds$k,
    geometry_eigenvalues = mds$eigvals,
    map_eigenvalues = map_fit$eigenvalues
  )
  if (!is.null(sab_info)) {
    extras$feature_correspondence <- list(
      type = "spatial_anchor_basis",
      anchors = sab_info$anchors,
      sigma = sab_info$sigma,
      nn_anchors = sab_info$nn_anchors,
      transform = sab_info$transform,
      transform_params = sab_info$transform_params,
      B = if (isTRUE(sab_info$store_B)) sab_info$B_list else NULL,
      anchor_loadings = map_fit$anchor_loadings,
      smoothness = if (is.null(sab_info$smoothness)) NULL else list(
        lambda = sab_info$smoothness$lambda,
        knn = sab_info$smoothness$knn,
        sigma = sab_info$smoothness$sigma,
        weight_mode = sab_info$smoothness$weight_mode,
        laplacian = if (isTRUE(sab_info$smoothness$store_laplacian)) sab_info$smoothness$laplacian else NULL
      )
    )
  }

  result <- new_alignment_result(
    scores = map_fit$scores,
    loadings = map_fit$loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "global_geo_align",
    extras = extras
  )

  result
}

#' @rdname global_geo_align
#' @method global_geo_align list
#' @export
global_geo_align.list <- function(
  data,
  y = NULL,
  correspondences = NULL,
  feature_correspondence = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = global_geo_align_control(),
  ...
) {
  if (!is.list(data) || length(data) < 2L) {
    stop("`data` must be a list of at least two matrices.", call. = FALSE)
  }
  if (!all(vapply(data, is.matrix, logical(1)))) {
    stop("All entries in `data` must be matrices.", call. = FALSE)
  }

  strata <- lapply(seq_along(data), function(i) {
    Xi <- as.matrix(data[[i]])
    list(
      x = Xi,
      design = data.frame(row_id = seq_len(nrow(Xi)))
    )
  })
  names(strata) <- names(data)
  if (is.null(names(strata)) || any(names(strata) == "")) {
    names(strata) <- paste0("domain", seq_along(strata))
  }
  class(strata) <- "hyperdesign"

  global_geo_align.hyperdesign(
    data = strata,
    y = y,
    correspondences = correspondences,
    feature_correspondence = feature_correspondence,
    preproc = preproc,
    ncomp = ncomp,
    control = control,
    ...
  )
}

#' @rdname global_geo_align
#' @method global_geo_align default
#' @export
global_geo_align.default <- function(data, ...) {
  stop("No applicable method for global_geo_align().", call. = FALSE)
}

.resolve_correspondences <- function(
  pdata,
  domain_names,
  domain_sizes,
  y_expr,
  correspondences,
  max_pairs_per_label,
  seed
) {
  corr_df <- NULL
  if (!is.null(correspondences)) {
    corr_df <- .normalize_correspondence_table(correspondences, domain_names, domain_sizes)
  } else if (!is.null(y_expr)) {
    y_sym <- rlang::ensym(y_expr)
    y_name <- rlang::as_name(y_sym)
    corr_df <- .build_label_correspondences(
      pdata = pdata,
      y_name = y_name,
      max_pairs_per_label = max_pairs_per_label,
      seed = seed
    )
  } else {
    corr_df <- data.frame(
      domain_i = integer(0),
      index_i = integer(0),
      domain_j = integer(0),
      index_j = integer(0)
    )
  }

  if (!nrow(corr_df)) {
    return(corr_df)
  }

  ds <- pmin(corr_df$domain_i, corr_df$domain_j)
  dt <- pmax(corr_df$domain_i, corr_df$domain_j)
  is <- ifelse(corr_df$domain_i <= corr_df$domain_j, corr_df$index_i, corr_df$index_j)
  jt <- ifelse(corr_df$domain_i <= corr_df$domain_j, corr_df$index_j, corr_df$index_i)

  out <- unique(data.frame(
    domain_i = as.integer(ds),
    index_i = as.integer(is),
    domain_j = as.integer(dt),
    index_j = as.integer(jt)
  ))
  out <- out[out$domain_i != out$domain_j, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.normalize_correspondence_table <- function(correspondences, domain_names, domain_sizes) {
  if (!is.data.frame(correspondences)) {
    stop("`correspondences` must be a data.frame.", call. = FALSE)
  }

  col_pick <- function(cands) {
    hit <- cands[cands %in% names(correspondences)]
    if (!length(hit)) return(NULL)
    hit[[1]]
  }

  c_di <- col_pick(c("domain_i", "s", "from_domain", "domain1", "i_domain"))
  c_ii <- col_pick(c("index_i", "i", "from_index", "sample_i", "idx_i"))
  c_dj <- col_pick(c("domain_j", "t", "to_domain", "domain2", "j_domain"))
  c_ij <- col_pick(c("index_j", "j", "to_index", "sample_j", "idx_j"))

  if (any(vapply(list(c_di, c_ii, c_dj, c_ij), is.null, logical(1)))) {
    stop(
      "`correspondences` must provide columns for domain/index pairs. ",
      "Accepted names include domain_i/index_i/domain_j/index_j or s/i/t/j.",
      call. = FALSE
    )
  }

  di <- correspondences[[c_di]]
  dj <- correspondences[[c_dj]]
  ii <- correspondences[[c_ii]]
  ij <- correspondences[[c_ij]]

  map_domain <- function(v) {
    if (is.character(v) || is.factor(v)) {
      idx <- match(as.character(v), domain_names)
    } else {
      idx <- as.integer(v)
    }
    idx
  }

  di_idx <- map_domain(di)
  dj_idx <- map_domain(dj)
  ii <- as.integer(ii)
  ij <- as.integer(ij)

  ok <- is.finite(di_idx) & is.finite(dj_idx) & is.finite(ii) & is.finite(ij)
  out <- data.frame(
    domain_i = as.integer(di_idx[ok]),
    index_i = as.integer(ii[ok]),
    domain_j = as.integer(dj_idx[ok]),
    index_j = as.integer(ij[ok])
  )
  out <- out[
    out$domain_i >= 1 & out$domain_i <= length(domain_names) &
      out$domain_j >= 1 & out$domain_j <= length(domain_names),
    ,
    drop = FALSE
  ]
  out <- out[
    out$index_i >= 1 & out$index_j >= 1 &
      out$index_i <= domain_sizes[out$domain_i] &
      out$index_j <= domain_sizes[out$domain_j],
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

.build_label_correspondences <- function(pdata, y_name, max_pairs_per_label, seed) {
  set.seed(seed)
  m <- length(pdata)
  rows <- vector("list", max(1L, m * (m - 1L) / 2L))
  ptr <- 0L

  labels <- lapply(seq_len(m), function(i) {
    dsg <- pdata[[i]]$design
    if (is.null(dsg) || !is.data.frame(dsg) || !y_name %in% names(dsg)) {
      stop("Label column `", y_name, "` missing in one or more domains.", call. = FALSE)
    }
    dsg[[y_name]]
  })

  for (i in seq_len(m - 1L)) {
    li <- labels[[i]]
    for (j in (i + 1L):m) {
      lj <- labels[[j]]
      shared <- intersect(unique(li[!is.na(li)]), unique(lj[!is.na(lj)]))
      if (!length(shared)) next

      di <- vector("list", length(shared))
      dj <- vector("list", length(shared))
      cnt <- 0L
      for (lab in shared) {
        idx_i <- which(li == lab)
        idx_j <- which(lj == lab)
        if (!length(idx_i) || !length(idx_j)) next
        k <- min(length(idx_i), length(idx_j), max_pairs_per_label)
        if (k < 1L) next
        if (length(idx_i) > k) idx_i <- sample(idx_i, k)
        if (length(idx_j) > k) idx_j <- sample(idx_j, k)
        if (length(idx_i) != length(idx_j)) {
          k <- min(length(idx_i), length(idx_j))
          idx_i <- idx_i[seq_len(k)]
          idx_j <- idx_j[seq_len(k)]
        }
        cnt <- cnt + 1L
        di[[cnt]] <- idx_i
        dj[[cnt]] <- idx_j
      }
      if (cnt > 0L) {
        ptr <- ptr + 1L
        ii <- unlist(di[seq_len(cnt)], use.names = FALSE)
        jj <- unlist(dj[seq_len(cnt)], use.names = FALSE)
        rows[[ptr]] <- data.frame(
          domain_i = i,
          index_i = ii,
          domain_j = j,
          index_j = jj
        )
      }
    }
  }

  if (ptr == 0L) {
    return(data.frame(
      domain_i = integer(0),
      index_i = integer(0),
      domain_j = integer(0),
      index_j = integer(0)
    ))
  }

  out <- do.call(rbind, rows[seq_len(ptr)])
  rownames(out) <- NULL
  out
}

.estimate_domain_scales <- function(pdata, corr_df, method, reg, anchor_max, seed) {
  m <- length(pdata)
  if (method == "none" || !nrow(corr_df)) {
    return(rep(1, m))
  }

  set.seed(seed)
  pair_keys <- unique(cbind(corr_df$domain_i, corr_df$domain_j))
  ratios <- list()
  for (k in seq_len(nrow(pair_keys))) {
    s <- pair_keys[k, 1]
    t <- pair_keys[k, 2]
    block <- corr_df[corr_df$domain_i == s & corr_df$domain_j == t, , drop = FALSE]
    if (nrow(block) < 2L) next
    if (nrow(block) > anchor_max) {
      block <- block[sample.int(nrow(block), anchor_max), , drop = FALSE]
    }
    Xs <- pdata[[s]]$x[block$index_i, , drop = FALSE]
    Xt <- pdata[[t]]$x[block$index_j, , drop = FALSE]
    Da <- as.matrix(stats::dist(Xs))
    Db <- as.matrix(stats::dist(Xt))
    if (!length(Da) || !length(Db)) next

    if (method == "pairwise_eta") {
      den <- sum(Db * Db)
      if (den <= 0) next
      r_st <- sum(Db * Da) / den
    } else {
      uu <- upper.tri(Da)
      rv <- Da[uu] / pmax(Db[uu], 1e-12)
      rv <- rv[is.finite(rv) & rv > 0]
      if (!length(rv)) next
      r_st <- stats::median(rv)
    }
    if (!is.finite(r_st) || r_st <= 0) next
    ratios[[length(ratios) + 1L]] <- c(s, t, r_st)
  }

  if (!length(ratios) || m == 1L) {
    return(rep(1, m))
  }

  ratio_mat <- do.call(rbind, ratios)
  A <- matrix(0, nrow(ratio_mat), m - 1L)
  b <- log(ratio_mat[, 3])
  for (r in seq_len(nrow(ratio_mat))) {
    s <- as.integer(ratio_mat[r, 1])
    t <- as.integer(ratio_mat[r, 2])
    if (s != 1L) A[r, s - 1L] <- A[r, s - 1L] + 1
    if (t != 1L) A[r, t - 1L] <- A[r, t - 1L] - 1
  }

  lhs <- crossprod(A) + reg * diag(m - 1L)
  rhs <- crossprod(A, b)
  x <- tryCatch(
    as.numeric(solve(lhs, rhs)),
    error = function(e) rep(0, m - 1L)
  )
  scales <- exp(c(0, x))
  scales / stats::median(scales)
}

.build_global_graph_edges <- function(
  pdata,
  domain_sizes,
  offsets,
  corr_df,
  scales,
  knn,
  metric,
  corr_edge_weight
) {
  edge_parts <- vector("list", length(pdata) + 1L)
  ptr <- 0L

  for (d in seq_along(pdata)) {
    X <- as.matrix(pdata[[d]]$x)
    n <- nrow(X)
    k_use <- min(knn + 1L, n)
    if (k_use <= 1L) next

    Xq <- X
    if (metric == "cosine") {
      rn <- sqrt(rowSums(Xq * Xq))
      rn[rn == 0] <- 1
      Xq <- Xq / rn
    }

    nn <- RANN::nn2(Xq, query = Xq, k = k_use)
    idx <- nn$nn.idx[, -1, drop = FALSE]
    dst <- nn$nn.dists[, -1, drop = FALSE]

    i_local <- rep.int(seq_len(n), times = ncol(idx))
    j_local <- as.vector(idx)
    w <- as.vector(dst)
    if (metric == "cosine") {
      w <- 0.5 * (w^2)
    }
    w <- pmax(w, .Machine$double.eps) * scales[d]

    ptr <- ptr + 1L
    edge_parts[[ptr]] <- data.frame(
      from = as.character(offsets[d] + i_local),
      to = as.character(offsets[d] + j_local),
      weight = w
    )
  }

  if (nrow(corr_df)) {
    from <- offsets[corr_df$domain_i] + corr_df$index_i
    to <- offsets[corr_df$domain_j] + corr_df$index_j
    ptr <- ptr + 1L
    edge_parts[[ptr]] <- data.frame(
      from = as.character(from),
      to = as.character(to),
      weight = rep(corr_edge_weight, length(from))
    )
  }

  if (ptr == 0L) {
    return(data.frame(from = character(0), to = character(0), weight = numeric(0)))
  }
  do.call(rbind, edge_parts[seq_len(ptr)])
}

.select_landmarks <- function(
  n_landmarks_total,
  landmarks_per_dataset,
  domain_sizes,
  offsets,
  corr_df,
  corr_frac,
  seed
) {
  set.seed(seed)
  total_n <- sum(domain_sizes)
  r_total <- min(as.integer(n_landmarks_total), total_n)
  if (r_total <= 0L) return(integer(0))
  if (r_total >= total_n) return(seq_len(total_n))

  corr_nodes <- integer(0)
  if (nrow(corr_df)) {
    corr_nodes <- unique(c(
      offsets[corr_df$domain_i] + corr_df$index_i,
      offsets[corr_df$domain_j] + corr_df$index_j
    ))
  }

  n_corr <- min(length(corr_nodes), floor(r_total * corr_frac))
  landmarks <- integer(0)
  if (n_corr > 0L) {
    landmarks <- sample(corr_nodes, n_corr)
  }

  remaining <- r_total - length(landmarks)
  if (remaining <= 0L) {
    return(sort(unique(landmarks))[seq_len(r_total)])
  }

  m <- length(domain_sizes)
  if (is.null(landmarks_per_dataset)) {
    raw <- remaining * domain_sizes / sum(domain_sizes)
    alloc <- floor(raw)
    rem <- remaining - sum(alloc)
    if (rem > 0L) {
      ord <- order(raw - alloc, decreasing = TRUE)
      alloc[ord[seq_len(rem)]] <- alloc[ord[seq_len(rem)]] + 1L
    }
  } else {
    alloc <- rep.int(min(landmarks_per_dataset, max(domain_sizes)), m)
    if (sum(alloc) > remaining) {
      raw <- alloc / sum(alloc) * remaining
      alloc <- floor(raw)
      rem <- remaining - sum(alloc)
      if (rem > 0L) {
        ord <- order(raw - alloc, decreasing = TRUE)
        alloc[ord[seq_len(rem)]] <- alloc[ord[seq_len(rem)]] + 1L
      }
    }
  }

  for (d in seq_len(m)) {
    if (alloc[d] <= 0L) next
    ids <- offsets[d] + seq_len(domain_sizes[d])
    ids <- setdiff(ids, landmarks)
    if (!length(ids)) next
    k <- min(alloc[d], length(ids))
    landmarks <- c(landmarks, sample(ids, k))
  }

  if (length(landmarks) < r_total) {
    pool <- setdiff(seq_len(total_n), landmarks)
    if (length(pool)) {
      landmarks <- c(landmarks, sample(pool, min(length(pool), r_total - length(landmarks))))
    }
  }

  sort(unique(landmarks))[seq_len(min(r_total, length(unique(landmarks))))]
}

.landmark_mds <- function(D_LL, k_embed, eig_tol) {
  r <- nrow(D_LL)
  if (r < 2L) {
    stop("Need at least two landmarks for MDS.", call. = FALSE)
  }

  S <- D_LL^2
  row_mean <- rowMeans(S)
  col_mean <- colMeans(S)
  total_mean <- mean(S)
  B <- -0.5 * (S - row_mean - matrix(col_mean, nrow = r, ncol = r, byrow = TRUE) + total_mean)
  B <- 0.5 * (B + t(B))

  k_target <- min(k_embed, r - 1L)
  if (k_target < 1L) {
    stop("`k_embed` too small for the landmark set.", call. = FALSE)
  }

  eigvals <- NULL
  eigvecs <- NULL
  if (requireNamespace("RSpectra", quietly = TRUE) && k_target < (r - 1L)) {
    eg <- RSpectra::eigs_sym(B, k = k_target, which = "LA")
    ord <- order(eg$values, decreasing = TRUE)
    eigvals <- eg$values[ord]
    eigvecs <- eg$vectors[, ord, drop = FALSE]
  } else {
    eg <- eigen(B, symmetric = TRUE)
    eigvals <- eg$values[seq_len(k_target)]
    eigvecs <- eg$vectors[, seq_len(k_target), drop = FALSE]
  }

  keep <- which(eigvals > eig_tol)
  if (!length(keep)) {
    stop("No positive landmark eigenvalues found.", call. = FALSE)
  }

  eigvals <- eigvals[keep]
  eigvecs <- eigvecs[, keep, drop = FALSE]
  Tmat <- diag(1 / sqrt(eigvals), nrow = length(eigvals)) %*% t(eigvecs)

  list(
    T = Tmat,
    row_mean = row_mean,
    total_mean = total_mean,
    eigvals = eigvals,
    k = length(eigvals)
  )
}

.embed_from_landmarks <- function(D_L_all, mds, block_size) {
  r <- nrow(D_L_all)
  n <- ncol(D_L_all)
  k <- mds$k
  Y <- matrix(0, nrow = n, ncol = k)

  mu_row <- mds$row_mean
  mu_total <- mds$total_mean
  Tmat <- mds$T

  for (start in seq.int(1L, n, by = block_size)) {
    end <- min(start + block_size - 1L, n)
    idx <- start:end
    S_blk <- D_L_all[, idx, drop = FALSE]^2
    mu_col <- colMeans(S_blk)
    B_blk <- -0.5 * (
      S_blk -
        matrix(mu_row, nrow = r, ncol = length(idx), byrow = FALSE) -
        matrix(mu_col, nrow = r, ncol = length(idx), byrow = TRUE) +
        mu_total
    )
    Y[idx, ] <- t(Tmat %*% B_blk)
  }

  Y
}

.solve_linear_maps_low_rank <- function(
  pdata,
  Y,
  offsets,
  ncomp,
  ridge_eps,
  ridge_relative,
  ginv_method,
  eig_tol,
  feature_penalty = NULL
) {
  m <- length(pdata)
  k <- ncol(Y)
  GinvC_list <- vector("list", m)
  R_small <- matrix(0, nrow = k, ncol = k)

  if (!is.null(feature_penalty)) {
    if (!is.list(feature_penalty) || length(feature_penalty) != m) {
      stop("`feature_penalty` must be NULL or a list with one matrix per domain.", call. = FALSE)
    }
  }

  for (d in seq_len(m)) {
    Xi <- as.matrix(pdata[[d]]$x)
    idx <- (offsets[d] + 1L):offsets[d + 1L]
    Yi <- Y[idx, , drop = FALSE]

    Ci <- crossprod(Xi, Yi)
    eps <- ridge_eps
    if (ridge_relative) {
      feat_energy <- sum(colSums(Xi^2)) / ncol(Xi)
      if (is.finite(feat_energy) && feat_energy > 0) {
        eps <- ridge_eps * feat_energy
      }
    }
    eps <- max(eps, .Machine$double.eps)

    GinvCi <- .solve_ginv_c(
      X = Xi,
      C = Ci,
      eps = eps,
      method = ginv_method,
      penalty = if (is.null(feature_penalty)) NULL else feature_penalty[[d]]
    )
    GinvC_list[[d]] <- GinvCi
    R_small <- R_small + crossprod(Ci, GinvCi)
  }

  R_small <- 0.5 * (R_small + t(R_small))

  d_out <- min(as.integer(ncomp), ncol(R_small))
  if (d_out < 1L) {
    stop("No components available after low-rank reduction.", call. = FALSE)
  }

  if (requireNamespace("RSpectra", quietly = TRUE) && d_out < (ncol(R_small) - 1L)) {
    eg <- RSpectra::eigs_sym(R_small, k = d_out, which = "LA")
    ord <- order(eg$values, decreasing = TRUE)
    vals <- eg$values[ord]
    U <- eg$vectors[, ord, drop = FALSE]
  } else {
    eg <- eigen(R_small, symmetric = TRUE)
    vals <- eg$values[seq_len(d_out)]
    U <- eg$vectors[, seq_len(d_out), drop = FALSE]
  }

  keep <- which(vals > eig_tol)
  if (!length(keep)) {
    stop("All map eigenvalues are numerically zero.", call. = FALSE)
  }
  vals <- vals[keep]
  U <- U[, keep, drop = FALSE]

  scale <- diag(1 / sqrt(vals), nrow = length(vals))
  A_list <- lapply(GinvC_list, function(GinvCi) GinvCi %*% U %*% scale)

  score_list <- lapply(seq_len(m), function(i) {
    Xi <- as.matrix(pdata[[i]]$x)
    Xi %*% A_list[[i]]
  })

  list(
    loadings = do.call(rbind, A_list),
    scores = do.call(rbind, score_list),
    eigenvalues = vals
  )
}

.solve_ginv_c <- function(X, C, eps, method = c("auto", "woodbury", "direct"), penalty = NULL) {
  method <- match.arg(method)
  n <- nrow(X)
  p <- ncol(X)
  has_penalty <- !is.null(penalty)
  use_woodbury <- (method == "woodbury" || (method == "auto" && p > n)) && !has_penalty

  if (use_woodbury) {
    XXt <- tcrossprod(X)
    M <- diag(n) + XXt / eps
    R <- chol(M)
    XC <- X %*% C
    tmp <- forwardsolve(t(R), XC, upper.tri = FALSE, transpose = FALSE)
    Minv_XC <- backsolve(R, tmp, upper.tri = TRUE, transpose = FALSE)
    return(C / eps - crossprod(X, Minv_XC) / (eps^2))
  }

  G <- crossprod(X) + diag(eps, p)
  if (has_penalty) {
    P <- as.matrix(penalty)
    if (!all(dim(P) == c(p, p))) {
      stop("Feature penalty matrix has incompatible dimensions.", call. = FALSE)
    }
    G <- G + P
  }
  G <- 0.5 * (G + t(G))
  R <- chol(G)
  tmp <- forwardsolve(t(R), C, upper.tri = FALSE, transpose = FALSE)
  backsolve(R, tmp, upper.tri = TRUE, transpose = FALSE)
}

.apply_spatial_anchor_basis <- function(pdata, domain_names, spec, seed) {
  if (!is.list(spec)) {
    stop("`feature_correspondence` must be a named list.", call. = FALSE)
  }
  enabled <- if (is.null(spec$enabled)) TRUE else isTRUE(spec$enabled)
  if (!enabled) {
    return(NULL)
  }

  voxel_coords <- spec$voxel_coords
  if (is.null(voxel_coords)) {
    voxel_coords <- spec$coords
  }
  if (is.null(voxel_coords)) {
    stop(
      "`feature_correspondence` requires `voxel_coords` (or `coords`) as a list of p_s x 3 matrices.",
      call. = FALSE
    )
  }

  coords_list <- .resolve_voxel_coords(
    voxel_coords = voxel_coords,
    pdata = pdata,
    domain_names = domain_names
  )

  transform <- if (is.null(spec$transform)) "identity" else as.character(spec$transform[[1L]])
  transform <- match.arg(transform, c("identity", "pca_normalized"))
  transform_fit <- lapply(coords_list, .fit_coordinate_transform, mode = transform)
  coords_tx <- lapply(transform_fit, `[[`, "coords")
  transform_params <- lapply(transform_fit, `[[`, "params")

  ref_idx <- .resolve_reference_dataset(
    ref = spec$reference_dataset,
    domain_names = domain_names
  )
  if (is.null(ref_idx)) {
    ref_coords <- do.call(rbind, coords_tx)
  } else {
    ref_coords <- coords_tx[[ref_idx]]
  }

  p_total <- nrow(ref_coords)
  k_default <- min(1500L, max(32L, as.integer(round(p_total / 4))))
  num_anchors <- if (is.null(spec$num_anchors)) k_default else as.integer(spec$num_anchors[[1L]])
  num_anchors <- max(2L, min(num_anchors, p_total))

  anchor_method <- if (is.null(spec$anchor_method)) "kmeans" else as.character(spec$anchor_method[[1L]])
  anchor_method <- match.arg(anchor_method, c("kmeans", "grid"))
  sab_seed <- if (is.null(spec$seed)) as.integer(seed) else as.integer(spec$seed[[1L]])
  anchors <- .select_anchor_points(
    coords = ref_coords,
    n_anchors = num_anchors,
    method = anchor_method,
    seed = sab_seed
  )

  nn_anchors <- if (is.null(spec$nn_anchors)) 6L else as.integer(spec$nn_anchors[[1L]])
  if (!is.finite(nn_anchors) || nn_anchors < 1L) {
    stop("`feature_correspondence$nn_anchors` must be >= 1.", call. = FALSE)
  }
  nn_anchors <- min(nn_anchors, nrow(anchors))

  sigma <- if (is.null(spec$sigma)) NULL else as.numeric(spec$sigma[[1L]])
  if (!is.null(sigma) && (!is.finite(sigma) || sigma <= 0)) {
    stop("`feature_correspondence$sigma` must be positive when provided.", call. = FALSE)
  }
  if (is.null(sigma)) {
    sigma <- .estimate_anchor_sigma(anchors)
  }
  sigma <- max(sigma, .Machine$double.eps)

  B_list <- lapply(coords_tx, function(ci) {
    .build_sab_matrix(
      voxel_coords = ci,
      anchors = anchors,
      sigma = sigma,
      nn_anchors = nn_anchors
    )
  })

  pdata_tilde <- lapply(seq_along(pdata), function(i) {
    Xi <- as.matrix(pdata[[i]]$x)
    list(
      x = as.matrix(Xi %*% B_list[[i]]),
      design = pdata[[i]]$design
    )
  })
  names(pdata_tilde) <- domain_names

  store_B <- if (is.null(spec$store_B)) TRUE else isTRUE(spec$store_B)
  smoothness <- .build_anchor_smoothness(anchors = anchors, smoothness_spec = spec$smoothness)

  list(
    pdata_tilde = pdata_tilde,
    anchors = anchors,
    sigma = sigma,
    nn_anchors = nn_anchors,
    transform = transform,
    transform_params = transform_params,
    B_list = B_list,
    store_B = store_B,
    smoothness = smoothness
  )
}

.resolve_voxel_coords <- function(voxel_coords, pdata, domain_names) {
  if (!is.list(voxel_coords) || !length(voxel_coords)) {
    stop("`feature_correspondence$voxel_coords` must be a non-empty list.", call. = FALSE)
  }

  if (!is.null(names(voxel_coords)) && all(names(voxel_coords) != "")) {
    missing_names <- setdiff(domain_names, names(voxel_coords))
    if (length(missing_names)) {
      stop(
        "Missing voxel coordinate entries for domain(s): ",
        paste(missing_names, collapse = ", "),
        call. = FALSE
      )
    }
    coords_list <- voxel_coords[domain_names]
  } else {
    if (length(voxel_coords) != length(pdata)) {
      stop(
        "`feature_correspondence$voxel_coords` must have one entry per domain.",
        call. = FALSE
      )
    }
    coords_list <- voxel_coords
    names(coords_list) <- domain_names
  }

  out <- lapply(seq_along(coords_list), function(i) {
    Ci <- as.matrix(coords_list[[i]])
    if (!is.numeric(Ci) || ncol(Ci) != 3L) {
      stop("Voxel coordinates must be numeric matrices with 3 columns.", call. = FALSE)
    }
    p_i <- ncol(pdata[[i]]$x)
    if (nrow(Ci) != p_i) {
      stop(
        "Voxel coordinate rows must match number of features for domain `",
        domain_names[[i]],
        "` (expected ",
        p_i,
        ", got ",
        nrow(Ci),
        ").",
        call. = FALSE
      )
    }
    Ci
  })
  names(out) <- domain_names
  out
}

.fit_coordinate_transform <- function(coords, mode = c("identity", "pca_normalized")) {
  mode <- match.arg(mode)
  if (mode == "identity") {
    return(list(coords = coords, params = list(mode = mode)))
  }

  center <- colMeans(coords)
  xc <- sweep(coords, 2L, center, FUN = "-")
  sv <- svd(xc)
  rot <- sv$v
  xr <- xc %*% rot
  scale <- apply(abs(xr), 2L, max)
  scale[!is.finite(scale) | scale <= .Machine$double.eps] <- 1
  xn <- sweep(xr, 2L, scale, FUN = "/")

  list(
    coords = xn,
    params = list(
      mode = mode,
      center = center,
      rotation = rot,
      scale = scale
    )
  )
}

.resolve_reference_dataset <- function(ref, domain_names) {
  if (is.null(ref)) {
    return(NULL)
  }
  if (is.numeric(ref) && length(ref) == 1L && is.finite(ref)) {
    idx <- as.integer(ref)
    if (idx < 1L || idx > length(domain_names)) {
      stop("`feature_correspondence$reference_dataset` index is out of bounds.", call. = FALSE)
    }
    return(idx)
  }
  if ((is.character(ref) || is.factor(ref)) && length(ref) == 1L) {
    idx <- match(as.character(ref), domain_names)
    if (is.na(idx)) {
      stop(
        "`feature_correspondence$reference_dataset` name not found in domains: ",
        as.character(ref),
        call. = FALSE
      )
    }
    return(as.integer(idx))
  }
  stop("`feature_correspondence$reference_dataset` must be a single index or domain name.", call. = FALSE)
}

.select_anchor_points <- function(coords, n_anchors, method = c("kmeans", "grid"), seed) {
  method <- match.arg(method)
  if (n_anchors >= nrow(coords)) {
    return(coords)
  }

  set.seed(seed)
  if (method == "grid") {
    mins <- apply(coords, 2L, min)
    maxs <- apply(coords, 2L, max)
    side <- max(2L, ceiling(n_anchors^(1 / 3)))
    gx <- seq(mins[1], maxs[1], length.out = side)
    gy <- seq(mins[2], maxs[2], length.out = side)
    gz <- seq(mins[3], maxs[3], length.out = side)
    grid <- as.matrix(expand.grid(gx, gy, gz))
    if (nrow(grid) > n_anchors) {
      grid <- grid[sample.int(nrow(grid), n_anchors), , drop = FALSE]
    }
    return(grid)
  }

  km <- stats::kmeans(coords, centers = n_anchors, iter.max = 50L, nstart = 3L)
  km$centers
}

.estimate_anchor_sigma <- function(anchors) {
  if (nrow(anchors) < 2L) {
    return(1)
  }
  nn <- RANN::nn2(data = anchors, query = anchors, k = 2L)
  d <- nn$nn.dists[, 2L]
  d <- d[is.finite(d) & d > 0]
  if (!length(d)) {
    return(1)
  }
  as.numeric(stats::median(d))
}

.build_sab_matrix <- function(voxel_coords, anchors, sigma, nn_anchors) {
  nn <- RANN::nn2(data = anchors, query = voxel_coords, k = nn_anchors)
  idx <- nn$nn.idx
  d <- nn$nn.dists
  w <- exp(-(d^2) / (2 * sigma^2))
  rs <- rowSums(w)
  rs[rs <= .Machine$double.eps] <- 1
  w <- w / rs

  p <- nrow(voxel_coords)
  K <- nrow(anchors)
  ii <- rep.int(seq_len(p), times = nn_anchors)
  jj <- as.vector(idx)
  xx <- as.vector(w)
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(p, K))
}

.build_anchor_smoothness <- function(anchors, smoothness_spec = NULL) {
  if (is.null(smoothness_spec)) {
    return(NULL)
  }
  if (!is.list(smoothness_spec)) {
    stop("`feature_correspondence$smoothness` must be a named list.", call. = FALSE)
  }
  enabled <- if (is.null(smoothness_spec$enabled)) TRUE else isTRUE(smoothness_spec$enabled)
  if (!enabled) {
    return(NULL)
  }

  lambda <- if (is.null(smoothness_spec$lambda)) 0 else as.numeric(smoothness_spec$lambda[[1L]])
  if (!is.finite(lambda) || lambda < 0) {
    stop("`feature_correspondence$smoothness$lambda` must be non-negative.", call. = FALSE)
  }
  if (lambda <= 0) {
    return(NULL)
  }

  knn <- if (is.null(smoothness_spec$knn)) 12L else as.integer(smoothness_spec$knn[[1L]])
  if (!is.finite(knn) || knn < 1L) {
    stop("`feature_correspondence$smoothness$knn` must be >= 1.", call. = FALSE)
  }
  knn <- min(knn, max(1L, nrow(anchors) - 1L))

  weight_mode <- if (is.null(smoothness_spec$weight_mode)) "heat" else as.character(smoothness_spec$weight_mode[[1L]])
  weight_mode <- match.arg(weight_mode, c("heat", "binary"))

  sigma <- if (is.null(smoothness_spec$sigma)) .estimate_anchor_sigma(anchors) else as.numeric(smoothness_spec$sigma[[1L]])
  if (!is.finite(sigma) || sigma <= 0) {
    stop("`feature_correspondence$smoothness$sigma` must be positive.", call. = FALSE)
  }

  S <- neighborweights::spatial_adjacency(
    anchors,
    sigma = sigma,
    weight_mode = weight_mode,
    nnk = knn,
    normalized = FALSE,
    stochastic = FALSE,
    handle_isolates = "self_loop"
  )
  S <- Matrix::Matrix(S, sparse = TRUE)
  S <- Matrix::forceSymmetric(S, uplo = "U")
  S <- as(S, "dgCMatrix")
  Matrix::diag(S) <- 0
  S <- Matrix::drop0(S)

  deg <- Matrix::rowSums(S)
  L <- Matrix::Diagonal(x = deg) - S
  store_laplacian <- if (is.null(smoothness_spec$store_laplacian)) FALSE else isTRUE(smoothness_spec$store_laplacian)

  list(
    lambda = lambda,
    knn = knn,
    sigma = sigma,
    weight_mode = weight_mode,
    laplacian = L,
    store_laplacian = store_laplacian
  )
}
