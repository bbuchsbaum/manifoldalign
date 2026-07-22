# RAPID-MA sparse relation graph foundation ---------------------------------

# These helpers implement the node-relation contract used by RAPID-MA. They
# deliberately stay internal until the top-level rapid_ma() API is available.
# All returned relations are nonnegative dgCMatrix objects. Production builders
# use bounded-width views and capped nearest-neighbour queries; they never form
# a node-by-node distance matrix.

.rapid_empty_relation <- function(n, dimnames = NULL) {
  Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    x = numeric(0),
    dims = c(as.integer(n), as.integer(n)),
    dimnames = dimnames
  )
}

.rapid_force_dgC <- function(A) {
  if (!inherits(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }
  if (methods::is(A, "symmetricMatrix") ||
      methods::is(A, "triangularMatrix") ||
      methods::is(A, "diagonalMatrix")) {
    A <- methods::as(A, "generalMatrix")
  }
  A <- Matrix::drop0(Matrix::Matrix(A, sparse = TRUE))
  sm <- Matrix::summary(A)
  if (!nrow(sm)) {
    return(.rapid_empty_relation(nrow(A), dimnames = dimnames(A)))
  }
  Matrix::sparseMatrix(
    i = as.integer(sm$i),
    j = as.integer(sm$j),
    x = as.numeric(sm$x),
    dims = dim(A),
    dimnames = dimnames(A)
  )
}

.rapid_int_scalar <- function(x, name, min_value = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      abs(x - round(x)) > 1e-8 || x < min_value) {
    stop("`", name, "` must be an integer >= ", min_value, ".", call. = FALSE)
  }
  as.integer(x)
}

.rapid_number_scalar <- function(x, name, lower = -Inf, strict = FALSE) {
  bad <- length(x) != 1L || !is.numeric(x) || !is.finite(x)
  if (!bad) {
    bad <- if (strict) x <= lower else x < lower
  }
  if (bad) {
    op <- if (strict) ">" else ">="
    stop("`", name, "` must be a finite number ", op, " ", lower, ".",
         call. = FALSE)
  }
  as.numeric(x)
}

.rapid_validate_numeric_matrix <- function(X, name, n_rows = NULL) {
  if (!(is.matrix(X) || inherits(X, "Matrix"))) {
    stop("`", name, "` must be a numeric matrix or Matrix object.", call. = FALSE)
  }
  if (length(dim(X)) != 2L || nrow(X) < 1L || ncol(X) < 1L) {
    stop("`", name, "` must have at least one row and one column.", call. = FALSE)
  }
  if (!is.null(n_rows) && nrow(X) != n_rows) {
    stop("`", name, "` must have ", n_rows, " rows.", call. = FALSE)
  }
  vals <- if (inherits(X, "Matrix")) X@x else X
  if (!is.numeric(vals)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }
  invisible(TRUE)
}

.rapid_hash_int <- function(key, seed = 1L, modulo = 2147483629) {
  codes <- utf8ToInt(enc2utf8(as.character(key)))
  h <- (as.double(seed) %% modulo + 104729) %% modulo
  if (length(codes)) {
    for (code in codes) {
      h <- (h * 131 + as.double(code) + 17) %% modulo
    }
  }
  h
}

.rapid_hash_bucket <- function(key, width, seed = 1L) {
  as.integer(.rapid_hash_int(key, seed = seed) %% width) + 1L
}

.rapid_hash_sign <- function(key, seed = 1L) {
  if ((.rapid_hash_int(paste0(key, ":sign"), seed = seed) %% 2) == 0) -1 else 1
}

.rapid_scale_dense <- function(Z, drop_constant = FALSE) {
  Z <- as.matrix(Z)
  if (!ncol(Z)) {
    return(list(
      view = matrix(numeric(0), nrow(Z), 0L),
      center = numeric(0),
      scale = numeric(0),
      active = logical(0)
    ))
  }
  center <- colMeans(Z)
  center[!is.finite(center)] <- 0
  Zc <- sweep(Z, 2L, center, "-")
  denom <- max(nrow(Zc) - 1L, 1L)
  scale <- sqrt(colSums(Zc * Zc) / denom)
  active <- is.finite(scale) & scale > 1e-12
  scale[!active] <- 1
  out <- sweep(Zc, 2L, scale, "/")
  out[!is.finite(out)] <- 0
  if (isTRUE(drop_constant)) {
    out <- out[, active, drop = FALSE]
  }
  list(view = out, center = center, scale = scale, active = active)
}

.rapid_feature_view <- function(X, max_dim = 32L, seed = 1L) {
  .rapid_validate_numeric_matrix(X, "X")
  vals <- if (inherits(X, "Matrix")) X@x else X
  if (any(!is.finite(vals))) {
    stop("`X` must contain only finite values.", call. = FALSE)
  }
  max_dim <- .rapid_int_scalar(max_dim, "feature_sketch_dim")
  p <- ncol(X)

  projected <- p > max_dim
  buckets <- signs <- NULL
  if (projected) {
    keys <- paste0("feature:", seq_len(p))
    buckets <- vapply(keys, .rapid_hash_bucket, integer(1),
                      width = max_dim, seed = seed)
    signs <- vapply(keys, .rapid_hash_sign, numeric(1), seed = seed)
    projection <- Matrix::sparseMatrix(
      i = seq_len(p),
      j = buckets,
      x = signs,
      dims = c(p, max_dim)
    )
    Z <- as.matrix(X %*% projection)
    colnames(Z) <- paste0("feature_hash_", seq_len(max_dim))
  } else {
    Z <- as.matrix(X)
    if (is.null(colnames(Z))) colnames(Z) <- paste0("feature_", seq_len(p))
  }

  scaled <- .rapid_scale_dense(Z, drop_constant = FALSE)
  list(
    view = scaled$view,
    metadata = list(
      original_dim = as.integer(p),
      output_dim = as.integer(ncol(scaled$view)),
      projected = projected,
      buckets = buckets,
      signs = signs,
      center = scaled$center,
      scale = scaled$scale,
      seed = as.integer(seed)
    )
  )
}

.rapid_validate_attributes <- function(attributes, n) {
  if (!is.data.frame(attributes)) {
    stop("`attributes` must be a data.frame.", call. = FALSE)
  }
  if (nrow(attributes) != n) {
    stop("`attributes` must have ", n, " rows.", call. = FALSE)
  }
  bad <- vapply(attributes, function(x) {
    !(is.numeric(x) || is.integer(x) || is.logical(x) ||
      is.factor(x) || is.character(x))
  }, logical(1))
  if (any(bad)) {
    stop(
      "Unsupported attribute columns: ",
      paste(names(attributes)[bad], collapse = ", "),
      ". Use numeric, integer, logical, factor, or character columns.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rapid_encode_attributes <- function(attributes, n, hash_dim = 32L, seed = 1L) {
  .rapid_validate_attributes(attributes, n)
  hash_dim <- .rapid_int_scalar(hash_dim, "attribute_hash_dim")

  column_names <- names(attributes)
  encoding_names <- vapply(seq_along(attributes), function(j) {
    name <- column_names[[j]]
    if (is.null(name) || !nzchar(name)) paste0("attribute_", j) else name
  }, character(1))

  numeric_hash <- matrix(0, n, hash_dim)
  categorical_hash <- matrix(0, n, hash_dim)
  missing_hash <- matrix(0, n, hash_dim)
  numeric_stats <- list()
  column_types <- character(length(attributes))
  names(column_types) <- names(attributes)
  missing_by_column <- integer(length(attributes))
  names(missing_by_column) <- names(attributes)
  missing_any <- logical(n)

  for (j in seq_along(attributes)) {
    name <- encoding_names[[j]]
    value <- attributes[[j]]
    missing <- is.na(value)
    if (is.numeric(value) && !is.factor(value) && !is.logical(value)) {
      missing <- missing | !is.finite(as.numeric(value))
    }
    missing_by_column[j] <- sum(missing)
    missing_any <- missing_any | missing

    if (is.numeric(value) && !is.factor(value) && !is.logical(value)) {
      column_types[j] <- "numeric"
      x <- as.numeric(value)
      observed <- x[!missing & is.finite(x)]
      center <- if (length(observed)) stats::median(observed) else 0
      spread <- if (length(observed) > 1L) {
        stats::mad(observed, center = center, constant = 1.4826)
      } else {
        1
      }
      if (!is.finite(spread) || spread < 1e-12) {
        spread <- if (length(observed) > 1L) stats::sd(observed) else 1
      }
      if (!is.finite(spread) || spread < 1e-12) spread <- 1
      z <- (x - center) / spread
      z[missing | !is.finite(z)] <- 0
      z <- pmax(pmin(z, 10), -10)
      key <- paste0("numeric:", name)
      bucket <- .rapid_hash_bucket(key, hash_dim, seed)
      sign <- .rapid_hash_sign(key, seed)
      numeric_hash[, bucket] <- numeric_hash[, bucket] + sign * z
      numeric_stats[[name]] <- list(center = center, scale = spread)
    } else {
      column_types[j] <- "categorical"
      x <- as.character(value)
      rows <- which(!missing)
      if (length(rows)) {
        values <- enc2utf8(x[rows])
        unique_values <- unique(values)
        keys <- paste0("categorical:", name, "=", unique_values)
        buckets <- vapply(keys, .rapid_hash_bucket, integer(1),
                          width = hash_dim, seed = seed)
        signs <- vapply(keys, .rapid_hash_sign, numeric(1), seed = seed)
        value_index <- match(values, unique_values)
        cells <- cbind(rows, buckets[value_index])
        categorical_hash[cells] <- categorical_hash[cells] + signs[value_index]
      }
    }

    if (any(missing)) {
      key <- paste0("missing:", name)
      bucket <- .rapid_hash_bucket(key, hash_dim, seed)
      sign <- .rapid_hash_sign(key, seed)
      rows <- which(missing)
      cells <- cbind(rows, rep.int(bucket, length(rows)))
      missing_hash[cells] <- missing_hash[cells] + sign
    }
  }

  colnames(numeric_hash) <- paste0("attribute_numeric_hash_", seq_len(hash_dim))
  colnames(categorical_hash) <- paste0("attribute_category_hash_", seq_len(hash_dim))
  colnames(missing_hash) <- paste0("attribute_missing_hash_", seq_len(hash_dim))
  raw <- cbind(numeric_hash, categorical_hash, missing_hash)
  scaled <- .rapid_scale_dense(raw, drop_constant = TRUE)

  list(
    view = scaled$view,
    metadata = list(
      column_names = column_names,
      encoding_names = encoding_names,
      column_types = column_types,
      missing_by_column = missing_by_column,
      missing_rows = as.integer(which(missing_any)),
      n_missing_rows = as.integer(sum(missing_any)),
      numeric_stats = numeric_stats,
      hash_dim = hash_dim,
      seed = as.integer(seed),
      center = scaled$center,
      scale = scaled$scale,
      active_columns = which(scaled$active),
      output_dim = as.integer(ncol(scaled$view))
    )
  )
}

.rapid_prepare_positions <- function(positions, n,
                                     mode = c("relative", "common", "shared", "ignore")) {
  mode <- match.arg(mode)
  if (identical(mode, "shared")) mode <- "common"
  .rapid_validate_numeric_matrix(positions, "positions", n_rows = n)
  P <- as.matrix(positions)
  valid <- rowSums(!is.finite(P)) == 0L
  missing_rows <- which(!valid)

  if (identical(mode, "ignore")) {
    return(list(
      relation_view = NULL,
      view = NULL,
      metadata = list(
        mode = mode,
        input_dim = as.integer(ncol(P)),
        missing_indicator = FALSE,
        valid_rows = as.integer(sum(valid)),
        missing_rows = as.integer(missing_rows),
        center = numeric(ncol(P)),
        scale = 1,
        translation_invariant = NA,
        rotation_invariant = NA,
        scale_invariant = NA
      )
    ))
  }

  center <- numeric(ncol(P))
  scale <- 1
  relation_view <- P
  if (identical(mode, "relative") && any(valid)) {
    center <- colMeans(P[valid, , drop = FALSE])
    centered <- sweep(P[valid, , drop = FALSE], 2L, center, "-")
    scale <- sqrt(mean(rowSums(centered * centered)))
    if (!is.finite(scale) || scale < 1e-12) scale <- 1
    relation_view[valid, ] <- centered / scale
  }

  encoder_view <- relation_view
  encoder_view[!is.finite(encoder_view)] <- 0
  if (length(missing_rows)) {
    encoder_view <- cbind(encoder_view, position_missing = as.numeric(!valid))
  }

  list(
    relation_view = relation_view,
    view = encoder_view,
    metadata = list(
      mode = mode,
      input_dim = as.integer(ncol(P)),
      missing_indicator = length(missing_rows) > 0L,
      valid_rows = as.integer(sum(valid)),
      missing_rows = as.integer(missing_rows),
      center = center,
      scale = scale,
      translation_invariant = TRUE,
      rotation_invariant = TRUE,
      scale_invariant = identical(mode, "relative")
    )
  )
}

.rapid_degree_counts <- function(A) {
  if (!length(A@x)) {
    return(list(
      out_degree = integer(nrow(A)),
      in_degree = integer(ncol(A))
    ))
  }
  sm <- Matrix::summary(A)
  list(
    out_degree = tabulate(as.integer(sm$i), nbins = nrow(A)),
    in_degree = tabulate(as.integer(sm$j), nbins = ncol(A))
  )
}

.rapid_mutual_relation <- function(A) {
  if (!length(A@x)) return(A)
  At <- Matrix::t(A)
  mask <- (A != 0) * (At != 0)
  .rapid_force_dgC(Matrix::drop0(((A + At) / 2) * mask))
}

.rapid_cap_out_degree <- function(A, degree_cap) {
  if (!length(A@x)) return(A)
  degree_cap <- .rapid_int_scalar(degree_cap, "degree_cap")
  sm <- Matrix::summary(A)
  counts <- tabulate(as.integer(sm$i), nbins = nrow(A))
  if (!length(counts) || max(counts) <= degree_cap) return(A)

  ord <- order(sm$i, -sm$x, sm$j)
  sm <- sm[ord, , drop = FALSE]
  row_rank <- ave(seq_len(nrow(sm)), sm$i, FUN = seq_along)
  sm <- sm[row_rank <= degree_cap, , drop = FALSE]
  Matrix::sparseMatrix(
    i = as.integer(sm$i),
    j = as.integer(sm$j),
    x = as.numeric(sm$x),
    dims = dim(A),
    dimnames = dimnames(A)
  )
}

.rapid_build_neighbor_relation <- function(view,
                                           k,
                                           degree_cap,
                                           radius = NULL,
                                           weight_mode = c("heat", "binary"),
                                           sigma = NULL,
                                           symmetrize = c("mutual", "none"),
                                           source = "feature") {
  weight_mode <- match.arg(weight_mode)
  symmetrize <- match.arg(symmetrize)
  view <- as.matrix(view)
  n <- nrow(view)
  if (!ncol(view)) {
    return(list(
      A = .rapid_empty_relation(n),
      metadata = list(
        source = source, k = 0L, radius = radius, sigma = NA_real_,
        weight_mode = weight_mode, symmetrize = symmetrize,
        valid_rows = 0L, missing_rows = as.integer(seq_len(n)),
        degree_cap = as.integer(degree_cap), nnz_directed = 0L
      )
    ))
  }

  degree_cap <- .rapid_int_scalar(degree_cap, "degree_cap")
  k <- .rapid_int_scalar(k, "k")
  if (!is.null(radius)) {
    radius <- .rapid_number_scalar(radius, "radius", lower = 0, strict = TRUE)
  }
  if (!is.null(sigma)) {
    sigma <- .rapid_number_scalar(sigma, "sigma", lower = 0, strict = TRUE)
  }

  valid <- rowSums(!is.finite(view)) == 0L
  valid_idx <- which(valid)
  nv <- length(valid_idx)
  if (nv <= 1L) {
    return(list(
      A = .rapid_empty_relation(n),
      metadata = list(
        source = source, k = 0L, radius = radius,
        sigma = if (is.null(sigma)) NA_real_ else sigma,
        weight_mode = weight_mode, symmetrize = symmetrize,
        valid_rows = as.integer(nv), missing_rows = as.integer(which(!valid)),
        degree_cap = degree_cap, nnz_directed = 0L
      )
    ))
  }

  limit <- if (is.null(radius)) min(k, degree_cap, nv - 1L) else min(degree_cap, nv - 1L)
  request_k <- min(nv, limit + 1L)
  data_valid <- view[valid_idx, , drop = FALSE]
  nn <- RANN::nn2(data = data_valid, query = data_valid, k = request_k)

  # Flatten the bounded ANN result row-wise. This avoids an O(n) interpreter
  # loop while retaining O(n * degree_cap) memory. Lexicographic ranks of the
  # node views provide a permutation-stable tie break for distinct rows.
  query <- rep(seq_len(nv), each = request_k)
  candidate <- as.integer(c(t(nn$nn.idx)))
  distance <- as.numeric(c(t(nn$nn.dists)))
  keep <- candidate != query & candidate >= 1L & candidate <= nv &
    is.finite(distance)
  if (!is.null(radius)) keep <- keep & distance <= radius
  query <- query[keep]
  candidate <- candidate[keep]
  distance <- distance[keep]

  if (length(query)) {
    tie_args <- lapply(seq_len(ncol(data_valid)), function(j) data_valid[, j])
    tie_args[[length(tie_args) + 1L]] <- seq_len(nv)
    tie_args$method <- "radix"
    tie_order <- do.call(order, tie_args)
    tie_rank <- integer(nv)
    tie_rank[tie_order] <- seq_len(nv)

    ord <- order(query, distance, tie_rank[candidate], candidate, method = "radix")
    query <- query[ord]
    candidate <- candidate[ord]
    distance <- distance[ord]
    within_query <- sequence(rle(query)$lengths)
    keep <- within_query <= limit
    query <- query[keep]
    candidate <- candidate[keep]
    distance <- distance[keep]
  }

  ii <- as.integer(valid_idx[query])
  jj <- as.integer(valid_idx[candidate])
  dd <- as.numeric(distance)
  if (!length(ii)) {
    A <- .rapid_empty_relation(n)
    sigma_use <- if (is.null(sigma)) NA_real_ else sigma
  } else {
    sigma_use <- sigma
    if (is.null(sigma_use)) {
      positive <- dd[is.finite(dd) & dd > 1e-12]
      sigma_use <- if (length(positive)) stats::median(positive) else 1
    }
    weights <- if (identical(weight_mode, "binary")) {
      rep(1, length(dd))
    } else {
      exp(pmax(-0.5 * (dd / sigma_use)^2, -700))
    }
    A <- Matrix::sparseMatrix(i = ii, j = jj, x = weights, dims = c(n, n))
    A <- .rapid_force_dgC(Matrix::drop0(A))
  }
  nnz_directed <- Matrix::nnzero(A)
  if (identical(symmetrize, "mutual")) A <- .rapid_mutual_relation(A)
  A <- .rapid_force_dgC(A)

  list(
    A = A,
    metadata = list(
      source = source,
      k = if (is.null(radius)) as.integer(limit) else NA_integer_,
      radius = radius,
      sigma = sigma_use,
      weight_mode = weight_mode,
      symmetrize = symmetrize,
      valid_rows = as.integer(nv),
      missing_rows = as.integer(which(!valid)),
      degree_cap = degree_cap,
      nnz_directed = as.integer(nnz_directed),
      deterministic = TRUE
    )
  )
}

.rapid_edges_to_sparse <- function(edges, n, name) {
  needed <- c("from", "to")
  missing <- setdiff(needed, names(edges))
  if (length(missing)) {
    stop("Edge relation `", name, "` is missing columns: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  from <- edges$from
  to <- edges$to
  if (!is.numeric(from) || !is.numeric(to) ||
      any(!is.finite(from)) || any(!is.finite(to)) ||
      any(abs(from - round(from)) > 1e-8) ||
      any(abs(to - round(to)) > 1e-8)) {
    stop("Edge indices for relation `", name, "` must be finite integers.",
         call. = FALSE)
  }
  from <- as.integer(from)
  to <- as.integer(to)
  if (any(from < 1L | from > n | to < 1L | to > n)) {
    stop("Edge indices for relation `", name, "` must lie in 1..", n, ".",
         call. = FALSE)
  }
  weight <- if ("weight" %in% names(edges)) edges$weight else rep(1, length(from))
  if (!is.numeric(weight) || any(!is.finite(weight)) || any(weight < 0)) {
    stop("Edge weights for relation `", name,
         "` must be finite and nonnegative.", call. = FALSE)
  }
  Matrix::sparseMatrix(
    i = from,
    j = to,
    x = as.numeric(weight),
    dims = c(n, n)
  )
}

.rapid_as_sparse_relation <- function(x,
                                      n,
                                      name,
                                      dense_max_n = 2000L,
                                      degree_cap = 64L,
                                      symmetrize = c("preserve", "mutual", "none")) {
  symmetrize <- match.arg(symmetrize)
  dense_max_n <- .rapid_int_scalar(dense_max_n, "dense_max_n")
  degree_cap <- .rapid_int_scalar(degree_cap, "degree_cap")
  input_type <- NULL
  input_was_dense <- FALSE

  if (is.data.frame(x)) {
    input_type <- "edge_data_frame"
    A <- .rapid_edges_to_sparse(x, n = n, name = name)
  } else if (inherits(x, "Matrix")) {
    input_type <- "sparse_matrix"
    A <- x
  } else if (is.matrix(x)) {
    input_type <- "dense_matrix"
    input_was_dense <- TRUE
    if (n > dense_max_n) {
      stop(
        "Dense relation `", name, "` has ", n, " nodes, exceeding ",
        "`dense_max_n = ", dense_max_n, "`. Pass a sparse Matrix or edge data.frame.",
        call. = FALSE
      )
    }
    A <- Matrix::Matrix(x, sparse = TRUE)
  } else {
    stop("Relation `", name,
         "` must be a sparse Matrix, a small dense matrix, or an edge data.frame.",
         call. = FALSE)
  }

  if (length(dim(A)) != 2L || !identical(as.integer(dim(A)), c(n, n))) {
    stop("Relation `", name, "` must have dimensions ", n, " x ", n, ".",
         call. = FALSE)
  }
  A <- .rapid_force_dgC(A)
  if (any(!is.finite(A@x)) || any(A@x < 0)) {
    stop("Relation `", name, "` must contain finite nonnegative weights.",
         call. = FALSE)
  }
  if (n > 0L) Matrix::diag(A) <- 0
  A <- .rapid_force_dgC(Matrix::drop0(A))
  symmetric_input <- isTRUE(isSymmetric(A, tol = 1e-12))
  nnz_before_cap <- Matrix::nnzero(A)
  A <- .rapid_cap_out_degree(A, degree_cap = degree_cap)

  symmetrize_used <- symmetrize
  if (identical(symmetrize, "preserve")) {
    symmetrize_used <- if (symmetric_input) "mutual" else "none"
  }
  if (identical(symmetrize_used, "mutual")) A <- .rapid_mutual_relation(A)
  A <- .rapid_force_dgC(A)

  list(
    A = A,
    metadata = list(
      source = "custom",
      input_type = input_type,
      input_was_dense = input_was_dense,
      symmetric_input = symmetric_input,
      symmetrize = symmetrize_used,
      degree_cap = degree_cap,
      dropped_for_degree_cap = as.integer(nnz_before_cap - Matrix::nnzero(A))
    )
  )
}

.rapid_normalize_relation <- function(A,
                                      normalization = c("symmetric", "random_walk")) {
  normalization <- match.arg(normalization)
  A <- .rapid_force_dgC(A)
  out_strength <- as.numeric(Matrix::rowSums(A))
  in_strength <- as.numeric(Matrix::colSums(A))

  if (identical(normalization, "random_walk")) {
    left <- ifelse(out_strength > 0, 1 / out_strength, 0)
    normalized <- Matrix::Diagonal(x = left) %*% A
  } else {
    left <- ifelse(out_strength > 0, 1 / sqrt(out_strength), 0)
    right <- ifelse(in_strength > 0, 1 / sqrt(in_strength), 0)
    normalized <- Matrix::Diagonal(x = left) %*% A %*%
      Matrix::Diagonal(x = right)
  }
  normalized <- .rapid_force_dgC(Matrix::drop0(normalized))
  if (any(!is.finite(normalized@x)) || any(normalized@x < 0)) {
    stop("Relation normalization produced invalid values.", call. = FALSE)
  }

  counts <- .rapid_degree_counts(A)
  isolates <- which((counts$out_degree + counts$in_degree) == 0L)
  list(
    A = normalized,
    metadata = list(
      normalization = normalization,
      nnz = as.integer(Matrix::nnzero(normalized)),
      density = as.numeric(Matrix::nnzero(normalized) / max(nrow(A)^2, 1)),
      max_out_degree = if (length(counts$out_degree)) {
        as.integer(max(counts$out_degree))
      } else {
        0L
      },
      mean_out_degree = if (length(counts$out_degree)) {
        mean(counts$out_degree)
      } else {
        0
      },
      isolates = as.integer(isolates),
      n_isolates = as.integer(length(isolates)),
      directed = !isTRUE(isSymmetric(A, tol = 1e-12)),
      out_strength = out_strength,
      in_strength = in_strength,
      volume = as.numeric(sum(abs(normalized@x)))
    )
  )
}

.rapid_relation_input_list <- function(relations) {
  if (is.null(relations)) return(list())
  if (is.data.frame(relations) || is.matrix(relations) || inherits(relations, "Matrix")) {
    return(list(topology = relations))
  }
  if (!is.list(relations)) {
    stop("`relations` must be NULL, a relation object, or a named list.",
         call. = FALSE)
  }
  if (!length(relations)) return(relations)
  relation_names <- names(relations)
  if (is.null(relation_names) || anyNA(relation_names) || any(!nzchar(relation_names))) {
    stop("Every custom relation must have a non-empty name.", call. = FALSE)
  }
  if (anyDuplicated(relation_names)) {
    stop("Custom relation names must be unique.", call. = FALSE)
  }
  relations
}

.rapid_metadata_value <- function(x, name, default = NULL) {
  if (!is.null(x[[name]])) x[[name]] else default
}

# Prepare sparse named relations for one RAPID-MA domain.
#
# Validates custom topology/relations, optionally builds capped feature,
# spatial, and mixed-attribute neighbour graphs, and returns normalized sparse
# matrices plus relation metadata and diagnostics. `spatial_mode = "relative"`
# centers coordinates and divides by their RMS radius, making Euclidean graph
# construction invariant to translation, rotation, and uniform scale.
# `spatial_mode = "common"` preserves coordinate units; kNN topology remains
# translation/rotation invariant, while a fixed radius or sigma is scale
# dependent. Radius graphs are deliberately capped: they retain the nearest
# `degree_cap` nodes inside the radius.
.rapid_prepare_relations <- function(
  X,
  relations = NULL,
  positions = NULL,
  attributes = NULL,
  build_feature = is.null(relations),
  feature_k = 15L,
  feature_sketch_dim = 32L,
  spatial_mode = c("relative", "common", "shared", "ignore"),
  spatial_k = 8L,
  spatial_radius = NULL,
  spatial_sigma = NULL,
  attribute_mode = c("relation", "shared", "ignore"),
  attribute_k = 15L,
  attribute_hash_dim = 32L,
  normalization = c("symmetric", "random_walk"),
  weight_mode = c("heat", "binary"),
  degree_cap = 64L,
  custom_symmetrize = c("preserve", "mutual", "none"),
  dense_max_n = 2000L,
  seed = 1L
) {
  .rapid_validate_numeric_matrix(X, "X")
  xvals <- if (inherits(X, "Matrix")) X@x else X
  if (any(!is.finite(xvals))) {
    stop("`X` must contain only finite values.", call. = FALSE)
  }
  n <- nrow(X)
  normalization <- match.arg(normalization)
  weight_mode <- match.arg(weight_mode)
  spatial_mode <- match.arg(spatial_mode)
  if (identical(spatial_mode, "shared")) spatial_mode <- "common"
  attribute_mode <- match.arg(attribute_mode)
  custom_symmetrize <- match.arg(custom_symmetrize)
  if (length(build_feature) != 1L || !is.logical(build_feature) || is.na(build_feature)) {
    stop("`build_feature` must be TRUE or FALSE.", call. = FALSE)
  }
  degree_cap <- .rapid_int_scalar(degree_cap, "degree_cap")
  if (degree_cap > 4096L) {
    stop("`degree_cap` must be <= 4096 to preserve bounded ANN memory.",
         call. = FALSE)
  }
  dense_max_n <- .rapid_int_scalar(dense_max_n, "dense_max_n")
  seed <- .rapid_int_scalar(seed, "seed", min_value = 0L)
  feature_k <- .rapid_int_scalar(feature_k, "feature_k")
  spatial_k <- .rapid_int_scalar(spatial_k, "spatial_k")
  attribute_k <- .rapid_int_scalar(attribute_k, "attribute_k")

  raw <- list()
  source_metadata <- list()
  skipped <- character(0)
  views <- list(feature = NULL, position = NULL, attribute = NULL)

  custom <- .rapid_relation_input_list(relations)
  if (length(custom)) {
    for (name in names(custom)) {
      prepared <- .rapid_as_sparse_relation(
        custom[[name]],
        n = n,
        name = name,
        dense_max_n = dense_max_n,
        degree_cap = degree_cap,
        symmetrize = custom_symmetrize
      )
      raw[[name]] <- prepared$A
      source_metadata[[name]] <- prepared$metadata
    }
  }

  if (isTRUE(build_feature)) {
    feature <- .rapid_feature_view(
      X,
      max_dim = feature_sketch_dim,
      seed = seed
    )
    views$feature <- feature
    if ("feature" %in% names(raw)) {
      skipped <- c(skipped, "feature:custom relation supplied")
    } else {
      built <- .rapid_build_neighbor_relation(
        feature$view,
        k = feature_k,
        degree_cap = degree_cap,
        weight_mode = weight_mode,
        symmetrize = "mutual",
        source = "feature"
      )
      raw$feature <- built$A
      source_metadata$feature <- utils::modifyList(
        built$metadata,
        list(view = feature$metadata)
      )
    }
  }

  if (!is.null(positions)) {
    position <- .rapid_prepare_positions(positions, n = n, mode = spatial_mode)
    if (!identical(spatial_mode, "ignore")) {
      scale_invariant <- identical(spatial_mode, "relative") ||
        (is.null(spatial_radius) && is.null(spatial_sigma))
      position$metadata$scale_invariant <- scale_invariant
    }
    views$position <- position
    if (!identical(spatial_mode, "ignore")) {
      if ("spatial" %in% names(raw)) {
        skipped <- c(skipped, "spatial:custom relation supplied")
      } else {
        built <- .rapid_build_neighbor_relation(
          position$relation_view,
          k = spatial_k,
          degree_cap = degree_cap,
          radius = spatial_radius,
          weight_mode = weight_mode,
          sigma = spatial_sigma,
          symmetrize = "mutual",
          source = "spatial"
        )
        raw$spatial <- built$A
        source_metadata$spatial <- utils::modifyList(
          built$metadata,
          list(coordinates = position$metadata)
        )
      }
    }
  }

  if (!is.null(attributes)) {
    .rapid_validate_attributes(attributes, n)
    if (!identical(attribute_mode, "ignore")) {
      attribute <- .rapid_encode_attributes(
        attributes,
        n = n,
        hash_dim = attribute_hash_dim,
        seed = seed
      )
      views$attribute <- attribute
      if (identical(attribute_mode, "relation")) {
        if ("attribute" %in% names(raw)) {
          skipped <- c(skipped, "attribute:custom relation supplied")
        } else if (!ncol(attribute$view)) {
          skipped <- c(skipped, "attribute:no informative encoded columns")
        } else {
          built <- .rapid_build_neighbor_relation(
            attribute$view,
            k = attribute_k,
            degree_cap = degree_cap,
            weight_mode = weight_mode,
            symmetrize = "mutual",
            source = "attribute"
          )
          raw$attribute <- built$A
          source_metadata$attribute <- utils::modifyList(
            built$metadata,
            list(encoding = attribute$metadata)
          )
        }
      }
    } else {
      views$attribute <- list(
        view = NULL,
        metadata = list(ignored = TRUE, n_columns = ncol(attributes))
      )
    }
  }

  normalized <- vector("list", length(raw))
  metadata <- vector("list", length(raw))
  names(normalized) <- names(metadata) <- names(raw)
  if (length(raw)) {
    for (name in names(raw)) {
      norm <- .rapid_normalize_relation(raw[[name]], normalization = normalization)
      normalized[[name]] <- norm$A
      metadata[[name]] <- utils::modifyList(source_metadata[[name]], norm$metadata)
    }
  }

  if (length(metadata)) {
    summary_rows <- lapply(names(metadata), function(name) {
      meta <- metadata[[name]]
      data.frame(
        relation = name,
        source = as.character(.rapid_metadata_value(meta, "source", "custom")),
        nnz = as.integer(.rapid_metadata_value(meta, "nnz", 0L)),
        max_out_degree = as.integer(.rapid_metadata_value(meta, "max_out_degree", 0L)),
        n_isolates = as.integer(.rapid_metadata_value(meta, "n_isolates", n)),
        directed = isTRUE(.rapid_metadata_value(meta, "directed", FALSE)),
        normalization = normalization,
        stringsAsFactors = FALSE
      )
    })
    summary <- do.call(rbind, summary_rows)
    rownames(summary) <- NULL
  } else {
    summary <- data.frame(
      relation = character(0), source = character(0), nnz = integer(0),
      max_out_degree = integer(0), n_isolates = integer(0),
      directed = logical(0), normalization = character(0),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      relations = normalized,
      metadata = metadata,
      diagnostics = list(
        summary = summary,
        relation_names = names(normalized),
        n_relations = as.integer(length(normalized)),
        total_nnz = as.integer(sum(vapply(normalized, Matrix::nnzero, numeric(1)))),
        max_out_degree = if (nrow(summary)) max(summary$max_out_degree) else 0L,
        skipped_builders = skipped,
        seed = seed
      ),
      views = views,
      n_nodes = as.integer(n),
      normalization = normalization
    ),
    class = "rapid_ma_relations"
  )
}
