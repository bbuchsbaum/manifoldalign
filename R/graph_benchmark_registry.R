.graph_benchmark_base_coords <- function(n, d, structure) {
  if (identical(structure, "ring")) {
    angles <- 2 * pi * seq_len(n) / n
    X <- cbind(cos(angles), sin(angles))
    if (d > 2L) {
      X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
    } else if (d == 1L) {
      X <- matrix(cos(angles), n, 1L)
    }
    return(X)
  }

  if (identical(structure, "grid")) {
    grid_size <- ceiling(sqrt(n))
    grid_coords <- expand.grid(
      x = seq_len(grid_size) / grid_size,
      y = seq_len(grid_size) / grid_size
    )[seq_len(n), ]
    X <- cbind(grid_coords$x, grid_coords$y)
    if (d > 2L) {
      X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
    } else if (d == 1L) {
      X <- matrix(grid_coords$x, n, 1L)
    }
    return(as.matrix(X))
  }

  if (identical(structure, "community")) {
    k <- 4L
    X <- matrix(0, n, max(2L, d))
    n_per <- ceiling(n / k)
    for (g in seq_len(k)) {
      idx <- seq.int((g - 1L) * n_per + 1L, min(g * n_per, n))
      center <- c(cos(2 * pi * g / k), sin(2 * pi * g / k))
      X[idx, 1:2] <- matrix(stats::rnorm(length(idx) * 2, sd = 0.2), length(idx), 2) +
        matrix(center, length(idx), 2, byrow = TRUE)
    }
    if (d > 2L) {
      X[, 3:d] <- matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L)
    }
    return(X[, seq_len(d), drop = FALSE])
  }

  matrix(stats::rnorm(n * d), n, d)
}

synthetic_graph_alignment_registry <- function(sizes = c(50L, 100L),
                                               structures = c("ring", "grid", "random", "community"),
                                               d = 3L,
                                               noise_sd = 0.05,
                                               permute_fraction = 1,
                                               n_anchors = 0L,
                                               n_reps = 3L,
                                               seed = 1L) {
  structures <- match.arg(structures, c("ring", "grid", "random", "community"), several.ok = TRUE)
  sizes <- as.integer(sizes)
  sizes <- sizes[is.finite(sizes) & sizes >= 3L]
  if (!length(sizes)) {
    stop("`sizes` must contain integers >= 3.", call. = FALSE)
  }

  rows <- list()
  ptr <- 1L
  for (structure in structures) {
    for (s in seq_along(sizes)) {
      n <- sizes[[s]]
      for (r in seq_len(as.integer(n_reps))) {
        generation_seed <- as.integer(seed) + match(structure, structures) * 10007L + s * 1009L + r * 17L
        rows[[ptr]] <- data.frame(
          scenario_family = "graph",
          scenario = paste(structure, n, r, sep = "_"),
          structure = structure,
          n = n,
          d = as.integer(d),
          noise_sd = as.numeric(noise_sd),
          permute_fraction = as.numeric(permute_fraction),
          n_anchors = as.integer(n_anchors),
          rep = as.integer(r),
          generation_seed = generation_seed,
          stringsAsFactors = FALSE
        )
        ptr <- ptr + 1L
      }
    }
  }

  do.call(rbind, rows)
}

synthetic_graph_alignment_case <- function(n,
                                           d = 3L,
                                           structure = c("ring", "grid", "random", "community"),
                                           noise_sd = 0.05,
                                           permute_fraction = 1,
                                           n_anchors = 0L,
                                           seed = 1L) {
  structure <- match.arg(structure)
  set.seed(as.integer(seed))

  X1 <- .graph_benchmark_base_coords(n = as.integer(n), d = as.integer(d), structure = structure)

  n_permute <- floor(as.integer(n) * as.numeric(permute_fraction))
  perm_idx <- if (n_permute >= 2L) sort(sample.int(as.integer(n), n_permute)) else integer(0)
  perm21 <- seq_len(as.integer(n))
  if (length(perm_idx)) {
    perm21[perm_idx] <- perm_idx[sample.int(length(perm_idx))]
  }

  X2 <- X1[perm21, , drop = FALSE] +
    matrix(stats::rnorm(as.integer(n) * as.integer(d), sd = noise_sd), as.integer(n), as.integer(d))
  map12 <- match(seq_len(as.integer(n)), perm21)

  a1 <- rep(NA_integer_, as.integer(n))
  a2 <- rep(NA_integer_, as.integer(n))
  anchor_idx <- integer(0)
  if (as.integer(n_anchors) > 0L) {
    available <- setdiff(seq_len(as.integer(n)), perm_idx)
    if (length(available) < as.integer(n_anchors)) {
      available <- seq_len(as.integer(n))
    }
    anchor_idx <- sort(sample(available, min(as.integer(n_anchors), length(available))))
    a1[anchor_idx] <- seq_along(anchor_idx)
    a2[map12[anchor_idx]] <- seq_along(anchor_idx)
  }

  hd <- as_hyperdesign(
    X_list = list(domain1 = X1, domain2 = X2),
    labels = list(a1, a2),
    label_name = "anchors"
  )

  correspondences <- NULL
  if (length(anchor_idx)) {
    correspondences <- data.frame(
      domain_i = rep(1L, length(anchor_idx)),
      index_i = anchor_idx,
      domain_j = rep(2L, length(anchor_idx)),
      index_j = map12[anchor_idx],
      stringsAsFactors = FALSE
    )
  }

  list(
    X1 = X1,
    X2 = X2,
    perm21 = perm21,
    map12 = map12,
    permuted_index = perm_idx,
    anchor_idx = anchor_idx,
    anchors1 = a1,
    anchors2 = a2,
    correspondences = correspondences,
    hd = hd,
    anchors_name = "anchors",
    benchmark_meta = list(
      scenario_family = "graph",
      structure = structure,
      n = as.integer(n),
      d = as.integer(d),
      noise_sd = as.numeric(noise_sd),
      permute_fraction = as.numeric(permute_fraction),
      n_anchors = as.integer(n_anchors),
      seed = as.integer(seed)
    )
  )
}
