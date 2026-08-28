# benchmark_global_geo_vs_baselines.R
# ------------------------------------------------------------
# Compare global_geo_align against multiple existing algorithms
# on the same multi-domain input and correspondence structure.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION") &&
      grepl("^Package:\\s*manifoldalign", paste(readLines("DESCRIPTION", n = 5L), collapse = "\n"))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(manifoldalign)
  library(multidesign)
  library(multivarious)
  library(adjoin)
  library(dplyr)
})

build_benchmark_hyperdesign <- function(raw_domains = NULL) {
  if (is.null(raw_domains)) {
    raw_domains <- manifoldalign::alignment_benchmark$domains
  }
  domains <- lapply(raw_domains, function(dom) {
    multidesign::multidesign(dom$x, dom$design)
  })
  hd <- multidesign::hyperdesign(domains)
  list(hd = hd, raw = raw_domains)
}

deep_copy_domains <- function(domains) {
  out <- lapply(domains, function(dom) {
    list(
      x = as.matrix(dom$x),
      design = as.data.frame(dom$design)
    )
  })
  names(out) <- names(domains)
  out
}

ensure_eval_columns <- function(domains) {
  out <- deep_copy_domains(domains)
  for (i in seq_along(out)) {
    di <- out[[i]]$design
    if (!("eval_id" %in% names(di)) && ("id" %in% names(di))) {
      di$eval_id <- di$id
    }
    if (!("condition_true" %in% names(di)) && ("condition" %in% names(di))) {
      di$condition_true <- di$condition
    }
    out[[i]]$design <- di
  }
  out
}

apply_supervision_budget <- function(domains, id_keep_frac = 1, condition_keep_frac = 1, seed = 1L) {
  set.seed(seed)
  out <- ensure_eval_columns(domains)

  id_pool <- unlist(lapply(out, function(dom) as.character(dom$design$eval_id)), use.names = FALSE)
  id_counts <- table(id_pool)
  # Keep supervision candidates that appear in at least two domains (pairwise overlap),
  # rather than requiring all-domain overlap.
  ids_shared <- names(id_counts[id_counts >= 2L])
  if (!length(ids_shared)) return(out)

  id_keep_n <- max(2L, min(length(ids_shared), floor(id_keep_frac * length(ids_shared))))
  cond_keep_n <- max(2L, min(length(ids_shared), floor(condition_keep_frac * length(ids_shared))))
  id_keep <- if (id_keep_n >= length(ids_shared)) ids_shared else sample(ids_shared, size = id_keep_n, replace = FALSE)
  cond_keep <- if (cond_keep_n >= length(ids_shared)) ids_shared else sample(ids_shared, size = cond_keep_n, replace = FALSE)

  for (i in seq_along(out)) {
    di <- out[[i]]$design
    eval_id <- as.character(di$eval_id)

    keep_id <- eval_id %in% id_keep
    id_train <- eval_id
    id_train[!keep_id] <- sprintf("private_id_%d_%s", i, eval_id[!keep_id])
    di$id <- id_train

    cond_true <- as.character(di$condition_true)
    keep_cond <- eval_id %in% cond_keep
    cond_train <- cond_true
    cond_train[!keep_cond] <- sprintf("private_cond_%d_%d", i, which(!keep_cond))
    di$condition_train <- factor(cond_train)

    out[[i]]$design <- di
  }

  out
}

mask_id_labels <- function(domains, keep_frac = 0.3, seed = 1L) {
  set.seed(seed)
  out <- ensure_eval_columns(domains)
  for (i in seq_along(out)) {
    n <- nrow(out[[i]]$x)
    keep_n <- max(2L, min(n, floor(keep_frac * n)))
    keep <- sample.int(n, size = keep_n, replace = FALSE)
    id_train <- as.character(out[[i]]$design$eval_id)
    drop_idx <- setdiff(seq_len(n), keep)
    if (length(drop_idx)) {
      id_train[drop_idx] <- sprintf("private_%d_%s", i, id_train[drop_idx])
    }
    out[[i]]$design$id <- id_train
  }
  out
}

noisify_condition_labels <- function(domains, flip_frac = 0.15, seed = 1L) {
  set.seed(seed)
  out <- ensure_eval_columns(domains)
  for (i in seq_along(out)) {
    cond <- as.character(out[[i]]$design$condition_true)
    n <- length(cond)
    flip_n <- floor(flip_frac * n)
    if (flip_n < 1L) next
    idx <- sample.int(n, size = flip_n, replace = FALSE)
    lv <- unique(cond)
    if (length(lv) < 2L) next
    cond[idx] <- vapply(cond[idx], function(cj) sample(lv[lv != cj], 1L), character(1))
    out[[i]]$design$condition <- factor(cond, levels = lv)
  }
  out
}

add_feature_shift_noise <- function(domains, noise_sd = 0.25, seed = 1L) {
  set.seed(seed)
  out <- ensure_eval_columns(domains)
  for (i in seq_along(out)) {
    X <- as.matrix(out[[i]]$x)
    d <- ncol(X)
    Q <- qr.Q(qr(matrix(rnorm(d * d), nrow = d)))
    scales <- runif(d, min = 0.75, max = 1.25)
    X_new <- X %*% Q %*% diag(scales, nrow = d) +
      matrix(rnorm(length(X), sd = noise_sd), nrow = nrow(X), ncol = d)
    out[[i]]$x <- scale(X_new, center = TRUE, scale = FALSE)
  }
  out
}

apply_rep_perturbation <- function(domains, noise_sd = 0.01, seed = 1L) {
  set.seed(seed)
  out <- deep_copy_domains(domains)
  for (i in seq_along(out)) {
    Xi <- as.matrix(out[[i]]$x)
    base_sd <- stats::sd(as.numeric(Xi))
    if (!is.finite(base_sd) || base_sd <= 0) base_sd <- 1
    Xi_pert <- Xi + matrix(rnorm(length(Xi), sd = noise_sd * base_sd), nrow = nrow(Xi), ncol = ncol(Xi))
    out[[i]]$x <- Xi_pert
  }
  out
}

build_holdout_pairs <- function(domains, id_train_col = "id", eval_col = "eval_id") {
  m <- length(domains)
  out <- vector("list", m * (m - 1L) / 2L)
  ptr <- 0L
  for (i in seq_len(m - 1L)) {
    eval_i <- as.character(domains[[i]]$design[[eval_col]])
    train_i <- as.character(domains[[i]]$design[[id_train_col]])
    for (j in (i + 1L):m) {
      eval_j <- as.character(domains[[j]]$design[[eval_col]])
      train_j <- as.character(domains[[j]]$design[[id_train_col]])

      eval_shared <- intersect(eval_i, eval_j)
      train_shared <- intersect(intersect(eval_shared, train_i), train_j)
      holdout_ids <- setdiff(eval_shared, train_shared)
      if (!length(holdout_ids)) next

      ptr <- ptr + 1L
      out[[ptr]] <- data.frame(
        domain_i = i,
        index_i = match(holdout_ids, eval_i),
        domain_j = j,
        index_j = match(holdout_ids, eval_j),
        stringsAsFactors = FALSE
      )
    }
  }

  if (ptr == 0L) {
    return(data.frame(
      domain_i = integer(0),
      index_i = integer(0),
      domain_j = integer(0),
      index_j = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out[seq_len(ptr)])
}

append_note <- function(note, addendum) {
  if (is.na(note) || !nzchar(note)) return(addendum)
  paste(note, addendum, sep = ";")
}

benchmark_uot_params <- function(scenario_name, hks = FALSE, fallback = FALSE) {
  eps <- 0.5
  rho1 <- 1.0
  rho2 <- 1.0
  lambda_feat <- if (hks) 0.25 else 0

  if (identical(scenario_name, "shifted_noisy")) {
    eps <- 1.0
    rho1 <- 2.0
    rho2 <- 2.0
    if (hks) lambda_feat <- 0.40
  }

  if (isTRUE(fallback)) {
    eps <- eps * 1.8
    rho1 <- rho1 * 1.5
    rho2 <- rho2 * 1.5
    if (hks) lambda_feat <- lambda_feat * 0.8
  }

  list(epsilon = eps, rho1 = rho1, rho2 = rho2, lambda_feat = lambda_feat)
}

is_uot_unstable <- function(scores) {
  if (!is.matrix(scores)) scores <- as.matrix(scores)
  if (!all(is.finite(scores))) return(TRUE)
  max_abs <- max(abs(scores))
  mean_energy <- mean(rowSums(scores^2))
  max_abs > 40 || mean_energy > 200
}

uot_metric_guardrails <- function(scenario_name) {
  if (identical(scenario_name, "shifted_noisy")) {
    return(list(identity_mse_max = 2.0, holdout_identity_mse_max = 2.0))
  }
  if (identical(scenario_name, "sparse_supervision")) {
    return(list(identity_mse_max = 1.0, holdout_identity_mse_max = 1.0))
  }
  list(identity_mse_max = 1.0, holdout_identity_mse_max = 1.0)
}

is_uot_metric_unstable <- function(metrics, scenario_name) {
  if (is.null(metrics) || !length(metrics)) return(FALSE)
  limits <- uot_metric_guardrails(scenario_name)
  bad_id <- is.finite(metrics$identity_mse) && metrics$identity_mse > limits$identity_mse_max
  bad_holdout <- is.finite(metrics$holdout_identity_mse) && metrics$holdout_identity_mse > limits$holdout_identity_mse_max
  bad_id || bad_holdout
}

subsample_domain_overlap <- function(domains, keep_fracs = c(1, 0.75, 0.65), seed = 1L) {
  set.seed(seed)
  out <- ensure_eval_columns(domains)
  if (length(keep_fracs) != length(out)) {
    keep_fracs <- rep(min(keep_fracs), length.out = length(out))
  }
  for (i in seq_along(out)) {
    n <- nrow(out[[i]]$x)
    keep_n <- max(20L, min(n, floor(keep_fracs[i] * n)))
    keep <- sort(sample.int(n, size = keep_n, replace = FALSE))
    out[[i]]$x <- out[[i]]$x[keep, , drop = FALSE]
    out[[i]]$design <- out[[i]]$design[keep, , drop = FALSE]
  }
  out
}

build_benchmark_scenarios <- function(base_domains, seed = 1L) {
  make_gaussian_class_domains <- function(
    seed = 1L,
    n_domains = 3L,
    n_classes = 4L,
    n_per_class = 40L,
    d = 12L,
    class_sep = 2.0,
    domain_noise_sd = 0.12
  ) {
    set.seed(seed)
    n <- as.integer(n_classes * n_per_class)
    class_levels <- sprintf("class_%02d", seq_len(n_classes))
    class_idx <- rep(seq_len(n_classes), each = n_per_class)
    class_labels <- class_levels[class_idx]
    eval_id <- sprintf("sample_%04d", seq_len(n))

    means <- matrix(stats::rnorm(n_classes * d), nrow = n_classes, ncol = d)
    means <- scale(means, center = TRUE, scale = FALSE)
    means <- means * class_sep

    latent <- matrix(0, nrow = n, ncol = d)
    for (cc in seq_len(n_classes)) {
      idx <- which(class_idx == cc)
      A <- matrix(stats::rnorm(d * d), nrow = d, ncol = d)
      Q <- qr.Q(qr(A))
      vals <- stats::runif(d, min = 0.25, max = 1.0)
      cov_half <- Q %*% diag(sqrt(vals), nrow = d)
      Z <- matrix(stats::rnorm(length(idx) * d), nrow = length(idx), ncol = d)
      latent[idx, ] <- Z %*% cov_half + matrix(means[cc, ], nrow = length(idx), ncol = d, byrow = TRUE)
    }
    latent <- scale(latent, center = TRUE, scale = FALSE)

    domains <- vector("list", n_domains)
    for (dd in seq_len(n_domains)) {
      R <- qr.Q(qr(matrix(stats::rnorm(d * d), nrow = d, ncol = d)))
      scales <- stats::runif(d, min = 0.8, max = 1.2)
      shift <- stats::rnorm(d, sd = 0.10)
      Xd <- latent %*% R %*% diag(scales, nrow = d) +
        matrix(shift, nrow = n, ncol = d, byrow = TRUE) +
        matrix(stats::rnorm(n * d, sd = domain_noise_sd), nrow = n, ncol = d)
      Xd <- scale(Xd, center = TRUE, scale = FALSE)
      design <- data.frame(
        id = eval_id,
        eval_id = eval_id,
        condition = factor(class_labels, levels = class_levels),
        condition_true = factor(class_labels, levels = class_levels),
        stringsAsFactors = FALSE
      )
      domains[[dd]] <- list(x = Xd, design = design)
      names(domains)[dd] <- sprintf("domain%d", dd)
    }
    domains
  }

  make_open_set_gaussian_domains <- function(
    seed = 1L,
    n_per_class = 32L,
    d = 12L,
    class_sep = 2.0,
    domain_noise_sd = 0.14
  ) {
    set.seed(seed)

    domain_class_map <- list(
      domain1 = c("class_A", "class_B", "class_C", "class_D"),
      domain2 = c("class_B", "class_C", "class_D", "class_E"),
      domain3 = c("class_C", "class_D", "class_E", "class_F")
    )

    class_levels <- unique(unlist(domain_class_map))
    n_classes <- length(class_levels)

    class_means <- matrix(stats::rnorm(n_classes * d), nrow = n_classes, ncol = d)
    class_means <- scale(class_means, center = TRUE, scale = FALSE) * class_sep
    rownames(class_means) <- class_levels

    class_latent <- setNames(vector("list", length(class_levels)), class_levels)
    for (cl in class_levels) {
      A <- matrix(stats::rnorm(d * d), nrow = d, ncol = d)
      Q <- qr.Q(qr(A))
      vals <- stats::runif(d, min = 0.20, max = 1.10)
      cov_half <- Q %*% diag(sqrt(vals), nrow = d)
      Z <- matrix(stats::rnorm(n_per_class * d), nrow = n_per_class, ncol = d)
      class_latent[[cl]] <- Z %*% cov_half + matrix(class_means[cl, ], nrow = n_per_class, ncol = d, byrow = TRUE)
    }

    class_counts <- table(unlist(domain_class_map))
    domains <- vector("list", length(domain_class_map))
    names(domains) <- names(domain_class_map)

    for (dd in seq_along(domain_class_map)) {
      domain_name <- names(domain_class_map)[dd]
      domain_classes <- domain_class_map[[dd]]

      R <- qr.Q(qr(matrix(stats::rnorm(d * d), nrow = d, ncol = d)))
      scales <- stats::runif(d, min = 0.80, max = 1.25)
      shift <- stats::rnorm(d, sd = 0.12)

      X_blocks <- vector("list", length(domain_classes))
      id_blocks <- vector("list", length(domain_classes))
      cond_blocks <- vector("list", length(domain_classes))

      for (ii in seq_along(domain_classes)) {
        cl <- domain_classes[[ii]]
        base <- class_latent[[cl]]
        Xi <- base %*% R %*% diag(scales, nrow = d) +
          matrix(shift, nrow = nrow(base), ncol = d, byrow = TRUE) +
          matrix(stats::rnorm(length(base), sd = domain_noise_sd), nrow = nrow(base), ncol = d)
        X_blocks[[ii]] <- Xi

        if (class_counts[[cl]] >= 2L) {
          id_blocks[[ii]] <- sprintf("shared_%s_%03d", cl, seq_len(n_per_class))
        } else {
          id_blocks[[ii]] <- sprintf("private_%s_%s_%03d", domain_name, cl, seq_len(n_per_class))
        }
        cond_blocks[[ii]] <- rep(cl, n_per_class)
      }

      Xd <- do.call(rbind, X_blocks)
      eval_id <- unlist(id_blocks, use.names = FALSE)
      cond <- unlist(cond_blocks, use.names = FALSE)

      ord <- sample.int(nrow(Xd))
      Xd <- Xd[ord, , drop = FALSE]
      eval_id <- eval_id[ord]
      cond <- cond[ord]

      Xd <- scale(Xd, center = TRUE, scale = FALSE)

      domains[[dd]] <- list(
        x = Xd,
        design = data.frame(
          id = eval_id,
          eval_id = eval_id,
          condition = factor(cond, levels = class_levels),
          condition_true = factor(cond, levels = class_levels),
          stringsAsFactors = FALSE
        )
      )
    }

    domains
  }

  list(
    clean = list(
      name = "clean",
      domains = ensure_eval_columns(base_domains),
      id_keep_frac = 0.8,
      condition_keep_frac = 0.8,
      perturb_sd = 0.01
    ),
    sparse_supervision = list(
      name = "sparse_supervision",
      domains = noisify_condition_labels(
        mask_id_labels(base_domains, keep_frac = 0.30, seed = seed + 11L),
        flip_frac = 0.15,
        seed = seed + 12L
      ),
      id_keep_frac = 0.30,
      condition_keep_frac = 0.30,
      perturb_sd = 0.015
    ),
    shifted_noisy = list(
      name = "shifted_noisy",
      domains = add_feature_shift_noise(base_domains, noise_sd = 0.30, seed = seed + 21L),
      id_keep_frac = 0.8,
      condition_keep_frac = 0.8,
      perturb_sd = 0.02
    ),
    gaussian_classes = list(
      name = "gaussian_classes",
      domains = make_gaussian_class_domains(
        seed = seed + 31L,
        n_domains = 3L,
        n_classes = 4L,
        n_per_class = 40L,
        d = 12L,
        class_sep = 2.2,
        domain_noise_sd = 0.10
      ),
      id_keep_frac = 0.7,
      condition_keep_frac = 0.35,
      perturb_sd = 0.010
    ),
    gaussian_overlap = list(
      name = "gaussian_overlap",
      domains = make_gaussian_class_domains(
        seed = seed + 32L,
        n_domains = 3L,
        n_classes = 4L,
        n_per_class = 40L,
        d = 12L,
        class_sep = 1.2,
        domain_noise_sd = 0.18
      ),
      id_keep_frac = 0.7,
      condition_keep_frac = 0.25,
      perturb_sd = 0.015
    ),
    open_set_overlap = list(
      name = "open_set_overlap",
      domains = make_open_set_gaussian_domains(
        seed = seed + 33L,
        n_per_class = 32L,
        d = 12L,
        class_sep = 1.6,
        domain_noise_sd = 0.12
      ),
      id_keep_frac = 0.85,
      condition_keep_frac = 0.65,
      perturb_sd = 0.012
    )
  )
}

build_id_correspondences <- function(domains, id_col = "id") {
  m <- length(domains)
  out <- vector("list", m * (m - 1) / 2)
  ptr <- 0L

  for (i in seq_len(m - 1L)) {
    ids_i <- domains[[i]]$design[[id_col]]
    ids_i <- ids_i[!is.na(ids_i)]
    for (j in (i + 1L):m) {
      ids_j <- domains[[j]]$design[[id_col]]
      ids_j <- ids_j[!is.na(ids_j)]
      common <- intersect(ids_i, ids_j)
      if (!length(common)) next
      ptr <- ptr + 1L
      out[[ptr]] <- data.frame(
        domain_i = i,
        index_i = match(common, domains[[i]]$design[[id_col]]),
        domain_j = j,
        index_j = match(common, domains[[j]]$design[[id_col]])
      )
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
  do.call(rbind, out[seq_len(ptr)])
}

split_scores_by_domain <- function(scores, domain_sizes) {
  ends <- cumsum(domain_sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  lapply(seq_along(domain_sizes), function(i) {
    scores[starts[i]:ends[i], , drop = FALSE]
  })
}

identity_pair_mse <- function(score_blocks, raw_domains) {
  m <- length(score_blocks)
  vals <- numeric(0)
  for (i in seq_len(m - 1L)) {
    id_col_i <- if ("eval_id" %in% names(raw_domains[[i]]$design)) "eval_id" else "id"
    id_i <- raw_domains[[i]]$design[[id_col_i]]
    for (j in (i + 1L):m) {
      id_col_j <- if ("eval_id" %in% names(raw_domains[[j]]$design)) "eval_id" else "id"
      id_j <- raw_domains[[j]]$design[[id_col_j]]
      common <- intersect(id_i, id_j)
      if (!length(common)) next
      ii <- match(common, id_i)
      jj <- match(common, id_j)
      dif <- score_blocks[[i]][ii, , drop = FALSE] - score_blocks[[j]][jj, , drop = FALSE]
      vals <- c(vals, rowSums(dif * dif))
    }
  }
  if (!length(vals)) return(NA_real_)
  mean(vals)
}

holdout_identity_mse <- function(score_blocks, holdout_pairs) {
  if (is.null(holdout_pairs) || !nrow(holdout_pairs)) return(NA_real_)
  vals <- numeric(nrow(holdout_pairs))
  for (r in seq_len(nrow(holdout_pairs))) {
    di <- holdout_pairs$domain_i[r]
    dj <- holdout_pairs$domain_j[r]
    ii <- holdout_pairs$index_i[r]
    jj <- holdout_pairs$index_j[r]
    dif <- score_blocks[[di]][ii, , drop = FALSE] - score_blocks[[dj]][jj, , drop = FALSE]
    vals[r] <- sum(dif * dif)
  }
  mean(vals)
}

geometry_spearman <- function(score_blocks, raw_domains) {
  vals <- vapply(seq_along(score_blocks), function(i) {
    x_orig <- as.matrix(raw_domains[[i]]$x)
    x_emb <- as.matrix(score_blocks[[i]])
    d_orig <- as.matrix(dist(x_orig))
    d_emb <- as.matrix(dist(x_emb))
    ut <- upper.tri(d_orig)
    stats::cor(d_orig[ut], d_emb[ut], method = "spearman")
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

label_knn_accuracy <- function(scores, labels, k = 5L) {
  y <- as.character(labels)
  n <- nrow(scores)
  if (length(y) != n) stop("labels length mismatch")
  if (n <= k) return(NA_real_)

  D <- as.matrix(dist(scores))
  diag(D) <- Inf

  pred <- character(n)
  for (i in seq_len(n)) {
    nn <- head(order(D[i, ]), k)
    votes <- table(y[nn])
    pred[i] <- names(votes)[which.max(votes)]
  }
  mean(pred == y)
}

holdout_condition_knn_accuracy <- function(scores, raw_domains, k = 5L) {
  truth <- unlist(lapply(raw_domains, function(dom) {
    design <- dom$design
    if ("condition_true" %in% names(design)) {
      as.character(design$condition_true)
    } else if ("condition" %in% names(design)) {
      as.character(design$condition)
    } else {
      rep(NA_character_, nrow(dom$x))
    }
  }), use.names = FALSE)

  train_lbl <- unlist(lapply(raw_domains, function(dom) {
    design <- dom$design
    if ("condition_train" %in% names(design)) {
      as.character(design$condition_train)
    } else if ("condition" %in% names(design)) {
      as.character(design$condition)
    } else {
      rep(NA_character_, nrow(dom$x))
    }
  }), use.names = FALSE)

  truth_levels <- unique(stats::na.omit(truth))
  is_train <- !is.na(train_lbl) & (train_lbl %in% truth_levels)
  is_holdout <- !is.na(truth) & !is_train

  train_idx <- which(is_train)
  hold_idx <- which(is_holdout)
  if (!length(train_idx) || !length(hold_idx)) return(NA_real_)

  Xtrain <- as.matrix(scores[train_idx, , drop = FALSE])
  Xhold <- as.matrix(scores[hold_idx, , drop = FALSE])
  ytrain <- train_lbl[train_idx]
  yhold <- truth[hold_idx]

  k_eff <- min(as.integer(k), nrow(Xtrain))
  if (k_eff < 1L) return(NA_real_)

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = Xtrain, query = Xhold, k = k_eff)$nn.idx
    if (!is.matrix(nn)) nn <- matrix(nn, ncol = 1L)
    pred <- apply(nn, 1, function(idx) {
      votes <- table(ytrain[idx])
      names(votes)[which.max(votes)]
    })
    return(mean(pred == yhold))
  }

  D <- as.matrix(stats::dist(rbind(Xhold, Xtrain)))
  D <- D[seq_len(nrow(Xhold)), nrow(Xhold) + seq_len(nrow(Xtrain)), drop = FALSE]
  pred <- apply(D, 1, function(row_d) {
    idx <- head(order(row_d), k_eff)
    votes <- table(ytrain[idx])
    names(votes)[which.max(votes)]
  })
  mean(pred == yhold)
}

evaluate_fit <- function(fit, raw_domains, holdout_pairs = NULL) {
  scores <- extract_scores(fit)
  domain_sizes <- vapply(raw_domains, function(d) nrow(d$x), integer(1))
  score_blocks <- split_scores_by_domain(scores, domain_sizes)
  pooled_labels <- unlist(lapply(raw_domains, function(dom) {
    design <- dom$design
    label_col <- if ("condition_true" %in% names(design)) "condition_true" else "condition"
    as.character(design[[label_col]])
  }), use.names = FALSE)

  list(
    ncomp = ncol(scores),
    identity_mse = identity_pair_mse(score_blocks, raw_domains),
    holdout_identity_mse = holdout_identity_mse(score_blocks, holdout_pairs),
    geometry_spearman = geometry_spearman(score_blocks, raw_domains),
    condition_knn_acc = label_knn_accuracy(scores, pooled_labels, k = 5L),
    holdout_condition_knn_acc = holdout_condition_knn_accuracy(scores, raw_domains, k = 5L)
  )
}

extract_scores <- function(fit) {
  if (!is.null(fit$s)) {
    return(as.matrix(fit$s))
  }
  if (inherits(fit, "grasp_multiset") && !is.null(fit$embeddings)) {
    return(do.call(rbind, lapply(fit$embeddings, as.matrix)))
  }
  if (inherits(fit, "grasp_multiset") && !is.null(fit$embeddings_raw) && !is.null(fit$permutations)) {
    anchor_idx <- if (is.numeric(fit$anchor)) as.integer(fit$anchor) else 1L
    n_anchor <- nrow(fit$embeddings_raw[[anchor_idx]])
    aligned <- lapply(seq_along(fit$embeddings_raw), function(s) {
      emb <- as.matrix(fit$embeddings_raw[[s]])
      if (s == anchor_idx) return(emb)
      perm <- as.integer(fit$permutations[[s]])
      if (length(perm) != n_anchor) return(emb)
      out <- matrix(0, nrow = n_anchor, ncol = ncol(emb))
      valid <- perm > 0L & perm <= nrow(emb)
      if (any(valid)) {
        out[valid, ] <- emb[perm[valid], , drop = FALSE]
      }
      out
    })
    return(do.call(rbind, aligned))
  }
  stop("Benchmark cannot extract sample scores from class: ", paste(class(fit), collapse = ", "))
}

pca_stack_scores <- function(raw_domains, ncomp) {
  blocks <- lapply(raw_domains, function(dom) {
    X <- scale(as.matrix(dom$x), center = TRUE, scale = FALSE)
    k <- min(ncomp, ncol(X), nrow(X))
    if (k <= 0L) {
      return(matrix(0, nrow = nrow(X), ncol = ncomp))
    }
    scores <- stats::prcomp(X, center = FALSE, scale. = FALSE)$x[, seq_len(k), drop = FALSE]
    out <- matrix(0, nrow(scores), ncomp)
    out[, seq_len(k)] <- scores
    out
  })
  do.call(rbind, blocks)
}

run_one <- function(method,
                    supervision,
                    hd,
                    raw_domains,
                    correspondences,
                    holdout_pairs,
                    ncomp,
                    scenario_name,
                    method_seed) {
  if (!supervision %in% c("id", "condition", "none")) {
    stop("Unknown supervision mode: ", supervision)
  }
  set.seed(as.integer(method_seed))

  is_open_set <- identical(scenario_name, "open_set_overlap")

  t0 <- proc.time()[["elapsed"]]
  note <- NA_character_
  precomputed_metrics <- NULL
  fit <- switch(method,
    global_geo = {
      gg_knn <- if (is_open_set) 24L else 12L
      gg_landmarks <- if (is_open_set) 160L else 120L
      gg_max_pairs <- if (is_open_set) 120L else 80L
      if (supervision == "id") {
        manifoldalign::global_geo_align(
          hd,
          correspondences = correspondences,
          preproc = NULL,
          ncomp = ncomp,
          control = manifoldalign::global_geo_align_control(
            knn = gg_knn,
            n_landmarks_total = gg_landmarks,
            k_embed = max(3L * ncomp, ncomp + 6L),
            scale_method = "none",
            ridge_eps = 1e-4,
            verbose = FALSE,
            seed = as.integer(method_seed)
          )
        )
      } else if (supervision == "condition") {
        manifoldalign::global_geo_align(
          hd,
          y = condition_train,
          preproc = NULL,
          ncomp = ncomp,
          control = manifoldalign::global_geo_align_control(
            knn = gg_knn,
            n_landmarks_total = gg_landmarks,
            k_embed = max(3L * ncomp, ncomp + 6L),
            scale_method = "none",
            ridge_eps = 1e-4,
            max_pairs_per_label = gg_max_pairs,
            verbose = FALSE,
            seed = as.integer(method_seed)
          )
        )
      } else {
        stop("global_geo requires supervision mode 'id' or 'condition'")
      }
    },
    kema = {
      if (supervision == "id") {
        manifoldalign::kema(
          hd, y = id,
          preproc = multivarious::center(),
          ncomp = ncomp, knn = 5, solver = "exact",
          sample_frac = 0.8, lambda = 1e-3
        )
      } else if (supervision == "condition") {
        manifoldalign::kema(
          hd, y = condition_train,
          preproc = multivarious::center(),
          ncomp = ncomp, knn = 5, solver = "exact",
          sample_frac = 0.8, lambda = 1e-3
        )
      } else {
        stop("kema requires supervision mode 'id' or 'condition'")
      }
    },
    lra = {
      if (supervision == "id") {
        manifoldalign::lowrank_align(
          hd, y = id,
          ncomp = ncomp,
          simfun = adjoin::binary_label_matrix,
          solver = "operator", sv_thresh = 1, lambda = 1e-2
        )
      } else if (supervision == "condition") {
        manifoldalign::lowrank_align(
          hd, y = condition_train,
          ncomp = ncomp,
          simfun = adjoin::binary_label_matrix,
          solver = "operator", sv_thresh = 1, lambda = 1e-2
        )
      } else {
        stop("lra requires supervision mode 'id' or 'condition'")
      }
    },
    ssma = {
      ssma_ctrl <- manifoldalign::ssma_align_control(
        knn = 12L,
        rank_per_domain = 64L,
        solver = "reduced",
        max_pairs_per_label = 80L,
        verbose = FALSE,
        seed = as.integer(method_seed)
      )
      if (supervision == "id") {
        manifoldalign::ssma_align(
          hd,
          correspondences = correspondences,
          preproc = multivarious::center(),
          ncomp = ncomp,
          control = ssma_ctrl
        )
      } else if (supervision == "condition") {
        manifoldalign::ssma_align(
          hd,
          y = condition_train,
          preproc = multivarious::center(),
          ncomp = ncomp,
          control = ssma_ctrl
        )
      } else {
        stop("ssma requires supervision mode 'id' or 'condition'")
      }
    },
    gpca = {
      if (supervision == "id") {
        manifoldalign::gpca_align(
          hd, y = id,
          preproc = multivarious::center(),
          ncomp = ncomp, u = 0.5, lambda = 0.1,
          control = manifoldalign::gpca_align_control(knn = 10, verbose = FALSE)
        )
      } else if (supervision == "condition") {
        manifoldalign::gpca_align(
          hd, y = condition_train,
          preproc = multivarious::center(),
          ncomp = ncomp, u = 0.5, lambda = 0.1,
          control = manifoldalign::gpca_align_control(knn = 10, verbose = FALSE)
        )
      } else {
        stop("gpca requires supervision mode 'id' or 'condition'")
      }
    },
    spectral_mnn = {
      if (supervision == "id") {
        manifoldalign::spectral_mnn_align(
          hd,
          correspondences = correspondences,
          preproc = multivarious::center(),
          ncomp = ncomp,
          ref_idx = 1L,
          control = manifoldalign::spectral_mnn_align_control(
            knn = 12,
            q_descriptors = max(12L, 2L * ncomp),
            mnn_k = 10L,
            max_anchors = 512L,
            dist_quantile = 0.3,
            refine_rounds = 0L,
            force_SO = FALSE,
            seed = as.integer(method_seed),
            verbose = FALSE
          )
        )
      } else if (supervision == "condition") {
        manifoldalign::spectral_mnn_align(
          hd,
          y = condition_train,
          preproc = multivarious::center(),
          ncomp = ncomp,
          ref_idx = 1L,
          control = manifoldalign::spectral_mnn_align_control(
            knn = 12,
            q_descriptors = max(12L, 2L * ncomp),
            mnn_k = 10L,
            max_anchors = 512L,
            max_pairs_per_label = 80L,
            dist_quantile = 0.3,
            refine_rounds = 1L,
            force_SO = TRUE,
            seed = as.integer(method_seed),
            verbose = FALSE
          )
        )
      } else {
        manifoldalign::spectral_mnn_align(
          hd,
          preproc = multivarious::center(),
          ncomp = ncomp,
          ref_idx = 1L,
          control = manifoldalign::spectral_mnn_align_control(
            knn = 12,
            q_descriptors = max(12L, 2L * ncomp),
            mnn_k = 10L,
            max_anchors = 512L,
            dist_quantile = 0.3,
            refine_rounds = 1L,
            force_SO = FALSE,
            seed = as.integer(method_seed),
            verbose = FALSE
          )
        )
      }
    },
    cone_multi = {
      if (supervision != "none") stop("cone_multi uses supervision='none'")
      cone_start <- if (is_open_set) min(as.integer(ncomp), 7L) else as.integer(ncomp)
      cone_grid <- seq.int(from = cone_start, to = 2L, by = -1L)
      cone_fit <- NULL
      cone_err <- NULL
      for (kc in cone_grid) {
        fit_try <- tryCatch(
          manifoldalign::cone_align_multiple(
            hd,
            ref_idx = 1L,
            preproc = multivarious::center(),
            ncomp = kc,
            sigma = 0.73,
            lambda = 0.1,
            solver = "linear",
            max_iter = 30,
            tol = 0.01,
            knn = 10
          ),
          error = function(e) e
        )
        if (!inherits(fit_try, "error")) {
          cone_fit <- fit_try
          break
        }
        cone_err <- fit_try
      }
      if (is.null(cone_fit)) {
        stop(cone_err$message)
      }
      cone_fit
    },
    grasp_multi = {
      if (supervision != "none") stop("grasp_multi uses supervision='none'")
      manifoldalign::grasp_multiset(
        hd,
        preproc = multivarious::center(),
        ncomp = ncomp,
        q_descriptors = max(20L, 2L * ncomp),
        sigma = 0.73,
        lambda = 0.1,
        solver = "linear",
        max_iter = 60,
        tol = 1e-5
      )
    },
    uot = {
      if (supervision != "none") stop("uot uses supervision='none'")

      # Anchor/template = domain 1 points; map each domain into anchor score space.
      anchor <- raw_domains[[1]]$x
      d_raw <- min(ncomp, ncol(anchor))
      anchor_raw <- stats::prcomp(anchor, center = FALSE, scale. = FALSE)$x[, seq_len(d_raw), drop = FALSE]
      anchor_scores <- matrix(0, nrow(anchor_raw), ncomp)
      anchor_scores[, seq_len(d_raw)] <- anchor_raw

      # Uniform masses (balanced across nodes); UOT handles mass drift via rho.
      beta0 <- rep(1, nrow(anchor))
      params <- benchmark_uot_params(scenario_name, hks = FALSE, fallback = FALSE)
      eps <- params$epsilon
      rho1 <- params$rho1
      rho2 <- params$rho2

      logsumexp <- function(x) {
        mx <- max(x)
        if (!is.finite(mx)) return(mx)
        mx + log(sum(exp(x - mx)))
      }

      map_to_anchor <- function(Xk) {
        alpha0 <- rep(1, nrow(Xk))
        fit_k <- manifoldalign::uot_fit_pair(
          X = Xk,
          Y = anchor,
          alpha = alpha0,
          beta = beta0,
          lambda_anat = 1,
          lambda_feat = 0,
          neighbor_mode = "dense",
          epsilon = eps,
          rho1 = rho1,
          rho2 = rho2,
          max_iter = 2000,
          tol = 1e-10
        )

        C <- fit_k$cost
        fbar <- as.numeric(fit_k$sol$fbar)
        gbar <- as.numeric(fit_k$sol$gbar)
        logb <- rep(0, nrow(anchor)) # log(beta0) with beta0 uniform

        out <- matrix(0, nrow(Xk), ncomp)
        inv_eps <- 1 / eps
        for (i in seq_len(nrow(Xk))) {
          z <- logb + (fbar[i] + gbar - C[i, ]) * inv_eps
          lz <- logsumexp(z)
          if (!is.finite(lz)) next
          w <- exp(z - lz)
          w[!is.finite(w)] <- 0
          s <- sum(w)
          if (!(s > 0)) next
          w <- w / s
          out[i, ] <- as.numeric(crossprod(w, anchor_scores))
        }
        out
      }

      blocks <- lapply(seq_along(raw_domains), function(k) {
        if (k == 1L) {
          anchor_scores
        } else {
          map_to_anchor(raw_domains[[k]]$x)
        }
      })
      scores_uot <- do.call(rbind, blocks)
      fit_uot <- list(s = scores_uot)
      met_uot <- evaluate_fit(fit_uot, raw_domains = raw_domains, holdout_pairs = holdout_pairs)
      score_unstable <- is_uot_unstable(scores_uot)
      metric_unstable <- is_uot_metric_unstable(met_uot, scenario_name)
      if (score_unstable || metric_unstable) {
        if (score_unstable) note <- append_note(note, "uot_guardrail_retune")
        if (metric_unstable) note <- append_note(note, "uot_metric_guardrail_retune")
        params <- benchmark_uot_params(scenario_name, hks = FALSE, fallback = TRUE)
        eps <- params$epsilon
        rho1 <- params$rho1
        rho2 <- params$rho2
        blocks <- lapply(seq_along(raw_domains), function(k) {
          if (k == 1L) {
            anchor_scores
          } else {
            map_to_anchor(raw_domains[[k]]$x)
          }
        })
        scores_uot <- do.call(rbind, blocks)
        fit_uot <- list(s = scores_uot)
        met_uot <- evaluate_fit(fit_uot, raw_domains = raw_domains, holdout_pairs = holdout_pairs)
        score_unstable <- is_uot_unstable(scores_uot)
        metric_unstable <- is_uot_metric_unstable(met_uot, scenario_name)
        if (score_unstable) note <- append_note(note, "uot_guardrail_unstable")
        if (metric_unstable) note <- append_note(note, "uot_metric_guardrail_unstable")
      }
      attr(fit_uot, "benchmark_metrics") <- met_uot
      fit_uot
    },
    uot_hks = {
      if (supervision != "none") stop("uot_hks uses supervision='none'")

      # Anchor/template = domain 1 points; map each domain into anchor score space.
      anchor <- raw_domains[[1]]$x
      d_raw <- min(ncomp, ncol(anchor))
      anchor_raw <- stats::prcomp(anchor, center = FALSE, scale. = FALSE)$x[, seq_len(d_raw), drop = FALSE]
      anchor_scores <- matrix(0, nrow(anchor_raw), ncomp)
      anchor_scores[, seq_len(d_raw)] <- anchor_raw

      # Structural augmentation: HKS descriptors per domain, used as OT features.
      knn_struct <- 12L
      q_desc <- max(16L, 2L * ncomp)
      k_embed <- max(3L * ncomp, ncomp + 10L)

      compute_desc <- function(X) {
        sigma_use <- tryCatch(manifoldalign::choose_sigma(as.matrix(X)), error = function(e) 0.5)
        gw <- adjoin::graph_weights(
          as.matrix(X),
          k = min(knn_struct, nrow(X) - 1L),
          weight_mode = "heat",
          sigma = sigma_use,
          neighbor_mode = "knn"
        )
        A <- adjoin::adjacency(gw)
        A <- (A + Matrix::t(A)) / 2
        manifoldalign::compute_hks_from_adjacency(
          A,
          k_embed = k_embed,
          q = q_desc,
          use_normalized_laplacian = TRUE,
          time_mode = "auto",
          normalize = "both"
        )
      }

      desc_list <- lapply(raw_domains, function(d) compute_desc(d$x))
      desc_anchor <- desc_list[[1]]

      beta0 <- rep(1, nrow(anchor))
      params <- benchmark_uot_params(scenario_name, hks = TRUE, fallback = FALSE)
      eps <- params$epsilon
      rho1 <- params$rho1
      rho2 <- params$rho2
      lambda_feat <- params$lambda_feat

      logsumexp <- function(x) {
        mx <- max(x)
        if (!is.finite(mx)) return(mx)
        mx + log(sum(exp(x - mx)))
      }

      map_to_anchor <- function(Xk, Fk) {
        alpha0 <- rep(1, nrow(Xk))
        fit_k <- manifoldalign::uot_fit_pair(
          X = Xk,
          Y = anchor,
          alpha = alpha0,
          beta = beta0,
          F = Fk,
          G = desc_anchor,
          lambda_anat = 1,
          lambda_feat = lambda_feat,
          neighbor_mode = "dense",
          epsilon = eps,
          rho1 = rho1,
          rho2 = rho2,
          max_iter = 2000,
          tol = 1e-10
        )

        C <- fit_k$cost
        fbar <- as.numeric(fit_k$sol$fbar)
        gbar <- as.numeric(fit_k$sol$gbar)
        logb <- rep(0, nrow(anchor)) # log(beta0) with beta0 uniform

        out <- matrix(0, nrow(Xk), ncomp)
        inv_eps <- 1 / eps
        for (i in seq_len(nrow(Xk))) {
          z <- logb + (fbar[i] + gbar - C[i, ]) * inv_eps
          lz <- logsumexp(z)
          if (!is.finite(lz)) next
          w <- exp(z - lz)
          w[!is.finite(w)] <- 0
          s <- sum(w)
          if (!(s > 0)) next
          w <- w / s
          out[i, ] <- as.numeric(crossprod(w, anchor_scores))
        }
        out
      }

      blocks <- lapply(seq_along(raw_domains), function(k) {
        if (k == 1L) {
          anchor_scores
        } else {
          map_to_anchor(raw_domains[[k]]$x, desc_list[[k]])
        }
      })
      scores_uot <- do.call(rbind, blocks)
      fit_uot <- list(s = scores_uot)
      met_uot <- evaluate_fit(fit_uot, raw_domains = raw_domains, holdout_pairs = holdout_pairs)
      score_unstable <- is_uot_unstable(scores_uot)
      metric_unstable <- is_uot_metric_unstable(met_uot, scenario_name)
      if (score_unstable || metric_unstable) {
        if (score_unstable) note <- append_note(note, "uot_hks_guardrail_retune")
        if (metric_unstable) note <- append_note(note, "uot_hks_metric_guardrail_retune")
        params <- benchmark_uot_params(scenario_name, hks = TRUE, fallback = TRUE)
        eps <- params$epsilon
        rho1 <- params$rho1
        rho2 <- params$rho2
        lambda_feat <- params$lambda_feat
        blocks <- lapply(seq_along(raw_domains), function(k) {
          if (k == 1L) {
            anchor_scores
          } else {
            map_to_anchor(raw_domains[[k]]$x, desc_list[[k]])
          }
        })
        scores_uot <- do.call(rbind, blocks)
        fit_uot <- list(s = scores_uot)
        met_uot <- evaluate_fit(fit_uot, raw_domains = raw_domains, holdout_pairs = holdout_pairs)
        score_unstable <- is_uot_unstable(scores_uot)
        metric_unstable <- is_uot_metric_unstable(met_uot, scenario_name)
        if (score_unstable) note <- append_note(note, "uot_hks_guardrail_unstable")
        if (metric_unstable) note <- append_note(note, "uot_hks_metric_guardrail_unstable")
      }
      attr(fit_uot, "benchmark_metrics") <- met_uot
      fit_uot
    },
    ot_procrustes = {
      if (supervision != "none") stop("ot_procrustes uses supervision='none'")

      # Alvarez-Melis et al. note that whitening is a common preprocessing step
      # for correspondence / feature alignment problems, and it can convert a
      # general linear relation into an (approximately) orthogonal one.
      whiten <- function(X, ridge = 1e-6) {
        X <- scale(as.matrix(X), center = TRUE, scale = FALSE)
        S <- stats::cov(X)
        eig <- eigen(S, symmetric = TRUE)
        vals <- pmax(eig$values, ridge)
        W <- eig$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(eig$vectors)
        X %*% W
      }

      mats <- lapply(raw_domains, function(d) {
        X <- as.matrix(d$x)
        if (!all(is.finite(X))) stop("ot_procrustes: domain contains non-finite values.", call. = FALSE)
        whiten(X)
      })

      d_feat <- ncol(mats[[1]])
      if (any(vapply(mats, ncol, integer(1)) != d_feat)) {
        stop("ot_procrustes requires matching feature dimensions across domains.", call. = FALSE)
      }

      res <- manifoldalign::align_many(
        mats,
        manifoldalign::ot_procrustes_aligner(),
        graph = "complete",
        consensus = "sync",
        k = d_feat,
        parallel = FALSE,
        epsilon0 = 1,
        epsilon_min = 0.05,
        decay = 0.8,
        n_init = 2L,
        seed = as.integer(method_seed),
        max_iter = 40,
        tol = 1e-6,
        sinkhorn_max_iter = 200,
        sinkhorn_tol = 1e-6,
        stabilized = FALSE,
        init = "identity",
        store_transport = FALSE
      )

      blocks <- lapply(res$embeddings, function(Z) {
        Z <- as.matrix(Z)
        if (ncol(Z) >= ncomp) return(Z[, seq_len(ncomp), drop = FALSE])
        out <- matrix(0, nrow(Z), ncomp)
        out[, seq_len(ncol(Z))] <- Z
        out
      })

      list(s = do.call(rbind, blocks))
    },
    multiscale = {
      if (supervision == "id") {
        manifoldalign::multiscale_manifold_align(
          hd,
          correspondences = correspondences,
          ncomp = ncomp,
          control = manifoldalign::multiscale_manifold_align_control(
            enabled = TRUE,
            backend = "spectral",
            rank_per_domain = 128L,
            knn = 12L,
            cross_edge_weight = 1.0,
            seed = as.integer(method_seed),
            verbose = FALSE
          ),
          verbose = FALSE
        )
      } else if (supervision == "condition") {
        manifoldalign::multiscale_manifold_align(
          hd,
          y = condition_train,
          ncomp = ncomp,
          control = manifoldalign::multiscale_manifold_align_control(
            enabled = TRUE,
            backend = "spectral",
            rank_per_domain = 128L,
            knn = 12L,
            max_pairs_per_label = 80L,
            cross_edge_weight = 1.0,
            seed = as.integer(method_seed),
            verbose = FALSE
          ),
          verbose = FALSE
        )
      } else {
        manifoldalign::multiscale_manifold_align(
          hd,
          ncomp = ncomp,
          control = manifoldalign::multiscale_manifold_align_control(
            enabled = TRUE,
            backend = "spectral",
            rank_per_domain = 128L,
            knn = 12L,
            mnn_k = 8L,
            max_unsup_pairs_per_pair = 300L,
            cross_edge_weight = 1.0,
            seed = as.integer(method_seed),
            verbose = FALSE
          ),
          verbose = FALSE
        )
      }
    },
    stop("Unknown method: ", method)
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  met <- attr(fit, "benchmark_metrics", exact = TRUE)
  if (is.null(met)) {
    met <- evaluate_fit(fit, raw_domains = raw_domains, holdout_pairs = holdout_pairs)
  }
  c(list(method = method, supervision = supervision, runtime_sec = elapsed, note = note), met)
}

run_benchmark <- function(n_reps = 3L, ncomp = 8L, seed = 20260208L) {
  set.seed(seed)
  base_domains <- manifoldalign::alignment_benchmark$domains
  scenarios <- build_benchmark_scenarios(base_domains, seed = seed)

  methods_supervision <- list(
    global_geo = c("id", "condition"),
    kema = c("id", "condition"),
    lra = c("id", "condition"),
    ssma = c("id", "condition"),
    gpca = c("id", "condition"),
    uot = "none",
    uot_hks = "none",
    ot_procrustes = "none",
    multiscale = c("id", "condition", "none"),
    spectral_mnn = c("id", "condition", "none"),
    cone_multi = "none",
    grasp_multi = "none"
  )
  rows <- list()

  for (rep_id in seq_len(n_reps)) {
    for (scenario_name in names(scenarios)) {
      scenario <- scenarios[[scenario_name]]
      scenario_seed <- as.integer(seed + rep_id * 1000L + match(scenario_name, names(scenarios)) * 100L)
      scenario_domains <- apply_supervision_budget(
        scenario$domains,
        id_keep_frac = scenario$id_keep_frac,
        condition_keep_frac = scenario$condition_keep_frac,
        seed = scenario_seed
      )
      scenario_domains <- apply_rep_perturbation(
        scenario_domains,
        noise_sd = scenario$perturb_sd,
        seed = scenario_seed + 1L
      )

      dat <- build_benchmark_hyperdesign(scenario_domains)
      hd <- dat$hd
      raw_domains <- dat$raw
      corr <- build_id_correspondences(raw_domains, id_col = "id")
      holdout_pairs <- build_holdout_pairs(raw_domains, id_train_col = "id", eval_col = "eval_id")

      for (m in names(methods_supervision)) {
        for (sup in methods_supervision[[m]]) {
          sup_idx <- match(sup, c("id", "condition", "none"))
          method_seed <- as.integer(scenario_seed + match(m, names(methods_supervision)) * 10L + sup_idx)
          message(sprintf("[rep %d/%d][%s] %s (%s)", rep_id, n_reps, scenario_name, m, sup))
          rows[[length(rows) + 1L]] <- tryCatch(
            {
              out <- run_one(
                method = m,
                supervision = sup,
                hd = hd,
                raw_domains = raw_domains,
                correspondences = corr,
                holdout_pairs = holdout_pairs,
                ncomp = ncomp,
                scenario_name = scenario_name,
                method_seed = method_seed
              )
              c(out, list(scenario = scenario_name, rep = rep_id, ok = TRUE, error = NA_character_))
            },
            error = function(e) {
              list(
                method = m,
                supervision = sup,
                runtime_sec = NA_real_,
                note = NA_character_,
                ncomp = NA_integer_,
                identity_mse = NA_real_,
                holdout_identity_mse = NA_real_,
                geometry_spearman = NA_real_,
                condition_knn_acc = NA_real_,
                holdout_condition_knn_acc = NA_real_,
                scenario = scenario_name,
                rep = rep_id,
                ok = FALSE,
                error = e$message
              )
            }
          )
        }
      }
    }
  }

  out <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(out) <- NULL
  out
}

summarize_benchmark <- function(results) {
  ok <- results[results$ok, , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  key <- paste(ok$scenario, ok$method, ok$supervision, sep = "||")
  split_ok <- split(ok, key)
  out <- lapply(names(split_ok), function(k) {
    df <- split_ok[[k]]
    scenario <- df$scenario[[1L]]
    method <- df$method[[1L]]
    supervision <- df$supervision[[1L]]
    data.frame(
      scenario = scenario,
      method = method,
      supervision = supervision,
      runtime_mean = mean(df$runtime_sec, na.rm = TRUE),
      runtime_sd = stats::sd(df$runtime_sec, na.rm = TRUE),
      identity_mse_mean = mean(df$identity_mse, na.rm = TRUE),
      identity_mse_sd = stats::sd(df$identity_mse, na.rm = TRUE),
      holdout_identity_mse_mean = mean(df$holdout_identity_mse, na.rm = TRUE),
      holdout_identity_mse_sd = stats::sd(df$holdout_identity_mse, na.rm = TRUE),
      geometry_spearman_mean = mean(df$geometry_spearman, na.rm = TRUE),
      geometry_spearman_sd = stats::sd(df$geometry_spearman, na.rm = TRUE),
      condition_knn_acc_mean = mean(df$condition_knn_acc, na.rm = TRUE),
      condition_knn_acc_sd = stats::sd(df$condition_knn_acc, na.rm = TRUE),
      holdout_condition_knn_acc_mean = mean(df$holdout_condition_knn_acc, na.rm = TRUE),
      holdout_condition_knn_acc_sd = stats::sd(df$holdout_condition_knn_acc, na.rm = TRUE)
    )
  })
  do.call(rbind, out)
}

if (sys.nframe() == 0L) {
  message("Running global_geo vs multi-method baseline benchmark (id + condition supervision modes)...")
  results <- run_benchmark(n_reps = 3L, ncomp = 8L, seed = 20260208L)
  summary_tbl <- summarize_benchmark(results)

  print(results)
  cat("\nSummary (mean +/- sd across successful reps):\n")
  print(summary_tbl)

  out_csv <- "dev-scripts/benchmarks/global_geo_vs_baselines_results.csv"
  out_summary_csv <- "dev-scripts/benchmarks/global_geo_vs_baselines_summary.csv"
  utils::write.csv(results, out_csv, row.names = FALSE)
  utils::write.csv(summary_tbl, out_summary_csv, row.names = FALSE)
  # Backward-compatible path used in previous runs.
  utils::write.csv(results, "dev-scripts/benchmarks/global_geo_vs_kema_lra_results.csv", row.names = FALSE)
  message("Saved: ", out_csv)
  message("Saved: ", out_summary_csv)
}
