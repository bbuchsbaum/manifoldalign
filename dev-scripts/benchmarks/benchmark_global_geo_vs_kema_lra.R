# benchmark_global_geo_vs_kema_lra.R
# ------------------------------------------------------------
# Compare global_geo_align vs KEMA vs lowrank_align on the same
# multi-domain input and correspondence structure.
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
  library(neighborweights)
  library(dplyr)
})

build_benchmark_hyperdesign <- function() {
  bench <- manifoldalign::alignment_benchmark
  domains <- lapply(bench$domains, function(dom) {
    multidesign::multidesign(dom$x, dom$design)
  })
  hd <- multidesign::hyperdesign(domains)
  list(hd = hd, raw = bench$domains, labels = bench$labels)
}

build_id_correspondences <- function(domains) {
  domain_names <- names(domains)
  m <- length(domains)
  out <- vector("list", m * (m - 1) / 2)
  ptr <- 0L

  for (i in seq_len(m - 1L)) {
    ids_i <- domains[[i]]$design$id
    for (j in (i + 1L):m) {
      ids_j <- domains[[j]]$design$id
      common <- intersect(ids_i, ids_j)
      if (!length(common)) next
      ptr <- ptr + 1L
      out[[ptr]] <- data.frame(
        domain_i = i,
        index_i = match(common, ids_i),
        domain_j = j,
        index_j = match(common, ids_j)
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
    id_i <- raw_domains[[i]]$design$id
    for (j in (i + 1L):m) {
      id_j <- raw_domains[[j]]$design$id
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

evaluate_fit <- function(fit, raw_domains, labels_condition) {
  domain_sizes <- vapply(raw_domains, function(d) nrow(d$x), integer(1))
  score_blocks <- split_scores_by_domain(as.matrix(fit$s), domain_sizes)
  pooled_labels <- rep(as.character(labels_condition), times = length(domain_sizes))

  list(
    ncomp = ncol(fit$s),
    identity_mse = identity_pair_mse(score_blocks, raw_domains),
    geometry_spearman = geometry_spearman(score_blocks, raw_domains),
    condition_knn_acc = label_knn_accuracy(as.matrix(fit$s), pooled_labels, k = 5L)
  )
}

run_one <- function(method, hd, raw_domains, correspondences, ncomp) {
  t0 <- proc.time()[["elapsed"]]
  fit <- switch(
    method,
    global_geo = manifoldalign::global_geo_align(
      hd,
      correspondences = correspondences,
      preproc = NULL,
      ncomp = ncomp,
      control = manifoldalign::global_geo_align_control(
        knn = 12,
        n_landmarks_total = 120,
        k_embed = max(3L * ncomp, ncomp + 6L),
        scale_method = "none",
        ridge_eps = 1e-4,
        verbose = FALSE,
        seed = 1
      )
    ),
    kema = manifoldalign::kema(
      hd,
      y = id,                 # use exact sample-id correspondences
      preproc = multivarious::center(),
      ncomp = ncomp,
      knn = 5,
      solver = "regression",
      sample_frac = 0.8,
      lambda = 1e-3
    ),
    lra = manifoldalign::lowrank_align(
      hd,
      y = id,                 # use exact sample-id correspondences
      ncomp = ncomp,
      simfun = neighborweights::binary_label_matrix,
      solver = "operator",
      sv_thresh = 1,
      lambda = 1e-2
    ),
    stop("Unknown method: ", method)
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  # condition labels from benchmark (shared 2-class latent structure)
  labels_condition <- raw_domains[[1]]$design$condition
  met <- evaluate_fit(fit, raw_domains = raw_domains, labels_condition = labels_condition)
  c(list(method = method, runtime_sec = elapsed), met)
}

run_benchmark <- function(n_reps = 3L, ncomp = 8L, seed = 20260208L) {
  set.seed(seed)
  dat <- build_benchmark_hyperdesign()
  hd <- dat$hd
  raw_domains <- dat$raw
  corr <- build_id_correspondences(raw_domains)

  methods <- c("global_geo", "kema", "lra")
  rows <- vector("list", length(methods) * n_reps)
  ptr <- 0L

  for (rep_id in seq_len(n_reps)) {
    for (m in methods) {
      message(sprintf("[rep %d/%d] %s", rep_id, n_reps, m))
      ptr <- ptr + 1L
      rows[[ptr]] <- tryCatch(
        {
          out <- run_one(
            method = m,
            hd = hd,
            raw_domains = raw_domains,
            correspondences = corr,
            ncomp = ncomp
          )
          c(out, list(rep = rep_id, ok = TRUE, error = NA_character_))
        },
        error = function(e) {
          list(
            method = m,
            runtime_sec = NA_real_,
            ncomp = NA_integer_,
            identity_mse = NA_real_,
            geometry_spearman = NA_real_,
            condition_knn_acc = NA_real_,
            rep = rep_id,
            ok = FALSE,
            error = e$message
          )
        }
      )
    }
  }

  out <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(out) <- NULL
  out
}

summarize_benchmark <- function(results) {
  ok <- results[results$ok, , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  split_ok <- split(ok, ok$method)
  out <- lapply(names(split_ok), function(m) {
    df <- split_ok[[m]]
    data.frame(
      method = m,
      runtime_mean = mean(df$runtime_sec, na.rm = TRUE),
      runtime_sd = stats::sd(df$runtime_sec, na.rm = TRUE),
      identity_mse_mean = mean(df$identity_mse, na.rm = TRUE),
      identity_mse_sd = stats::sd(df$identity_mse, na.rm = TRUE),
      geometry_spearman_mean = mean(df$geometry_spearman, na.rm = TRUE),
      geometry_spearman_sd = stats::sd(df$geometry_spearman, na.rm = TRUE),
      condition_knn_acc_mean = mean(df$condition_knn_acc, na.rm = TRUE),
      condition_knn_acc_sd = stats::sd(df$condition_knn_acc, na.rm = TRUE)
    )
  })
  do.call(rbind, out)
}

if (sys.nframe() == 0L) {
  message("Running global_geo vs kema vs lra benchmark...")
  results <- run_benchmark(n_reps = 3L, ncomp = 8L, seed = 20260208L)
  summary_tbl <- summarize_benchmark(results)

  print(results)
  cat("\nSummary (mean +/- sd across successful reps):\n")
  print(summary_tbl)

  out_csv <- "dev-scripts/benchmarks/global_geo_vs_kema_lra_results.csv"
  utils::write.csv(results, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
}
