#!/usr/bin/env Rscript

# Derive a small, real-data CI fixture from the governed So2Sat validation and
# testing artifacts. The fixture is class-balanced to exercise schemas and
# leakage boundaries; it is not a representative accuracy benchmark. Official
# accuracy must use the full frozen split through run_real_alignment_benchmarks.R.

parse_args <- function(args) {
  data_dir <- Sys.getenv("MANIFOLDALIGN_REAL_DATA_DIR", unset = "")
  output <- file.path(
    "inst", "extdata", "alignment-real", "so2sat_lcz42_ci_features.csv"
  )
  per_class <- 4L
  chunk_size <- 128L
  for (arg in args) {
    if (startsWith(arg, "--data-dir=")) {
      data_dir <- sub("^--data-dir=", "", arg)
    } else if (startsWith(arg, "--output=")) {
      output <- sub("^--output=", "", arg)
    } else if (startsWith(arg, "--per-class=")) {
      value <- sub("^--per-class=", "", arg)
      per_class <- if (identical(value, "all")) Inf else as.integer(value)
    } else if (startsWith(arg, "--chunk-size=")) {
      chunk_size <- as.integer(sub("^--chunk-size=", "", arg))
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  if (!nzchar(data_dir) || !dir.exists(data_dir)) {
    stop("Provide --data-dir containing verified So2Sat v4 artifacts.",
      call. = FALSE
    )
  }
  if (is.na(per_class) || per_class < 1L) {
    stop("--per-class must be a positive integer or `all`.", call. = FALSE)
  }
  if (is.na(chunk_size) || chunk_size < 1L) {
    stop("--chunk-size must be a positive integer.", call. = FALSE)
  }
  list(
    data_dir = normalizePath(data_dir),
    output = normalizePath(output, mustWork = FALSE),
    per_class = per_class,
    chunk_size = chunk_size
  )
}

gunzip_to <- function(source, destination) {
  input <- gzfile(source, open = "rb")
  output <- file(destination, open = "wb")
  on.exit(close(input), add = TRUE)
  on.exit(close(output), add = TRUE)
  repeat {
    block <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(block)) break
    writeBin(block, output)
  }
  destination
}

as_sample_major <- function(x, n_samples, feature_count = NULL) {
  x <- as.matrix(x)
  if (nrow(x) == n_samples) return(x)
  if (ncol(x) == n_samples) return(t(x))
  if (!is.null(feature_count) && ncol(x) == feature_count) return(x)
  stop("HDF5 array does not have a recognizable sample dimension.",
    call. = FALSE
  )
}

decode_city <- function(x, n_samples) {
  x <- as.vector(x)
  if (length(x) != n_samples) {
    return(rep("unknown", n_samples))
  }
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "unknown"
  x
}

patch_features <- function(patch) {
  values <- as.numeric(patch)
  values <- values[is.finite(values)]
  if (!length(values)) return(rep(NA_real_, 10L))
  # Exact means/extrema use every value. Quantiles use a fixed, evenly spaced
  # grid so full-split extraction remains linear-time and reproducible rather
  # than sorting every 8-10k-value patch.
  quantile_values <- if (length(values) > 512L) {
    values[unique(as.integer(round(seq(
      1, length(values), length.out = 512L
    ))))]
  } else {
    values
  }
  quantile <- as.numeric(stats::quantile(
    quantile_values, probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
    names = FALSE, type = 8
  ))
  centered <- values - mean(values)
  scale <- stats::sd(values)
  if (!is.finite(scale) || scale < 1e-12) scale <- 1
  standardized <- centered / scale
  c(
    mean = mean(values), sd = stats::sd(values), min = min(values),
    q10 = quantile[[1L]], q25 = quantile[[2L]], median = quantile[[3L]],
    q75 = quantile[[4L]], q90 = quantile[[5L]], max = max(values),
    abs_z = mean(abs(standardized))
  )
}

select_balanced <- function(labels, per_class) {
  if (is.infinite(per_class)) return(seq_along(labels))
  classes <- sort(unique(labels))
  unlist(lapply(classes, function(label) {
    head(which(labels == label), per_class)
  }), use.names = FALSE)
}

sensor_feature_matrix <- function(dataset, index, n_samples, chunk_size) {
  if (length(dataset$dims) != 4L || tail(dataset$dims, 1L) != n_samples) {
    stop("Expected hdf5r sensor arrays with samples in the final R dimension.",
      call. = FALSE
    )
  }
  feature_names <- names(patch_features(1:10))
  out <- matrix(
    NA_real_, nrow = length(index), ncol = length(feature_names),
    dimnames = list(NULL, feature_names)
  )
  chunk <- split(
    seq_along(index),
    ceiling(seq_along(index) / chunk_size)
  )
  for (positions in chunk) {
    selected <- index[positions]
    block <- dataset[, , , selected]
    if (length(selected) == 1L) {
      block <- array(block, dim = c(dim(block), 1L))
    }
    for (j in seq_along(positions)) {
      out[positions[[j]], ] <- patch_features(block[, , , j])
    }
  }
  out
}

read_split <- function(data_dir, split, public_split, per_class, chunk_size) {
  compressed <- file.path(data_dir, "v4", paste0(split, ".h5.gz"))
  geo_path <- file.path(data_dir, "v4", paste0(split, "_geo.h5"))
  if (!file.exists(compressed) || !file.exists(geo_path)) {
    stop("Missing verified So2Sat artifacts for split ", split, ".",
      call. = FALSE
    )
  }
  uncompressed <- tempfile(paste0("so2sat-", split, "-"), fileext = ".h5")
  on.exit(unlink(uncompressed), add = TRUE)
  gunzip_to(compressed, uncompressed)

  data_h5 <- hdf5r::H5File$new(uncompressed, mode = "r")
  geo_h5 <- hdf5r::H5File$new(geo_path, mode = "r")
  on.exit(data_h5$close_all(), add = TRUE)
  on.exit(geo_h5$close_all(), add = TRUE)

  labels_raw <- data_h5[["label"]]$read()
  label_matrix <- as.matrix(labels_raw)
  if (ncol(label_matrix) != 17L && nrow(label_matrix) == 17L) {
    label_matrix <- t(label_matrix)
  }
  if (ncol(label_matrix) != 17L) {
    stop("Expected 17 one-hot So2Sat classes.", call. = FALSE)
  }
  labels <- max.col(label_matrix, ties.method = "first")
  index <- sort(select_balanced(labels, per_class))
  n_samples <- nrow(label_matrix)

  tfw <- as_sample_major(geo_h5[["tfw"]]$read(), n_samples)
  epsg <- as.vector(geo_h5[["epsg"]]$read())
  if (length(epsg) != n_samples) epsg <- rep(NA_integer_, n_samples)
  city <- if ("city" %in% names(geo_h5)) {
    decode_city(geo_h5[["city"]]$read(), n_samples)
  } else {
    rep("unknown", n_samples)
  }

  sen1 <- data_h5[["sen1"]]
  sen2 <- data_h5[["sen2"]]
  f1 <- sensor_feature_matrix(sen1, index, n_samples, chunk_size)
  f2 <- sensor_feature_matrix(sen2, index, n_samples, chunk_size)
  rows <- data.frame(
    split = public_split,
    sample_id = paste0(split, "-", index),
    source_index = index,
    city = city[index],
    epsg = as.integer(epsg[index]),
    x = as.numeric(tfw[index, ncol(tfw) - 1L]),
    y = as.numeric(tfw[index, ncol(tfw)]),
    label = as.integer(labels[index]),
    stringsAsFactors = FALSE
  )
  for (k in seq_len(ncol(f1))) {
    rows[[paste0("s1_", colnames(f1)[[k]])]] <- f1[, k]
    rows[[paste0("s2_", colnames(f2)[[k]])]] <- f2[, k]
  }
  rows
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("Building the So2Sat CI fixture requires the `hdf5r` package.",
    call. = FALSE
  )
}
validation <- read_split(
  args$data_dir, "validation", "validation_west", args$per_class,
  args$chunk_size
)
testing <- read_split(
  args$data_dir, "testing", "testing_east", args$per_class,
  args$chunk_size
)
fixture <- rbind(validation, testing)
fixture <- fixture[order(
  match(fixture$split, c("validation_west", "testing_east")),
  fixture$label, fixture$source_index
), , drop = FALSE]
rownames(fixture) <- NULL
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(fixture, args$output, row.names = FALSE, na = "")
message(
  "Wrote ", nrow(fixture), " real So2Sat rows (per-class=",
  if (is.infinite(args$per_class)) "all" else args$per_class, ") to ",
  normalizePath(args$output)
)
