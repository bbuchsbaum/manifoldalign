# Internal validation helpers for leakage-safe RAPID-MA benchmarks.

.rapid_validation_int <- function(value, name, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      abs(value - round(value)) > sqrt(.Machine$double.eps) ||
      value < minimum) {
    stop("`", name, "` must be one integer >= ", minimum, ".",
         call. = FALSE)
  }
  as.integer(value)
}

.rapid_validation_positions <- function(position) {
  position <- as.matrix(position)
  if (!is.numeric(position) || ncol(position) < 2L || nrow(position) < 1L) {
    stop("`position` must be a non-empty numeric matrix with at least two columns.",
         call. = FALSE)
  }
  position <- position[, 1:2, drop = FALSE]
  storage.mode(position) <- "double"
  if (any(!is.finite(position))) {
    stop("`position` must contain only finite coordinates.", call. = FALSE)
  }
  position
}

.rapid_validation_grid <- function(grid) {
  if (!is.numeric(grid) || length(grid) != 2L || any(!is.finite(grid)) ||
      any(abs(grid - round(grid)) > sqrt(.Machine$double.eps)) ||
      any(grid < 1L)) {
    stop("`grid` must contain two positive integers.", call. = FALSE)
  }
  as.integer(grid)
}

.rapid_spatial_axis_bin <- function(value, bins) {
  span <- range(value)
  if (span[[2L]] <= span[[1L]]) return(rep.int(1L, length(value)))
  normalized <- (value - span[[1L]]) / (span[[2L]] - span[[1L]])
  pmin.int(bins, as.integer(floor(normalized * bins)) + 1L)
}

# Assign equal-extent spatial tiles. The partition is invariant to row order
# and positive affine changes of either coordinate.
.rapid_spatial_blocks <- function(position, grid = c(4L, 4L)) {
  position <- .rapid_validation_positions(position)
  grid <- .rapid_validation_grid(grid)
  x_bin <- .rapid_spatial_axis_bin(position[, 1L], grid[[1L]])
  y_bin <- .rapid_spatial_axis_bin(position[, 2L], grid[[2L]])
  as.integer(x_bin + (y_bin - 1L) * grid[[1L]])
}

.rapid_validation_sample_by_class <- function(rows, labels, budget) {
  if (!length(rows) || budget == 0L) return(integer(0))
  observed <- as.character(labels[rows])
  keep <- !is.na(observed) & nzchar(observed)
  rows <- rows[keep]
  observed <- observed[keep]
  classes <- sort(unique(observed), method = "radix")
  selected <- lapply(classes, function(class) {
    candidate <- rows[observed == class]
    sample(candidate, min(length(candidate), budget), replace = FALSE)
  })
  as.integer(unlist(selected, use.names = FALSE))
}

# Build a strict inductive split. Test-block labels are deliberately never
# indexed while sampling fit or calibration state.
.rapid_blocked_inductive_split <- function(
  labels,
  position,
  test_block,
  calibration_block,
  grid = c(4L, 4L),
  labels_per_class = 4L,
  calibration_per_class = 4L,
  max_fit_rows = Inf,
  seed = 1L
) {
  position <- .rapid_validation_positions(position)
  labels <- as.character(labels)
  if (length(labels) != nrow(position)) {
    stop("`labels` must match the number of position rows.", call. = FALSE)
  }
  grid <- .rapid_validation_grid(grid)
  test_block <- .rapid_validation_int(test_block, "test_block")
  calibration_block <- .rapid_validation_int(
    calibration_block, "calibration_block"
  )
  if (test_block == calibration_block) {
    stop("Test and calibration blocks must be different.", call. = FALSE)
  }
  labels_per_class <- .rapid_validation_int(
    labels_per_class, "labels_per_class", minimum = 0L
  )
  calibration_per_class <- .rapid_validation_int(
    calibration_per_class, "calibration_per_class", minimum = 0L
  )
  seed <- .rapid_validation_int(seed, "seed", minimum = 0L)
  if (!is.numeric(max_fit_rows) || length(max_fit_rows) != 1L ||
      is.na(max_fit_rows) || max_fit_rows <= 0 ||
      (!is.infinite(max_fit_rows) &&
       abs(max_fit_rows - round(max_fit_rows)) > sqrt(.Machine$double.eps))) {
    stop("`max_fit_rows` must be a positive integer or Inf.", call. = FALSE)
  }

  block_id <- .rapid_spatial_blocks(position, grid)
  present <- unique(block_id)
  if (!test_block %in% present) {
    stop("`test_block` is not present in the spatial partition.", call. = FALSE)
  }
  if (!calibration_block %in% present) {
    stop("`calibration_block` is not present in the spatial partition.",
         call. = FALSE)
  }
  test_rows <- which(block_id == test_block)
  calibration_candidates <- which(block_id == calibration_block)
  fit_candidates <- which(
    block_id != test_block & block_id != calibration_block
  )

  set.seed(seed)
  labeled_rows <- .rapid_validation_sample_by_class(
    fit_candidates, labels, labels_per_class
  )
  calibration_rows <- .rapid_validation_sample_by_class(
    calibration_candidates, labels, calibration_per_class
  )
  if (is.infinite(max_fit_rows)) {
    fit_rows <- fit_candidates
  } else {
    max_fit_rows <- as.integer(max_fit_rows)
    if (length(labeled_rows) > max_fit_rows) {
      stop("`max_fit_rows` is smaller than the required labeled-row budget.",
           call. = FALSE)
    }
    remainder <- setdiff(fit_candidates, labeled_rows)
    available <- max_fit_rows - length(labeled_rows)
    sampled <- if (length(remainder) > available) {
      sample(remainder, available, replace = FALSE)
    } else {
      remainder
    }
    fit_rows <- sort(c(labeled_rows, sampled), method = "radix")
  }

  list(
    fit_rows = as.integer(fit_rows),
    labeled_rows = as.integer(sort(labeled_rows, method = "radix")),
    calibration_rows = as.integer(sort(calibration_rows, method = "radix")),
    test_rows = as.integer(test_rows),
    block_id = block_id,
    grid = grid,
    test_block = test_block,
    calibration_block = calibration_block,
    provenance = list(
      seed = seed,
      labels_per_class = labels_per_class,
      calibration_per_class = calibration_per_class,
      max_fit_rows = max_fit_rows,
      test_labels_read_during_split = FALSE,
      fit_excludes_calibration = !any(fit_rows %in% calibration_rows),
      fit_excludes_test = !any(fit_rows %in% test_rows)
    )
  )
}

.rapid_validation_probability <- function(probability, labels = NULL) {
  probability <- as.matrix(probability)
  if (!is.numeric(probability) || nrow(probability) < 1L ||
      ncol(probability) < 2L || any(!is.finite(probability)) ||
      any(probability < 0)) {
    stop("`probability` must be a finite nonnegative matrix with >=2 columns.",
         call. = FALSE)
  }
  if (is.null(colnames(probability)) || any(!nzchar(colnames(probability))) ||
      anyDuplicated(colnames(probability))) {
    stop("`probability` must have unique non-empty class column names.",
         call. = FALSE)
  }
  total <- rowSums(probability)
  if (any(total <= 0)) {
    stop("Every probability row must have positive mass.", call. = FALSE)
  }
  probability <- probability / total
  if (!is.null(labels)) {
    labels <- as.character(labels)
    if (length(labels) != nrow(probability) || anyNA(labels) ||
        any(!labels %in% colnames(probability))) {
      stop("`labels` must name one represented class per probability row.",
           call. = FALSE)
    }
  }
  probability
}

.rapid_temperature_apply <- function(probability, temperature) {
  probability <- .rapid_validation_probability(probability)
  if (!is.numeric(temperature) || length(temperature) != 1L ||
      !is.finite(temperature) || temperature <= 0) {
    stop("`temperature` must be one finite positive number.", call. = FALSE)
  }
  logits <- log(pmax(probability, 1e-15)) / temperature
  logits <- sweep(logits, 1L, apply(logits, 1L, max), "-")
  output <- exp(logits)
  output <- output / rowSums(output)
  dimnames(output) <- dimnames(probability)
  output
}

.rapid_probability_nll <- function(probability, labels) {
  probability <- .rapid_validation_probability(probability, labels)
  index <- match(as.character(labels), colnames(probability))
  -mean(log(pmax(probability[cbind(seq_len(nrow(probability)), index)], 1e-15)))
}

.rapid_temperature_fit <- function(
  probability,
  labels,
  interval = c(0.05, 20)
) {
  probability <- .rapid_validation_probability(probability, labels)
  if (!is.numeric(interval) || length(interval) != 2L ||
      any(!is.finite(interval)) || interval[[1L]] <= 0 ||
      interval[[2L]] <= interval[[1L]] ||
      interval[[1L]] > 1 || interval[[2L]] < 1) {
    stop("`interval` must be finite, positive, increasing, and contain 1.",
         call. = FALSE)
  }
  objective <- function(temperature) {
    .rapid_probability_nll(
      .rapid_temperature_apply(probability, temperature), labels
    )
  }
  baseline <- objective(1)
  optimized <- stats::optimize(objective, interval = interval)
  temperature <- optimized$minimum
  after <- optimized$objective
  if (!is.finite(after) || after > baseline) {
    temperature <- 1
    after <- baseline
  }
  structure(
    list(
      temperature = as.numeric(temperature),
      nll_before = baseline,
      nll_after = after,
      n = nrow(probability),
      class_levels = colnames(probability),
      provenance = list(training_only = TRUE, test_labels_read = FALSE)
    ),
    class = "rapid_ma_temperature"
  )
}

.rapid_inductive_interpolate <- function(
  training,
  values,
  query,
  k = 8L,
  zero_tolerance = sqrt(.Machine$double.eps)
) {
  training <- as.matrix(training)
  values <- as.matrix(values)
  query <- as.matrix(query)
  if (!is.numeric(training) || !is.numeric(values) || !is.numeric(query) ||
      nrow(training) < 1L || ncol(training) < 1L ||
      nrow(values) != nrow(training) || ncol(query) != ncol(training) ||
      any(!is.finite(training)) || any(!is.finite(values)) ||
      any(!is.finite(query))) {
    stop("Training, values, and query must be finite conformable matrices.",
         call. = FALSE)
  }
  k <- .rapid_validation_int(k, "k")
  neighbours <- .rapid_oos_weights(
    training, query, min(k, nrow(training)), zero_tolerance
  )
  value <- as.matrix(neighbours$interpolation %*% values)
  colnames(value) <- colnames(values)
  list(
    value = value,
    neighbours = neighbours$index,
    distance = neighbours$distance,
    weight = neighbours$weight,
    exact_reapplications =
      rowSums(neighbours$distance <= zero_tolerance) > 0L
  )
}
