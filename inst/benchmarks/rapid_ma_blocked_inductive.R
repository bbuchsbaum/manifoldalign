#!/usr/bin/env Rscript

# Leakage-safe real-scene blocked-inductive benchmark for RAPID-MA.
#
# Scene: Chikusei airborne hyperspectral imagery (real reflectance and ground
# truth). The target modality is a frozen ten-band Gaussian approximation to
# Sentinel-2 VNIR response functions, following the simulated HS-MS setting in
# Hong et al. (2019). Test tiles are absent from fitting, preprocessing,
# relation construction, and calibration. This is real-scene evidence, not a
# genuinely independent-sensor benchmark.

args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(
  VECLIB_MAXIMUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)
scene_dir <- if (length(args) >= 1L) args[[1L]] else {
  file.path(Sys.getenv("RAPID_MA_DATA_DIR", tempdir()), "Chikusei_ENVI")
}
output_prefix <- if (length(args) >= 2L) args[[2L]] else {
  file.path("docs", "rapid-ma-blocked-inductive")
}
profile <- Sys.getenv("RAPID_MA_PROFILE", "full")
if (!profile %in% c("quick", "full", "exhaustive")) {
  stop("RAPID_MA_PROFILE must be quick, full, or exhaustive.")
}

if (!"manifoldalign" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(manifoldalign)
  }
}

scene_files <- c(
  high = "HyperspecVNIR_Chikusei_20140729.bsq",
  high_header = "HyperspecVNIR_Chikusei_20140729.hdr",
  truth = "HyperspecVNIR_Chikusei_20140729_Ground_Truth.bsq",
  truth_header = "HyperspecVNIR_Chikusei_20140729_Ground_Truth.hdr",
  readme = "Readme.txt"
)
scene_paths <- file.path(scene_dir, scene_files)
if (!all(file.exists(scene_paths))) {
  stop(
    "Chikusei ENVI files are unavailable at ", scene_dir, ". Run ",
    "RAPID_MA_ACCEPT_CHIKUSEI_LICENSE=yes Rscript ",
    "inst/benchmarks/download_rapid_ma_scene.R <data-dir>."
  )
}

class_names <- c(
  "water", "bare_soil_school", "bare_soil_park", "bare_soil_farmland",
  "natural_plants", "weeds_farmland", "forest", "grass",
  "rice_grown", "rice_first_stage", "row_crops", "plastic_house",
  "manmade_non_dark", "manmade_dark", "manmade_blue", "manmade_red",
  "manmade_grass", "asphalt", "paved_ground"
)
wavelength_um <- c(
  0.36259, 0.36775, 0.37290, 0.37807, 0.38323, 0.38839, 0.39355,
  0.39871, 0.40387, 0.40903, 0.41419, 0.41936, 0.42452, 0.42968,
  0.43484, 0.44000, 0.44516, 0.45032, 0.45548, 0.46064, 0.46580,
  0.47096, 0.47612, 0.48129, 0.48645, 0.49161, 0.49677, 0.50193,
  0.50709, 0.51225, 0.51741, 0.52257, 0.52773, 0.53289, 0.53806,
  0.54321, 0.54838, 0.55354, 0.55870, 0.56386, 0.56902, 0.57418,
  0.57934, 0.58450, 0.58966, 0.59483, 0.59999, 0.60514, 0.61031,
  0.61547, 0.62063, 0.62579, 0.63095, 0.63611, 0.64127, 0.64643,
  0.65159, 0.65675, 0.66192, 0.66707, 0.67224, 0.67740, 0.68256,
  0.68772, 0.69288, 0.69804, 0.70320, 0.70836, 0.71352, 0.71868,
  0.72385, 0.72901, 0.73417, 0.73933, 0.74449, 0.74965, 0.75481,
  0.75997, 0.76513, 0.77029, 0.77545, 0.78061, 0.78578, 0.79094,
  0.79610, 0.80126, 0.80642, 0.81158, 0.81674, 0.82190, 0.82706,
  0.83223, 0.83738, 0.84254, 0.84771, 0.85287, 0.85803, 0.86319,
  0.86835, 0.87351, 0.87867, 0.88383, 0.88899, 0.89416, 0.89931,
  0.90448, 0.90964, 0.91480, 0.91996, 0.92512, 0.93028, 0.93544,
  0.94060, 0.94576, 0.95092, 0.95609, 0.96125, 0.96641, 0.97157,
  0.97673, 0.98189, 0.98705, 0.99221, 0.99737, 1.00253, 1.00769,
  1.01286, 1.01802
)
bad_bands <- c(65:74, 77:80, 84:93, 105:123)
good_bands <- setdiff(seq_along(wavelength_um), bad_bands)

sentinel_center_um <- c(
  coastal = 0.443, blue = 0.490, green = 0.560, red = 0.665,
  red_edge_1 = 0.705, red_edge_2 = 0.740, red_edge_3 = 0.783,
  nir_wide = 0.842, nir_narrow = 0.865, water_vapor = 0.945
)
sentinel_fwhm_um <- c(
  0.020, 0.065, 0.035, 0.030, 0.015,
  0.015, 0.020, 0.115, 0.020, 0.020
)

gaussian_response <- function(wavelength, center, fwhm) {
  sigma <- fwhm / (2 * sqrt(2 * log(2)))
  weight <- exp(-0.5 * ((wavelength - center) / sigma)^2)
  weight / sum(weight)
}

response <- vapply(seq_along(sentinel_center_um), function(j) {
  gaussian_response(
    wavelength_um[good_bands],
    sentinel_center_um[[j]],
    sentinel_fwhm_um[[j]]
  )
}, numeric(length(good_bands)))
colnames(response) <- names(sentinel_center_um)

read_chikusei_scene <- function(scene_dir, cache_path) {
  loader_version <- 2L
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$metadata$loader_version, loader_version)) {
      return(cached)
    }
  }

  samples <- 2335L
  lines <- 2517L
  pixels <- samples * lines
  truth_path <- file.path(
    scene_dir, "HyperspecVNIR_Chikusei_20140729_Ground_Truth.bsq"
  )
  high_path <- file.path(
    scene_dir, "HyperspecVNIR_Chikusei_20140729.bsq"
  )
  truth_connection <- file(truth_path, "rb")
  truth <- readBin(
    truth_connection, integer(), n = pixels,
    size = 1L, signed = FALSE, endian = "little"
  )
  close(truth_connection)
  labeled_index <- which(truth > 0L)
  row <- ((labeled_index - 1L) %/% samples) + 1L
  column <- ((labeled_index - 1L) %% samples) + 1L
  position <- cbind(
    easting_m = 408041.5 + (column - 1L) * 2.5,
    northing_m = 4020265.5 - (row - 1L) * 2.5
  )

  high <- matrix(
    NA_real_, length(labeled_index), length(good_bands),
    dimnames = list(NULL, paste0("hs_", good_bands))
  )
  high_connection <- file(high_path, "rb")
  replacement_count <- 0L
  for (j in seq_along(good_bands)) {
    band <- good_bands[[j]]
    seek(
      high_connection,
      where = as.double(band - 1L) * as.double(pixels) * 2,
      origin = "start",
      rw = "read"
    )
    full_band <- readBin(
      high_connection, integer(), n = pixels,
      size = 2L, signed = TRUE, endian = "little"
    )
    value <- as.numeric(full_band[labeled_index])
    invalid <- !is.finite(value) | value < 0 | value == 15000
    if (any(invalid)) {
      replacement <- stats::median(value[!invalid])
      value[invalid] <- replacement
      replacement_count <- replacement_count + sum(invalid)
    }
    high[, j] <- value / 10000
    if (j %% 10L == 0L || j == length(good_bands)) {
      message("loaded Chikusei band ", j, "/", length(good_bands))
    }
  }
  close(high_connection)
  low <- high %*% response
  colnames(low) <- paste0("s2_", colnames(response))

  nearest_band <- function(center) which.min(abs(wavelength_um[good_bands] - center))
  high_red <- high[, nearest_band(0.665)]
  high_nir <- high[, nearest_band(0.842)]
  high_attribute <- cbind(
    brightness = rowMeans(high),
    spectral_sd = apply(high, 1L, stats::sd),
    ndvi = (high_nir - high_red) / pmax(high_nir + high_red, 1e-8)
  )
  low_attribute <- cbind(
    brightness = rowMeans(low),
    spectral_sd = apply(low, 1L, stats::sd),
    ndvi = (low[, "s2_nir_wide"] - low[, "s2_red"]) /
      pmax(low[, "s2_nir_wide"] + low[, "s2_red"], 1e-8)
  )

  scene <- list(
    high = high,
    low = low,
    labels = class_names[truth[labeled_index]],
    position = position,
    attributes = list(high = high_attribute, low = low_attribute),
    raster_index = labeled_index,
    metadata = list(
      loader_version = loader_version,
      samples = samples,
      lines = lines,
      bands_original = length(wavelength_um),
      bands_retained = length(good_bands),
      labeled_pixels = length(labeled_index),
      invalid_values_replaced = replacement_count,
      target_bands = ncol(low),
      target_response = "frozen Gaussian approximation to Sentinel-2 VNIR",
      real_scene = TRUE,
      independent_target_sensor = FALSE,
      archive_sha256 =
        "192825c63e764a16e9324470e005992f343c24317178e8a5546147aadde838dc",
      license = "CC BY 3.0",
      citation = paste(
        "N. Yokoya and A. Iwasaki, Airborne hyperspectral data over",
        "Chikusei, SAL-2016-05-27, 2016."
      )
    )
  )
  saveRDS(scene, cache_path, compress = FALSE)
  scene
}

cache_path <- Sys.getenv(
  "RAPID_MA_CHIKUSEI_CACHE",
  file.path(dirname(scene_dir), "rapid-ma-chikusei-labelled-v2.rds")
)
scene <- read_chikusei_scene(scene_dir, cache_path)

fit_standardizer <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 1e-12] <- 1
  list(center = center, scale = scale)
}

apply_standardizer <- function(x, standardizer) {
  value <- sweep(as.matrix(x), 2L, standardizer$center, "-")
  sweep(value, 2L, standardizer$scale, "/")
}

fit_centroid <- function(embedding, labels) {
  labels <- as.character(labels)
  keep <- !is.na(labels) & nzchar(labels)
  embedding <- as.matrix(embedding)[keep, , drop = FALSE]
  labels <- labels[keep]
  classes <- sort(unique(labels), method = "radix")
  centroids <- do.call(rbind, lapply(classes, function(class) {
    colMeans(embedding[labels == class, , drop = FALSE])
  }))
  rownames(centroids) <- classes
  training_distance <- vapply(seq_along(classes), function(j) {
    rowSums((
      embedding - matrix(
        centroids[j, ], nrow(embedding), ncol(embedding), byrow = TRUE
      )
    )^2)
  }, numeric(nrow(embedding)))
  own <- training_distance[cbind(seq_along(labels), match(labels, classes))]
  distance_scale <- stats::median(own[is.finite(own) & own > 0])
  if (!is.finite(distance_scale) || distance_scale <= 1e-12) {
    distance_scale <- stats::median(
      training_distance[is.finite(training_distance) & training_distance > 0]
    )
  }
  if (!is.finite(distance_scale) || distance_scale <= 1e-12) distance_scale <- 1
  list(centroids = centroids, scale = distance_scale, classes = classes)
}

predict_centroid <- function(model, embedding) {
  embedding <- as.matrix(embedding)
  distance <- vapply(seq_len(nrow(model$centroids)), function(j) {
    rowSums((
      embedding - matrix(
        model$centroids[j, ], nrow(embedding), ncol(embedding), byrow = TRUE
      )
    )^2)
  }, numeric(nrow(embedding)))
  distance <- as.matrix(distance)
  probability <- exp(pmax(-distance / model$scale, -700))
  probability <- probability / pmax(rowSums(probability), 1e-15)
  colnames(probability) <- model$classes
  probability
}

structured_oos_weights <- c(feature = 0.5, position = 0.25, attribute = 0.25)

structured_target_predict <- function(
  training_x,
  training_position,
  training_attribute,
  training_probability,
  query_x,
  query_position,
  query_attribute,
  k = 8L
) {
  position_state <- manifoldalign:::.rapid_prepare_positions(
    training_position, nrow(training_x), mode = "relative"
  )
  attribute_state <- manifoldalign:::.rapid_encode_attributes(
    training_attribute, nrow(training_x), hash_dim = 32L, seed = 1L
  )
  training_views <- list(
    feature = training_x,
    position = position_state$view,
    attribute = attribute_state$view
  )
  query_views <- list(
    feature = query_x,
    position = manifoldalign:::.rapid_apply_oos_position_view(
      query_position, position_state$metadata, nrow(query_x)
    ),
    attribute = manifoldalign:::.rapid_apply_oos_attribute_view(
      query_attribute, attribute_state$metadata, nrow(query_x)
    )
  )
  probability <- matrix(
    0, nrow(query_x), ncol(training_probability),
    dimnames = list(NULL, colnames(training_probability))
  )
  for (view in names(structured_oos_weights)) {
    neighbours <- manifoldalign:::.rapid_oos_weights(
      training_views[[view]], query_views[[view]], k,
      sqrt(.Machine$double.eps)
    )
    probability <- probability + structured_oos_weights[[view]] * as.matrix(
      neighbours$interpolation %*% training_probability
    )
  }
  probability / pmax(rowSums(probability), 1e-15)
}

lema_predict_low <- function(fit, new_low) {
  low <- sweep(as.matrix(new_low), 2L, fit$preprocessing$low_center, "-")
  low <- sweep(low, 2L, fit$preprocessing$low_scale, "/")
  high_width <- length(fit$preprocessing$high_center)
  joint <- rbind(
    matrix(0, high_width, nrow(low)),
    t(low)
  )
  embedding <- fit$theta %*% joint
  score <- t(fit$P %*% embedding)
  score <- sweep(score, 1L, apply(score, 1L, max), "-")
  probability <- exp(score)
  probability <- probability / rowSums(probability)
  colnames(probability) <- fit$class_levels
  probability
}

classification_metrics <- function(truth, prediction, probability) {
  truth <- as.character(truth)
  prediction <- as.character(prediction)
  classes <- sort(unique(truth), method = "radix")
  precision <- recall <- f1 <- iou <- stats::setNames(
    numeric(length(classes)), classes
  )
  support <- table(factor(truth, levels = classes))
  for (class in classes) {
    tp <- sum(truth == class & prediction == class)
    fp <- sum(truth != class & prediction == class)
    fn <- sum(truth == class & prediction != class)
    precision[[class]] <- if (tp + fp) tp / (tp + fp) else 0
    recall[[class]] <- if (tp + fn) tp / (tp + fn) else 0
    f1[[class]] <- if (precision[[class]] + recall[[class]] > 0) {
      2 * precision[[class]] * recall[[class]] /
        (precision[[class]] + recall[[class]])
    } else 0
    iou[[class]] <- if (tp + fp + fn) tp / (tp + fp + fn) else 0
  }
  confidence <- apply(probability, 1L, max)
  correct <- prediction == truth
  bins <- pmin(10L, floor(confidence * 10L) + 1L)
  ece <- sum(vapply(seq_len(10L), function(bin) {
    rows <- which(bins == bin)
    if (!length(rows)) return(0)
    length(rows) / length(truth) *
      abs(mean(correct[rows]) - mean(confidence[rows]))
  }, numeric(1)))
  rare_class <- names(which.min(support))[[1L]]
  c(
    oa = mean(correct),
    macro_f1 = mean(f1),
    miou = mean(iou),
    rare_recall = recall[[rare_class]],
    ece = ece
  )
}

make_fixture <- function(scene, split) {
  target_rows <- split$fit_rows
  source_rows <- split$labeled_rows
  paired_target <- match(source_rows, target_rows)
  target_labels <- rep(NA_character_, length(target_rows))
  target_labels[paired_target] <- scene$labels[source_rows]
  list(
    data = list(
      source = scene$high[source_rows, , drop = FALSE],
      target = scene$low[target_rows, , drop = FALSE]
    ),
    truth = list(
      source = scene$labels[source_rows],
      target = scene$labels[target_rows]
    ),
    labels = list(
      source = scene$labels[source_rows],
      target = target_labels
    ),
    positions = list(
      source = scene$position[source_rows, , drop = FALSE],
      target = scene$position[target_rows, , drop = FALSE]
    ),
    attributes = list(
      source = as.data.frame(
        scene$attributes$high[source_rows, , drop = FALSE]
      ),
      target = as.data.frame(
        scene$attributes$low[target_rows, , drop = FALSE]
      )
    ),
    correspondence = paired_target,
    metadata = list(
      split = split$provenance,
      test_labels_read_during_fit = FALSE,
      target_rows = target_rows,
      source_rows = source_rows
    )
  )
}

fit_and_predict <- function(
  method,
  fixture,
  calibration_x,
  test_x,
  calibration_position,
  test_position,
  calibration_attribute,
  test_attribute,
  seed
) {
  target_x <- fixture$data$target
  target_standardizer <- fit_standardizer(target_x)
  target_scaled <- apply_standardizer(target_x, target_standardizer)
  calibration_scaled <- apply_standardizer(calibration_x, target_standardizer)
  test_scaled <- apply_standardizer(test_x, target_standardizer)

  if (method %in% c("target_centroid", "target_structured")) {
    observed <- !is.na(fixture$labels$target)
    model <- fit_centroid(target_scaled[observed, , drop = FALSE],
                          fixture$labels$target[observed])
    if (identical(method, "target_structured")) {
      training_probability <- predict_centroid(model, target_scaled)
      return(list(
        calibration = structured_target_predict(
          target_scaled, fixture$positions$target,
          fixture$attributes$target, training_probability,
          calibration_scaled, calibration_position, calibration_attribute
        ),
        test = structured_target_predict(
          target_scaled, fixture$positions$target,
          fixture$attributes$target, training_probability,
          test_scaled, test_position, test_attribute
        ),
        fit = list(
          centroid = model,
          target_standardizer = target_standardizer,
          structured_oos_weights = structured_oos_weights
        ),
        classes = model$classes,
        oos = paste0(
          "target-only bounded feature/position/attribute smoothing ",
          "(0.50/0.25/0.25)"
        )
      ))
    }
    return(list(
      calibration = predict_centroid(model, calibration_scaled),
      test = predict_centroid(model, test_scaled),
      fit = NULL,
      classes = model$classes,
      oos = "direct scaled target features"
    ))
  }

  benchmark_method <- if (identical(method, "rapid_ma_structured")) {
    "rapid_ma_semisup"
  } else {
    method
  }
  comparison <- benchmark_rapid_ma(
    fixture,
    methods = benchmark_method,
    ncomp = 8L,
    seed = seed,
    keep_fits = TRUE
  )
  if (!identical(comparison$results$status, "ok")) {
    stop(comparison$results$error)
  }
  fit <- comparison$fits[[benchmark_method]]
  if (method %in% c("rapid_ma_semisup", "rapid_ma_structured")) {
    structured <- identical(method, "rapid_ma_structured")
    weights <- if (structured) structured_oos_weights else c(feature = 1)
    return(list(
      calibration = oos_predict(
        fit, calibration_x, side = "target", type = "probabilities", k = 8L,
        new_positions = if (structured) calibration_position else NULL,
        new_attributes = if (structured) calibration_attribute else NULL,
        view_weights = weights
      ),
      test = oos_predict(
        fit, test_x, side = "target", type = "probabilities", k = 8L,
        new_positions = if (structured) test_position else NULL,
        new_attributes = if (structured) test_attribute else NULL,
        view_weights = weights
      ),
      fit = fit,
      classes = fit$prototypes$class_levels,
      oos = if (structured) {
        paste0(
          "bounded feature/position/attribute interpolation of relation-teacher ",
          "state (0.50/0.25/0.25)"
        )
      } else {
        "bounded target-feature interpolation of relation-teacher state"
      }
    ))
  }
  if (identical(method, "lema")) {
    return(list(
      calibration = lema_predict_low(fit, calibration_x),
      test = lema_predict_low(fit, test_x),
      fit = fit,
      classes = fit$class_levels,
      oos = "direct learned low-modality linear projection"
    ))
  }

  sizes <- vapply(fixture$data, nrow, integer(1))
  embedding <- manifoldalign:::.rapid_benchmark_split_scores(
    as.matrix(fit$s), sizes
  )
  observed_embedding <- do.call(rbind, lapply(seq_along(embedding), function(m) {
    embedding[[m]][!is.na(fixture$labels[[m]]), , drop = FALSE]
  }))
  observed_labels <- unlist(lapply(fixture$labels, function(label) {
    label[!is.na(label)]
  }), use.names = FALSE)
  model <- fit_centroid(observed_embedding, observed_labels)
  calibration_embedding <- manifoldalign:::.rapid_inductive_interpolate(
    target_scaled, embedding[[2L]], calibration_scaled, k = 8L
  )$value
  test_embedding <- manifoldalign:::.rapid_inductive_interpolate(
    target_scaled, embedding[[2L]], test_scaled, k = 8L
  )$value
  list(
    calibration = predict_centroid(model, calibration_embedding),
    test = predict_centroid(model, test_embedding),
    fit = fit,
    classes = model$classes,
    oos = "bounded target-feature interpolation of training embedding"
  )
}

sample_process_rss <- function(pid) {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  tryCatch(
    as.numeric(ps::ps_memory_info(ps::ps_handle(pid))[["rss"]]) / 1024^2,
    error = function(...) NA_real_
  )
}

run_method <- function(method, fixture, split, scene, seed) {
  calibration_x <- scene$low[split$calibration_rows, , drop = FALSE]
  test_x <- scene$low[split$test_rows, , drop = FALSE]
  prediction_state <- fit_and_predict(
    method,
    fixture,
    calibration_x,
    test_x,
    scene$position[split$calibration_rows, , drop = FALSE],
    scene$position[split$test_rows, , drop = FALSE],
    as.data.frame(scene$attributes$low[split$calibration_rows, , drop = FALSE]),
    as.data.frame(scene$attributes$low[split$test_rows, , drop = FALSE]),
    seed
  )

  # Calibration labels are accessed only after fitting and both OOS prediction
  # matrices have been generated. Test labels are accessed only below.
  calibration_truth <- scene$labels[split$calibration_rows]
  represented_calibration <- calibration_truth %in%
    colnames(prediction_state$calibration)
  if (sum(represented_calibration) >= 2L) {
    calibration <- manifoldalign:::.rapid_temperature_fit(
      prediction_state$calibration[represented_calibration, , drop = FALSE],
      calibration_truth[represented_calibration]
    )
  } else {
    calibration <- list(
      temperature = 1,
      nll_before = NA_real_,
      nll_after = NA_real_,
      n = sum(represented_calibration),
      provenance = list(training_only = TRUE, test_labels_read = FALSE)
    )
  }
  calibrated_test <- manifoldalign:::.rapid_temperature_apply(
    prediction_state$test, calibration$temperature
  )
  uncalibrated_prediction <- colnames(prediction_state$test)[
    max.col(prediction_state$test, ties.method = "first")
  ]
  prediction <- colnames(calibrated_test)[
    max.col(calibrated_test, ties.method = "first")
  ]

  test_truth <- scene$labels[split$test_rows]
  overall <- classification_metrics(test_truth, prediction, calibrated_test)
  uncalibrated <- classification_metrics(
    test_truth, uncalibrated_prediction, prediction_state$test
  )
  seen <- test_truth %in% prediction_state$classes
  seen_metrics <- if (any(seen)) {
    classification_metrics(
      test_truth[seen], prediction[seen], calibrated_test[seen, , drop = FALSE]
    )
  } else {
    stats::setNames(rep(NA_real_, 5L), names(overall))
  }
  list(
    row = data.frame(
      method = method,
      status = "ok",
      oos = prediction_state$oos,
      n_fit = length(split$fit_rows),
      n_source = length(split$labeled_rows),
      n_calibration = length(split$calibration_rows),
      n_test = length(split$test_rows),
      n_classes_train = length(prediction_state$classes),
      n_classes_test = length(unique(test_truth)),
      unseen_fraction = mean(!seen),
      oa = overall[["oa"]],
      macro_f1 = overall[["macro_f1"]],
      miou = overall[["miou"]],
      rare_recall = overall[["rare_recall"]],
      ece = overall[["ece"]],
      ece_uncalibrated = uncalibrated[["ece"]],
      seen_oa = seen_metrics[["oa"]],
      seen_macro_f1 = seen_metrics[["macro_f1"]],
      seen_miou = seen_metrics[["miou"]],
      temperature = calibration$temperature,
      calibration_n = calibration$n,
      calibration_nll_before = calibration$nll_before,
      calibration_nll_after = calibration$nll_after,
      test_labels_read_during_fit = FALSE,
      uses_source_alignment = method %in% c(
        "rapid_ma_semisup", "rapid_ma_structured", "lema", "ssma", "kema"
      ),
      uses_oos_structure = method %in% c(
        "target_structured", "rapid_ma_structured"
      ),
      error = NA_character_,
      stringsAsFactors = FALSE
    ),
    retained_mb = as.numeric(utils::object.size(prediction_state$fit)) / 1024^2
  )
}

failure_row <- function(method, message) {
  data.frame(
    method = method, status = "failed", oos = NA_character_,
    n_fit = NA_integer_, n_source = NA_integer_,
    n_calibration = NA_integer_, n_test = NA_integer_,
    n_classes_train = NA_integer_, n_classes_test = NA_integer_,
    unseen_fraction = NA_real_, oa = NA_real_, macro_f1 = NA_real_,
    miou = NA_real_, rare_recall = NA_real_, ece = NA_real_,
    ece_uncalibrated = NA_real_, seen_oa = NA_real_,
    seen_macro_f1 = NA_real_, seen_miou = NA_real_,
    temperature = NA_real_, calibration_n = NA_integer_,
    calibration_nll_before = NA_real_, calibration_nll_after = NA_real_,
    test_labels_read_during_fit = FALSE,
    uses_source_alignment = method %in% c(
      "rapid_ma_semisup", "rapid_ma_structured", "lema", "ssma", "kema"
    ),
    uses_oos_structure = method %in% c(
      "target_structured", "rapid_ma_structured"
    ),
    error = message,
    stringsAsFactors = FALSE
  )
}

run_isolated_method <- function(method, fixture, split, scene, seed, timeout_sec) {
  start <- proc.time()[["elapsed"]]
  job <- parallel::mcparallel(
    run_method(method, fixture, split, scene, seed),
    name = paste(method, split$test_block, seed, sep = "-"),
    mc.set.seed = FALSE,
    silent = TRUE
  )
  initial <- sample_process_rss(job$pid)
  peak <- initial
  value <- NULL
  repeat {
    collected <- parallel::mccollect(job, wait = FALSE)
    current <- sample_process_rss(job$pid)
    if (is.finite(current)) peak <- max(peak, current, na.rm = TRUE)
    if (!is.null(collected)) {
      value <- collected[[1L]]
      break
    }
    elapsed <- proc.time()[["elapsed"]] - start
    if (elapsed >= timeout_sec) {
      parallel:::mckill(job, signal = 15L)
      parallel::mccollect(job, wait = TRUE, timeout = 1)
      row <- failure_row(
        method,
        paste0("method exceeded isolated ", timeout_sec, " second timeout")
      )
      row$runtime_sec <- elapsed
      row$retained_mb <- NA_real_
      row$sampled_peak_rss_mb <- peak
      return(row)
    }
    Sys.sleep(0.02)
  }
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(value, "try-error") || !is.list(value) || is.null(value$row)) {
    row <- failure_row(
      method,
      if (inherits(value, "try-error")) as.character(value) else
        "isolated method returned no result"
    )
    row$runtime_sec <- elapsed
    row$retained_mb <- NA_real_
  } else {
    row <- value$row
    row$runtime_sec <- elapsed
    row$retained_mb <- value$retained_mb
  }
  row$sampled_peak_rss_mb <- peak
  row
}

grid <- c(4L, 4L)
block_id <- manifoldalign:::.rapid_spatial_blocks(scene$position, grid)
present_blocks <- sort(unique(block_id), method = "radix")
default_count <- switch(profile, quick = 2L, full = 5L, exhaustive = length(present_blocks))
test_blocks <- utils::head(present_blocks, default_count)
if (nzchar(Sys.getenv("RAPID_MA_BLOCKS"))) {
  test_blocks <- as.integer(strsplit(
    Sys.getenv("RAPID_MA_BLOCKS"), ",", fixed = TRUE
  )[[1L]])
}
if (!length(test_blocks) || any(!test_blocks %in% present_blocks)) {
  stop("RAPID_MA_BLOCKS must name present spatial block IDs.")
}

max_fit_rows <- as.integer(Sys.getenv("RAPID_MA_MAX_FIT_ROWS", "1500"))
labels_per_class <- as.integer(Sys.getenv("RAPID_MA_LABELS_PER_CLASS", "4"))
calibration_per_class <- as.integer(
  Sys.getenv("RAPID_MA_CALIBRATION_PER_CLASS", "4")
)
timeout_sec <- as.numeric(Sys.getenv("RAPID_MA_METHOD_TIMEOUT", "120"))
available_methods <- c(
  "target_centroid", "target_structured",
  "rapid_ma_semisup", "rapid_ma_structured",
  "lema", "ssma", "kema"
)
methods <- available_methods
if (nzchar(Sys.getenv("RAPID_MA_METHODS"))) {
  methods <- strsplit(Sys.getenv("RAPID_MA_METHODS"), ",", fixed = TRUE)[[1L]]
  if (!length(methods) || any(!methods %in% available_methods) ||
      anyDuplicated(methods)) {
    stop(
      "RAPID_MA_METHODS must be a unique comma-separated subset of: ",
      paste(available_methods, collapse = ", "), "."
    )
  }
}

rows <- list()
index <- 0L
for (fold in seq_along(test_blocks)) {
  test_block <- test_blocks[[fold]]
  test_position <- match(test_block, present_blocks)
  calibration_block <- present_blocks[[
    (test_position %% length(present_blocks)) + 1L
  ]]
  seed <- 5100L + test_block
  split <- manifoldalign:::.rapid_blocked_inductive_split(
    scene$labels,
    scene$position,
    test_block = test_block,
    calibration_block = calibration_block,
    grid = grid,
    labels_per_class = labels_per_class,
    calibration_per_class = calibration_per_class,
    max_fit_rows = max_fit_rows,
    seed = seed
  )
  fixture <- make_fixture(scene, split)
  method_rows <- lapply(methods, function(method) {
    suppressWarnings(run_isolated_method(
      method, fixture, split, scene, seed, timeout_sec
    ))
  })
  block <- do.call(rbind, method_rows)
  rownames(block) <- NULL
  block$fold <- fold
  block$test_block <- test_block
  block$calibration_block <- calibration_block
  block$seed <- seed
  block$labels_per_class <- labels_per_class
  block$calibration_per_class <- calibration_per_class
  block$max_fit_rows <- max_fit_rows
  index <- index + 1L
  rows[[index]] <- block
  message(
    "completed fold=", fold,
    " test_block=", test_block,
    " failures=", sum(block$status != "ok")
  )
}

raw <- do.call(rbind, rows)
rownames(raw) <- NULL
metric_names <- c(
  "runtime_sec", "retained_mb", "sampled_peak_rss_mb", "unseen_fraction",
  "oa", "macro_f1", "miou", "rare_recall", "ece", "ece_uncalibrated",
  "seen_oa", "seen_macro_f1", "seen_miou", "temperature",
  "calibration_n", "calibration_nll_before", "calibration_nll_after"
)
groups <- split(seq_len(nrow(raw)), raw$method)
summary_rows <- lapply(groups, function(row_index) {
  block <- raw[row_index, , drop = FALSE]
  output <- data.frame(
    method = block$method[[1L]],
    n_folds = nrow(block),
    n_ok = sum(block$status == "ok"),
    n_failed = sum(block$status != "ok"),
    stringsAsFactors = FALSE
  )
  for (metric in metric_names) {
    value <- block[[metric]][block$status == "ok"]
    output[[paste0(metric, "_mean")]] <- if (length(value)) {
      mean(value, na.rm = TRUE)
    } else NA_real_
    output[[paste0(metric, "_sd")]] <- if (sum(is.finite(value)) > 1L) {
      stats::sd(value, na.rm = TRUE)
    } else NA_real_
  }
  output
})
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL
summary <- summary[order(-summary$seen_macro_f1_mean, summary$method), ]

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
raw_path <- paste0(output_prefix, "-raw.csv")
summary_path <- paste0(output_prefix, "-summary.csv")
metadata_path <- paste0(output_prefix, "-metadata.txt")
utils::write.csv(raw, raw_path, row.names = FALSE)
utils::write.csv(summary, summary_path, row.names = FALSE)

source_files <- c(
  "R/rapid_ma_validation.R", "R/rapid_ma_benchmark.R", "R/rapid_ma.R",
  "R/rapid_ma_semisup.R", "inst/benchmarks/lema_reference.R",
  "inst/benchmarks/rapid_ma_blocked_inductive.R"
)
source_hash <- unname(tools::md5sum(source_files[file.exists(source_files)]))
metadata <- c(
  "benchmark=RAPID-MA Chikusei real-scene blocked-inductive comparison",
  paste0("profile=", profile),
  paste0("test_blocks=", paste(test_blocks, collapse = ",")),
  paste0("grid=", paste(grid, collapse = "x")),
  paste0("max_fit_rows=", max_fit_rows),
  paste0("labels_per_class=", labels_per_class),
  paste0("calibration_per_class=", calibration_per_class),
  paste0("method_timeout_sec=", timeout_sec),
  paste0("methods=", paste(methods, collapse = ",")),
  "split=strict: test and calibration tiles absent from fit and preprocessing",
  "calibration=temperature fitted only on a separate calibration tile",
  paste0(
    "structured_oos_weights=",
    paste(names(structured_oos_weights), structured_oos_weights,
          sep = ":", collapse = ",")
  ),
  "structured_oos_weights_frozen_before_heldout_block_evaluation=true",
  "target_modality=frozen ten-band Gaussian approximation to Sentinel-2 VNIR",
  "real_scene=true; independent_target_sensor=false",
  "test_labels_read_during_fit=false",
  "license=CC BY 3.0",
  paste0("dataset_archive_sha256=", scene$metadata$archive_sha256),
  paste0("R.version=", R.version.string),
  paste0("platform=", R.version$platform),
  paste0(names(source_hash), " md5=", source_hash)
)
writeLines(metadata, metadata_path)

print(summary[, c(
  "method", "n_ok", "oa_mean", "macro_f1_mean", "miou_mean",
  "seen_oa_mean", "seen_macro_f1_mean", "seen_miou_mean",
  "ece_mean", "runtime_sec_mean"
)], row.names = FALSE)
