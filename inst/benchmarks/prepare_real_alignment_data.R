#!/usr/bin/env Rscript

# Acquire and verify the governed external inputs used by the real-data
# alignment benchmark. Network fixtures are small and redistributable under
# Apache-2.0. So2Sat is large (about 18.2 GB compressed); it is downloaded only
# after an explicit license acknowledgement and is never added to the package.

find_package_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(
          path, "inst", "benchmarks", "real_alignment_artifacts.csv"
        ))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the manifoldalign package root.", call. = FALSE)
    }
    path <- parent
  }
}

sha256_file <- function(path) {
  sha256sum <- Sys.which("sha256sum")
  shasum <- Sys.which("shasum")
  if (nzchar(sha256sum)) {
    out <- system2(sha256sum, path, stdout = TRUE, stderr = TRUE)
  } else if (nzchar(shasum)) {
    out <- system2(shasum, c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  } else {
    stop("A `sha256sum` or `shasum` executable is required.", call. = FALSE)
  }
  if (!length(out)) stop("SHA-256 command returned no output.", call. = FALSE)
  tolower(strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]])
}

download_verified <- function(url, path, bytes, sha256) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path) || file.info(path)$size != bytes ||
      !identical(sha256_file(path), sha256)) {
    tmp <- paste0(path, ".partial")
    curl <- Sys.which("curl")
    if (nzchar(curl)) {
      status <- system2(curl, c(
        "--location", "--fail", "--retry", "5", "--retry-delay", "2",
        "--continue-at", "-", "--progress-bar", "--output", tmp, url
      ))
      if (!identical(status, 0L)) {
        stop(
          "Resumable download failed for ", basename(path),
          "; partial data remain at ", tmp, ".", call. = FALSE
        )
      }
    } else {
      old_timeout <- getOption("timeout")
      on.exit(options(timeout = old_timeout), add = TRUE)
      options(timeout = max(7200, old_timeout))
      utils::download.file(url, tmp, mode = "wb", quiet = FALSE)
    }
    observed_bytes <- unname(file.info(tmp)$size)
    observed_sha <- sha256_file(tmp)
    if (!identical(as.numeric(observed_bytes), as.numeric(bytes)) ||
        !identical(observed_sha, sha256)) {
      if (is.finite(observed_bytes) && observed_bytes >= bytes) {
        unlink(tmp)
      }
      stop(
        "Artifact verification failed for ", basename(path),
        ": expected ", bytes, " bytes / ", sha256,
        ", observed ", observed_bytes, " bytes / ", observed_sha, ".",
        call. = FALSE
      )
    }
    if (!file.rename(tmp, path)) {
      stop("Could not move verified artifact into the cache.", call. = FALSE)
    }
  }
  invisible(path)
}

parse_args <- function(args, root) {
  output <- Sys.getenv(
    "MANIFOLDALIGN_REAL_DATA_DIR",
    unset = file.path(root, ".cache", "alignment-real")
  )
  download_so2sat <- FALSE
  so2sat_artifacts <- "all"
  extract_networks_to <- ""
  for (arg in args) {
    if (identical(arg, "--download-so2sat")) {
      download_so2sat <- TRUE
    } else if (startsWith(arg, "--output-dir=")) {
      output <- sub("^--output-dir=", "", arg)
    } else if (startsWith(arg, "--so2sat-artifacts=")) {
      so2sat_artifacts <- strsplit(
        sub("^--so2sat-artifacts=", "", arg), ",", fixed = TRUE
      )[[1L]]
    } else if (startsWith(arg, "--extract-network-fixtures-to=")) {
      extract_networks_to <- sub(
        "^--extract-network-fixtures-to=", "", arg
      )
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  list(
    output = normalizePath(output, mustWork = FALSE),
    download_so2sat = download_so2sat,
    so2sat_artifacts = so2sat_artifacts,
    extract_networks_to = extract_networks_to
  )
}

root <- find_package_root()
args <- parse_args(commandArgs(trailingOnly = TRUE), root)
artifacts <- utils::read.csv(
  file.path(root, "inst", "benchmarks", "real_alignment_artifacts.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)

network_source <- artifacts[artifacts$artifact_id == "multinet_source", ]
network_archive <- file.path(args$output, network_source$relative_path)
download_verified(
  network_source$source_url, network_archive, network_source$bytes,
  network_source$sha256
)

if (nzchar(args$extract_networks_to)) {
  destination <- normalizePath(args$extract_networks_to, mustWork = FALSE)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  wanted <- c(
    "multinet/inst/extdata/aucs.mpx",
    "multinet/inst/extdata/monastery.mpx"
  )
  scratch <- tempfile("manifoldalign-multinet-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  utils::untar(network_archive, files = wanted, exdir = scratch)
  for (source in wanted) {
    target <- file.path(destination, basename(source))
    if (!file.copy(file.path(scratch, source), target, overwrite = TRUE)) {
      stop("Could not extract ", basename(source), ".", call. = FALSE)
    }
  }
}

status <- data.frame(
  dataset_id = "multinet",
  artifact_id = "multinet_source",
  status = "available",
  status_code = "verified",
  path = normalizePath(network_archive),
  stringsAsFactors = FALSE
)

so2sat <- artifacts[
  artifacts$dataset_id == "so2sat_lcz42" &
    artifacts$redistribution == "cache only",
  , drop = FALSE
]
if (isTRUE(args$download_so2sat)) {
  accepted <- tolower(Sys.getenv("MANIFOLDALIGN_ACCEPT_SO2SAT_LICENSE", ""))
  if (!accepted %in% c("1", "true", "yes")) {
    stop(
      "Set MANIFOLDALIGN_ACCEPT_SO2SAT_LICENSE=yes after reviewing CC BY 4.0 ",
      "and the Copernicus data terms.", call. = FALSE
    )
  }
  if (!identical(args$so2sat_artifacts, "all")) {
    unknown <- setdiff(args$so2sat_artifacts, so2sat$artifact_id)
    if (length(unknown)) {
      stop(
        "Unknown So2Sat artifact(s): ", paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    so2sat <- so2sat[so2sat$artifact_id %in% args$so2sat_artifacts, ]
  }
  for (i in seq_len(nrow(so2sat))) {
    path <- file.path(args$output, so2sat$relative_path[[i]])
    download_verified(
      so2sat$source_url[[i]], path, so2sat$bytes[[i]], so2sat$sha256[[i]]
    )
    status <- rbind(status, data.frame(
      dataset_id = "so2sat_lcz42",
      artifact_id = so2sat$artifact_id[[i]],
      status = "available",
      status_code = "verified",
      path = normalizePath(path),
      stringsAsFactors = FALSE
    ))
  }
} else {
  status <- rbind(status, data.frame(
    dataset_id = "so2sat_lcz42",
    artifact_id = so2sat$artifact_id,
    status = "dataset_unavailable",
    status_code = "download_not_requested",
    path = file.path(args$output, so2sat$relative_path),
    stringsAsFactors = FALSE
  ))
}

utils::write.csv(
  status, file.path(args$output, "dataset-status.csv"), row.names = FALSE,
  na = ""
)
message("Real-data acquisition status written to ", args$output)
