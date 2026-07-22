#!/usr/bin/env Rscript

# Reproducible downloader for the Chikusei real-scene benchmark.
# The archive is not bundled with manifoldalign. It is provided by the Space
# Application Laboratory, University of Tokyo under CC BY 3.0.

args <- commandArgs(trailingOnly = TRUE)
destination <- if (length(args)) args[[1L]] else {
  Sys.getenv(
    "RAPID_MA_DATA_DIR",
    file.path(tempdir(), "rapid-ma-data")
  )
}
accepted <- tolower(Sys.getenv("RAPID_MA_ACCEPT_CHIKUSEI_LICENSE", ""))
if (!accepted %in% c("yes", "true", "1")) {
  stop(
    "Set RAPID_MA_ACCEPT_CHIKUSEI_LICENSE=yes after reviewing ",
    "https://naotoyokoya.com/Download.html and CC BY 3.0."
  )
}

url <- paste0(
  "https://www.sal.t.u-tokyo.ac.jp/hyperdata/",
  "Hyperspec_Chikusei_ENVI.zip"
)
expected_bytes <- 1066530806
expected_sha256 <-
  "192825c63e764a16e9324470e005992f343c24317178e8a5546147aadde838dc"
archive <- file.path(destination, "Hyperspec_Chikusei_ENVI.zip")
scene_dir <- file.path(destination, "Chikusei_ENVI")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)

valid_archive <- file.exists(archive) &&
  identical(as.numeric(file.info(archive)$size), expected_bytes)
if (!valid_archive) {
  message("Downloading 1.0 GB Chikusei ENVI archive from the official host...")
  utils::download.file(url, archive, mode = "wb", method = "libcurl")
}
if (!identical(as.numeric(file.info(archive)$size), expected_bytes)) {
  stop("Downloaded archive has an unexpected byte length: ", archive)
}

shasum <- Sys.which("shasum")
if (nzchar(shasum)) {
  observed <- system2(shasum, c("-a", "256", shQuote(archive)), stdout = TRUE)
  observed <- strsplit(observed[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!identical(observed, expected_sha256)) {
    stop("Downloaded archive failed SHA-256 verification: ", archive)
  }
}

required <- file.path(
  scene_dir,
  c(
    "HyperspecVNIR_Chikusei_20140729.bsq",
    "HyperspecVNIR_Chikusei_20140729.hdr",
    "HyperspecVNIR_Chikusei_20140729_Ground_Truth.bsq",
    "HyperspecVNIR_Chikusei_20140729_Ground_Truth.hdr",
    "Readme.txt"
  )
)
if (!all(file.exists(required))) {
  utils::unzip(archive, exdir = destination)
}
if (!all(file.exists(required))) {
  stop("The Chikusei archive did not contain the expected ENVI files.")
}

message("Chikusei scene ready at: ", normalizePath(scene_dir))
message(
  "Required acknowledgement: The authors gratefully acknowledge Space ",
  "Application Laboratory, Department of Advanced Interdisciplinary Studies, ",
  "the University of Tokyo for providing the hyperspectral data."
)
message(
  "Citation: N. Yokoya and A. Iwasaki, Airborne hyperspectral data over ",
  "Chikusei, SAL-2016-05-27, 2016."
)
