## Centralized benchmark/heavy test toggle helper
##
## Preferred control: options(manifoldalign.run_benchmarks = TRUE)
## Backward-compat: legacy env vars still enable when set to truthy values.

.truthy <- function(x) {
  val <- tolower(as.character(x %||% ""))
  val %in% c("1", "true", "t", "yes", "y")
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || identical(a, "")) b else a

manifoldalign_benchmarks_enabled <- function() {
  # Preferred option
  if (isTRUE(getOption("manifoldalign.run_benchmarks", FALSE))) return(TRUE)

  # Legacy env vars (kept for compatibility)
  legacy_env <- c(
    Sys.getenv("WITH_PARROT_PERF_TESTS", ""),
    Sys.getenv("WITH_SCALABLE_TESTS", "")
  )
  any(vapply(legacy_env, .truthy, logical(1)))
}

skip_if_benchmarks_disabled <- function(msg = "Benchmark tests disabled; enable with options(manifoldalign.run_benchmarks = TRUE)") {
  if (!manifoldalign_benchmarks_enabled()) testthat::skip(msg)
}

