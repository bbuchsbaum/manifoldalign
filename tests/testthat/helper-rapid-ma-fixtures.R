make_tiny_rapid_benchmark_fixture <- function(seed = 17L) {
  rapid_ma_benchmark_fixture(
    n_common = 36L,
    extra_target = 6L,
    source_features = 10L,
    target_features = 7L,
    labels_per_class = 2L,
    seed = seed
  )
}

rapid_benchmark_held_out <- function(fixture) {
  is.na(fixture$labels[[2L]])
}
