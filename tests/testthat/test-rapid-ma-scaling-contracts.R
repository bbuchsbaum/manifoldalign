test_that("RAPID retained state and sparse supports remain bounded", {
  fits <- lapply(c(120L, 240L), function(n) {
    fixture <- rapid_ma_benchmark_fixture(
      n_common = n,
      extra_target = as.integer(n / 5L),
      labels_per_class = 3L,
      seed = 71L
    )
    rapid_ma(
      fixture$data,
      labels = fixture$labels,
      positions = fixture$positions,
      attributes = fixture$attributes,
      ncomp = 6L,
      control = list(max_iter = 1L, min_iter = 1L, seed = 71L)
    )
  })

  sizes <- vapply(fits, function(x) as.numeric(utils::object.size(x)), numeric(1))
  expect_lt(sizes[[2L]] / sizes[[1L]], 3)
  for (fit in fits) {
    for (m in seq_along(fit$domain_sizes)) {
      n <- fit$domain_sizes[[m]]
      q <- fit$control$uot$q
      expect_lte(Matrix::nnzero(fit$couplings[[m]]),
                 q * n + nrow(fit$prototypes$embedding))
      for (relation in fit$relations[[m]]$relations) {
        expect_lte(Matrix::nnzero(relation),
                   fit$control$relation$degree_cap * n)
      }
    }
  }
})

test_that("scalability fixture size does not change the frozen label budget", {
  small <- rapid_ma_benchmark_fixture(n_common = 120L, labels_per_class = 4L)
  large <- rapid_ma_benchmark_fixture(n_common = 1200L, labels_per_class = 4L)

  expect_equal(sum(!is.na(small$labels[[1L]])), 12L)
  expect_equal(sum(!is.na(large$labels[[1L]])), 12L)
  expect_equal(sum(!is.na(small$labels[[2L]])), 12L)
  expect_equal(sum(!is.na(large$labels[[2L]])), 12L)
})

test_that("multi-view matcher retains only a bounded linear candidate graph", {
  matchings <- lapply(c(120L, 240L), function(n) {
    fixture <- rapid_ma_benchmark_fixture(
      n_common = n,
      extra_target = as.integer(n / 5L),
      labels_per_class = 3L,
      seed = 73L
    )
    fit <- rapid_ma(
      fixture$data,
      labels = fixture$labels,
      positions = fixture$positions,
      attributes = fixture$attributes,
      ncomp = 6L,
      control = list(max_iter = 1L, min_iter = 1L, seed = 73L)
    )
    rapid_ma_match(
      fit,
      assignment = "independent",
      candidate_cap = 16L,
      view_candidate_k = 6L
    )
  })

  for (matching in matchings) {
    expect_lte(
      matching$candidate_edges,
      nrow(matching$matches) * matching$candidate_cap
    )
    expect_equal(matching$coverage, 1)
    expect_false(matching$dense_pairwise_allocated)
  }
  retained <- vapply(
    matchings, function(x) as.numeric(utils::object.size(x)), numeric(1)
  )
  expect_lt(retained[[2L]] / retained[[1L]], 3)
})
