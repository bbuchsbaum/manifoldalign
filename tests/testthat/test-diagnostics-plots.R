test_that("plot helpers return ggplot objects and summarize_bad_edges works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  diagnostics <- list(
    edge_residuals = data.frame(i = c(1, 1, 2), j = c(2, 3, 3), residual = c(0.1, 0.3, 0.2)),
    edge_weights   = data.frame(i = c(1, 1, 2), j = c(2, 3, 3), weight = c(0.9, 0.8, 0.85)),
    cycle_residuals = data.frame(i = 1, j = 2, k = 3, residual = 0.25)
  )

  p1 <- plot_edge_residuals_heatmap(diagnostics)
  p2 <- plot_cycle_consistency(diagnostics, aggregate = "min")
  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))

  bad <- summarize_bad_edges(diagnostics, top = 2)
  expect_true(is.data.frame(bad))
  expect_equal(nrow(bad), 2)
  expect_true(all(c("i", "j", "residual", "weight") %in% names(bad)))
  # ensure sorted by residual desc
  expect_gte(bad$residual[1], bad$residual[2])
})

