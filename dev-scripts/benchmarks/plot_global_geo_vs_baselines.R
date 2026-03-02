#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
results_csv <- if (length(args) >= 1L) args[[1L]] else "dev-scripts/benchmarks/global_geo_vs_baselines_results.csv"
summary_csv <- if (length(args) >= 2L) args[[2L]] else "dev-scripts/benchmarks/global_geo_vs_baselines_summary.csv"
out_pdf <- if (length(args) >= 3L) args[[3L]] else "dev-scripts/benchmarks/global_geo_vs_baselines_plots.pdf"

if (!file.exists(results_csv)) stop("Missing results CSV: ", results_csv, call. = FALSE)
if (!file.exists(summary_csv)) stop("Missing summary CSV: ", summary_csv, call. = FALSE)

results <- read.csv(results_csv, stringsAsFactors = FALSE)
summary_tbl <- read.csv(summary_csv, stringsAsFactors = FALSE)

scenario_pref <- c("clean", "sparse_supervision", "shifted_noisy", "gaussian_classes", "gaussian_overlap")
scenario_levels <- c(scenario_pref[scenario_pref %in% unique(summary_tbl$scenario)],
                     setdiff(sort(unique(summary_tbl$scenario)), scenario_pref))
supervision_levels <- c("id", "condition", "none")
metric_dirs <- c(
  identity_mse_mean = "min",
  holdout_identity_mse_mean = "min",
  geometry_spearman_mean = "max",
  condition_knn_acc_mean = "max",
  holdout_condition_knn_acc_mean = "max",
  runtime_mean = "min"
)
metric_dirs <- metric_dirs[names(metric_dirs) %in% names(summary_tbl)]

summary_tbl <- summary_tbl %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    supervision = factor(supervision, levels = supervision_levels)
  )

results <- results %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    supervision = factor(supervision, levels = supervision_levels)
  )

method_levels <- summary_tbl %>%
  group_by(method) %>%
  summarise(overall_identity = mean(identity_mse_mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(overall_identity) %>%
  pull(method)

summary_tbl$method <- factor(summary_tbl$method, levels = method_levels)
results$method <- factor(results$method, levels = method_levels)

method_cols <- setNames(
  grDevices::hcl.colors(length(method_levels), palette = "Dark 3"),
  method_levels
)

base_theme <- theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.position = "none"
  )

plot_metric <- function(metric_mean, metric_sd, title, subtitle, log_scale = FALSE, ylab = NULL) {
  pdat <- summary_tbl %>%
    mutate(
      value = .data[[metric_mean]],
      spread = .data[[metric_sd]],
      ymin = value - spread,
      ymax = value + spread
    )

  if (log_scale) {
    eps <- 1e-8
    pdat <- pdat %>%
      mutate(
        value = pmax(value, eps),
        ymin = pmax(ymin, eps),
        ymax = pmax(ymax, eps)
      )
  }

  p <- ggplot(pdat, aes(x = method, y = value, color = method)) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, alpha = 0.55, linewidth = 0.45) +
    geom_point(size = 2.3, alpha = 0.95) +
    facet_grid(supervision ~ scenario, scales = "free_y", drop = TRUE) +
    scale_color_manual(values = method_cols) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = if (is.null(ylab)) metric_mean else ylab
    ) +
    base_theme

  if (log_scale) {
    p <- p + scale_y_log10()
  }
  p
}

rank_long <- summary_tbl %>%
  mutate(setting = paste(as.character(scenario), as.character(supervision), sep = " | ")) %>%
  select(method, setting, all_of(names(metric_dirs))) %>%
  pivot_longer(cols = all_of(names(metric_dirs)), names_to = "metric", values_to = "value") %>%
  group_by(setting, metric) %>%
  mutate(
    rank = if (metric_dirs[metric[1]] == "min") {
      rank(value, ties.method = "min")
    } else {
      rank(-value, ties.method = "min")
    }
  ) %>%
  ungroup() %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "identity_mse_mean",
        "holdout_identity_mse_mean",
        "geometry_spearman_mean",
        "condition_knn_acc_mean",
        "holdout_condition_knn_acc_mean",
        "runtime_mean"
      ),
      labels = c(
        "Identity MSE (lower is better)",
        "Holdout Identity MSE (lower is better)",
        "Geometry Spearman (higher is better)",
        "Condition kNN Accuracy (higher is better)",
        "Holdout Condition kNN Accuracy (higher is better)",
        "Runtime (lower is better)"
      )
    ),
    setting = factor(
      setting,
      levels = as.vector(outer(scenario_levels, supervision_levels, paste, sep = " | "))
    )
  )

p_title <- ggplot() +
  annotate("text", x = 0.02, y = 0.92, hjust = 0, label = "Global Geo Benchmark: Full Results", size = 10, fontface = "bold") +
  annotate("text", x = 0.02, y = 0.80, hjust = 0, label = sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), size = 4) +
  annotate(
    "text",
    x = 0.02,
    y = 0.72,
    hjust = 0,
    label = sprintf(
      "Input rows: %d per-rep results, %d summary rows",
      nrow(results), nrow(summary_tbl)
    ),
    size = 4
  ) +
  annotate(
    "text",
    x = 0.02,
    y = 0.64,
    hjust = 0,
    label = sprintf("Methods: %s", paste(levels(summary_tbl$method), collapse = ", ")),
    size = 3.5
  ) +
  annotate(
    "text",
    x = 0.02,
    y = 0.56,
    hjust = 0,
    label = sprintf("Scenarios: %s", paste(scenario_levels, collapse = ", ")),
    size = 3.5
  ) +
  annotate(
    "text",
    x = 0.02,
    y = 0.48,
    hjust = 0,
    label = sprintf("Supervision modes: %s", paste(supervision_levels, collapse = ", ")),
    size = 3.5
  ) +
  annotate(
    "text",
    x = 0.02,
    y = 0.36,
    hjust = 0,
    label = "Panels use mean ± sd across reps from summary CSV.\nRuntime and MSE panels are shown on log scale where appropriate.",
    size = 3.8
  ) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(20, 20, 20, 20))

p_identity <- plot_metric(
  "identity_mse_mean", "identity_mse_sd",
  title = "Identity Alignment Error",
  subtitle = "Mean ± sd across reps; lower is better (log scale).",
  log_scale = TRUE,
  ylab = "Identity MSE"
)

p_holdout <- plot_metric(
  "holdout_identity_mse_mean", "holdout_identity_mse_sd",
  title = "Holdout Identity Alignment Error",
  subtitle = "Mean ± sd across reps; lower is better (log scale).",
  log_scale = TRUE,
  ylab = "Holdout Identity MSE"
)

p_geo <- plot_metric(
  "geometry_spearman_mean", "geometry_spearman_sd",
  title = "Within-Domain Geometry Preservation",
  subtitle = "Mean ± sd Spearman correlation of pairwise distances; higher is better.",
  log_scale = FALSE,
  ylab = "Geometry Spearman"
)

p_cond <- plot_metric(
  "condition_knn_acc_mean", "condition_knn_acc_sd",
  title = "Condition kNN Accuracy",
  subtitle = "Mean ± sd pooled-domain classification accuracy; higher is better.",
  log_scale = FALSE,
  ylab = "Condition kNN Accuracy"
)

p_holdout_cond <- NULL
if ("holdout_condition_knn_acc_mean" %in% names(summary_tbl) &&
    "holdout_condition_knn_acc_sd" %in% names(summary_tbl)) {
  p_holdout_cond <- plot_metric(
    "holdout_condition_knn_acc_mean", "holdout_condition_knn_acc_sd",
    title = "Holdout Condition kNN Accuracy",
    subtitle = "Classification on held-out condition labels only; higher is better.",
    log_scale = FALSE,
    ylab = "Holdout Condition kNN Accuracy"
  )
}

p_runtime <- ggplot(results, aes(x = method, y = runtime_sec, color = method)) +
  geom_boxplot(outlier.alpha = 0.25, alpha = 0.20, width = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.40, size = 1.2) +
  facet_grid(supervision ~ scenario, scales = "free_y", drop = TRUE) +
  scale_color_manual(values = method_cols) +
  scale_y_log10() +
  labs(
    title = "Runtime Distribution Across Repetitions",
    subtitle = "Per-run runtime (seconds), log scale.",
    x = NULL,
    y = "Runtime (sec)"
  ) +
  base_theme

p_rank <- ggplot(rank_long, aes(x = setting, y = method, fill = rank)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = rank), size = 3.0, color = "white", fontface = "bold") +
  facet_wrap(~ metric, ncol = 1) +
  scale_fill_gradient(low = "#0b7285", high = "#d9480f") +
  labs(
    title = "Method Ranks by Scenario/Supervision",
    subtitle = "Rank 1 is best for each metric and setting.",
    x = "Scenario | Supervision",
    y = NULL
  ) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(out_pdf, width = 18, height = 11, onefile = TRUE)
print(p_title)
print(p_identity)
print(p_holdout)
print(p_geo)
print(p_cond)
if (!is.null(p_holdout_cond)) print(p_holdout_cond)
print(p_runtime)
print(p_rank)
grDevices::dev.off()

cat("Saved plot report to:", out_pdf, "\n")
