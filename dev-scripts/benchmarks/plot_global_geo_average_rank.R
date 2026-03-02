#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
summary_csv <- if (length(args) >= 1L) args[[1L]] else "dev-scripts/benchmarks/global_geo_vs_baselines_summary.csv"
out_pdf <- if (length(args) >= 2L) args[[2L]] else "dev-scripts/benchmarks/global_geo_vs_baselines_average_rank.pdf"
out_csv <- if (length(args) >= 3L) args[[3L]] else "dev-scripts/benchmarks/global_geo_vs_baselines_average_rank.csv"

if (!file.exists(summary_csv)) {
  stop("Missing summary CSV: ", summary_csv, call. = FALSE)
}

sumry <- read.csv(summary_csv, stringsAsFactors = FALSE)

metric_dirs <- c(
  identity_mse_mean = "min",
  holdout_identity_mse_mean = "min",
  geometry_spearman_mean = "max",
  condition_knn_acc_mean = "max",
  holdout_condition_knn_acc_mean = "max",
  runtime_mean = "min"
)
metric_dirs <- metric_dirs[names(metric_dirs) %in% names(sumry)]

metric_labels <- c(
  identity_mse_mean = "Identity MSE",
  holdout_identity_mse_mean = "Holdout Identity MSE",
  geometry_spearman_mean = "Geometry Spearman",
  condition_knn_acc_mean = "Condition kNN Accuracy",
  holdout_condition_knn_acc_mean = "Holdout Condition kNN Accuracy",
  runtime_mean = "Runtime"
)

long <- sumry %>%
  mutate(setting = paste(scenario, supervision, sep = " | ")) %>%
  select(method, scenario, supervision, setting, all_of(names(metric_dirs))) %>%
  pivot_longer(cols = all_of(names(metric_dirs)), names_to = "metric", values_to = "value")

ranked <- long %>%
  group_by(setting, metric) %>%
  mutate(
    rank = if (metric_dirs[first(metric)] == "min") {
      rank(value, ties.method = "average")
    } else {
      rank(-value, ties.method = "average")
    }
  ) %>%
  ungroup()

avg_rank_by_method <- ranked %>%
  group_by(method) %>%
  summarise(
    avg_rank = mean(rank, na.rm = TRUE),
    n_settings = n_distinct(setting),
    n_metric_rows = n(),
    .groups = "drop"
  ) %>%
  arrange(avg_rank)

avg_rank_by_metric <- ranked %>%
  group_by(method, metric) %>%
  summarise(avg_rank = mean(rank, na.rm = TRUE), .groups = "drop") %>%
  mutate(metric = factor(metric, levels = names(metric_dirs), labels = unname(metric_labels[names(metric_dirs)])))

method_levels <- avg_rank_by_method$method
avg_rank_by_method$method <- factor(avg_rank_by_method$method, levels = method_levels)
avg_rank_by_metric$method <- factor(avg_rank_by_metric$method, levels = method_levels)

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

palette_vals <- setNames(grDevices::hcl.colors(length(method_levels), "Dark 3"), method_levels)

p1 <- ggplot(avg_rank_by_method, aes(x = method, y = avg_rank, fill = method)) +
  geom_col(width = 0.72, alpha = 0.95, color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", avg_rank)), vjust = -0.35, size = 3.3) +
  scale_fill_manual(values = palette_vals) +
  labs(
    title = "Average Rank Across All Scenarios and Metrics",
    subtitle = "Ranks computed within each Scenario|Supervision|Metric block; lower is better.",
    x = NULL,
    y = "Average rank"
  ) +
  theme_minimal(base_size = 13, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

p2 <- ggplot(avg_rank_by_metric, aes(x = metric, y = method, fill = avg_rank)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", avg_rank)), color = "white", size = 3.2, fontface = "bold") +
  scale_fill_gradient(low = "#0b7285", high = "#d9480f") +
  labs(
    title = "Average Rank by Metric",
    subtitle = "Lower values are better.",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

grDevices::pdf(out_pdf, width = 14, height = 8.5, onefile = TRUE)
print(p1)
print(p2)
grDevices::dev.off()

utils::write.csv(avg_rank_by_method, out_csv, row.names = FALSE)

cat("Saved:", out_pdf, "\n")
cat("Saved:", out_csv, "\n")
