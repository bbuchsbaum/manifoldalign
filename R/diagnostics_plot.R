#' Plot cycle-consistency residuals as a heatmap
#'
#' Aggregates triangle (i,j,k) residuals into a per-pair matrix and draws a
#' heatmap. Useful for spotting inconsistent domains or weak edges after
#' synchronization.
#'
#' @param diagnostics List returned in `align_many(...)[["diagnostics"]]`
#' @param N Optional number of domains; inferred from diagnostics when NULL
#' @param aggregate How to aggregate multiple j-paths per (i,k): "min",
#'   "mean", or "max"
#' @param na_color Fill color for missing cells (no triangles)
#' @return A ggplot object
#' @export
plot_cycle_consistency <- function(diagnostics,
                                   N = NULL,
                                   aggregate = c("min", "mean", "max"),
                                   na_color = "#f0f0f0") {
  if (is.null(diagnostics$cycle_residuals)) {
    stop("No cycle_residuals found in diagnostics; dataset graph may lack triangles.")
  }
  tri <- diagnostics$cycle_residuals
  aggregate <- match.arg(aggregate)

  if (is.null(N)) {
    if (!is.null(diagnostics$edge_residuals)) {
      N <- max(c(diagnostics$edge_residuals$i, diagnostics$edge_residuals$j))
    } else {
      N <- max(c(tri$i, tri$j, tri$k))
    }
  }

  # Build symmetric (i,k) pairs by folding triangles
  tmp <- rbind(
    data.frame(i = tri$i, k = tri$k, value = tri$residual),
    data.frame(i = tri$k, k = tri$i, value = tri$residual)
  )

  agg_fun <- switch(aggregate, min = min, mean = mean, max = max)
  agg_df <- stats::aggregate(value ~ i + k, data = tmp, FUN = agg_fun)

  grid <- expand.grid(i = seq_len(N), k = seq_len(N))
  full_df <- merge(grid, agg_df, by = c("i", "k"), all.x = TRUE)

  ggplot2::ggplot(full_df, ggplot2::aes(x = k, y = i, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = na_color) +
    ggplot2::labs(
      title = "Cycle-consistency residuals",
      x = "k (target domain)",
      y = "i (source domain)",
      fill = paste0(aggregate, " residual")
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()
}

#' Plot edge residuals as a heatmap
#'
#' Uses diagnostics$edge_residuals (i,j,residual) to build an N x N matrix of
#' symmetric residuals (min over directions when both i->j and j->i are given),
#' then draws a heatmap.
#'
#' @param diagnostics List returned by align_many(...)[["diagnostics"]]
#' @param N Optional number of domains; inferred when NULL
#' @param na_color Fill color for missing cells
#' @return ggplot object
#' @export
plot_edge_residuals_heatmap <- function(diagnostics, N = NULL, na_color = "#f0f0f0") {
  if (is.null(diagnostics$edge_residuals)) stop("diagnostics$edge_residuals missing")
  ed <- diagnostics$edge_residuals
  if (is.null(N)) N <- max(c(ed$i, ed$j))
  # symmetrize via min
  tmp <- rbind(data.frame(i = ed$i, j = ed$j, v = ed$residual), data.frame(i = ed$j, j = ed$i, v = ed$residual))
  agg <- stats::aggregate(v ~ i + j, data = tmp, FUN = min)
  grid <- expand.grid(i = seq_len(N), j = seq_len(N))
  full_df <- merge(grid, agg, by = c("i", "j"), all.x = TRUE)
  ggplot2::ggplot(full_df, ggplot2::aes(x = j, y = i, fill = v)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = na_color) +
    ggplot2::labs(title = "Edge residuals", x = "j", y = "i", fill = "residual") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()
}

#' Summarize the worst edges by residual
#'
#' Returns a data.frame of the top-K edges with largest residuals; includes
#' edge weights when available.
#'
#' @param diagnostics List returned by align_many(...)[["diagnostics"]]
#' @param top Number of worst edges to return
#' @return data.frame with columns i, j, residual, weight
#' @export
summarize_bad_edges <- function(diagnostics, top = 10) {
  if (is.null(diagnostics$edge_residuals)) stop("diagnostics$edge_residuals missing")
  ed <- diagnostics$edge_residuals
  wt <- diagnostics$edge_weights
  if (is.null(wt)) {
    weight <- rep(NA_real_, nrow(ed))
  } else {
    weight <- rep(NA_real_, nrow(ed))
    # match by (i,j)
    key_ed <- paste(ed$i, ed$j, sep = ":"); key_w <- paste(wt$i, wt$j, sep = ":")
    m <- match(key_ed, key_w)
    weight <- ifelse(is.na(m), NA_real_, wt$weight[m])
  }
  df <- data.frame(i = ed$i, j = ed$j, residual = ed$residual, weight = weight)
  df <- df[order(-df$residual), , drop = FALSE]
  head(df, top)
}
