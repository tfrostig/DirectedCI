#' @title Plotting Functions for DirectedCI
#' @description Visualization utilities for acceptance regions and confidence intervals.

#' Plot Conditional Acceptance Region
#'
#' Creates a plot showing the acceptance region envelope as a function of theta.
#'
#' @param ct Numeric > 0. Critical value.
#' @param r Numeric >= 1. Inflation factor for direction-preferring AR. Default 1.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param theta_min Numeric. Minimum theta for plot. Default -3.
#' @param theta_max Numeric. Maximum theta for plot. Default 3.
#' @param n_theta Integer. Number of theta values. Default 400.
#' @param type Character. Type of AR to plot: "dp" (direction preferring),
#'   "mp" (modified pratt), or "shortest". Default "dp".
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plot_conditional_ar(ct = qnorm(0.975), r = 1.3, type = "dp")
#' }
#'
#' @export
plot_conditional_ar <- function(ct,
                                 r = 1,
                                 alpha = 0.05,
                                 theta_min = -3,
                                 theta_max = 3,
                                 n_theta = 400,
                                 type = c("dp", "mp", "shortest")) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  type <- match.arg(type)
  thetas <- seq(theta_min, theta_max, length.out = n_theta)

  df <- data.frame(
    theta = thetas,
    ar_min = NA_real_,
    ar_max = NA_real_,
    sar_min = NA_real_,
    sar_max = NA_real_
  )

  for (i in seq_along(thetas)) {
    theta_i <- thetas[i]

    # Inflated/modified acceptance region
    if (type == "dp") {
      ar <- dp_conditional_ar(theta_i, ct = ct, r = r, alpha = alpha)$ar
    } else if (type == "mp") {
      ar <- mp_conditional_ar(theta_i, ct = ct, r = r, alpha = alpha)$ar
    } else {
      ar <- shortest_conditional_ar(theta_i, ct = ct, alpha = alpha)$ar
    }

    if (!all(is.na(ar))) {
      vals <- ar[!is.na(ar)]
      df$ar_min[i] <- min(vals)
      df$ar_max[i] <- max(vals)
    }

    # Shortest acceptance region (for reference)
    sar <- shortest_conditional_ar(theta_i, ct = ct, alpha = alpha)$ar
    if (!all(is.na(sar))) {
      vals <- sar[!is.na(sar)]
      df$sar_min[i] <- min(vals)
      df$sar_max[i] <- max(vals)
    }
  }

  type_label <- switch(type,
    "dp" = "Direction Preferring",
    "mp" = "Modified Pratt",
    "shortest" = "Shortest"
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = theta)) +
    # Main AR envelope
    ggplot2::geom_line(
      ggplot2::aes(y = ar_min),
      color = "blue",
      linewidth = 1.0,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = ar_max),
      color = "blue",
      linewidth = 1.0,
      na.rm = TRUE
    )

  # Add shortest AR for comparison if not already plotting shortest
  if (type != "shortest") {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(y = sar_min),
        color = "red",
        linewidth = 0.8,
        linetype = "dashed",
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = sar_max),
        color = "red",
        linewidth = 0.8,
        linetype = "dashed",
        na.rm = TRUE
      )
  }

  p <- p +
    ggplot2::geom_hline(
      yintercept = c(-ct, ct),
      linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::labs(
      x = expression(theta),
      y = "Acceptance region",
      title = paste0(type_label, " Conditional AR (r = ", r, ")"),
      caption = if (type != "shortest") "Blue: selected AR, Red dashed: shortest AR" else NULL
    ) +
    ggplot2::theme_bw()

  p
}
