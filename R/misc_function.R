
library(ggplot2)

plot_conditional_ar_outer <- function(ct,
                                      r = 1,
                                      alpha = 0.05,
                                      theta_min = -3,
                                      theta_max = 3,
                                      n_theta = 400) {
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

    ## Inflated acceptance region
    ar <- pp_conditional_ar(theta_i,
                            ct = ct,
                            r = r,
                            alpha = alpha)
    if (!all(is.na(ar))) {
      vals <- ar[!is.na(ar)]
      df$ar_min[i] <- min(vals)
      df$ar_max[i] <- max(vals)
    }

    ## Shortest acceptance region
    sar <- Shortest.AR(theta_i, ct = ct, alpha = alpha)$A
    if (!all(is.na(sar))) {
      vals <- sar[!is.na(sar)]
      df$sar_min[i] <- min(vals)
      df$sar_max[i] <- max(vals)
    }
  }



  min_val <- min(pp_conditional_ar(theta_1_minus_r - 10^-6, ct, r, alpha),
                 na.rm = TRUE)

  ggplot(df, aes(x = theta)) +
    # Inflated region envelope (blue)
    geom_line(
      aes(y = ar_min),
      color = "blue",
      linewidth = 1.0,
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = ar_max),
      color = "blue",
      linewidth = 1.0,
      na.rm = TRUE
    ) +

    # Shortest AR envelope (red)
    geom_line(
      aes(y = sar_min),
      color = "red",
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = sar_max),
      color = "red",
      linewidth = 0.8,
      na.rm = TRUE
    ) +

    # Truncation boundaries
    geom_hline(
      yintercept = c(-ct, ct),
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(yintercept = min_val,
               linetype = "dotted",
               color = "black") +
    labs(x = expression(theta),
         y = "Acceptance region",
         title = "Acceptance Region Envelope with Shortest AR (in red)") +
    theme_bw()
}
