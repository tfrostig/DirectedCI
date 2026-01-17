#' Coverage Simulation for DirectedCI
#'
#' This script checks empirical coverage of all CI methods across a range of theta values.

# Load the package (assuming we're in the package directory)
devtools::load_all(".")

library(ggplot2)

#' Check if theta is contained in a CI
#'
#' @param theta True parameter value
#' @param ci Numeric vector of length 2: c(lower, upper)
#' @return Logical: TRUE if theta is in the CI
in_ci <- function(theta, ci) {
  lower <- ci[1]
  upper <- ci[2]

  !is.na(lower) && !is.na(upper) && theta >= lower && theta <= upper
}

#' Run coverage simulation for a single theta value
#'
#' @param theta True parameter value
#' @param K Number of simulations
#' @param alpha Significance level
#' @param r_marginal Relative length factor for marginal DP (>= 1)
#' @param r_conditional Inflation factor for conditional methods (>= 1)
#' @param ct Critical value for conditional CIs
#' @param seed Random seed for reproducibility
#' @return Data frame with coverage results
simulate_coverage_single_theta <- function(theta,
                                            K = 1000,
                                            alpha = 0.05,
                                            r_marginal = 1.5,
                                            r_conditional = 1.3,
                                            ct = qnorm(0.975),
                                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate K observations from N(theta, 1)
  y_samples <- rnorm(K, mean = theta, sd = 1)

  # Initialize coverage counters
  coverage <- list(
    # Marginal
    marginal_dp_pos = 0,
    marginal_dp_neg = 0,
    marginal_mp = 0,
    marginal_shortest = 0,
    # Conditional
    conditional_dp_pos = 0,
    conditional_dp_neg = 0,
    conditional_mp = 0,
    conditional_shortest = 0
  )

  # Count how many observations are significant (for conditional CI denominator)
  n_significant <- 0

  for (i in seq_len(K)) {
    y <- y_samples[i]

    # =========================================================================
    # MARGINAL CIs (always computed)
    # =========================================================================

    # Direction Preferring - Positive
    tryCatch({
      ci <- dp_marginal_ci(y, r_l = r_marginal, alpha = alpha)
      if (in_ci(theta, ci)) coverage$marginal_dp_pos <- coverage$marginal_dp_pos + 1
    }, error = function(e) {
      message(sprintf("Error in dp_marginal_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
    })

    # Direction Preferring - Negative
    tryCatch({
      ci <- dn_marginal_ci(y, r_l = r_marginal, alpha = alpha)
      if (in_ci(theta, ci)) coverage$marginal_dp_neg <- coverage$marginal_dp_neg + 1
    }, error = function(e) {
      message(sprintf("Error in dn_marginal_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
    })

    # Modified Pratt (symmetric)
    tryCatch({
      ci <- mp_marginal_ci(y, r = r_conditional, alpha = alpha)
      if (in_ci(theta, ci)) coverage$marginal_mp <- coverage$marginal_mp + 1
    }, error = function(e) {
      message(sprintf("Error in mp_marginal_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
    })

    # Shortest
    tryCatch({
      ci <- shortest_marginal_ci(y, alpha = alpha)
      if (in_ci(theta, ci)) coverage$marginal_shortest <- coverage$marginal_shortest + 1
    }, error = function(e) {
      message(sprintf("Error in shortest_marginal_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
    })

    # =========================================================================
    # CONDITIONAL CIs (only when |y| > ct)
    # =========================================================================

    if (abs(y) > ct) {
      n_significant <- n_significant + 1

      # Direction Preferring - Positive
      tryCatch({
        ci <- dp_conditional_ci(y, ct = ct, r = r_conditional, alpha = alpha)
        if (in_ci(theta, ci)) coverage$conditional_dp_pos <- coverage$conditional_dp_pos + 1
      }, error = function(e) {
        message(sprintf("Error in dp_conditional_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
      })

      # Direction Preferring - Negative
      tryCatch({
        ci <- dn_conditional_ci(y, ct = ct, r = r_conditional, alpha = alpha)
        if (in_ci(theta, ci)) coverage$conditional_dp_neg <- coverage$conditional_dp_neg + 1
      }, error = function(e) {
        message(sprintf("Error in dn_conditional_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
      })

      # Modified Pratt (symmetric)
      tryCatch({
        ci <- mp_conditional_ci(y, ct = ct, r = r_conditional, alpha = alpha)
        if (in_ci(theta, ci)) coverage$conditional_mp <- coverage$conditional_mp + 1
      }, error = function(e) {
        message(sprintf("Error in mp_conditional_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
      })

      # Shortest
      tryCatch({
        ci <- shortest_conditional_ci(y, ct = ct, alpha = alpha)
        if (in_ci(theta, ci)) coverage$conditional_shortest <- coverage$conditional_shortest + 1
      }, error = function(e) {
        message(sprintf("Error in shortest_conditional_ci (theta=%.2f, y=%.2f): %s", theta, y, e$message))
      })
    }
  }

  # Compute coverage proportions
  data.frame(
    theta = theta,
    K = K,
    n_significant = n_significant,
    # Marginal coverage
    marginal_dp_pos = coverage$marginal_dp_pos / K,
    marginal_dp_neg = coverage$marginal_dp_neg / K,
    marginal_mp = coverage$marginal_mp / K,
    marginal_shortest = coverage$marginal_shortest / K,
    # Conditional coverage (denominator is n_significant)
    conditional_dp_pos = if (n_significant > 0) coverage$conditional_dp_pos / n_significant else NA_real_,
    conditional_dp_neg = if (n_significant > 0) coverage$conditional_dp_neg / n_significant else NA_real_,
    conditional_mp = if (n_significant > 0) coverage$conditional_mp / n_significant else NA_real_,
    conditional_shortest = if (n_significant > 0) coverage$conditional_shortest / n_significant else NA_real_
  )
}

#' Run full coverage simulation across theta values
#'
#' @param theta_seq Sequence of theta values
#' @param K Number of simulations per theta
#' @param alpha Significance level
#' @param r_marginal Relative length for marginal DP
#' @param r_conditional Inflation factor for conditional methods
#' @param ct Critical value
#' @param seed Base seed (will be offset by theta index)
#' @return Data frame with all coverage results
run_coverage_simulation <- function(theta_seq = seq(-3, 3, by = 0.1),
                                     K = 1000,
                                     alpha = 0.05,
                                     r_marginal = 1.5,
                                     r_conditional = 1.3,
                                     ct = qnorm(0.975),
                                     seed = 12345) {

  cat("Running coverage simulation...\n")
  cat(sprintf("  theta range: [%.1f, %.1f], step = %.2f\n",
              min(theta_seq), max(theta_seq), diff(theta_seq)[1]))
  cat(sprintf("  K = %d simulations per theta\n", K))
  cat(sprintf("  alpha = %.3f, r_marginal = %.2f, r_conditional = %.2f\n",
              alpha, r_marginal, r_conditional))

  results <- vector("list", length(theta_seq))

  for (i in seq_along(theta_seq)) {
    theta <- theta_seq[i]
    if (i %% 10 == 1) {
      cat(sprintf("  Processing theta = %.2f (%d/%d)\n", theta, i, length(theta_seq)))
    }

    results[[i]] <- simulate_coverage_single_theta(
      theta = theta,
      K = K,
      alpha = alpha,
      r_marginal = r_marginal,
      r_conditional = r_conditional,
      ct = ct,
      seed = seed + i
    )
  }

  do.call(rbind, results)
}

#' Plot coverage results
#'
#' @param results Data frame from run_coverage_simulation
#' @param alpha Nominal coverage level
#' @return List of ggplot objects
plot_coverage <- function(results, alpha = 0.05) {
  nominal_coverage <- 1 - alpha

  # Reshape marginal data to long format using base R
  marginal_cols <- c("marginal_shortest", "marginal_dp_pos", "marginal_dp_neg", "marginal_mp")
  marginal_long <- data.frame(
    theta = rep(results$theta, length(marginal_cols)),
    method = rep(gsub("marginal_", "", marginal_cols), each = nrow(results)),
    coverage = unlist(results[marginal_cols])
  )
  marginal_long$method <- factor(marginal_long$method,
                                  levels = c("shortest", "dp_pos", "dp_neg", "mp"))

  # Reshape conditional data to long format
  conditional_cols <- c("conditional_shortest", "conditional_dp_pos", "conditional_dp_neg", "conditional_mp")
  conditional_long <- data.frame(
    theta = rep(results$theta, length(conditional_cols)),
    method = rep(gsub("conditional_", "", conditional_cols), each = nrow(results)),
    coverage = unlist(results[conditional_cols])
  )
  conditional_long$method <- factor(conditional_long$method,
                                     levels = c("shortest", "dp_pos", "dp_neg", "mp"))

  # Plot marginal coverage
  p_marginal <- ggplot(marginal_long, aes(x = theta, y = coverage, color = method)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = nominal_coverage, linetype = "dashed", color = "black") +
    labs(
      title = "Marginal CI Coverage",
      x = expression(theta),
      y = "Empirical Coverage",
      color = "Method"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylim(c(min(0.85, min(marginal_long$coverage, na.rm = TRUE) - 0.02), 1))

  # Plot conditional coverage
  p_conditional <- ggplot(conditional_long, aes(x = theta, y = coverage, color = method)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = nominal_coverage, linetype = "dashed", color = "black") +
    labs(
      title = "Conditional CI Coverage",
      x = expression(theta),
      y = "Empirical Coverage",
      color = "Method"
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "bottom") +
    ylim(c(min(0.85, min(conditional_long$coverage, na.rm = TRUE) - 0.02), 1))

  # Plot number of significant observations
  K <- results$K[1]
  p_n_sig <- ggplot(results, aes(x = theta, y = n_significant / K)) +
    geom_line(linewidth = 0.8, color = "darkblue") +
    geom_hline(yintercept = alpha, linetype = "dashed", color = "red") +
    labs(
      title = "Proportion of Significant Observations",
      x = expression(theta),
      y = "P(|Y| > ct)",
      caption = "Red dashed line = alpha (expected under H0)"
    ) +
    theme_bw()

  list(
    marginal = p_marginal,
    conditional = p_conditional,
    n_significant = p_n_sig
  )
}

# =============================================================================
# RUN SIMULATION
# =============================================================================

# Parameters
theta_seq <- seq(-3, 3, by = 0.5)
K <- 1000  # Number of simulations per theta
alpha <- 0.05
r_marginal <- 1.5      # Relative length for marginal DP (must be >= 1)
r_conditional <- 1.3   # Inflation factor for conditional
ct <- qnorm(1 - alpha / 2)

# Run simulation
results <- run_coverage_simulation(
  theta_seq = theta_seq,
  K = K,
  alpha = alpha,
  r_marginal = r_marginal,
  r_conditional = r_conditional,
  ct = ct,
  seed = 12345
)

# Save results
saveRDS(results, "simulations/coverage_results.rds")

# Create plots
plots <- plot_coverage(results, alpha = alpha)

# Display plots
print(plots$marginal)
print(plots$conditional)
print(plots$n_significant)

# Save plots
ggsave("simulations/coverage_marginal.png", plots$marginal, width = 10, height = 6)
ggsave("simulations/coverage_conditional.png", plots$conditional, width = 10, height = 6)
ggsave("simulations/coverage_n_significant.png", plots$n_significant, width = 10, height = 6)

# Print summary
cat("\n\n=== COVERAGE SUMMARY ===\n")
cat("\nMarginal CI Coverage (should be >= ", 1 - alpha, "):\n", sep = "")
cat(sprintf("  Shortest:    %.3f (min) - %.3f (max)\n",
            min(results$marginal_shortest), max(results$marginal_shortest)))
cat(sprintf("  DP Positive: %.3f (min) - %.3f (max)\n",
            min(results$marginal_dp_pos), max(results$marginal_dp_pos)))
cat(sprintf("  DP Negative: %.3f (min) - %.3f (max)\n",
            min(results$marginal_dp_neg), max(results$marginal_dp_neg)))
cat(sprintf("  MP:          %.3f (min) - %.3f (max)\n",
            min(results$marginal_mp), max(results$marginal_mp)))

cat("\nConditional CI Coverage (should be >= ", 1 - alpha, "):\n", sep = "")
cat(sprintf("  Shortest:    %.3f (min) - %.3f (max)\n",
            min(results$conditional_shortest, na.rm = TRUE),
            max(results$conditional_shortest, na.rm = TRUE)))
cat(sprintf("  DP Positive: %.3f (min) - %.3f (max)\n",
            min(results$conditional_dp_pos, na.rm = TRUE),
            max(results$conditional_dp_pos, na.rm = TRUE)))
cat(sprintf("  DP Negative: %.3f (min) - %.3f (max)\n",
            min(results$conditional_dp_neg, na.rm = TRUE),
            max(results$conditional_dp_neg, na.rm = TRUE)))
cat(sprintf("  MP:          %.3f (min) - %.3f (max)\n",
            min(results$conditional_mp, na.rm = TRUE),
            max(results$conditional_mp, na.rm = TRUE)))
