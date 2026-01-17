#' Diagnostic Script: Find where each CI function fails
#'
#' Tests all CI functions over a grid of y values to identify failure regions.

devtools::load_all(".")

# Parameters
alpha <- 0.05
ct <- qnorm(1 - alpha / 2)  # ~1.96
r_marginal <- 1.5
r_conditional <- 1.3

# Grid of y values
y_grid <- seq(-5, 5, by = 0.1)

# Initialize results
results <- data.frame(
  y = y_grid,
  # Marginal
  marginal_shortest = NA,
  marginal_dp_pos = NA,
  marginal_dp_neg = NA,
  marginal_mp = NA,
  # Conditional (only for |y| > ct)
  conditional_shortest = NA,
  conditional_dp_pos = NA,
  conditional_dp_neg = NA,
  conditional_mp = NA
)

cat("Testing CI functions over y grid...\n")
cat(sprintf("  y range: [%.1f, %.1f], step = 0.1\n", min(y_grid), max(y_grid)))
cat(sprintf("  alpha = %.3f, ct = %.3f\n", alpha, ct))
cat(sprintf("  r_marginal = %.2f, r_conditional = %.2f\n\n", r_marginal, r_conditional))

for (i in seq_along(y_grid)) {
  y <- y_grid[i]

  # =========================================================================
  # MARGINAL CIs
  # =========================================================================

  # Shortest marginal
  results$marginal_shortest[i] <- tryCatch({
    ci <- shortest_marginal_ci(y, alpha = alpha)
    "OK"
  }, error = function(e) paste("ERROR:", e$message))

  # DP positive marginal
  results$marginal_dp_pos[i] <- tryCatch({
    ci <- dp_marginal_ci(y, r_l = r_marginal, alpha = alpha)
    "OK"
  }, error = function(e) paste("ERROR:", e$message))

  # DP negative marginal
  results$marginal_dp_neg[i] <- tryCatch({
    ci <- dn_marginal_ci(y, r_l = r_marginal, alpha = alpha)
    "OK"
  }, error = function(e) paste("ERROR:", e$message))

  # MP marginal (symmetric)
  results$marginal_mp[i] <- tryCatch({
    ci <- mp_marginal_ci(y, r = r_conditional, alpha = alpha)
    "OK"
  }, error = function(e) paste("ERROR:", e$message))

  # =========================================================================
  # CONDITIONAL CIs (only when |y| > ct)
  # =========================================================================

  if (abs(y) > ct) {
    # Shortest conditional
    results$conditional_shortest[i] <- tryCatch({
      ci <- shortest_conditional_ci(y, ct = ct, alpha = alpha)
      "OK"
    }, error = function(e) paste("ERROR:", e$message))

    # DP positive conditional
    results$conditional_dp_pos[i] <- tryCatch({
      ci <- dp_conditional_ci(y, ct = ct, r = r_conditional, alpha = alpha)
      "OK"
    }, error = function(e) paste("ERROR:", e$message))

    # DP negative conditional
    results$conditional_dp_neg[i] <- tryCatch({
      ci <- dn_conditional_ci(y, ct = ct, r = r_conditional, alpha = alpha)
      "OK"
    }, error = function(e) paste("ERROR:", e$message))

    # MP conditional (symmetric)
    results$conditional_mp[i] <- tryCatch({
      ci <- mp_conditional_ci(y, ct = ct, r = r_conditional, alpha = alpha)
      "OK"
    }, error = function(e) paste("ERROR:", e$message))
  } else {
    # Not significant - mark as N/A
    results$conditional_shortest[i] <- "N/A (|y| <= ct)"
    results$conditional_dp_pos[i] <- "N/A (|y| <= ct)"
    results$conditional_dp_neg[i] <- "N/A (|y| <= ct)"
    results$conditional_mp[i] <- "N/A (|y| <= ct)"
  }
}

# =========================================================================
# REPORT RESULTS
# =========================================================================

cat("=== MARGINAL CI FAILURES ===\n\n")

for (col in c("marginal_shortest", "marginal_dp_pos", "marginal_dp_neg", "marginal_mp")) {
  failures <- results[results[[col]] != "OK", c("y", col)]
  cat(sprintf("%s:\n", col))
  if (nrow(failures) == 0) {
    cat("  No failures\n")
  } else {
    cat(sprintf("  %d failures\n", nrow(failures)))
    cat(sprintf("  y range with failures: [%.1f, %.1f]\n",
                min(failures$y), max(failures$y)))
    # Show first few errors
    cat("  First few errors:\n")
    for (j in 1:min(3, nrow(failures))) {
      cat(sprintf("    y = %.1f: %s\n", failures$y[j], failures[[col]][j]))
    }
  }
  cat("\n")
}

cat("\n=== CONDITIONAL CI FAILURES ===\n\n")

for (col in c("conditional_shortest", "conditional_dp_pos", "conditional_dp_neg", "conditional_mp")) {
  # Exclude N/A entries
  valid <- results[!grepl("N/A", results[[col]]), ]
  failures <- valid[valid[[col]] != "OK", c("y", col)]

  cat(sprintf("%s:\n", col))
  cat(sprintf("  Tested %d y values (|y| > ct)\n", nrow(valid)))
  if (nrow(failures) == 0) {
    cat("  No failures\n")
  } else {
    cat(sprintf("  %d failures\n", nrow(failures)))
    cat(sprintf("  y range with failures: [%.1f, %.1f]\n",
                min(failures$y), max(failures$y)))
    # Show first few errors
    cat("  First few errors:\n")
    for (j in 1:min(3, nrow(failures))) {
      cat(sprintf("    y = %.1f: %s\n", failures$y[j], failures[[col]][j]))
    }
  }
  cat("\n")
}

# Save full results
write.csv(results, "simulations/ci_diagnostic_results.csv", row.names = FALSE)
cat("\nFull results saved to simulations/ci_diagnostic_results.csv\n")
