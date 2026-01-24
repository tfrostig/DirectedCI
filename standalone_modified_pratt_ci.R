# =============================================================================
# Standalone Conditional Modified Pratt Confidence Interval
# =============================================================================
#
# This script provides a self-contained implementation of the Conditional
# Modified Pratt Confidence Interval for normally distributed estimators
# conditioned upon exceeding (in absolute value) a threshold.
#
# Based on: "Selection Adjusted Confidence Intervals with More Power to
# Determine the Sign" by Asaf Weinstein, William Fithian, and Yoav Benjamini
# (JASA 2012)
#
# Parameters:
#   x     - The observed value (must satisfy |x| > ct)
#   r     - The allowed expansion of acceptance region length relative to the
#           shortest CI (default: 1.2, meaning 20% expansion allowed)
#   ct    - Critical threshold for selection (default: qnorm(0.975) = 1.96)
#   alpha - Significance level (default: 0.05 for 95% CI)
#
# =============================================================================


# -----------------------------------------------------------------------------
# Helper Function: Conditional Shortest Acceptance Region
# -----------------------------------------------------------------------------
# Computes the shortest acceptance region for a given theta, conditional on
# the observation exceeding the threshold ct in absolute value.
#
# Returns a list with:
#   A - Vector of acceptance region bounds (ll, ul, lr, ur)
#   l - Length of the acceptance region

Shortest.AR <- function(theta, ct, alpha) {
  Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)

  # Compute critical theta values
  f <- function(theta) (pnorm(ct + theta) - pnorm(ct - theta)) - (1 - alpha) * Q(theta)
  theta1 <- uniroot(f, c(0, ct + qnorm(1 - alpha)))$root

  f <- function(theta) 2 * pnorm(theta - ct) - 1 - (1 - alpha) * Q(theta)
  theta2 <- uniroot(f, c(theta1, ct + qnorm(1 - alpha/2)))$root

  # Handle negative theta by symmetry
  is.neg <- 0
  if (theta < 0) is.neg <- 1
  theta <- abs(theta)

  # Compute acceptance region bounds based on theta region
  if (theta == 0) {
    ll <- -qnorm(1 - alpha/2 * Q(0))
    ul <- -ct
    lr <- ct
    ur <- qnorm(1 - alpha/2 * Q(0))
  }

  if (0 < theta && theta < theta1) {
    ll <- theta - qnorm(1 - alpha/2 * Q(theta))
    ul <- -ct
    lr <- ct
    ur <- theta + qnorm(1 - alpha/2 * Q(theta))
  }

  if (theta1 <= theta && theta < theta2) {
    ll <- -ct
    ul <- -ct
    lr <- ct
    ur <- theta + qnorm(pnorm(ct - theta) + (1 - alpha) * Q(theta))
  }

  if (theta2 <= theta) {
    ll <- -ct
    ul <- -ct
    lr <- theta - qnorm(0.5 * (1 + (1 - alpha) * Q(theta)))
    ur <- theta + qnorm(0.5 * (1 + (1 - alpha) * Q(theta)))
  }

  # Compute total length
  l <- (ul - ll) + (ur - lr)

  if (ll == -ct) {
    ll <- NA
    ul <- NA
  }

  A <- c(ll, ul, lr, ur)
  if (is.neg == 1) A <- c(-ur, -lr, -ul, -ll)

  return(list(A = A, l = l))
}


# -----------------------------------------------------------------------------
# Helper Function: theta_1 finder
# -----------------------------------------------------------------------------
# Finds the critical theta_1 value used in the Modified Pratt construction

theta_1_finder <- function(ct, alpha) {
  Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)
  f <- function(theta) (pnorm(ct + theta) - pnorm(ct - theta)) - (1 - alpha) * Q(theta)
  return(uniroot(f, c(0, ct + qnorm(1 - alpha)))$root)
}


# -----------------------------------------------------------------------------
# Main Function: Conditional Modified Pratt Confidence Interval
# -----------------------------------------------------------------------------
#
# Computes the (1-alpha)*100% Conditional Modified Pratt CI for an observed
# value x that has been selected because |x| > ct.
#
# Arguments:
#   x       - Observed value (must satisfy |x| > ct)
#   r       - Length expansion factor (default 1.2 = 20% expansion)
#   ct      - Selection threshold (default qnorm(0.975) ≈ 1.96)
#   alpha   - Significance level (default 0.05)
#   epsilon - Small perturbation for boundary cases (default 0)
#
# Returns:
#   Named vector with 'lower' and 'upper' CI bounds

MP.CI <- function(x, r = 1.2, ct = qnorm(0.975), alpha = 0.05, epsilon = 0) {

  # Validate input
  if (abs(x) <= ct) {
    warning("Observed value |x| must exceed threshold ct. Returning NA.")
    return(c(lower = NA, upper = NA))
  }

  # Selection probability function
  Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)

  # Compute critical quantities
  theta1 <- theta_1_finder(ct, alpha)

  f <- function(theta) {
    pnorm(ct + r * Shortest.AR(theta, ct, alpha)$l - theta) -
      pnorm(ct - theta) - (1 - alpha) * Q(theta)
  }
  thetatilde1 <- uniroot(f, c(0, theta1))$root

  zalphahalf <- qnorm(1 - 0.5 * alpha * Q(0))
  lzero <- Shortest.AR(0, ct, alpha)$l

  f <- function(x) {
    1 - pnorm(x) + 1 - pnorm(2 * ct + r * lzero - x) - 2 * alpha * (1 - pnorm(ct))
  }
  xtilde1 <- uniroot(f, c(ct, zalphahalf))$root
  xtilde2 <- uniroot(f, c(zalphahalf, ct + r * 2 * qnorm(1 - alpha/2)))$root

  ltilde1 <- Shortest.AR(thetatilde1, ct, alpha)$l

  # Handle negative x by symmetry
  is.neg <- 0
  if (x < 0) is.neg <- 1
  x <- abs(x)

  # Obtain lower bound of CI based on x region
  if (ct < x && x < xtilde1) {
    f <- function(theta) {
      1 - pnorm(x - theta) +
        1 - pnorm(theta - x + 2 * ct + r * Shortest.AR(theta, ct, alpha)$l) -
        alpha * Q(theta)
    }
    lower <- uniroot(f, c(-thetatilde1 - 1e-3, 0))$root
  }

  if (xtilde1 <= x && x < zalphahalf) {
    lower <- 0
  }

  if (zalphahalf <= x && x < xtilde2) {
    lower <- 0 + epsilon
  }

  if (xtilde2 <= x && x < ct + r * ltilde1) {
    f <- function(theta) {
      1 - pnorm(x - theta) +
        1 - pnorm(theta - x + 2 * ct + r * Shortest.AR(theta, ct, alpha)$l) -
        alpha * Q(theta)
    }
    lower <- uniroot(f, c(0, thetatilde1))$root
  }

  if (x >= ct + r * ltilde1) {
    f <- function(theta) {
      pnorm(x - theta) -
        pnorm(x - r * Shortest.AR(theta, ct, alpha)$l - theta) -
        (1 - alpha) * Q(theta)
    }
    m <- optimize(f, c(thetatilde1, x), maximum = TRUE)$maximum
    lower <- uniroot(f, c(thetatilde1, m))$root
  }

  # Obtain upper bound of CI
  f <- function(theta) {
    pnorm(x + r * Shortest.AR(theta, ct, alpha)$l - theta) -
      pnorm(x - theta) - (1 - alpha) * Q(theta)
  }
  m <- optimize(f, c(x, x + r * 2 * qnorm(1 - alpha/2)), maximum = TRUE)$maximum
  upper <- uniroot(f, c(0, m))$root

  CI <- c(lower = lower, upper = upper)

  # Apply symmetry for negative observations
  if (is.neg == 1) {
    CI <- c(lower = -upper, upper = -lower)
  }

  return(CI)
}


# =============================================================================
# Example Usage
# =============================================================================

if (interactive() || !exists("SOURCED_AS_LIBRARY")) {

  cat("\n")
  cat("=======================================================\n")
  cat("  Conditional Modified Pratt Confidence Interval\n")
  cat("=======================================================\n\n")

  # Default parameters
  ct <- qnorm(0.975)  # Selection threshold (~1.96)
  alpha <- 0.05       # 95% CI
  r <- 1.2            # 20% length expansion

  cat("Parameters:\n")
  cat(sprintf("  Selection threshold (ct): %.4f\n", ct))
  cat(sprintf("  Significance level (alpha): %.2f\n", alpha))
  cat(sprintf("  Length expansion factor (r): %.1f\n", r))
  cat("\n")

  # Example: Single observation
  cat("Example 1: Single observation\n")
  cat("-----------------------------\n")
  x <- 2.5
  ci <- MP.CI(x, r = r, ct = ct, alpha = alpha)
  cat(sprintf("  Observed x = %.2f\n", x))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci["lower"], ci["upper"]))
  cat(sprintf("  CI width: %.4f\n", ci["upper"] - ci["lower"]))
  cat("\n")

  # Example: Negative observation
  cat("Example 2: Negative observation (uses symmetry)\n")
  cat("------------------------------------------------\n")
  x <- -2.5
  ci <- MP.CI(x, r = r, ct = ct, alpha = alpha)
  cat(sprintf("  Observed x = %.2f\n", x))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci["lower"], ci["upper"]))
  cat("\n")

  # Example: Observation close to threshold
  cat("Example 3: Observation close to threshold\n")
  cat("------------------------------------------\n")
  x <- 2.0
  ci <- MP.CI(x, r = r, ct = ct, alpha = alpha)
  cat(sprintf("  Observed x = %.2f\n", x))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci["lower"], ci["upper"]))
  cat("\n")

  # Example: Large observation
  cat("Example 4: Large observation\n")
  cat("----------------------------\n")
  x <- 4.0
  ci <- MP.CI(x, r = r, ct = ct, alpha = alpha)
  cat(sprintf("  Observed x = %.2f\n", x))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ci["lower"], ci["upper"]))
  cat("\n")

  # Example: Multiple observations
  cat("Example 5: Multiple observations\n")
  cat("---------------------------------\n")
  observations <- c(2.1, 2.5, 3.0, 3.5, 4.0)
  results <- data.frame(
    x = observations,
    lower = NA,
    upper = NA,
    width = NA
  )

  for (i in seq_along(observations)) {
    ci <- MP.CI(observations[i], r = r, ct = ct, alpha = alpha)
    results$lower[i] <- ci["lower"]
    results$upper[i] <- ci["upper"]
    results$width[i] <- ci["upper"] - ci["lower"]
  }

  print(results, row.names = FALSE)
  cat("\n")

  # Custom function for easy batch processing
  cat("=======================================================\n")
  cat("  Batch Processing Function\n")
  cat("=======================================================\n\n")

  compute_mp_ci_batch <- function(x_values, r = 1.2, ct = qnorm(0.975), alpha = 0.05) {
    results <- lapply(x_values, function(x) {
      ci <- MP.CI(x, r = r, ct = ct, alpha = alpha)
      data.frame(
        x = x,
        lower = ci["lower"],
        upper = ci["upper"],
        width = ci["upper"] - ci["lower"],
        includes_zero = (ci["lower"] <= 0 && ci["upper"] >= 0)
      )
    })
    do.call(rbind, results)
  }

  cat("Usage: compute_mp_ci_batch(c(2.1, 2.5, 3.0))\n\n")

  # Demonstrate batch function
  batch_results <- compute_mp_ci_batch(c(2.1, 2.5, 3.0, -2.5, -3.0))
  print(batch_results, row.names = FALSE)

  cat("\n=======================================================\n")
}
