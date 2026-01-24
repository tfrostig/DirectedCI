# =============================================================================
# Standalone Non-Conditional Modified Pratt Confidence Interval
# =============================================================================
#
# This script provides self-contained implementations of Non-Conditional
# (Marginal) Modified Pratt Confidence Intervals for normally distributed
# estimators. Unlike conditional CIs, these do NOT condition on selection.
#
# Three variants are provided:
#   1. ci_sym_nc()       - Symmetric non-central CI (same expansion on both sides)
#   2. ci_ns_nc()        - Non-symmetric non-central CI (different expansion per side)
#   3. equivariant_pp_ci() - Equivariant Pratt-type CI
#
# =============================================================================


# -----------------------------------------------------------------------------
# Function 1: Symmetric Non-Central Marginal CI
# -----------------------------------------------------------------------------
#
# Computes a non-conditional CI with symmetric length expansion.
# The CI length is expanded by factor r relative to the standard CI.
#
# Arguments:
#   y       - Observed value
#   r       - Length expansion factor (r >= 1, where r=1 gives standard CI)
#   alpha   - Significance level (default 0.05 for 95% CI)
#   epsilon - Small value to make interval open at 0 when needed (default 0)
#
# Returns:
#   Vector c(lower, upper) with CI bounds

ci_sym_nc <- function(y, r, alpha = 0.05, epsilon = 0) {

  # Validate input
  if (r < 1) {
    stop("Expansion factor r must be >= 1")
  }

  # Length function for given beta
  len <- function(beta) qnorm(1 - (alpha - beta)) - qnorm(beta)

  # Find beta that achieves desired length expansion
  find_beta <- function(beta) len(beta) - r * 2 * qnorm(1 - alpha / 2)
  beta <- uniroot(find_beta, c(alpha / 2, alpha), tol = .Machine$double.eps)$root

  # Compute the two quantile values
  z_long  <- qnorm(1 - (alpha - beta))  # Longer tail
  z_short <- qnorm(1 - beta)            # Shorter tail

  # Handle negative y by symmetry
  if (y < 0) {
    return(sort(-ci_sym_nc(abs(y), r, alpha, epsilon = epsilon)))
  }

  # Compute CI based on y region
  if (y >= 0 && y < z_short) {
    ci <- c(y - z_short, y + z_short)
  }
  if (y >= z_short && y < qnorm(1 - alpha / 2)) {
    ci <- c(0, y + z_short)  # Closed at 0
  }
  if (y >= qnorm(1 - alpha / 2) && y < z_long) {
    ci <- c(0 + epsilon, y + z_short)  # Open at 0
  }
  if (y >= z_long) {
    ci <- c(y - z_long, y + z_short)
  }

  return(ci)
}


# -----------------------------------------------------------------------------
# Helper Function: Non-Symmetric Acceptance Region
# -----------------------------------------------------------------------------
# Computes the acceptance region for the non-symmetric non-central CI

ar_ns_nc <- function(theta, r_l, r_u, switch_val = 0, alpha = 0.05) {
  len <- function(beta) qnorm(1 - (alpha - beta)) - qnorm(beta)

  if (theta <= switch_val) {
    if (theta < -qnorm(1 - alpha / 2) || r_l == 1) {
      beta <- alpha / 2
    } else {
      find_beta <- function(beta) len(beta) - r_l * 2 * qnorm(1 - alpha / 2)
      beta <- uniroot(find_beta, c(0, alpha / 2), tol = .Machine$double.eps)$root
    }
  }

  if (theta > switch_val) {
    if (theta > qnorm(1 - alpha / 2) || r_u == 1) {
      beta <- alpha / 2
    } else {
      find_beta <- function(beta) len(beta) - r_u * 2 * qnorm(1 - alpha / 2)
      beta <- alpha - uniroot(find_beta, c(0, alpha / 2), tol = .Machine$double.eps)$root
    }
  }

  if (theta == switch_val && r_l == r_u) {
    beta <- alpha / 2
  }

  return(list(
    ar = c(theta + qnorm(beta), theta + qnorm(1 - (alpha - beta))),
    beta = beta
  ))
}


# -----------------------------------------------------------------------------
# Function 2: Non-Symmetric Non-Central Marginal CI
# -----------------------------------------------------------------------------
#
# Computes a non-conditional CI with different expansion factors for
# the left (negative) and right (positive) sides.
#
# Arguments:
#   y          - Observed value
#   r_l        - Left-side length expansion factor (for negative theta)
#   r_u        - Right-side length expansion factor (for positive theta, default 1)
#   switch_val - Value at which expansion switches from r_l to r_u (default 0)
#   alpha      - Significance level (default 0.05)
#   epsilon    - Small value to make interval open at 0 (default 1e-15)
#
# Returns:
#   Vector c(lower, upper) with CI bounds

ci_ns_nc <- function(y, r_l, r_u = 1, switch_val = 0, alpha = 0.05, epsilon = 1e-15) {

  ar_func <- function(theta) ar_ns_nc(theta, r_l, r_u, switch_val = 0, alpha)

  # Compute critical values
  c0 <- ar_func(qnorm(alpha / 2) + epsilon)$ar[1]
  c1 <- ar_func(switch_val - epsilon)$ar[1]
  c2 <- ar_func(switch_val + epsilon)$ar[1]
  c3 <- ar_func(switch_val - epsilon)$ar[2]
  c4 <- ar_func(switch_val + epsilon)$ar[2]

  l_switch_val_beta <- ar_func(switch_val - epsilon)$beta
  u_switch_val_beta <- ar_func(switch_val + epsilon)$beta

  # Determine CI based on y region
  if (y < c0) {
    return(c(y - qnorm(1 - alpha / 2), y - qnorm(alpha / 2)))
  }
  if (y >= c0 && y < c1) {
    return(c(y - qnorm(1 - alpha / 2), y + qnorm(1 - l_switch_val_beta)))
  }
  if (y >= c1 && y < c2) {
    return(c(y - qnorm(1 - alpha / 2), 0))  # Upper bound closed at 0
  }
  if (y >= c2 && y <= switch_val) {
    return(c(y - qnorm(1 - alpha / 2), y + qnorm(1 - u_switch_val_beta)))
  }
  if (y >= switch_val && y < c3) {
    return(c(y - qnorm(1 - (alpha - l_switch_val_beta)), y + qnorm(1 - alpha / 2)))
  }
  if (y >= c3 && y < c4) {
    return(c(0 + epsilon, y + qnorm(1 - alpha / 2)))  # Lower bound open at 0
  }
  if (y >= c4) {
    return(c(y - qnorm(1 - alpha / 2), y - qnorm(alpha / 2)))
  }
}


# -----------------------------------------------------------------------------
# Function 3: Equivariant Pratt-type CI
# -----------------------------------------------------------------------------
#
# Computes an equivariant (location-shift invariant) CI with length expansion.
# This is the simplest form - the CI is always of the form y + c(a, b).
#
# Arguments:
#   y     - Observed value
#   r     - Length expansion factor (r >= 1)
#   alpha - Significance level (default 0.05)
#
# Returns:
#   Vector c(lower, upper) with CI bounds

equivariant_pp_ci <- function(y, r, alpha = 0.05) {

  ci_shortest <- 2 * qnorm(1 - alpha / 2)

  if (r > 1) {
    find_beta <- function(beta) {
      (qnorm(1 - (alpha - beta)) - qnorm(beta)) - r * ci_shortest
    }
    beta <- uniroot(find_beta, c(alpha / 2 - 1e-9, alpha))$root
  } else if (r == 1) {
    beta <- alpha / 2
  } else {
    stop("Expansion factor r must be >= 1")
  }

  return(y + c(qnorm(beta), qnorm(1 - (alpha - beta))))
}


# -----------------------------------------------------------------------------
# Helper Function: Standard (Conventional) CI
# -----------------------------------------------------------------------------
# For comparison - the standard symmetric CI

standard_ci <- function(y, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  return(c(y - z, y + z))
}


# -----------------------------------------------------------------------------
# Helper Function: Symmetric Acceptance Region
# -----------------------------------------------------------------------------

ar_sym <- function(theta, r, alpha = 0.05) {
  len <- function(beta) qnorm(1 - (alpha - beta)) - qnorm(beta)
  find_beta <- function(beta) len(beta) - r * 2 * qnorm(1 - alpha / 2)
  beta <- uniroot(find_beta, c(alpha / 2, alpha), tol = .Machine$double.eps)$root

  if (theta == 0) {
    beta <- alpha / 2
  }

  if (theta <= 0) {
    ar <- c(theta - qnorm(1 - (alpha - beta)), theta + qnorm(1 - beta))
  }
  if (theta > 0) {
    ar <- c(theta - qnorm(1 - beta), theta + qnorm(1 - (alpha - beta)))
  }

  return(list(ar = ar, beta = beta))
}


# =============================================================================
# Batch Processing Functions
# =============================================================================

compute_sym_nc_batch <- function(y_values, r = 1.2, alpha = 0.05) {
  results <- lapply(y_values, function(y) {
    ci <- ci_sym_nc(y, r = r, alpha = alpha)
    data.frame(
      y = y,
      lower = ci[1],
      upper = ci[2],
      width = ci[2] - ci[1],
      includes_zero = (ci[1] <= 0 && ci[2] >= 0)
    )
  })
  do.call(rbind, results)
}

compute_ns_nc_batch <- function(y_values, r_l = 1.2, r_u = 1, alpha = 0.05) {
  results <- lapply(y_values, function(y) {
    ci <- ci_ns_nc(y, r_l = r_l, r_u = r_u, alpha = alpha)
    data.frame(
      y = y,
      lower = ci[1],
      upper = ci[2],
      width = ci[2] - ci[1],
      includes_zero = (ci[1] <= 0 && ci[2] >= 0)
    )
  })
  do.call(rbind, results)
}

compute_equivariant_batch <- function(y_values, r = 1.2, alpha = 0.05) {
  results <- lapply(y_values, function(y) {
    ci <- equivariant_pp_ci(y, r = r, alpha = alpha)
    data.frame(
      y = y,
      lower = ci[1],
      upper = ci[2],
      width = ci[2] - ci[1],
      includes_zero = (ci[1] <= 0 && ci[2] >= 0)
    )
  })
  do.call(rbind, results)
}


# =============================================================================
# Example Usage
# =============================================================================

if (interactive() || !exists("SOURCED_AS_LIBRARY")) {

  cat("\n")
  cat("=======================================================\n")
  cat("  Non-Conditional Modified Pratt Confidence Intervals\n")
  cat("=======================================================\n\n")

  alpha <- 0.05
  r <- 1.2  # 20% length expansion

  cat("Parameters:\n")
  cat(sprintf("  Significance level (alpha): %.2f\n", alpha))
  cat(sprintf("  Length expansion factor (r): %.1f\n", r))
  cat(sprintf("  Standard 95%% CI half-width: %.4f\n", qnorm(1 - alpha/2)))
  cat("\n")

  # -------------------------------------------------------------------------
  cat("=======================================================\n")
  cat("  1. Symmetric Non-Central CI: ci_sym_nc()\n")
  cat("=======================================================\n\n")

  cat("This CI expands symmetrically by factor r.\n\n")

  y_values <- c(-2.5, -1.5, 0, 1.5, 2.5)
  cat("Example observations:\n")
  for (y in y_values) {
    ci <- ci_sym_nc(y, r = r, alpha = alpha)
    cat(sprintf("  y = %5.2f  =>  95%% CI: [%7.4f, %7.4f]  width: %.4f\n",
                y, ci[1], ci[2], ci[2] - ci[1]))
  }
  cat("\n")

  # Compare with standard CI
  cat("Comparison with standard CI (r=1):\n")
  y <- 2.0
  ci_expanded <- ci_sym_nc(y, r = r, alpha = alpha)
  ci_standard <- standard_ci(y, alpha = alpha)
  cat(sprintf("  y = %.2f:\n", y))
  cat(sprintf("    Standard CI:  [%.4f, %.4f]  width: %.4f\n",
              ci_standard[1], ci_standard[2], ci_standard[2] - ci_standard[1]))
  cat(sprintf("    Expanded CI:  [%.4f, %.4f]  width: %.4f\n",
              ci_expanded[1], ci_expanded[2], ci_expanded[2] - ci_expanded[1]))
  cat(sprintf("    Expansion ratio: %.2f\n",
              (ci_expanded[2] - ci_expanded[1]) / (ci_standard[2] - ci_standard[1])))
  cat("\n")

  # -------------------------------------------------------------------------
  cat("=======================================================\n")
  cat("  2. Non-Symmetric Non-Central CI: ci_ns_nc()\n")
  cat("=======================================================\n\n")

  cat("This CI allows different expansion factors for left (r_l) and right (r_u).\n")
  cat("Useful when you want more power to detect effects in one direction.\n\n")

  r_l <- 1.3  # 30% expansion on left (negative) side
  r_u <- 1.0  # No expansion on right side

  cat(sprintf("  r_l (left expansion):  %.1f\n", r_l))
  cat(sprintf("  r_u (right expansion): %.1f\n\n", r_u))

  y_values <- c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5)
  cat("Example observations:\n")
  for (y in y_values) {
    ci <- ci_ns_nc(y, r_l = r_l, r_u = r_u, alpha = alpha)
    cat(sprintf("  y = %5.2f  =>  95%% CI: [%7.4f, %7.4f]  width: %.4f\n",
                y, ci[1], ci[2], ci[2] - ci[1]))
  }
  cat("\n")

  # -------------------------------------------------------------------------
  cat("=======================================================\n")
  cat("  3. Equivariant Pratt-type CI: equivariant_pp_ci()\n")
  cat("=======================================================\n\n")

  cat("This is the simplest form - CI is always y + c(a, b) where a, b are constants.\n")
  cat("The CI is location-shift invariant (equivariant).\n\n")

  y_values <- c(-2.5, -1.5, 0, 1.5, 2.5)
  cat("Example observations:\n")
  for (y in y_values) {
    ci <- suppressMessages(equivariant_pp_ci(y, r = r, alpha = alpha))
    cat(sprintf("  y = %5.2f  =>  95%% CI: [%7.4f, %7.4f]  width: %.4f\n",
                y, ci[1], ci[2], ci[2] - ci[1]))
  }
  cat("\n")

  # Note: Width is constant for equivariant CI
  cat("Note: Equivariant CI has constant width regardless of y.\n")
  cat(sprintf("      Width = %.4f (expanded by factor r = %.1f from %.4f)\n",
              ci[2] - ci[1], r, 2 * qnorm(1 - alpha/2)))
  cat("\n")

  # -------------------------------------------------------------------------
  cat("=======================================================\n")
  cat("  Comparison of All Three Methods\n")
  cat("=======================================================\n\n")

  y <- 1.8
  cat(sprintf("For y = %.2f with r = %.1f:\n\n", y, r))

  ci1 <- ci_sym_nc(y, r = r, alpha = alpha)
  ci2 <- ci_ns_nc(y, r_l = r, r_u = 1, alpha = alpha)
  ci3 <- suppressMessages(equivariant_pp_ci(y, r = r, alpha = alpha))
  ci0 <- standard_ci(y, alpha = alpha)

  cat(sprintf("  Standard CI:       [%7.4f, %7.4f]  width: %.4f\n",
              ci0[1], ci0[2], ci0[2] - ci0[1]))
  cat(sprintf("  Symmetric NC:      [%7.4f, %7.4f]  width: %.4f\n",
              ci1[1], ci1[2], ci1[2] - ci1[1]))
  cat(sprintf("  Non-Symmetric NC:  [%7.4f, %7.4f]  width: %.4f\n",
              ci2[1], ci2[2], ci2[2] - ci2[1]))
  cat(sprintf("  Equivariant:       [%7.4f, %7.4f]  width: %.4f\n",
              ci3[1], ci3[2], ci3[2] - ci3[1]))
  cat("\n")

  # -------------------------------------------------------------------------
  cat("=======================================================\n")
  cat("  Batch Processing Examples\n")
  cat("=======================================================\n\n")

  observations <- c(-2.0, -1.0, 0, 1.0, 2.0)

  cat("Symmetric NC batch (r = 1.2):\n")
  batch_sym <- compute_sym_nc_batch(observations, r = 1.2)
  print(batch_sym, row.names = FALSE)
  cat("\n")

  cat("Non-Symmetric NC batch (r_l = 1.3, r_u = 1.0):\n")
  batch_ns <- compute_ns_nc_batch(observations, r_l = 1.3, r_u = 1.0)
  print(batch_ns, row.names = FALSE)
  cat("\n")

  cat("=======================================================\n")
  cat("  Usage Summary\n")
  cat("=======================================================\n\n")

  cat("# Symmetric non-central CI:\n")
  cat("ci <- ci_sym_nc(y = 2.0, r = 1.2, alpha = 0.05)\n\n")

  cat("# Non-symmetric non-central CI:\n")
  cat("ci <- ci_ns_nc(y = 2.0, r_l = 1.3, r_u = 1.0, alpha = 0.05)\n\n")

  cat("# Equivariant Pratt CI:\n")
  cat("ci <- equivariant_pp_ci(y = 2.0, r = 1.2, alpha = 0.05)\n\n")

  cat("# Standard CI (for comparison):\n")
  cat("ci <- standard_ci(y = 2.0, alpha = 0.05)\n\n")

  cat("=======================================================\n")
}
