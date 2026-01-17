#' @title User-Facing Wrappers for Conditional Confidence Intervals
#' @description User-friendly functions that handle standardization and
#'   transform results back to the original scale.
#'   All conditional CIs require |y| > ct (significance condition).

#' Direction Preferring Conditional Confidence Interval
#'
#' Computes a conditional direction-preferring confidence interval for an
#' estimator Y ~ N(mu, sigma^2), conditional on significance (|Y/sd| > ct).
#'
#' @param y Numeric. Observed estimator value.
#' @param sd Numeric > 0. Standard deviation of the estimator.
#' @param ct Numeric > 0. Critical value on Z-scale. Default qnorm(0.975).
#' @param r Numeric >= 1. Inflation factor for acceptance region. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param direction Character. Either "positive" or "negative". Default "positive".
#' @param mean Numeric. Mean of the estimator (for centering). Default 0.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur) on original scale.
#'   For single interval, two values will be NA.
#'
#' @details
#' The conditional CI is valid only when the observed statistic is significant,
#' i.e., when |y/sd| > ct. The function will stop with an error if this
#' condition is not met.
#'
#' @examples
#' # Significant positive result
#' direction_preferring_conditional_ci(y = 3.0, sd = 1.0, ct = qnorm(0.975),
#'                                      r = 1.3, direction = "positive")
#'
#' @export
direction_preferring_conditional_ci <- function(y,
                                                 sd,
                                                 ct = qnorm(0.975),
                                                 r = 1.3,
                                                 alpha = 0.05,
                                                 direction = c("positive", "negative"),
                                                 mean = 0) {
  # Validation
  validate_y(y)
  validate_sd(sd)
  validate_ct(ct)
  validate_r_conditional(r)
  validate_alpha(alpha)
  direction <- validate_direction(direction)

  # Standardize to Z-scale
  z <- (y - mean) / sd

  # Check significance condition
  if (abs(z) <= ct) {
    stop("Conditional CI requires |y/sd| > ct (significance condition not met).")
  }

  # Compute CI on Z-scale
  if (direction == "positive") {
    ci_z <- dp_conditional_ci(y = z, ct = ct, r = r, alpha = alpha)
  } else {
    ci_z <- dn_conditional_ci(y = z, ct = ct, r = r, alpha = alpha)
  }

  # Transform back to original scale
  ci_original <- ci_z
  ci_original[3:4] <- mean + sd * ci_z[3:4]

  ci_original
}

#' Modified Pratt Conditional Confidence Interval
#'
#' Computes a conditional Modified Pratt confidence interval for an
#' estimator Y ~ N(mu, sigma^2), conditional on significance.
#'
#' @param y Numeric. Observed estimator value.
#' @param sd Numeric > 0. Standard deviation of the estimator.
#' @param ct Numeric > 0. Critical value on Z-scale. Default qnorm(0.975).
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param direction Character. Either "positive" or "negative". Default "positive".
#' @param mean Numeric. Mean of the estimator (for centering). Default 0.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur) on original scale.
#'
#' @examples
#' modified_pratt_conditional_ci(y = 3.0, sd = 1.0, ct = qnorm(0.975),
#'                                r = 1.3, direction = "positive")
#'
#' @export
modified_pratt_conditional_ci <- function(y,
                                           sd,
                                           ct = qnorm(0.975),
                                           r = 1.3,
                                           alpha = 0.05,
                                           direction = c("positive", "negative"),
                                           mean = 0) {
  # Validation
  validate_y(y)
  validate_sd(sd)
  validate_ct(ct)
  validate_r_conditional(r)
  validate_alpha(alpha)
  direction <- validate_direction(direction)

  # Standardize to Z-scale
  z <- (y - mean) / sd

  # Check significance condition
  if (abs(z) <= ct) {
    stop("Conditional CI requires |y/sd| > ct (significance condition not met).")
  }

  # Compute CI on Z-scale
  if (direction == "positive") {
    ci_z <- mp_conditional_ci(y = z, ct = ct, r = r, alpha = alpha)
  } else {
    # Use reflection for negative direction
    ci_z_pos <- mp_conditional_ci(y = -z, ct = ct, r = r, alpha = alpha)
    ci_z <- c(NA_real_, NA_real_, -ci_z_pos[4], -ci_z_pos[3])
  }

  # Transform back to original scale
  ci_original <- ci_z
  ci_original[3:4] <- mean + sd * ci_z[3:4]

  ci_original
}

#' Shortest Conditional Confidence Interval
#'
#' Computes the shortest (symmetric) conditional confidence interval for an
#' estimator Y ~ N(mu, sigma^2), conditional on significance.
#'
#' @param y Numeric. Observed estimator value.
#' @param sd Numeric > 0. Standard deviation of the estimator.
#' @param ct Numeric > 0. Critical value on Z-scale. Default qnorm(0.975).
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param mean Numeric. Mean of the estimator (for centering). Default 0.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur) on original scale.
#'
#' @details
#' The shortest CI is symmetric and does not have a direction parameter.
#' It provides the minimum-length CI with the specified coverage.
#'
#' @examples
#' shortest_conditional_ci_wrapper(y = 3.0, sd = 1.0, ct = qnorm(0.975), alpha = 0.05)
#'
#' @export
shortest_conditional_ci_wrapper <- function(y,
                                             sd,
                                             ct = qnorm(0.975),
                                             alpha = 0.05,
                                             mean = 0) {
  # Validation
  validate_y(y)
  validate_sd(sd)
  validate_ct(ct)
  validate_alpha(alpha)

  # Standardize to Z-scale
  z <- (y - mean) / sd

  # Check significance condition
  if (abs(z) <= ct) {
    stop("Conditional CI requires |y/sd| > ct (significance condition not met).")
  }

  # Compute CI on Z-scale
  ci_z <- shortest_conditional_ci(y = z, ct = ct, alpha = alpha)

  # Transform back to original scale
  ci_original <- ci_z
  ci_original[3:4] <- mean + sd * ci_z[3:4]

  ci_original
}
