#' @title User-Facing Wrappers for Marginal Confidence Intervals
#' @description User-friendly functions that handle standardization and
#'   transform results back to the original scale.

#' Direction Preferring Marginal Confidence Interval
#'
#' Computes a marginal (unconditional) direction-preferring confidence interval
#' for an estimator Y ~ N(mu, sigma^2).
#'
#' @param y Numeric. Observed estimator value.
#' @param sd Numeric > 0. Standard deviation of the estimator.
#' @param r Numeric >= 1. Inflation factor for acceptance region length.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param direction Character. Either "positive" or "negative". Default "positive".
#' @param mean Numeric. Mean of the estimator (for centering). Default 0.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur) on original scale.
#'   For single interval, two values will be NA.
#'
#' @examples
#' direction_preferring_marginal_ci(y = 3.0, sd = 1.5, r = 0.5, direction = "positive")
#' direction_preferring_marginal_ci(y = -2.0, sd = 1.0, r = 0.5, direction = "negative")
#'
#' @export
direction_preferring_marginal_ci <- function(y,
                                              sd,
                                              r,
                                              alpha = 0.05,
                                              direction = c("positive", "negative"),
                                              mean = 0) {
  # Validation
  validate_y(y)
  validate_sd(sd)
  validate_r_marginal(r)
  validate_alpha(alpha)
  direction <- validate_direction(direction)

  # Standardize to Z-scale
  z <- (y - mean) / sd

  # Compute CI on Z-scale
  if (direction == "positive") {
    ci_z <- dp_marginal_ci(y = z, r_l = r, alpha = alpha)
  } else {
    ci_z <- dn_marginal_ci(y = z, r_l = r, alpha = alpha)
  }

  # Transform back to original scale
  ci_original <- ci_z
  ci_original[3:4] <- mean + sd * ci_z[3:4]

  ci_original
}

#' Modified Pratt Marginal Confidence Interval
#'
#' Computes a marginal (unconditional) Modified Pratt confidence interval
#' for an estimator Y ~ N(mu, sigma^2).
#'
#' @param y Numeric. Observed estimator value.
#' @param sd Numeric > 0. Standard deviation of the estimator.
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param direction Character. Either "positive" or "negative". Default "positive".
#' @param mean Numeric. Mean of the estimator (for centering). Default 0.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur) on original scale.
#'
#' @examples
#' modified_pratt_marginal_ci(y = 2.5, sd = 1.0, r = 1.3, direction = "positive")
#'
#' @export
modified_pratt_marginal_ci <- function(y,
                                        sd,
                                        r = 1.3,
                                        alpha = 0.05,
                                        direction = c("positive", "negative"),
                                        mean = 0) {
  # Validation
  validate_y(y)
  validate_sd(sd)
  validate_r_conditional(r)  # MP uses inflation factor >= 1
  validate_alpha(alpha)
  direction <- validate_direction(direction)

  # Standardize to Z-scale
  z <- (y - mean) / sd

  # Compute CI on Z-scale
  if (direction == "positive") {
    ci_z <- mp_marginal_ci(y = z, r = r, alpha = alpha)
  } else {
    # Use reflection for negative direction
    ci_z_pos <- mp_marginal_ci(y = -z, r = r, alpha = alpha)
    ci_z <- c(NA_real_, NA_real_, -ci_z_pos[4], -ci_z_pos[3])
  }

  # Transform back to original scale
  ci_original <- ci_z
  ci_original[3:4] <- mean + sd * ci_z[3:4]

  ci_original
}

#' Shortest Marginal Confidence Interval
#'
#' Computes the shortest (symmetric) marginal confidence interval
#' for an estimator Y ~ N(mu, sigma^2).
#'
#' @param y Numeric. Observed estimator value.
#' @param sd Numeric > 0. Standard deviation of the estimator.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param mean Numeric. Mean of the estimator (for centering). Default 0.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur) on original scale.
#'   For the symmetric CI, returns c(NA, NA, lower, upper).
#'
#' @details
#' The shortest CI is symmetric: y +/- z_{1-alpha/2} * sd.
#' No direction parameter is needed as it is symmetric.
#'
#' @examples
#' shortest_marginal_ci_wrapper(y = 2.0, sd = 1.0, alpha = 0.05)
#'
#' @export
shortest_marginal_ci_wrapper <- function(y,
                                          sd,
                                          alpha = 0.05,
                                          mean = 0) {
  # Validation
  validate_y(y)
  validate_sd(sd)
  validate_alpha(alpha)

  # Standardize to Z-scale
  z <- (y - mean) / sd

  # Compute CI on Z-scale
  ci_z <- shortest_marginal_ci(y = z, alpha = alpha)

  # Transform back to original scale
  ci_original <- ci_z
  ci_original[3:4] <- mean + sd * ci_z[3:4]

  ci_original
}
