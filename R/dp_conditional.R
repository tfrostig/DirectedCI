#' @title Direction Preferring Conditional Confidence Intervals and Acceptance Regions
#' @description Direction-preferring CI and AR for conditional inference.
#'   Conditional on |Y| > ct, favoring a preferred direction.

#' Direction Preferring Conditional Acceptance Region (Positive Direction)
#'
#' Computes the direction-preferring acceptance region for testing a normal mean
#' conditional on significance, favoring positive values.
#'
#' @param theta Numeric. Parameter value (on Z-scale).
#' @param ct Numeric > 0. Critical value. Default qnorm(0.975).
#' @param r Numeric >= 1. Inflation factor for AR length. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return List with components:
#'   \item{ar}{Numeric vector of length 4: c(ll, ul, lr, ur).}
#'   \item{length}{Numeric. Total length of the acceptance region.}
#'
#' @examples
#' dp_conditional_ar(theta = 0, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#'
#' @export
dp_conditional_ar <- function(theta, ct = qnorm(0.975), r = 1.3, alpha = 0.05) {
  shortest <- shortest_conditional_ar(theta, ct, alpha)

  # If r = 1, return shortest AR
  if (r == 1) {
    return(shortest)
  }

  theta_1_minus_shortest <- -theta_1_finder(ct, alpha)
  theta_1_minus_r <- -find_tilde_1(ct, r, alpha)

  # Distance from inflated shortest
  f <- function(lb) {
    diff_from_shortest_given_lb(
      target_len = r * shortest$length,
      lb = lb,
      theta = theta,
      ct = ct,
      alpha = alpha
    )
  }

  # Determine regime and compute AR
  if (theta <= theta_1_minus_shortest) {
    ar <- shortest$ar
  } else if (theta > theta_1_minus_shortest && theta < theta_1_minus_r) {
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, theta, ct),
      upper = min(shortest$ar, na.rm = TRUE),
      tol = EPS
    )
    ar <- cond_ar_given_lb(lb$root, theta, ct, alpha)
  } else if (theta >= theta_1_minus_r && theta <= 0) {
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, theta, ct),
      upper = min(shortest$ar, na.rm = TRUE),
      tol = EPS
    )
    ar <- cond_ar_given_lb(lb$root, theta, ct, alpha)
  } else {
    # theta > 0
    ar <- shortest$ar
  }

  list(ar = ar, length = calculate_ar_length(ar))
}

#' Direction Preferring Conditional CI (Positive Direction)
#'
#' Computes the direction-preferring CI for a normal mean conditional on
#' significance, favoring positive values.
#'
#' @param y Numeric. Observed test statistic (on Z-scale). Must satisfy |y| > ct.
#' @param ct Numeric > 0. Critical value. Default qnorm(0.975).
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur).
#'
#' @examples
#' dp_conditional_ci(y = 2.5, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#'
#' @export
dp_conditional_ci <- function(y, ct = qnorm(0.975), r = 1.3, alpha = 0.05) {
  # Regime change points
  theta_1_minus_shortest <- -theta_1_finder(ct, alpha)
  theta_1_minus_r <- -find_tilde_1(ct, r, alpha)

  # AR wrapper that returns vector
  ar_func <- function(theta) dp_conditional_ar(theta, ct, r, alpha)$ar

  # Important thresholds
  neg_threshold_1 <- min(ar_func(theta_1_minus_r), na.rm = TRUE)
  neg_threshold_2 <- min(ar_func(0), na.rm = TRUE)
  neg_threshold_3 <- min(shortest_conditional_ar(0, ct, alpha)$ar, na.rm = TRUE)
  pos_threshold_1 <- max(ar_func(0), na.rm = TRUE)
  pos_threshold_2 <- max(shortest_conditional_ar(0, ct, alpha)$ar, na.rm = TRUE)

  # Functions for root finding
  lb_function <- function(theta) {
    ar <- ar_func(theta)
    max(ar, na.rm = TRUE) - y
  }
  ub_function <- function(theta) {
    ar <- ar_func(theta)
    min(ar, na.rm = TRUE) - y
  }

  # Determine regime and compute CI bounds
  if (y < neg_threshold_1) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  } else if (y >= neg_threshold_1 && y < neg_threshold_2) {
    lower_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root

    upper_interval <- c(theta_1_minus_r + EPS, 0 - EPS)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  } else if (y >= neg_threshold_2 && y < neg_threshold_3) {
    lower_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root
    upper_bound <- 0  # closed
  } else if (y >= neg_threshold_3 && y <= -ct) {
    lower_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  } else if (y >= ct && y < pos_threshold_1) {
    lower_interval <- c(0 - EPS, theta_1_minus_r + EPS)
    lower_bound <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  } else if (y >= pos_threshold_1 && y < pos_threshold_2) {
    lower_bound <- 0  # open

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  } else {
    # y >= pos_threshold_2
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  }

  c(NA_real_, NA_real_, lower_bound, upper_bound)
}

#' Direction Preferring Conditional CI (Negative Direction)
#'
#' Computes the direction-preferring CI for a normal mean conditional on
#' significance, favoring negative values. Uses reflection of dp_conditional_ci.
#'
#' @param y Numeric. Observed test statistic (on Z-scale). Must satisfy |y| > ct.
#' @param ct Numeric > 0. Critical value. Default qnorm(0.975).
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur).
#'
#' @export
dn_conditional_ci <- function(y, ct = qnorm(0.975), r = 1.3, alpha = 0.05) {
  ci_neg <- dp_conditional_ci(-y, ct, r, alpha)
  c(NA_real_, NA_real_, -ci_neg[4], -ci_neg[3])
}
