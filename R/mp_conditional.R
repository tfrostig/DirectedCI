#' @title Modified Pratt Conditional Confidence Intervals and Acceptance Regions
#' @description Modified Pratt CI and AR for conditional inference.
#'   Conditional on |Y| > ct, with asymmetric error allocation.

#' Modified Pratt Conditional Acceptance Region
#'
#' Computes the Modified Pratt acceptance region for testing a normal mean
#' conditional on significance. The AR is inflated symmetrically around theta.
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
#' mp_conditional_ar(theta = 0, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#' mp_conditional_ar(theta = -1, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#'
#' @export
mp_conditional_ar <- function(theta, ct = qnorm(0.975), r = 1.3, alpha = 0.05) {
  if (r == 1) {
    return(shortest_conditional_ar(theta, ct, alpha))
  }

  # Work with absolute value and apply symmetry at end
  tmp_theta <- if (theta <= 0) theta else -theta
  shortest <- shortest_conditional_ar(tmp_theta, ct, alpha)

  theta_1_plus_r <- find_tilde_1(ct, r, alpha)
  theta_1_minus_r <- -theta_1_plus_r

  # Distance from inflated shortest
  f <- function(lb) {
    diff_from_shortest_given_lb(
      target_len = r * shortest$length,
      lb = lb,
      theta = tmp_theta,
      ct = ct,
      alpha = alpha
    )
  }

  # Determine regime and compute AR
  if (tmp_theta < theta_1_minus_r) {
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, tmp_theta, ct),
      upper = min(shortest$ar, na.rm = TRUE) + EPS,
      tol = EPS
    )
    ar <- cond_ar_given_lb(lb$root, tmp_theta, ct, alpha)
  } else if (tmp_theta >= theta_1_minus_r && tmp_theta < 0) {
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, tmp_theta, ct),
      upper = min(shortest$ar, na.rm = TRUE),
      tol = EPS
    )
    ar <- cond_ar_given_lb(lb$root, tmp_theta, ct, alpha)
  } else {
    # tmp_theta == 0
    ar <- shortest$ar
  }

  # Apply symmetry for positive theta
  if (theta > 0) {
    ar <- sort(-ar, na.last = TRUE)
  }

  list(ar = ar, length = calculate_ar_length(ar))
}

#' Modified Pratt Conditional CI
#'
#' Computes the Modified Pratt CI for a normal mean conditional on significance.
#'
#' @param y Numeric. Observed test statistic (on Z-scale). Must satisfy |y| > ct.
#' @param ct Numeric > 0. Critical value. Default qnorm(0.975).
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 2: c(lower, upper).
#'
#' @examples
#' mp_conditional_ci(y = 2.5, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#' mp_conditional_ci(y = -2.5, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#'
#' @export
mp_conditional_ci <- function(y, ct = qnorm(0.975), r = 1.3, alpha = 0.05) {
  # Use symmetry: for positive y, compute for -y and reflect
  if (sign(y) == 1) {
    ci <- mp_conditional_ci(-y, ct, r, alpha)
    return(c(-ci[2], -ci[1]))
  }

  # Regime change points
  theta_1_minus_r <- -find_tilde_1(ct, r, alpha)

  # AR wrapper
  ar_func <- function(theta) mp_conditional_ar(theta, ct, r, alpha)$ar

  # Important thresholds
  neg_threshold_1 <- min(ar_func(theta_1_minus_r), na.rm = TRUE)
  neg_threshold_2 <- min(ar_func(0), na.rm = TRUE)
  neg_threshold_3 <- min(ar_func(EPS), na.rm = TRUE)

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

    upper_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = FALSE)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  } else if (y >= neg_threshold_1 && y < neg_threshold_2) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root
    upper_bound <- 0  # open (not containing 0)
  } else if (y >= neg_threshold_2 && y < neg_threshold_3) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root
    upper_bound <- 0  # closed (containing 0)
  } else {
    # y >= neg_threshold_3 && y <= -ct
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound <- uniroot(ub_function, upper_interval)$root
  }

  c(lower_bound, upper_bound)
}
