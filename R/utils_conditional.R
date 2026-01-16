#' @title Conditional Inference Utilities
#' @description Helper functions for conditional CI and AR computations.

#' Find theta_1: value where P(X <= -ct | theta) = alpha/2
#'
#' @param ct Numeric. Critical value.
#' @param alpha Numeric. Significance level.
#'
#' @return Numeric. The theta_1 value.
#' @keywords internal
theta_1_finder <- function(ct, alpha) {
  f <- function(theta) {
    cdf_tn_outer(-ct, theta, ct) - alpha / 2
  }
  uniroot(f, c(0, ct + qnorm(1 - alpha)))$root
}

#' Find tilde_1: value where min(AR(r)) > ct
#'
#' @param ct Numeric. Critical value.
#' @param r Numeric. Inflation factor.
#' @param alpha Numeric. Significance level.
#'
#' @return Numeric. The theta_tilde_1 value.
#' @keywords internal
find_tilde_1 <- function(ct, r, alpha) {
  theta1 <- theta_1_finder(ct, alpha)
  f <- function(theta) {
    cdf_tn_outer(ct, theta, ct) +
      1 - cdf_tn_outer(ct + r * shortest_conditional_ar(theta, ct, alpha)$length, theta, ct) -
      alpha
  }
  uniroot(f, c(0, theta1), tol = .Machine$double.eps)$root
}

#' Find upper bound given lower bound
#'
#' @param lb Numeric. Lower bound.
#' @param theta Numeric. Parameter value.
#' @param ct Numeric. Critical value.
#' @param alpha Numeric. Significance level.
#'
#' @return Numeric. Upper bound.
#' @keywords internal
ub_finder <- function(lb, theta, ct, alpha) {
  pl <- alpha - cdf_tn_outer(lb, theta, ct)
  quantile_tn_outer(1 - max(pl, 0), theta, ct)
}

#' Find lower bound given upper bound
#'
#' @param ub Numeric. Upper bound.
#' @param theta Numeric. Parameter value.
#' @param ct Numeric. Critical value.
#' @param alpha Numeric. Significance level.
#'
#' @return Numeric. Lower bound.
#' @keywords internal
lb_finder <- function(ub, theta, ct, alpha) {
  pr <- 1 - cdf_tn_outer(ub, theta, ct)
  pl <- max(alpha - pr, 0)
  quantile_tn_outer(pl, theta, ct)
}

#' Construct AR given lower and upper bounds
#'
#' @param lb Numeric. Lower bound.
#' @param ub Numeric. Upper bound.
#' @param ct Numeric. Critical value.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur).
#' @keywords internal
cond_ar_given_lb_ub <- function(lb, ub, ct) {
  if (lb >= ct) {
    return(c(NA_real_, NA_real_, lb, ub))
  } else if (ub <= -ct) {
    return(c(lb, ub, NA_real_, NA_real_))
  } else if (lb <= -ct && ub >= ct) {
    return(c(lb, -ct, ct, ub))
  }
  # Fallback
  c(NA_real_, NA_real_, NA_real_, NA_real_)
}

#' Construct AR given lower bound
#'
#' @param lb Numeric. Lower bound.
#' @param theta Numeric. Parameter value.
#' @param ct Numeric. Critical value.
#' @param alpha Numeric. Significance level.
#'
#' @return Numeric vector of length 4.
#' @keywords internal
cond_ar_given_lb <- function(lb, theta, ct, alpha) {
  ub <- ub_finder(lb, theta, ct, alpha)
  cond_ar_given_lb_ub(lb, ub, ct)
}

#' Construct AR given upper bound
#'
#' @param ub Numeric. Upper bound.
#' @param theta Numeric. Parameter value.
#' @param ct Numeric. Critical value.
#' @param alpha Numeric. Significance level.
#'
#' @return Numeric vector of length 4.
#' @keywords internal
cond_ar_given_ub <- function(ub, theta, ct, alpha) {
  lb <- lb_finder(ub, theta, ct, alpha)
  cond_ar_given_lb_ub(lb, ub, ct)
}

#' Calculate AR length
#'
#' @param ar Numeric vector of length 4.
#'
#' @return Numeric. Total length.
#' @keywords internal
calculate_ar_length <- function(ar) {
  na_ar <- is.na(ar)

  if (all(na_ar)) {
    return(NA_real_)
  } else if (any(na_ar)) {
    # Single interval
    return(max(ar, na.rm = TRUE) - min(ar, na.rm = TRUE))
  } else {
    # Two intervals: [ll, ul] and [lr, ur]
    return((ar[2] - ar[1]) + (ar[4] - ar[3]))
  }
}

#' Difference from shortest AR given lower bound
#'
#' @param target_len Numeric. Target length.
#' @param lb Numeric. Lower bound.
#' @param theta Numeric. Parameter value.
#' @param ct Numeric. Critical value.
#' @param alpha Numeric. Significance level.
#'
#' @return Numeric. Difference from target.
#' @keywords internal
diff_from_shortest_given_lb <- function(target_len, lb, theta, ct, alpha) {
  ar <- cond_ar_given_lb(lb, theta, ct, alpha)
  len_ar <- calculate_ar_length(ar)
  len_ar - target_len
}

#' Bounds for CI root finding
#'
#' @param y Numeric. Observed statistic.
#' @param ct Numeric. Critical value.
#' @param alpha Numeric. Significance level.
#' @param lower_ci_bound Logical. Finding lower (TRUE) or upper (FALSE) CI bound.
#'
#' @return Numeric vector of length 2: search interval.
#' @keywords internal
bounds_for_ci <- function(y, ct, alpha, lower_ci_bound = TRUE) {
  large_bound <- c(-1 / EPS, 1 / EPS)

  if (!lower_ci_bound) {
    f <- function(theta) quantile_tn_outer(EPS, theta, ct) - y
    theta_lower_bound <- uniroot(f, large_bound)$root

    f <- function(theta) quantile_tn_outer(alpha, theta, ct) - y
    theta_upper_bound <- uniroot(f, large_bound)$root
  } else {
    f <- function(theta) quantile_tn_outer(1 - EPS, theta, ct) - y
    theta_lower_bound <- uniroot(f, large_bound)$root

    f <- function(theta) quantile_tn_outer(1 - alpha, theta, ct) - y
    theta_upper_bound <- uniroot(f, large_bound)$root
  }

  c(theta_lower_bound, theta_upper_bound)
}

#' Invert acceptance regions to get CI
#'
#' @param y Numeric. Observed statistic.
#' @param ars Matrix. Acceptance regions.
#'
#' @return Numeric vector. CI bounds.
#' @keywords internal
invert_ars <- function(y, ars) {
  range(ars[ars[, 2] <= y & ars[, 3] >= y, 1])
}
