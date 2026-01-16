EPS <-  sqrt(.Machine$double.eps)

invert_ars <- function(y, ars) {
  return(range(ars[ars[, 2] <= y & ars[, 3] >= y, 1]))
}

ub_finder <- function(lb, theta, ct, alpha) {
  ## find the upper bound given the lower bound and coverage rate
  pl <- alpha - cdf_tn_outer(lb, theta, ct)
  quantile_tn_outer(1 - max(pl, 0), theta, ct)
}

lb_finder <- function(ub, theta, ct, alpha) {
  ## find the lower bound given the upper bound and coverage rate
  pr <- 1 - cdf_tn_outer(ub, theta, ct)
  pl <- max(alpha - pr, 0)
  quantile_tn_outer(pl, theta, ct)
}

theta_1_finder <- function(ct, alpha) {
  ## theta 1 is defined as the value of theta such that P(X <= -ct) = alpha / 2
  f <- function(theta)
    cdf_tn_outer(-ct, theta, ct) - alpha / 2
  return(uniroot(f, c(0, ct + qnorm(1 - alpha)))$root)
}

cond_ar_given_lb <- function(lb, theta, ct, alpha) {
  ## find the upper bound given the lower bound and return the AR
  ub <- ub_finder(lb, theta, ct, alpha)
  cond_ar_given_lb_ub(lb, ub, ct)
}

cond_ar_given_ub <- function(ub, theta, ct, alpha) {
  ## find the lower bound given the upper bound and return the AR
  lb <- lb_finder(ub, theta, ct, alpha)
  cond_ar_given_lb_ub(lb, ub, ct)
}

cond_ar_given_lb_ub <- function(lb, ub, ct) {
  ## utility function to construct AR given lb and ub
  if (lb >= ct) {
    return(c(NA, NA, lb, ub))
  } else if (ub <= -ct) {
    return(c(lb, ub, NA, NA))
  } else if (lb <= -ct && ub >= ct) {
    return(c(lb, -ct, ct, ub))
  }
}


diff_from_shortest_given_lb <- function(target_len, lb, theta, ct, alpha) {
  ar     <- cond_ar_given_lb(lb, theta, ct, alpha)
  len_ar <- calculate_ar_length(ar)
  return(len_ar - target_len)
}


calculate_ar_length <- function(ar) {
  na_ar <- is.na(ar)

  if (all(na_ar)) {
    # All values are NA
    return(NA)

  } else if (any(na_ar)) {
    # Mixed NA and non-NA values → two regions
    return(max(ar, na.rm = TRUE) - min(ar, na.rm = TRUE))

  } else if (all(!na_ar)) {
    # No NA values → single region
    return((ar[2] - ar[1]) + (ar[4] - ar[3]))
  }
}

find_tilde_1 <- function(ct, r, alpha) {
  # For a given value of r > 1, find the value of theta such that
  # the min(AR(r) > ct).
  # The function assumes that the acceptance region extends to the positive side
  theta1          <- theta_1_finder(ct, alpha)
  # find tilde_1
  f <- function(theta) {
    cdf_tn_outer(ct, theta, ct) +
      1 - cdf_tn_outer(ct + r *  Shortest.AR(theta, ct, alpha)$l, theta, ct) -
      alpha
  }
  theta_tilde_1 <- uniroot(f, c(0, theta1), tol = .Machine$double.eps)$root
  return(theta_tilde_1)
}


bounds_for_ci <- function(y, ct, alpha, lower_ci_bound = TRUE) {

  large_bound <- c(-1 / EPS, 1 / EPS)

  if (!lower_ci_bound) {
    f <- function(theta)
      quantile_tn_outer(EPS, theta, ct) - y
    theta_lower_bound <- uniroot(f, large_bound)$root
    f <- function(theta)
      quantile_tn_outer(alpha, theta, ct) - y
    theta_upper_bound <- uniroot(f, large_bound)$root
  }
  else {
    f <- function(theta)
      quantile_tn_outer(1 - EPS, theta, ct) - y
    theta_lower_bound <- uniroot(f, large_bound)$root
    f <- function(theta)
      quantile_tn_outer(1 - alpha, theta, ct) - y
    theta_upper_bound <- uniroot(f, large_bound)$root
  }

  return(c(theta_lower_bound, theta_upper_bound))
}

