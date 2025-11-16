postivie_inflated_conditional_acceptance_region <- function(theta, r, ct, alpha) {
  # Compute shortest AR
  shortest_ar <- Shortest.AR(theta = theta,
                             ct = ct,
                             alpha = alpha)

  # Lower/Upper bound for uniroot
  minimal_lower_bound <- min(shortest_ar$A, na.rm = TRUE)
  alpha_eff <- alpha -  sqrt(.Machine$double.xmin)
  maximal_lower_bound <- quantile_tn_outer(alpha_eff, theta, ct)

  if (r == 1) {
    return(shortest_ar$A)
  }
  # Function to solve for lower bound
  f <- function(lb) {
    diff_from_shortest_given_lb(
      target_len = r * shortest_ar$l,
      lb = lb,
      theta = theta,
      ct = ct,
      alpha = alpha
    )
  }

  if (minimal_lower_bound < -ct & maximal_lower_bound > ct) {
    if (sign(f(ct)) != sign(f(maximal_lower_bound))) {
      selected_minimal_lower_bound <- ct
      selected_maximal_lower_bound <- maximal_lower_bound
    }
    if (sign(f(-ct)) != sign(f(minimal_lower_bound))) {
      selected_maximal_lower_bound <- -ct
      selected_minimal_lower_bound <- minimal_lower_bound
    }
  } else {
    selected_minimal_lower_bound <- minimal_lower_bound
    selected_maximal_lower_bound <- maximal_lower_bound
  }

  # Solve for lb
  lower_bound <- uniroot(f,
                         lower = selected_minimal_lower_bound,
                         upper = selected_maximal_lower_bound,
                         tol = .Machine$double.xmin)$root


  # Construct final acceptance region
  ar <- cond_ar_given_lb(lower_bound, theta, ct, alpha)
  return(ar)
}

cond_ar_given_lb <- function(lb, theta, ct, alpha) {
  ub <- ub_finder(lb, theta, ct, alpha)
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

negative_inflated_conditional_acceptance_region <- function(theta, r, ct, alpha, two_regions = FALSE) {
  # Flip the sign of theta
  ar_pos <- positive_inflated_conditional_acceptance_region(
    theta = -theta,
    r = r,
    ct = ct,
    alpha = alpha,
    two_regions = two_regions
  )

  # Reverse and negate the bounds
  ar_neg <- -rev(ar_pos)

  return(ar_neg)
}

## Identifying the values of theta for which there is regime change
## largest theta for which the upper bound ar is -ct
quantile_ar_theta_1_finder <- function(ct, alpha, p = alpha / 2) {
  ## find theta for which the symmetric (in terms of alpha) ar lower-bound is ct
  f <- function(theta)
    cdf_tn_outer(-ct, theta, ct) - p
  return(uniroot(f, c(0, ct + qnorm(1 - alpha)))$root)
}

generate_lower_bound_function_diff <- function(r, ct, alpha) {
  f <- function(theta) {
    ar <- postivie_inflated_conditional_acceptance_region(
      theta = theta,
      r = r,
      ct = ct,
      alpha = alpha
    )
    min(ar, na.rm = TRUE) - ct
  }
  return(f)
}

## finding minimal theta to find that the lower-bound is at -ct (shortest)
ct <- qnorm(0.975)
alpha <- 0.05

f <- generate_lower_bound_function_diff(1, ct = ct, alpha = alpha)
min_theta_for_shortest <- uniroot(f, lower = 0, upper = theta_1_finder(ct, alpha))$root

f <- generate_lower_bound_function_diff(r = 1.3, ct = ct, alpha = alpha)

min_possible_theta  <- quantile_ar_theta_1_finder(ct    = qnorm(0.975),
                                                  alpha = 0.05,
                                                  p     = 0.05)
min_theta_for_nsnsc <- uniroot(f, lower = min_possible_theta, upper = min_theta_for_shortest)

theta <- min_possible_theta
theta <- min_theta_for_shortest
ct    <- qnorm(0.975)
alpha <- 0.05
r     <- 1.3

## searching for theta for which the lower-bound is at
