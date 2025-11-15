invert_ars <- function(y, ars) {
  return(range(ars[ars[, 2] <= y & ars[, 3] >= y, 1]))
}

ub_finder <- function(lb, theta, ct, alpha) {
  pl <- alpha - cdf_tn_outer(lb, theta, ct)
  quantile_tn_outer(1 - max(pl, 0), theta, ct)
}

lb_finder <- function(ub, theta, ct, alpha) {
  pr <- 1 - cdf_tn_outer(ub, theta, ct)
  pl <- max(alpha - pr, 0)
  quantile_tn_outer(pl, theta, ct)
}

theta_1_finder <- function(ct, alpha) {
  f <- function(theta) cdf_tn_outer(-ct, theta, ct) - alpha / 2
  return(uniroot(f, c(0, ct + qnorm(1 - alpha)))$root)
}

# for a given extension over the shortest AR length, for which theta does the AR lower bound move right of the truncation
find_tilde_1 <- function(ct, r, alpha) {
  theta1 <- theta_1_finder(ct, alpha)
  if (r == 1) {
    return(theta1)
  }
  # find tilde_1
  f <- function(theta) {
    cdf_tn_outer(ct, theta, ct) +
      1 - cdf_tn_outer(ct + r * Shortest.AR(theta, ct, alpha)$l, theta, ct) -
      alpha
  }
  theta_tilde_1 <- uniroot(f, c(0, theta1), tol = .Machine$double.eps)$root
  return(theta_tilde_1)
}

l_1_finder_left_of_crit <- function(theta, r, ct, alpha) {
  short_ar <- Shortest.AR(theta, ct, alpha)
  if (r == 1) {
    return(range(short_ar$A, na.rm = TRUE))
  }
  # function to solve
  f <- function(u) {
    cdf_tn_outer(u - 2 * ct - r * short_ar$l, theta, ct) +
      1 - cdf_tn_outer(u, theta, ct) - alpha
  }
  upper_bound <- uniroot(
    f,
    c(ct, max(short_ar$A, na.rm = TRUE)),
    tol = .Machine$double.eps
  )$root
  lower_bound <- lb_finder(upper_bound, theta, ct, alpha)
  ar <- c(lower_bound, upper_bound)
  return(ar)
}

l_1_finder_right_of_crit <- function(theta, r, ct, alpha) {
  short_ar <- Shortest.AR(theta, ct, alpha)
  if (r == 1) {
    return(range(short_ar$A, na.rm = TRUE))
  }
  # function to solve
  f <- function(l) {
    cdf_tn_outer(l, theta, ct) +
      1 - cdf_tn_outer(l + 2 * ct + r * short_ar$l, theta, ct) -
      alpha
  }
  # possible solutions
  lower_bound <- uniroot(
    f,
    c(min(short_ar$A, na.rm = TRUE), -ct)
  )$root
  upper_bound <- ub_finder(lower_bound, theta, ct, alpha)
  ar <- c(lower_bound, upper_bound)
  return(ar)
}

l_2_finder_right_of_crit <- function(theta, r, ct, alpha) {
  short_ar <- Shortest.AR(theta, ct, alpha)
  if (r == 1) {
    ar <- range(short_ar$A, na.rm = TRUE)
  }
  if (r > 1) {
    f <- function(l) {
      cdf_tn_outer(l, theta, ct) +
        1 - cdf_tn_outer(l + r * short_ar$l, theta, ct) -
        alpha
    }
    lower_uniroot_interval <- NULL
    upper_uniroot_interval <- quantile_tn_outer(1 - alpha / 2, theta, ct)
    # finding the lower bound
    lq <- quantile_tn_outer(alpha / 2, theta, ct)
    if (sign(f(ct)) != sign(f(upper_uniroot_interval))) {
      lower_uniroot_interval <- c(lower_uniroot_interval, ct)
    }
    if (sign(f(lq)) != sign(f(upper_uniroot_interval))) {
      lower_uniroot_interval <- c(lower_uniroot_interval, lq)
    }
    lower_uniroot_interval <- max(lower_uniroot_interval)
    if (is.infinite(lower_uniroot_interval)) {
      lower_bound <- optimize(function(x) abs(f(x)), c(-3, 3))$minimum
    } else {
      lower_bound <- uniroot(
        f,
        c(lower_uniroot_interval, upper_uniroot_interval),
        tol = .Machine$double.eps
      )$root
    }
    # finding the AR
    upper_bound <- ub_finder(lower_bound, theta, ct, alpha)
    ar <- c(lower_bound, upper_bound)
  }
  return(ar)
}

l_2_finder_left_of_crit <- function(theta, r, ct, alpha) {
  short_ar <- Shortest.AR(theta, ct, alpha)
  if (r == 1) {
    ar <- range(short_ar$A, na.rm = TRUE)
  }
  if (r > 1) {
    f <- function(u) {
      cdf_tn_outer(u - r * short_ar$l, theta, ct) +
        1 - cdf_tn_outer(u, theta, ct) -
        alpha
    }
    lower_uniroot_interval <- quantile_tn_outer(alpha / 2, theta, ct)
    upper_uniroot_interval <- NULL
    # finding the lower bound
    uq <- quantile_tn_outer(1 - alpha / 2, theta, ct)
    if (sign(f(-ct)) != sign(f(lower_uniroot_interval))) {
      upper_uniroot_interval <- c(upper_uniroot_interval, -ct)
    }
    if (sign(f(uq)) != sign(f(lower_uniroot_interval))) {
      upper_uniroot_interval <- c(upper_uniroot_interval, uq)
    }
    upper_uniroot_interval <- min(upper_uniroot_interval)
    if (is.infinite(upper_uniroot_interval)) {
      upper_bound <- optimize(function(x) abs(f(x)), c(-3, 3))$minimum
    } else {
      upper_bound <- uniroot(
        f,
        c(lower_uniroot_interval, upper_uniroot_interval),
        tol = .Machine$double.eps
      )$root
    }
    lower_bound <- lb_finder(upper_bound, theta, ct, alpha)
    ar <- c(lower_bound, upper_bound)
  }
  return(ar)
}
