mp_conditional_ci <- function(y, ct, r, alpha) {
  ## assume negative value
  if (sign(y) == 1) {
    ci <- mp_conditional_ci(-y, ct, r, alpha)
    return(rev(-ci))
  }
  ## theta values at regime changes
  theta_1_minus_r        <- -find_tilde_1(ct, r, alpha)

  neg_threshold_1        <- min(mp_conditional_ar(theta_1_minus_r, ct, r, alpha), na.rm = TRUE)
  neg_threshold_2        <- min(mp_conditional_ar(0, ct, r, alpha), na.rm = TRUE)
  neg_threshold_3        <- min(mp_conditional_ar(EPS, ct, r, alpha), na.rm = TRUE)

  lb_function <- function(theta) {
    ar <- mp_conditional_ar(theta, ct, r, alpha)
    max(ar, na.rm = TRUE) - y
  }
  ub_function <- function(theta) {
    ar <- mp_conditional_ar(theta, ct, r, alpha)
    min(ar, na.rm = TRUE) - y
  }

  if (y < neg_threshold_1) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = FALSE)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }
  if (y >= neg_threshold_1 & y < neg_threshold_2) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- c(theta_1_minus_r + EPS, 0 - EPS)
    upper_bound    <- 0 ## open (not containing 0)
  }
  if (y >= neg_threshold_2 & y < neg_threshold_3) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_bound    <- 0 ## closed (containing 0)
  }
  if (y >= neg_threshold_3 & y <= -ct) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }

  return(c(lower_bound, upper_bound))

}
