
pp_conditional_ar <- function(theta, ct, r, alpha) {
  shortest_ar <- Shortest.AR(theta, ct, alpha)
  ## regimes for acceptance region
  if (r == 1) {
    return(shortest_ar$A)
  }

  theta_1_minus_shortest <- -theta_1_finder(ct, alpha)
  theta_1_minus_r        <- -find_tilde_1(ct, r, alpha)

  ## distance from inflated shortest
  f <- function(lb) {
    diff_from_shortest_given_lb(
      target_len = r * shortest_ar$l,
      lb = lb,
      theta = theta,
      ct = ct,
      alpha = alpha
    )
  }

  ## regimes for acceptance region
  if (theta <= theta_1_minus_shortest) {
    ar <- shortest_ar$A
  }
  if (theta > theta_1_minus_shortest & theta < theta_1_minus_r) {
    ## two sided acceptance region
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, theta, ct),
      upper = min(shortest_ar$A, na.rm = TRUE),
      tol   = EPS
    )
    ar <- cond_ar_given_lb(lb$root, theta, ct, alpha)

  }
  if (theta >= theta_1_minus_r & theta <= 0) {
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, theta, ct),
      upper = min(shortest_ar$A, na.rm = TRUE),
      tol   = EPS
    )
    ar <- cond_ar_given_lb(lb$root, theta, ct, alpha)

  }
  if (theta > 0) {
    ar <- shortest_ar$A
  }
  return(ar)
}
