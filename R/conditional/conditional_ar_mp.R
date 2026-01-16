
mp_conditional_ar <- function(theta, ct = 1.96, r = 1.3, alpha = 0.05) {
  if (r == 1) {
    stop("Modified Pratt AR is identical to shortest AR when r = 1")
  }

  tmp_theta <- if (theta <= 0) theta else -theta
  shortest_ar <- Shortest.AR(tmp_theta, ct, alpha)

  theta_1_plus_r   <- find_tilde_1(ct, r, alpha)
  theta_1_minus_r  <- -theta_1_plus_r

  ## distance from inflated shortest
  f <- function(lb) {
    diff_from_shortest_given_lb(
      target_len = r * shortest_ar$l,
      lb = lb,
      theta = tmp_theta,
      ct = ct,
      alpha = alpha
    )
  }

  if (tmp_theta < theta_1_minus_r) {
    ## two sided acceptance region
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, tmp_theta, ct),
      upper = min(shortest_ar$A, na.rm = TRUE) + EPS,
      tol   = EPS
    )
    ar <- cond_ar_given_lb(lb$root, tmp_theta, ct, alpha)

  }
  if (tmp_theta >= theta_1_minus_r & tmp_theta < 0) {
    lb <- uniroot(
      f,
      lower = quantile_tn_outer(.Machine$double.eps, tmp_theta, ct),
      upper = min(shortest_ar$A, na.rm = TRUE),
      tol   = EPS
    )
    ar <- cond_ar_given_lb(lb$root, tmp_theta, ct, alpha)

  }
  if (tmp_theta == 0) {
    ar <- shortest_ar$A
  }

  # leveraging symmetry
  if (theta > 0) {
    ar <- sort(-ar, na.last = TRUE)
  }

  return(ar)
}



theta_seq <- seq(0, 3, 0.01)

# Store results
ar_lower <- numeric(length(theta_seq))
ar_upper <- numeric(length(theta_seq))

for (i in seq_along(theta_seq)) {
  ar <- mp_conditional_ar(theta_seq[i], ct = 1, r = 1.4, alpha = 0.05)
  ar_lower[i] <- min(ar, na.rm=TRUE)
  ar_upper[i] <- max(ar, na.rm=TRUE)
}

# Plot
plot(theta_seq, ar_upper, type = "l", col = "blue", lwd = 2,
     xlab = expression(theta), ylab = "Acceptance Region",
     main = "Modified Pratt Conditional AR (ct=1, r=1.4)",
     ylim = range(c(ar_lower, ar_upper), na.rm = TRUE))
lines(theta_seq, ar_lower, col = "red", lwd = 2)
legend("topright", legend = c("Upper", "Lower"),
       col = c("blue", "red"), lwd = 2)
