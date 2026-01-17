#' @title Shortest Conditional Confidence Intervals and Acceptance Regions
#' @description Shortest (symmetric) CI and AR for conditional inference.
#'   Conditional on |Y| > ct, where ct is the critical value.

#' Shortest Conditional Acceptance Region
#'
#' Computes the shortest acceptance region for testing a normal mean
#' conditional on significance (|Y| > ct).
#'
#' @param theta Numeric. Parameter value (on Z-scale).
#' @param ct Numeric > 0. Critical value. Default qnorm(0.975).
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return List with components:
#'   \item{ar}{Numeric vector of length 4: c(ll, ul, lr, ur).
#'     The acceptance region is [ll, ul] U [lr, ur]. NA indicates absence.}
#'   \item{length}{Numeric. Total length of the acceptance region.}
#'
#' @details
#' The acceptance region may consist of one or two disjoint intervals
#' depending on the value of theta. For theta near 0, the AR has two
#' components (one on each side of [-ct, ct]).
#'
#' @examples
#' shortest_conditional_ar(theta = 0, ct = qnorm(0.975), alpha = 0.05)
#' shortest_conditional_ar(theta = 2, ct = qnorm(0.975), alpha = 0.05)
#'
#' @export
shortest_conditional_ar <- function(theta, ct = qnorm(0.975), alpha = 0.05) {
  Q <- function(th) 1 - pnorm(ct + th) + 1 - pnorm(ct - th)

  # Compute regime boundaries
  f1 <- function(th) (pnorm(ct + th) - pnorm(ct - th)) - (1 - alpha) * Q(th)
  theta1 <- uniroot(f1, c(0, ct + qnorm(1 - alpha)))$root

  f2 <- function(th) 2 * pnorm(th - ct) - 1 - (1 - alpha) * Q(th)
  theta2 <- uniroot(f2, c(theta1, ct + qnorm(1 - alpha / 2)))$root

  # Handle negative theta via symmetry
  is_neg <- theta < 0
  theta_abs <- abs(theta)

  if (theta_abs == 0) {
    ll <- -qnorm(1 - alpha / 2 * Q(0))
    ul <- -ct
    lr <- ct
    ur <- qnorm(1 - alpha / 2 * Q(0))
  } else if (theta_abs < theta1) {
    ll <- theta_abs - qnorm(1 - alpha / 2 * Q(theta_abs))
    ul <- -ct
    lr <- ct
    ur <- theta_abs + qnorm(1 - alpha / 2 * Q(theta_abs))
  } else if (theta_abs < theta2) {
    ll <- NA_real_
    ul <- NA_real_
    lr <- ct
    ur <- theta_abs + qnorm(pnorm(ct - theta_abs) + (1 - alpha) * Q(theta_abs))
  } else {
    # theta_abs >= theta2
    ll <- NA_real_
    ul <- NA_real_
    lr <- theta_abs - qnorm(0.5 * (1 + (1 - alpha) * Q(theta_abs)))
    ur <- theta_abs + qnorm(0.5 * (1 + (1 - alpha) * Q(theta_abs)))
  }

  # Calculate length
  len <- 0
  if (!is.na(ll) && !is.na(ul)) len <- len + (ul - ll)
  if (!is.na(lr) && !is.na(ur)) len <- len + (ur - lr)

  # Apply reflection for negative theta
  if (is_neg) {
    ar <- c(-ur, -lr, -ul, -ll)
    # Handle NAs properly
    if (is.na(ll)) ar <- c(ar[1], ar[2], NA_real_, NA_real_)
  } else {
    ar <- c(ll, ul, lr, ur)
  }

  list(ar = ar, length = len)
}

# Alias for backwards compatibility
Shortest.AR <- shortest_conditional_ar

#' Shortest Conditional Confidence Interval
#'
#' Computes the shortest confidence interval for a normal mean
#' conditional on significance (|Y| > ct).
#'
#' @param y Numeric. Observed test statistic (on Z-scale). Must satisfy |y| > ct.
#' @param ct Numeric > 0. Critical value. Default qnorm(0.975).
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 2: c(lower, upper).
#'
#' @examples
#' shortest_conditional_ci(y = 2.5, ct = qnorm(0.975), alpha = 0.05)
#'
#' @export
shortest_conditional_ci <- function(y, ct = qnorm(0.975), alpha = 0.05) {
  stopifnot(length(y) == 1, length(ct) == 1, length(alpha) == 1)
  stopifnot(is.finite(y), is.finite(ct), is.finite(alpha))
  stopifnot(alpha > 0, alpha < 1)

  Q <- function(theta) 2 - pnorm(ct + theta) - pnorm(ct - theta)

  # Robust uniroot with interval expansion
  safe_uniroot <- function(fun, lower, upper, expand = 1.5, max_expand = 60, ...) {
    fl <- fun(lower)
    fu <- fun(upper)

    k <- 0
    while (is.finite(fl) && is.finite(fu) && fl * fu > 0 && k < max_expand) {
      mid <- 0.5 * (lower + upper)
      width <- (upper - lower) * expand
      lower <- mid - 0.5 * width
      upper <- mid + 0.5 * width
      fl <- fun(lower)
      fu <- fun(upper)
      k <- k + 1
    }

    if (!is.finite(fl) || !is.finite(fu) || fl * fu > 0) {
      stop("safe_uniroot: failed to bracket a root.")
    }

    uniroot(fun, c(lower, upper), ...)$root
  }

  # Compute regime boundaries
  f1 <- function(theta) {
    (pnorm(ct + theta) - pnorm(ct - theta)) - (1 - alpha) * Q(theta)
  }
  theta1 <- safe_uniroot(f1, 0, ct + qnorm(1 - alpha))

  f2 <- function(theta) {
    2 * pnorm(theta - ct) - 1 - (1 - alpha) * Q(theta)
  }
  theta2 <- safe_uniroot(f2, theta1, ct + qnorm(1 - alpha / 2))

  x1 <- theta1 + qnorm(
    pnorm(ct - theta1) +
      (1 - alpha) * (1 - pnorm(ct + theta1) + 1 - pnorm(ct - theta1))
  )
  x2 <- 2 * theta2 - ct

  fR <- function(theta) {
    denom_arg <- pnorm(ct - theta) +
      (1 - alpha) * (1 - pnorm(ct + theta) + 1 - pnorm(ct - theta))
    denom_arg <- pmin(pmax(denom_arg, .Machine$double.eps), 1 - .Machine$double.eps)
    denom <- dnorm(qnorm(denom_arg))
    if (!is.finite(denom) || denom <= .Machine$double.eps) return(Inf)
    1 + (-dnorm(ct - theta) + (1 - alpha) * (-dnorm(ct + theta) + dnorm(ct - theta))) / denom
  }
  R <- safe_uniroot(fR, theta1, theta2)

  # Handle negative y via symmetry
  is_neg <- y < 0
  y_abs <- abs(y)

  # Find lower bound
  if (ct < y_abs && y_abs < x1) {
    fl <- function(theta) 2 * (1 - pnorm(y_abs - theta)) - alpha * Q(theta)
    lower <- safe_uniroot(fl, -theta1 - 0.1, theta1)
  } else if (x1 <= y_abs && y_abs < x2) {
    fl <- function(theta) pnorm(y_abs - theta) - pnorm(ct - theta) - (1 - alpha) * Q(theta)
    lower <- safe_uniroot(fl, R, theta2)
  } else if (y_abs >= x2) {
    fl <- function(theta) 2 * pnorm(y_abs - theta) - 1 - (1 - alpha) * Q(theta)
    lower <- safe_uniroot(fl, theta2, y_abs)
  } else {
    stop("y is not in any expected region; check that |y| > ct.")
  }

  # Find upper bound
  fu <- function(theta) 2 * pnorm(theta - y_abs) - 1 - (1 - alpha) * Q(theta)
  upper <- safe_uniroot(fu, theta2, y_abs + 2 * qnorm(1 - alpha / 2))

  # Apply reflection for negative y
  if (is_neg) {
    c(-upper, -lower)
  } else {
    c(lower, upper)
  }
}

# Alias for backwards compatibility
Shortest.CI <- function(x, ct, alpha) {
  ci <- shortest_conditional_ci(x, ct, alpha)
  list(lower = ci[1], upper = ci[2])
}
