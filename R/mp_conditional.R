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
#' Based on "Selection Adjusted Confidence Intervals with More Power to
#' Determine the Sign" by Weinstein, Fithian, and Benjamini (JASA 2012).
#'
#' @param y Numeric. Observed test statistic (on Z-scale). Must satisfy |y| > ct.
#' @param ct Numeric > 0. Critical value. Default qnorm(0.975).
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param epsilon Numeric. Small value for open intervals at 0. Default 0.
#'
#' @return Numeric vector of length 2: c(lower, upper).
#'
#' @examples
#' mp_conditional_ci(y = 2.5, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#' mp_conditional_ci(y = -2.5, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#'
#' @export
mp_conditional_ci <- function(y, ct = qnorm(0.975), r = 1.3, alpha = 0.05, epsilon = 0) {
  # Validate input
  if (abs(y) <= ct) {
    warning("Observed value |y| must exceed threshold ct. Returning NA.")
    return(c(NA_real_, NA_real_))
  }

  # Selection probability function
  Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)

  # Compute critical quantities using the paper's formulas
  # theta1: regime boundary for shortest AR
  f1 <- function(theta) (pnorm(ct + theta) - pnorm(ct - theta)) - (1 - alpha) * Q(theta)
  theta1 <- uniroot(f1, c(0, ct + qnorm(1 - alpha)))$root

  # thetatilde1: where inflated AR first crosses ct
  f2 <- function(theta) {
    pnorm(ct + r * shortest_conditional_ar(theta, ct, alpha)$length - theta) -
      pnorm(ct - theta) - (1 - alpha) * Q(theta)
  }
  thetatilde1 <- uniroot(f2, c(0, theta1))$root

  # zalphahalf: quantile for theta=0
  zalphahalf <- qnorm(1 - 0.5 * alpha * Q(0))

  # lzero: shortest AR length at theta=0
  lzero <- shortest_conditional_ar(0, ct, alpha)$length

  # xtilde1, xtilde2: transition points for x
  f3 <- function(x) {
    1 - pnorm(x) + 1 - pnorm(2 * ct + r * lzero - x) - 2 * alpha * (1 - pnorm(ct))
  }
  xtilde1 <- uniroot(f3, c(ct, zalphahalf))$root
  xtilde2 <- uniroot(f3, c(zalphahalf, ct + r * 2 * qnorm(1 - alpha / 2)))$root

  # ltilde1: shortest AR length at thetatilde1
  ltilde1 <- shortest_conditional_ar(thetatilde1, ct, alpha)$length

  # Handle negative y by symmetry
  is_neg <- y < 0
  x <- abs(y)

  # Obtain lower bound of CI based on x region (5 regimes)
  if (ct < x && x < xtilde1) {
    # Regime 1: just above threshold
    f_lower <- function(theta) {
      1 - pnorm(x - theta) +
        1 - pnorm(theta - x + 2 * ct + r * shortest_conditional_ar(theta, ct, alpha)$length) -
        alpha * Q(theta)
    }
    lower <- uniroot(f_lower, c(-thetatilde1 - 1e-3, 0))$root
  } else if (xtilde1 <= x && x < zalphahalf) {
    # Regime 2: lower bound is 0 (closed)
    lower <- 0
  } else if (zalphahalf <= x && x < xtilde2) {
    # Regime 3: lower bound is 0 (open)
    lower <- 0 + epsilon
  } else if (xtilde2 <= x && x < ct + r * ltilde1) {
    # Regime 4: transition region
    f_lower <- function(theta) {
      1 - pnorm(x - theta) +
        1 - pnorm(theta - x + 2 * ct + r * shortest_conditional_ar(theta, ct, alpha)$length) -
        alpha * Q(theta)
    }
    lower <- uniroot(f_lower, c(0, thetatilde1))$root
  } else {
    # Regime 5: large x (x >= ct + r * ltilde1)
    f_lower <- function(theta) {
      pnorm(x - theta) -
        pnorm(x - r * shortest_conditional_ar(theta, ct, alpha)$length - theta) -
        (1 - alpha) * Q(theta)
    }
    m <- optimize(f_lower, c(thetatilde1, x), maximum = TRUE)$maximum
    lower <- uniroot(f_lower, c(thetatilde1, m))$root
  }

  # Obtain upper bound of CI (always via root finding)
  f_upper <- function(theta) {
    pnorm(x + r * shortest_conditional_ar(theta, ct, alpha)$length - theta) -
      pnorm(x - theta) - (1 - alpha) * Q(theta)
  }
  m <- optimize(f_upper, c(x, x + r * 2 * qnorm(1 - alpha / 2)), maximum = TRUE)$maximum
  upper <- uniroot(f_upper, c(0, m))$root

  # Apply symmetry for negative observations
  if (is_neg) {
    c(-upper, -lower)
  } else {
    c(lower, upper)
  }
}
