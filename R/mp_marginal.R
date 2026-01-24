#' @title Modified Pratt Marginal Confidence Intervals and Acceptance Regions
#' @description Modified Pratt CI and AR for unconditional inference.
#'   The MP method allows asymmetric allocation of type I error between tails.

#' Modified Pratt Marginal Acceptance Region
#'
#' Computes the Modified Pratt acceptance region for testing a normal mean.
#' The AR length is inflated by factor r in the non-preferred direction.
#'
#' @param theta Numeric. Parameter value (on Z-scale).
#' @param r Numeric >= 1. Inflation factor for AR length. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return List with components:
#'   \item{ar}{Numeric vector of length 4: c(ll, ul, lr, ur).
#'     For single interval, two values will be NA.}
#'   \item{length}{Numeric. Total length of the acceptance region.}
#'
#' @examples
#' mp_marginal_ar(theta = 0, r = 1.3, alpha = 0.05)
#' mp_marginal_ar(theta = -1, r = 1.3, alpha = 0.05)
#'
#' @export
mp_marginal_ar <- function(theta, r = 1.3, alpha = 0.05) {
  f <- function(beta) {
    qnorm(1 - (alpha - beta)) - qnorm(beta) - r * (2 * qnorm(1 - alpha / 2))
  }

  eps <- .Machine$double.eps

  if (theta < 0) {
    beta_star <- uniroot(
      f,
      c(eps, alpha / 2 - eps),
      tol = eps
    )$root
    lower <- theta + qnorm(beta_star)
    upper <- theta + qnorm(1 - (alpha - beta_star))
  } else if (theta == 0) {
    lower <- theta + qnorm(alpha / 2)
    upper <- theta - qnorm(alpha / 2)
  } else {
    beta_left <- uniroot(
      f,
      c(eps, alpha / 2 - eps),
      tol = eps
    )$root
    beta_star <- alpha - beta_left
    lower <- theta + qnorm(beta_star)
    upper <- theta + qnorm(1 - (alpha - beta_star))
  }

  list(
    ar = c(NA_real_, NA_real_, lower, upper),
    length = upper - lower
  )
}

#' Modified Pratt Marginal Confidence Interval
#'
#' Computes the Modified Pratt confidence interval for a normal mean.
#' Prefers positive direction by default. The CI truncates at 0 for
#' certain y values to give more power to determine the sign.
#'
#' @param y Numeric. Observed test statistic (on Z-scale).
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#' @param epsilon Numeric. Small value for open intervals at 0. Default 0.
#'
#' @return Numeric vector of length 2: c(lower, upper).
#'
#' @examples
#' mp_marginal_ci(y = 1.5, r = 1.3, alpha = 0.05)
#' mp_marginal_ci(y = -1.5, r = 1.3, alpha = 0.05)
#'
#' @export
mp_marginal_ci <- function(y, r = 1.3, alpha = 0.05, epsilon = 0) {
  # Handle negative y by symmetry: compute for |y| and reflect
  if (y < 0) {
    ci <- mp_marginal_ci(-y, r, alpha, epsilon)
    return(sort(-ci))
  }

  # Find beta that achieves desired length expansion
  # Length function: len(beta) = qnorm(1 - (alpha - beta)) - qnorm(beta)
  # We want: len(beta) = r * 2 * qnorm(1 - alpha/2)
  f <- function(beta) {
    qnorm(1 - (alpha - beta)) - qnorm(beta) - r * 2 * qnorm(1 - alpha / 2)
  }

  eps <- .Machine$double.eps
  beta <- uniroot(f, c(alpha / 2, alpha - eps), tol = eps)$root

  # Compute the two quantile values
  z_long <- qnorm(1 - (alpha - beta))   # Longer tail
  z_short <- qnorm(1 - beta)            # Shorter tail
  z_alpha <- qnorm(1 - alpha / 2)       # Standard threshold

  # Compute CI based on y region (for positive y)
  # Four regimes with truncation at 0 for middle regimes
  if (y >= 0 && y < z_short) {
    # Small y: full symmetric interval
    ci <- c(y - z_short, y + z_short)
  } else if (y >= z_short && y < z_alpha) {
    # Medium-small y: truncate lower bound to 0 (closed)
    ci <- c(0, y + z_short)
  } else if (y >= z_alpha && y < z_long) {
    # Medium-large y: truncate lower bound to 0 (open)
    ci <- c(0 + epsilon, y + z_short)
  } else {
    # Large y: full asymmetric interval
    ci <- c(y - z_long, y + z_short)
  }

  ci
}
