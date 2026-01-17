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
#' Prefers positive direction by default.
#'
#' @param y Numeric. Observed test statistic (on Z-scale).
#' @param r Numeric >= 1. Inflation factor. Default 1.3.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 2: c(lower, upper).
#'
#' @examples
#' mp_marginal_ci(y = 1.5, r = 1.3, alpha = 0.05)
#' mp_marginal_ci(y = -1.5, r = 1.3, alpha = 0.05)
#'
#' @export
mp_marginal_ci <- function(y, r = 1.3, alpha = 0.05) {
  # Use symmetry: for positive y, compute for -y and reflect
  if (y > 0) {
    ci <- mp_marginal_ci(-y, r, alpha)
    # Reflect: swap and negate
    return(c(-ci[2], -ci[1]))
  }

  f <- function(beta) {
    qnorm(1 - (alpha - beta)) - qnorm(beta) - r * (2 * qnorm(1 - alpha / 2))
  }

  eps <- .Machine$double.eps
  beta_star <- uniroot(f, c(eps, alpha / 2 - eps), tol = eps)$root

  neg_threshold_1 <- -qnorm(1 - alpha / 2)
  neg_threshold_2 <- -qnorm(1 - (alpha - beta_star))

  if (y < neg_threshold_1) {
    lower <- y - qnorm(1 - (alpha - beta_star))
    upper <- y - qnorm(beta_star)
  } else if (y >= neg_threshold_1 && y < neg_threshold_2) {
    lower <- y - qnorm(1 - (alpha - beta_star))
    upper <- 0  # open interval
  } else {
    # y >= neg_threshold_2 && y <= 0
    lower <- y + qnorm(alpha - beta_star)
    upper <- y - qnorm(alpha - beta_star)
  }

  c(lower, upper)
}
