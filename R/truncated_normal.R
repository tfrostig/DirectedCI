#' @title Truncated Normal Distribution Functions (Outer Truncation)
#' @description Functions for the normal distribution truncated outside an interval.
#'   These are used in conditional inference where observations fall outside [-ct, ct].

#' CDF of N(mu, sigma^2) truncated outside (-b, b)
#'
#' @param y Numeric. Point at which to evaluate the CDF.
#' @param mu Numeric. Mean of the underlying normal distribution.
#' @param b Numeric. Truncation bound (symmetric: truncated outside (-b, b)).
#' @param sigma Numeric. Standard deviation (default 1).
#' @param a Numeric. Lower truncation bound (default -b for symmetric truncation).
#'
#' @return Numeric. CDF value P(Y <= y) for the truncated distribution.
#'
#' @export
cdf_tn_outer <- function(y, mu, b, sigma = 1, a = -b) {
  p_below_a <- pnorm((a - mu) / sigma)
  p_above_b <- pnorm((b - mu) / sigma, lower.tail = FALSE)
  one_minus_W <- p_below_a + p_above_b

  if (y <= a) {
    return(pnorm((y - mu) / sigma) / one_minus_W)
  }
  if (y > a & y <= b) {
    return(p_below_a / one_minus_W)
  }
  if (y > b) {
    p_above_y <- pnorm((y - mu) / sigma, lower.tail = FALSE)
    return(1 - p_above_y / one_minus_W)
  }
}

#' Quantile function for N(mu, 1) truncated outside (-b, b)
#'
#' @param p Numeric in [0, 1]. Probability.
#' @param mu Numeric. Mean of the underlying normal distribution.
#' @param b Numeric. Truncation bound.
#'
#' @return Numeric. Quantile value.
#'
#' @export
quantile_tn_outer <- function(p, mu, b) {
  W <- pnorm(b - mu) - pnorm(-b - mu)
  qc <- 1 - W
  if (p == 1) {
    return(Inf)
  }
  if (p > cdf_tn_outer(b, mu, b)) {
    return(qnorm(p * qc + 1 - qc) + mu)
  } else {
    return(qnorm(p * qc) + mu)
  }
}

#' Density of N(mu, 1) truncated outside (-b, b)
#'
#' @param x Numeric vector. Points at which to evaluate the density.
#' @param mu Numeric. Mean of the underlying normal distribution.
#' @param b Numeric. Truncation bound.
#'
#' @return Numeric vector. Density values.
#'
#' @export
density_tn_outer <- function(x, mu, b) {
  W <- 1 - (pnorm(b - mu) - pnorm(-b - mu))
  dens <- dnorm(x - mu) / W
  dens[abs(x) < b] <- 0
  return(dens)
}
