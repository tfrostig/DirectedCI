#' @title Shortest Marginal Confidence Intervals and Acceptance Regions
#' @description Shortest (symmetric) CI and AR for unconditional inference.

#' Shortest Marginal Acceptance Region
#'
#' Computes the shortest (symmetric) acceptance region for testing a normal mean.
#'
#' @param theta Numeric. Parameter value (on Z-scale).
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return List with components:
#'   \item{ar}{Numeric vector of length 4: c(ll, ul, lr, ur). For symmetric AR,
#'     this is c(NA, NA, lower, upper).}
#'   \item{length}{Numeric. Total length of the acceptance region.}
#'
#' @examples
#' shortest_marginal_ar(theta = 0, alpha = 0.05)
#' shortest_marginal_ar(theta = 1.5, alpha = 0.05)
#'
#' @export
shortest_marginal_ar <- function(theta, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  lower <- theta - z
  upper <- theta + z

  list(
    ar = c(NA_real_, NA_real_, lower, upper),
    length = upper - lower
  )
}

#' Shortest Marginal Confidence Interval
#'
#' Computes the shortest (symmetric) confidence interval for a normal mean.
#' This is the standard symmetric CI: y +/- z_{1-alpha/2}.
#'
#' @param y Numeric. Observed test statistic (on Z-scale).
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 2: c(lower, upper).
#'
#' @examples
#' shortest_marginal_ci(y = 1.5, alpha = 0.05)
#'
#' @export
shortest_marginal_ci <- function(y, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  c(y - z, y + z)
}
