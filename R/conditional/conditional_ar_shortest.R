#' Shortest Acceptance Region for Normal Mean Testing
#'
#' Computes the shortest acceptance region for testing a normal mean with known
#' variance using a two-sided test.
#'
#' @param theta Numeric. The true parameter value (mean) under consideration.
#' @param ct Numeric. Critical value for the test. Default is \code{qnorm(0.975)}
#'   for a two-sided test at the 0.05 level.
#' @param alpha Numeric. Significance level. Default is 0.05.
#'
#' @return A list with two components:
#'   \item{A}{Numeric vector of length 4 giving the acceptance region bounds
#'     \code{c(ll, ul, lr, ur)} where the acceptance region is
#'     \code{[ll, ul] U [lr, ur]}. If there is no left interval, \code{ll} and
#'     \code{ul} are \code{NA}.}
#'   \item{l}{Numeric. The total length of the acceptance region.}
#'
#' @details
#' This function implements the shortest acceptance region approach for testing
#' \eqn{H_0: \theta = \theta_0} against \eqn{H_1: \theta \ne \theta_0} when
#' observations follow a normal distribution with unit variance. The acceptance
#' region may consist of one or two disjoint intervals depending on the value
#' of \code{theta}.
#'
#' @examples
#' # Compute shortest AR for theta = 0
#' shortest_ar(theta = 0, ct = qnorm(0.975), alpha = 0.05)
#'
#' # Compute for theta = 1
#' shortest_ar(theta = 1, ct = qnorm(0.975), alpha = 0.05)
#'
#' @export
shortest_ar <- function(theta, ct = qnorm(0.975), alpha = 0.05) {
  Q <- function(theta) 1 - pnorm(ct + theta) + 1 - pnorm(ct - theta)
  # compute useful quantities
  f <- function(theta) (pnorm(ct + theta) - pnorm(ct - theta)) - (1 - alpha) * Q(theta)
  theta1 <- uniroot(f, c(0, ct + qnorm(1 - alpha)))$root

  f <- function(theta) 2 * pnorm(theta - ct) - 1 - (1 - alpha) * Q(theta)
  theta2 <- uniroot(f, c(theta1, ct + qnorm(1 - alpha / 2)))$root

  # compute ends of the AR
  is.neg <- 0
  if (theta < 0) is.neg <- 1
  theta <- abs(theta)

  if (theta == 0) {
    ll <- -qnorm(1 - alpha / 2 * Q(0))
    ul <- -ct
    lr <- ct
    ur <- qnorm(1 - alpha / 2 * Q(0))
  }

  if (0 < theta && theta < theta1) {
    ll <- theta - qnorm(1 - alpha / 2 * Q(theta))
    ul <- -ct
    lr <- ct
    ur <- theta + qnorm(1 - alpha / 2 * Q(theta))
  }

  if (theta1 <= theta && theta < theta2) {
    ll <- -ct
    ul <- -ct
    lr <- ct
    ur <- theta + qnorm(pnorm(ct - theta) + (1 - alpha) * Q(theta))
  }

  if (theta2 <= theta) {
    ll <- -ct
    ul <- -ct
    lr <- theta - qnorm(.5 * (1 + (1 - alpha) * Q(theta)))
    ur <- theta + qnorm(.5 * (1 + (1 - alpha) * Q(theta)))
  }

  l <- (ul - ll) + (ur - lr)
  if (ll == -ct) {
    ll <- NA
    ul <- NA
  }
  A <- c(ll, ul, lr, ur)
  if (is.neg == 1) A <- c(-ur, -lr, -ul, -ll)
  v <- list(A = A, l = l)
  return(v)
}

# Create alias for backwards compatibility
Shortest.AR <- shortest_ar

