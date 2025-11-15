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


#' Non-Standard Non-Central MCP Acceptance Region
#'
#' Computes acceptance regions for the non-standard non-central multiple
#' comparison procedure.
#'
#' @param theta Numeric. The parameter value.
#' @param r_l Numeric. Left extension ratio relative to the shortest AR length.
#' @param r_u Numeric. Right extension ratio relative to the shortest AR length.
#' @param switch_val Numeric. Switching value for the procedure. Default is 0.
#' @param ct Numeric. Critical value. Default is \code{qnorm(0.975)}.
#' @param alpha Numeric. Significance level. Default is 0.05.
#'
#' @return Numeric vector of length 2 giving the acceptance region
#'   \code{c(lower, upper)}.
#'
#' @details
#' This function computes acceptance regions that extend the shortest acceptance
#' region by specified ratios \code{r_l} and \code{r_u} on the left and right
#' sides respectively. The procedure handles different cases based on where
#' \code{theta} falls relative to critical thresholds.
#'
#' @seealso \code{\link{shortest_ar}}, \code{\link{nsnc_mcp_ci}}
#'
#' @examples
#' nsnc_mcp_ar(theta = 0.5, r_l = 1.2, r_u = 1.2)
#' nsnc_mcp_ar(theta = 2, r_l = 1.5, r_u = 1.3, alpha = 0.1)
#'
#' @export
nsnc_mcp_ar <- function(theta, r_l, r_u,
                        switch_val = 0,
                        ct = qnorm(0.975),
                        alpha = 0.05) {
  # calculate required parameters
  theta_tilde_1_up   <- find_tilde_1(ct, r_u, alpha)
  theta_tilde_1_down <- -find_tilde_1(ct, r_l, alpha)

  if (theta <= theta_tilde_1_down) {
    ar <- l_2_finder_left_of_crit(theta, r_l, ct, alpha)
  }
  if (theta > theta_tilde_1_down & theta < switch_val) {
    ar <- l_1_finder_left_of_crit(theta, r_l, ct, alpha)
  }
  if (theta == switch_val) {
    ar <- range(Shortest.AR(theta, ct, alpha)$A, na.rm = TRUE)
  }
  if (theta < theta_tilde_1_up & theta > switch_val) {
    ar <- l_1_finder_right_of_crit(theta, r_u, ct, alpha)
  }
  if (theta >= theta_tilde_1_up) {
    ar <- l_2_finder_right_of_crit(theta, r_u, ct, alpha)
  }
  return(ar)
}


#' Non-Standard Non-Central MCP Confidence Interval
#'
#' Constructs confidence intervals using the non-standard non-central multiple
#' comparison procedure by inverting acceptance regions.
#'
#' @param y Numeric. Observed test statistic or data value.
#' @param r_l Numeric. Left extension ratio relative to the shortest AR length.
#' @param r_u Numeric. Right extension ratio relative to the shortest AR length.
#' @param switch_val Numeric. Switching value for the procedure. Default is 0.
#' @param ct Numeric. Critical value. Default is \code{qnorm(0.975)}.
#' @param alpha Numeric. Significance level. Default is 0.05.
#' @param epsilon Numeric. Small value for handling open/closed interval
#'   boundaries. Default is 0.
#'
#' @return Numeric vector of length 2 giving the confidence interval
#'   \code{c(lower_limit, upper_limit)}.
#'
#' @details
#' This function inverts the acceptance regions computed by
#' \code{\link{nsnc_mcp_ar}} to construct confidence intervals. The procedure
#' handles multiple cases based on the observed value \code{y} and uses root
#' finding to determine the confidence limits.
#'
#' When \code{r_l = r_u}, the procedure reverts to a symmetric extension of
#' the shortest acceptance region. The \code{epsilon} parameter allows for
#' fine control over whether interval boundaries are open or closed at the
#' switching value.
#'
#' @seealso \code{\link{nsnc_mcp_ar}}, \code{\link{shortest_ar}}
#'
#' @examples
#' # Construct CI for observed value y = 2
#' nsnc_mcp_ci(y = 2, r_l = 1.2, r_u = 1.2)
#'
#' # Asymmetric extension
#' nsnc_mcp_ci(y = 1.5, r_l = 1.5, r_u = 1.3)
#'
#' # Different significance level
#' nsnc_mcp_ci(y = 2.5, r_l = 1.1, r_u = 1.1, alpha = 0.1)
#'
#' @export
nsnc_mcp_ci <- function(y, r_l, r_u,
                        switch_val = 0,
                        ct = qnorm(0.975),
                        alpha = 0.05,
                        epsilon = 0) {
  eps <- 10^-11
  # calculate required parameters
  theta_tilde_1_up   <- find_tilde_1(ct, r_u, alpha)
  theta_tilde_1_down <- -find_tilde_1(ct, r_l, alpha)

  shortest_ar <- Shortest.AR(switch_val, ct, alpha)
  if (r_l == r_u) {
    null_ar <- shortest_ar$A
  } else {
    null_ar <- nsnc_mcp_ar(
      0,
      r_l = r_l, r_u = r_u,
      switch_val = 0,
      ct = ct,
      alpha = alpha
    )
  }

  # inversion functions
  find_lower_l1_left   <- function(theta, r) l_1_finder_left_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l1_left   <- function(theta, r) l_1_finder_left_of_crit(theta, r, ct, alpha)[1] - y

  find_lower_l1_right  <- function(theta, r) l_1_finder_right_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l1_right  <- function(theta, r) l_1_finder_right_of_crit(theta, r, ct, alpha)[1] - y

  find_lower_l2_left   <- function(theta, r) l_2_finder_left_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l2_left   <- function(theta, r) l_2_finder_left_of_crit(theta, r, ct, alpha)[1] - y

  find_lower_l2_right  <- function(theta, r) l_2_finder_right_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l2_right  <- function(theta, r) l_2_finder_right_of_crit(theta, r, ct, alpha)[1] - y

  # critical values
  c1 <- l_2_finder_left_of_crit(theta_tilde_1_down - eps, r_l, ct, alpha)[1]
  c2 <- l_1_finder_left_of_crit(switch_val - eps, r_l, ct, alpha)[1]
  c3 <- min(shortest_ar$A, na.rm = TRUE)
  c4 <- l_1_finder_right_of_crit(switch_val + eps, r_u, ct, alpha)[1]

  c5 <- l_1_finder_left_of_crit(switch_val - eps, r_l, ct, alpha)[2]
  c6 <- max(shortest_ar$A, na.rm = TRUE)
  c7 <- l_1_finder_right_of_crit(switch_val + eps, r_u, ct, alpha)[2]
  c8 <- l_2_finder_right_of_crit(theta_tilde_1_up + eps, r_u, ct, alpha)[2]

  # finding the CI
  if (y < c1) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, theta_tilde_1_down - eps), r = r_l)$root
    ul <- uniroot(find_upper_l2_left, c(y, theta_tilde_1_down - eps), r = r_l)$root
  }
  if (y >= c1 & y < c2) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, y), r = r_l)$root
    ul <- uniroot(find_upper_l1_left, c(theta_tilde_1_down + eps, switch_val), r = r_l)$root
  }
  if (y >= c2 & y < c3) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, -3 * y), r = r_l)$root
    ul <- switch_val - epsilon
  }
  if (y >= c3 & y < c4) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, -3 * y), r = r_l)$root
    ul <- switch_val
  }
  if (y >= c4 & y <= -ct) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, theta_tilde_1_down), r = r_l)$root
    ul <- uniroot(find_upper_l1_right,
                  c(-theta_tilde_1_up - eps, -theta_tilde_1_down),
                  r = r_u)$root
  }
  if (y >= ct & y < c5) {
    ll <- uniroot(find_lower_l1_left,
                  c(theta_tilde_1_down + eps, theta_tilde_1_up),
                  r = r_l)$root
    ul <- uniroot(find_upper_l2_right,
                  c(theta_tilde_1_up + eps, 3 * y),
                  r = r_u)$root
  }
  if (y >= c5 & y < c6) {
    ll <- switch_val
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root
  }
  if (y >= c6 & y < c7) {
    ll <- switch_val + epsilon
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root
  }
  if (y >= c7 & y < c8) {
    ll <- uniroot(find_lower_l1_right,
                  c(-theta_tilde_1_up, theta_tilde_1_up - eps),
                  r = r_u)$root
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root
  }
  if (y >= c8) {
    ll <- uniroot(find_lower_l2_right,
                  c(theta_tilde_1_up + eps, y),
                  r = r_u)$root
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root
  }

  return(c(ll, ul))
}
