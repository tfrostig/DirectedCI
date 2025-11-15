#' Non-Standard Non-Central MCP Confidence Interval
#'
#' Constructs confidence intervals using the non-standard non-central multiple
#' comparison procedure by inverting acceptance regions.
#'
#' @param y Numeric scalar. Observed test statistic or data value.
#' @param r_l Numeric scalar. Left extension ratio relative to the shortest AR length.
#' @param r_u Numeric scalar. Right extension ratio relative to the shortest AR length.
#' @param switch_val Numeric. Switching value for the procedure. Default is 0.
#' @param ct Numeric. Critical value (assumed > 0). Default is \code{qnorm(0.975)}.
#' @param alpha Numeric in (0, 1). Significance level. Default is 0.05.
#' @param epsilon Numeric. Small value for handling open/closed interval
#'   boundaries at the switching value. Default is 0.
#'
#' @return Numeric vector of length 2 giving the confidence interval
#'   \code{c(lower_limit, upper_limit)}.
#'
#' @details
#' This function inverts the acceptance regions computed by
#' \code{\link{nsnc_mcp_ar}} to construct confidence intervals. The procedure
#' partitions the real line according to the observed value \code{y}, and
#' uses root-finding to determine the corresponding confidence limits.
#'
#' When \code{r_l = r_u}, the procedure reverts to a symmetric extension of
#' the shortest acceptance region. The \code{epsilon} parameter allows for
#' fine control over whether interval boundaries are open or closed at the
#' switching value.
#'
#' @seealso \code{\link{nsnc_mcp_ar}}, \code{\link{Shortest.AR}}
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

  # ---- input checks ----
  if (!is.numeric(y) || length(y) != 1L || !is.finite(y)) {
    stop("`y` must be a finite numeric scalar.")
  }
  if (!is.numeric(r_l) || length(r_l) != 1L ||
      !is.numeric(r_u) || length(r_u) != 1L) {
    stop("`r_l` and `r_u` must be numeric scalars.")
  }
  if (r_l <= 0 || r_u <= 0) {
    stop("`r_l` and `r_u` must be positive.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).")
  }
  if (!is.numeric(ct) || length(ct) != 1L || ct <= 0) {
    stop("`ct` must be a positive numeric scalar.")
  }

  internal_eps <- 10^-11

  # ---- calculate required parameters ----
  theta_tilde_1_up   <- find_tilde_1(ct, r_u, alpha)
  theta_tilde_1_down <- -find_tilde_1(ct, r_l, alpha)

  shortest_ar <- Shortest.AR(switch_val, ct, alpha)

  # ---- inversion functions ----
  find_lower_l1_left   <- function(theta, r) l_1_finder_left_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l1_left   <- function(theta, r) l_1_finder_left_of_crit(theta, r, ct, alpha)[1] - y

  find_lower_l1_right  <- function(theta, r) l_1_finder_right_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l1_right  <- function(theta, r) l_1_finder_right_of_crit(theta, r, ct, alpha)[1] - y

  find_lower_l2_left   <- function(theta, r) l_2_finder_left_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l2_left   <- function(theta, r) l_2_finder_left_of_crit(theta, r, ct, alpha)[1] - y

  find_lower_l2_right  <- function(theta, r) l_2_finder_right_of_crit(theta, r, ct, alpha)[2] - y
  find_upper_l2_right  <- function(theta, r) l_2_finder_right_of_crit(theta, r, ct, alpha)[1] - y

  # ---- critical values ----
  c1 <- l_2_finder_left_of_crit(theta_tilde_1_down - internal_eps, r_l, ct, alpha)[1]
  c2 <- l_1_finder_left_of_crit(switch_val - internal_eps, r_l, ct, alpha)[1]
  c3 <- min(shortest_ar$A, na.rm = TRUE)
  c4 <- l_1_finder_right_of_crit(switch_val + internal_eps, r_u, ct, alpha)[1]

  c5 <- l_1_finder_left_of_crit(switch_val - internal_eps, r_l, ct, alpha)[2]
  c6 <- max(shortest_ar$A, na.rm = TRUE)
  c7 <- l_1_finder_right_of_crit(switch_val + internal_eps, r_u, ct, alpha)[2]
  c8 <- l_2_finder_right_of_crit(theta_tilde_1_up + internal_eps, r_u, ct, alpha)[2]

  # ---- find CI by region ----
  if (y < c1) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, theta_tilde_1_down - internal_eps), r = r_l)$root
    ul <- uniroot(find_upper_l2_left, c(y, theta_tilde_1_down - internal_eps), r = r_l)$root

  } else if (y >= c1 && y < c2) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, y), r = r_l)$root
    ul <- uniroot(find_upper_l1_left, c(theta_tilde_1_down + internal_eps, switch_val), r = r_l)$root

  } else if (y >= c2 && y < c3) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, -3 * y), r = r_l)$root
    ul <- switch_val - epsilon

  } else if (y >= c3 && y < c4) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, -3 * y), r = r_l)$root
    ul <- switch_val

  } else if (y >= c4 && y <= -ct) {
    ll <- uniroot(find_lower_l2_left, c(3 * y, theta_tilde_1_down), r = r_l)$root
    ul <- uniroot(find_upper_l1_right,
                  c(-theta_tilde_1_up - internal_eps, -theta_tilde_1_down),
                  r = r_u)$root

  } else if (y >= ct && y < c5) {
    ll <- uniroot(find_lower_l1_left,
                  c(theta_tilde_1_down + internal_eps, theta_tilde_1_up),
                  r = r_l)$root
    ul <- uniroot(find_upper_l2_right,
                  c(theta_tilde_1_up + internal_eps, 3 * y),
                  r = r_u)$root

  } else if (y >= c5 && y < c6) {
    ll <- switch_val
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root

  } else if (y >= c6 && y < c7) {
    ll <- switch_val + epsilon
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root

  } else if (y >= c7 && y < c8) {
    ll <- uniroot(find_lower_l1_right,
                  c(-theta_tilde_1_up, theta_tilde_1_up - internal_eps),
                  r = r_u)$root
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root

  } else if (y >= c8) {
    ll <- uniroot(find_lower_l2_right,
                  c(theta_tilde_1_up + internal_eps, y),
                  r = r_u)$root
    ul <- uniroot(find_upper_l2_right, c(y, 3 * y), r = r_u)$root

  } else {
    stop("`y` did not fall into any region; check critical values and logic.")
  }

  c(ll, ul)
}
