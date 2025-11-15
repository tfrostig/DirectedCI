
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
nsnc_mcp_ar <- function(theta, r_l, r_u, switch_val=0, ct=qnorm(0.975), alpha=0.05) {
  # calculate required parameters
  theta_tilde_1_up   <- find_tilde_1(ct, r_u, alpha)
  theta_tilde_1_down <- -find_tilde_1(ct, r_l, alpha)


  if (theta <= theta_tilde_1_down) {
    ar         <- range(Shortest.AR(theta, ct, alpha)$A, na.rm = TRUE)
  }
  if (theta > theta_tilde_1_down & theta <= switch_val) {
    ar         <- l_1_finder_left_of_crit(theta, r_l, ct, alpha)
  }
  if (theta <= theta_tilde_1_up & theta > switch_val) {
    ar         <- l_1_finder_right_of_crit(theta, r_u, ct, alpha)
  }
  if (theta >= theta_tilde_1_up) {
    ar         <- range(Shortest.AR(theta, ct, alpha)$A, na.rm = TRUE)
  }

  return(ar)

}
