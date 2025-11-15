#' Acceptance Region for a Normal Mean with Asymmetric Tail Allocation
#'
#' @description
#' Compute an acceptance region for a normal mean parameter \eqn{\theta} under
#' a two-sided test with total type I error \eqn{\alpha}, where the type I error
#' is split asymmetrically between the lower and upper tails.
#'
#' The asymmetry is governed by relative length parameters \code{r_l} and
#' \code{r_u} for the lower and upper parts of the confidence interval and a
#' switching value \code{switch_val} that separates the two regimes.
#'
#' @param theta Numeric scalar. Parameter value (on the Z-scale) for which the
#'   acceptance region is computed.
#' @param r_l Numeric scalar in \eqn{[0, 1]}. Relative length factor assigned to
#'   the lower side of the interval when \code{theta <= switch_val}.
#' @param r_u Numeric scalar in \eqn{[0, 1]}. Relative length factor assigned to
#'   the upper side of the interval when \code{theta >= switch_val}.
#' @param switch_val Numeric scalar. Switching point that determines whether the
#'   lower- or upper-side allocation rule is used.
#' @param alpha Numeric scalar in \eqn{(0, 1)}. Total type I error for the
#'   two-sided procedure.
#'
#' @details
#' The procedure partitions \eqn{\alpha} into \eqn{\beta} and
#' \eqn{\alpha - \beta}, which correspond to the lower and upper tail
#' probabilities of a normal acceptance region for \eqn{\theta}. The function
#' solves for \eqn{\beta} such that the length of the resulting (Z-scale)
#' interval matches a target length determined by \code{r_l}, \code{r_u}, and
#' the usual symmetric \eqn{1 - \alpha} normal confidence interval.
#'
#' When \code{theta} is below \code{switch_val}, the allocation is governed by
#' \code{r_l}. When \code{theta} is above \code{switch_val}, it is governed by
#' \code{r_u}. When \code{theta == switch_val} and \code{r_l == r_u}, the
#' procedure reverts to a symmetric split \eqn{\beta = \alpha / 2}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{ar}{Numeric vector of length 2 giving the acceptance region
#'     boundaries on the Z-scale for the given \code{theta}.}
#'   \item{beta}{Numeric scalar \eqn{\beta} giving the lower-tail type I error
#'     allocation used to define the region.}
#' }
#'
#' @examples
#' ar_ns_nc(0, r_l = 0.5, r_u = 1, switch_val = 0, alpha = 0.05)
#'
#' @export
ar_ns_nc <- function(theta, r_l, r_u, switch_val = 0, alpha = 0.05) {
  interval_length <- function(beta) {
    qnorm(1 - (alpha - beta)) - qnorm(beta)
  }

  # small offset to avoid qnorm(0) and qnorm(1)
  eps <- .Machine$double.eps^0.25
  z_alpha2 <- qnorm(1 - alpha / 2)
  target_len_l <- r_l * 2 * z_alpha2
  target_len_u <- r_u * 2 * z_alpha2

  if (theta < switch_val) {
    # Lower-side allocation
    if (theta < -z_alpha2 || r_l == 1) {
      beta <- alpha / 2
    } else {
      find_beta <- function(beta) interval_length(beta) - target_len_l
      beta <- uniroot(find_beta,
                      c(eps, alpha / 2 - eps),
                      tol = .Machine$double.eps)$root
    }
  } else if (theta > switch_val) {
    # Upper-side allocation
    if (theta > z_alpha2 || r_u == 1) {
      beta <- alpha / 2
    } else {
      find_beta <- function(beta) interval_length(beta) - target_len_u
      beta_left <- uniroot(find_beta,
                           c(eps, alpha / 2 - eps),
                           tol = .Machine$double.eps)$root
      beta <- alpha - beta_left
    }
  } else {  # theta == switch_val
    if (r_l == r_u) {
      beta <- alpha / 2
    } else {
      # Default: use lower-side allocation rule at the switch
      if (theta < -z_alpha2 || r_l == 1) {
        beta <- alpha / 2
      } else {
        find_beta <- function(beta) interval_length(beta) - target_len_l
        beta <- uniroot(find_beta,
                        c(eps, alpha / 2 - eps),
                        tol = .Machine$double.eps)$root
      }
    }
  }

  list(
    ar   = c(theta + qnorm(beta),
             theta + qnorm(1 - (alpha - beta))),
    beta = beta
  )
}


#' Direction-Preferring Confidence Interval for a Normal Mean
#'
#' @description
#' Construct a (potentially asymmetric) confidence interval for a normal mean
#' based on a direction-preferring acceptance region defined by
#' \code{\link{ar_ns_nc}}.
#'
#' The interval is derived by inverting the acceptance region as a function of
#' the observed statistic \code{y} on the Z-scale. The resulting interval can
#' be one-sided, two-sided, or truncated at 0 depending on the configuration of
#' \code{r_l}, \code{r_u}, and \code{switch_val}.
#'
#' @param y Numeric scalar. Observed value of the test statistic on the Z-scale.
#' @param r_l Numeric scalar in \eqn{[0, 1]}. Relative length factor for the
#'   lower side of the interval; see \code{\link{ar_ns_nc}}.
#' @param r_u Numeric scalar in \eqn{[0, 1]}. Relative length factor for the
#'   upper side of the interval. The current implementation assumes
#'   \code{r_u = 1} in its piecewise structure; other values are not fully
#'   supported and may require additional derivations.
#' @param switch_val Numeric scalar. Switching point for the allocation rule;
#'   see \code{\link{ar_ns_nc}}.
#' @param alpha Numeric scalar in \eqn{(0, 1)}. Total type I error for the
#'   two-sided procedure.
#' @param epsilon Numeric scalar. Small perturbation used to probe the
#'   acceptance region around \code{switch_val} and the lower boundary of the
#'   parameter space.
#'
#' @details
#' The function inverts the acceptance region induced by
#' \code{\link{ar_ns_nc}} by identifying critical values of \code{y} that
#' correspond to boundaries of the acceptance region as functions of
#' \eqn{\theta}. The resulting mapping yields a piecewise expression for the
#' confidence interval as a function of \code{y}.
#'
#' The current implementation is tailored to the case \code{r_u = 1}, i.e., the
#' upper side behaves like a standard symmetric interval beyond the switch
#' point. Using this function with \code{r_u != 1} is not recommended unless
#' the piecewise structure has been re-derived accordingly.
#'
#' @return
#' A numeric vector of length 2 giving the lower and upper bounds of the
#' confidence interval for the true mean on the Z-scale.
#'
#' @examples
#' ci_ns_nc(y = 1.2, r_l = 0.5, r_u = 1, switch_val = 0, alpha = 0.05)
#'
#' @export
ci_ns_nc <- function(y,
                     r_l,
                     r_u = 1,
                     switch_val = 0,
                     alpha = 0.05,
                     epsilon = 1e-15) {

  ar_func <- function(theta) {
    ar_ns_nc(theta, r_l = r_l, r_u = r_u,
             switch_val = switch_val, alpha = alpha)
  }

  z_alpha2 <- qnorm(1 - alpha / 2)

  # Boundary where the lower bound coincides with the usual symmetric interval
  c0 <- ar_func(qnorm(alpha / 2) + epsilon)$ar[1]

  # Behaviour just below and above the switch value
  left  <- ar_func(switch_val - epsilon)
  right <- ar_func(switch_val + epsilon)

  c1 <- left$ar[1]
  c3 <- left$ar[2]
  l_switch_val_beta <- left$beta

  c2 <- right$ar[1]
  c4 <- right$ar[2]
  u_switch_val_beta <- right$beta

  if (y < c0) {
    return(c(y - z_alpha2, y - qnorm(alpha / 2)))
  }
  if (y >= c0 && y < c1) {
    return(c(y - z_alpha2, y + qnorm(1 - l_switch_val_beta)))
  }
  if (y >= c1 && y < c2) {
    return(c(y - z_alpha2, 0)) # (..., 0]
  }
  if (y >= c2 && y <= switch_val) {
    return(c(y - z_alpha2, y + qnorm(1 - u_switch_val_beta)))
  }
  if (y >= switch_val && y < c3) {
    return(c(y - qnorm(1 - (alpha - l_switch_val_beta)), y + z_alpha2))
  }
  if (y >= c3 && y < c4) {
    return(c(epsilon, y + z_alpha2)) # (0, ...)
  }
  if (y >= c4) {
    return(c(y - z_alpha2, y - qnorm(alpha / 2)))
  }

  # Fallback (should not be reached)
  c(NA_real_, NA_real_)
}
