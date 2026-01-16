#' @title Direction Preferring Marginal Confidence Intervals and Acceptance Regions
#' @description Direction-preferring CI and AR for unconditional inference.
#'   These methods favor inference in a preferred direction (positive or negative).

#' Direction Preferring Marginal Acceptance Region
#'
#' Computes an acceptance region for a normal mean with asymmetric tail allocation.
#' The asymmetry is governed by relative length parameters r_l and r_u.
#'
#' @param theta Numeric. Parameter value (on Z-scale).
#' @param r_l Numeric in [0, 1]. Relative length factor for lower side.
#' @param r_u Numeric in [0, 1]. Relative length factor for upper side. Default 1.
#' @param switch_val Numeric. Switching point between allocation rules. Default 0.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return List with components:
#'   \item{ar}{Numeric vector of length 4: c(ll, ul, lr, ur).
#'     For single interval, two values will be NA.}
#'   \item{length}{Numeric. Total length of the acceptance region.}
#'   \item{beta}{Numeric. Lower-tail type I error allocation.}
#'
#' @examples
#' dp_marginal_ar(theta = 0, r_l = 0.5, alpha = 0.05)
#' dp_marginal_ar(theta = 1, r_l = 0.5, r_u = 1, alpha = 0.05)
#'
#' @export
dp_marginal_ar <- function(theta, r_l, r_u = 1, switch_val = 0, alpha = 0.05) {
  interval_length <- function(beta) {
    qnorm(1 - (alpha - beta)) - qnorm(beta)
  }

  eps <- .Machine$double.eps
  z_alpha2 <- qnorm(1 - alpha / 2)
  target_len_l <- r_l * 2 * z_alpha2
  target_len_u <- r_u * 2 * z_alpha2

  if (theta < switch_val) {
    # Lower-side allocation
    if (theta < -z_alpha2 || r_l == 1) {
      beta <- alpha / 2
    } else {
      find_beta <- function(beta) interval_length(beta) - target_len_l
      beta <- uniroot(find_beta, c(eps, alpha / 2), tol = eps)$root
    }
  } else if (theta > switch_val) {
    # Upper-side allocation
    if (theta > z_alpha2 || r_u == 1) {
      beta <- alpha / 2
    } else {
      find_beta <- function(beta) interval_length(beta) - target_len_u
      beta_left <- uniroot(find_beta, c(eps, alpha / 2 - eps), tol = eps)$root
      beta <- alpha - beta_left
    }
  } else {
    # theta == switch_val
    if (r_l == r_u) {
      beta <- alpha / 2
    } else {
      if (theta < -z_alpha2 || r_l == 1) {
        beta <- alpha / 2
      } else {
        find_beta <- function(beta) interval_length(beta) - target_len_l
        beta <- uniroot(find_beta, c(eps, alpha / 2 - eps), tol = eps)$root
      }
    }
  }

  lower <- theta + qnorm(beta)
  upper <- theta + qnorm(1 - (alpha - beta))

  list(
    ar = c(NA_real_, NA_real_, lower, upper),
    length = upper - lower,
    beta = beta
  )
}

#' Direction Preferring Marginal Confidence Interval (Positive Direction)
#'
#' Computes a direction-preferring CI for a normal mean favoring positive values.
#' Internal function operating on Z-scale.
#'
#' @param y Numeric. Observed test statistic (on Z-scale).
#' @param r_l Numeric in [0, 1]. Relative length factor for lower side.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur).
#'   For single interval, returns c(NA, NA, lower, upper).
#'
#' @examples
#' pp_marginal_ci(y = 2.5, r_l = 0.5, alpha = 0.05)
#'
#' @export
pp_marginal_ci <- function(y, r_l, alpha = 0.05) {
  r_u <- 1
  switch_val <- 0
  epsilon <- 1e-15

  ar_func <- function(theta) {
    dp_marginal_ar(theta, r_l = r_l, r_u = r_u,
                   switch_val = switch_val, alpha = alpha)
  }

  z_alpha2 <- qnorm(1 - alpha / 2)

  # Boundary where the lower bound coincides with the usual symmetric interval
  c0 <- ar_func(qnorm(alpha / 2) + epsilon)$ar[3]

  # Behaviour just below and above the switch value
  left <- ar_func(switch_val - epsilon)
  right <- ar_func(switch_val + epsilon)

  c1 <- left$ar[3]
  c3 <- left$ar[4]
  l_switch_val_beta <- left$beta

  c2 <- right$ar[3]
  c4 <- right$ar[4]
  u_switch_val_beta <- right$beta

  if (y < c0) {
    ci <- c(y - z_alpha2, y - qnorm(alpha / 2))
  } else if (y >= c0 && y < c1) {
    ci <- c(y - z_alpha2, y + qnorm(1 - l_switch_val_beta))
  } else if (y >= c1 && y < c2) {
    ci <- c(y - z_alpha2, 0)  # (..., 0]
  } else if (y >= c2 && y <= switch_val) {
    ci <- c(y - z_alpha2, y + qnorm(1 - u_switch_val_beta))
  } else if (y > switch_val && y < c3) {
    ci <- c(y - qnorm(1 - (alpha - l_switch_val_beta)), y + z_alpha2)
  } else if (y >= c3 && y < c4) {
    ci <- c(epsilon, y + z_alpha2)  # (0, ...)
  } else {
    # y >= c4
    ci <- c(y - z_alpha2, y - qnorm(alpha / 2))
  }

  c(NA_real_, NA_real_, ci[1], ci[2])
}

#' Direction Preferring Marginal Confidence Interval (Negative Direction)
#'
#' Computes a direction-preferring CI for a normal mean favoring negative values.
#' Uses reflection of pp_marginal_ci.
#'
#' @param y Numeric. Observed test statistic (on Z-scale).
#' @param r_l Numeric in [0, 1]. Relative length factor.
#' @param alpha Numeric in (0, 1). Significance level. Default 0.05.
#'
#' @return Numeric vector of length 4: c(ll, ul, lr, ur).
#'
#' @export
np_marginal_ci <- function(y, r_l, alpha = 0.05) {
  ci_pos <- pp_marginal_ci(-y, r_l, alpha)
  c(NA_real_, NA_real_, -ci_pos[4], -ci_pos[3])
}
