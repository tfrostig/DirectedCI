#' @title Shared Utilities for DirectedCI
#' @description Common constants and helper functions used across the package.

#' Small epsilon for numerical stability
#' @export
EPS <- sqrt(.Machine$double.eps)

#' Validate common CI/AR parameters
#'
#' @param y Observed statistic (for CI functions)
#' @param theta Parameter value (for AR functions)
#' @param alpha Significance level
#' @param r Inflation/relative length parameter
#' @param sd Standard deviation (for user wrappers)
#' @param direction Direction preference
#'
#' @return NULL (invisibly). Stops with error if validation fails
#' @keywords internal
validate_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).")
  }
  invisible(NULL)
}

validate_y <- function(y) {
  if (!is.numeric(y) || length(y) != 1 || !is.finite(y)) {
    stop("y must be a single finite numeric value.")
  }
  invisible(NULL)
}

validate_theta <- function(theta) {
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) {
    stop("theta must be a single finite numeric value.")
  }
  invisible(NULL)
}

validate_sd <- function(sd) {
  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0) {
    stop("sd must be a single positive numeric value.")
  }
  invisible(NULL)
}

validate_ct <- function(ct) {
  if (!is.numeric(ct) || length(ct) != 1 || !is.finite(ct) || ct <= 0) {
    stop("ct must be a single positive numeric value.")
  }
  invisible(NULL)
}

validate_r_conditional <- function(r) {
  if (!is.numeric(r) || length(r) != 1 || r < 1) {
    stop("r (inflation factor) must be a numeric value >= 1.")
  }
  invisible(NULL)
}

validate_r_marginal <- function(r) {
  if (!is.numeric(r) || length(r) != 1 || r < 1) {
    stop("r (inflation factor) must be a numeric value >= 1.")
  }
  invisible(NULL)
}

validate_direction <- function(direction) {
  direction <- match.arg(direction, choices = c("positive", "negative"))
  direction
}
