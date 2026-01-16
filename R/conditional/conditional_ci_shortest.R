Shortest.CI <- function(x, ct, alpha) {
  stopifnot(length(x) == 1, length(ct) == 1, length(alpha) == 1)
  stopifnot(is.finite(x), is.finite(ct), is.finite(alpha))
  stopifnot(alpha > 0, alpha < 1)

  # helper: Q(theta)
  Q <- function(theta) {
    2 - pnorm(ct + theta) - pnorm(ct - theta)
  }

  # helper: robust uniroot with interval expansion
  safe_uniroot <- function(fun, lower, upper, expand = 1.5, max_expand = 60, ...) {
    fl <- fun(lower)
    fu <- fun(upper)

    k <- 0
    while (is.finite(fl) && is.finite(fu) && fl * fu > 0 && k < max_expand) {
      mid <- 0.5 * (lower + upper)
      width <- (upper - lower) * expand
      lower <- mid - 0.5 * width
      upper <- mid + 0.5 * width
      fl <- fun(lower)
      fu <- fun(upper)
      k <- k + 1
    }

    if (!is.finite(fl) || !is.finite(fu) || fl * fu > 0) {
      stop("safe_uniroot: failed to bracket a root for uniroot().")
    }

    uniroot(fun, c(lower, upper), ...)$root
  }

  # --- compute useful quantities ---
  f1 <- function(theta) {
    (pnorm(ct + theta) - pnorm(ct - theta)) - (1 - alpha) * Q(theta)
  }
  theta1 <- safe_uniroot(f1, 0, ct + qnorm(1 - alpha))

  f2 <- function(theta) {
    2 * pnorm(theta - ct) - 1 - (1 - alpha) * Q(theta)
  }
  theta2 <- safe_uniroot(f2, theta1, ct + qnorm(1 - alpha / 2))

  # not used later, kept as in your original
  x0 <- qnorm(1 - 0.5 * alpha * Q(0))

  x1 <- theta1 + qnorm(
    pnorm(ct - theta1) +
      (1 - alpha) * (1 - pnorm(ct + theta1) + 1 - pnorm(ct - theta1))
  )

  x2 <- 2 * theta2 - ct

  fR <- function(theta) {
    denom_arg <- pnorm(ct - theta) +
      (1 - alpha) * (1 - pnorm(ct + theta) + 1 - pnorm(ct - theta))
    denom_arg <- pmin(pmax(denom_arg, .Machine$double.eps), 1 - .Machine$double.eps)

    denom <- dnorm(qnorm(denom_arg))
    if (!is.finite(denom) || denom <= .Machine$double.eps) return(Inf)

    1 + (
      -dnorm(ct - theta) +
        (1 - alpha) * (-dnorm(ct + theta) + dnorm(ct - theta))
    ) / denom
  }

  R <- safe_uniroot(fR, theta1, theta2)

  # --- obtain CI ends ---
  is.neg <- (x < 0)
  x <- abs(x)

  # lower end
  lower <- NA_real_
  if (ct < x && x < x1) {
    fl <- function(theta) 2 * (1 - pnorm(x - theta)) - alpha * Q(theta)
    lower <- safe_uniroot(fl, -theta1 - 0.1, theta1)
  } else if (x1 < x && x < x2) {
    fl <- function(theta) pnorm(x - theta) - pnorm(ct - theta) - (1 - alpha) * Q(theta)
    lower <- safe_uniroot(fl, R, theta2)
  } else if (x >= x2) {
    fl <- function(theta) 2 * pnorm(x - theta) - 1 - (1 - alpha) * Q(theta)
    lower <- safe_uniroot(fl, theta2, x)
  } else {
    stop("x is not in any expected region; check ct/x relationship or assumptions.")
  }

  # upper end
  fu <- function(theta) 2 * pnorm(theta - x) - 1 - (1 - alpha) * Q(theta)
  upper <- safe_uniroot(fu, theta2, x + 2 * qnorm(1 - alpha / 2))

  CI <- list(lower = lower, upper = upper)
  if (is.neg) CI <- list(lower = -upper, upper = -lower)
  CI
}
