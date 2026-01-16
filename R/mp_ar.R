mp_ar <- function(theta, r, alpha) {
  f <- function(beta) {
    qnorm(1 - (alpha - beta)) - qnorm(beta) - r * (2 * qnorm(1 - alpha / 2))
  }

  if (theta < 0) {
    beta_star <- uniroot(
      f,
      c(.Machine$double.eps, alpha / 2 - .Machine$double.eps),
      tol = .Machine$double.eps
    )$root
    return(c(theta + qnorm(beta_star), theta + qnorm(1 - (alpha - beta_star))))
  }

  if (theta == 0) {
    return(c(theta + qnorm(alpha / 2), theta - qnorm(alpha / 2)))
  }

  if (theta > 0) {
    beta_left <- uniroot(
      f,
      c(.Machine$double.eps, alpha / 2 - .Machine$double.eps),
      tol = .Machine$double.eps
    )$root
    beta_star <- alpha - beta_left
    return(c(theta + qnorm(beta_star), theta + qnorm(1 - (alpha - beta_star))))
  }

}
