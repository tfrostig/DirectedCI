mp_ci <- function(y, r, alpha) {
  if (y > 0) {
    return(sort(-mp_ci(-y, r, alpha)))
  }
  f <- function(beta) {
    qnorm(1 - (alpha - beta)) - qnorm(beta) - r * (2 * qnorm(1 - alpha / 2))
  }

  # important values
  beta_star <- uniroot(f, c(0, alpha / 2))$root
  neg_threshold_1 <- -qnorm(1 - alpha / 2)
  neg_threshold_2 <- -qnorm(1 - (alpha - beta_star))


  if (y < neg_threshold_1) {
    return(c(y - qnorm(1 - (alpha - beta_star)), y - qnorm(beta_star)))
  }
  if (y >= neg_threshold_1 & y < neg_threshold_2) {
    return(c(y - qnorm(1 - (alpha - beta_star)),  0)) ## open
  }
  if (y >= neg_threshold_2 & y <= 0) {
    return(c(y + qnorm(alpha - beta_star), y - qnorm(alpha - beta_star)))
  }

}
