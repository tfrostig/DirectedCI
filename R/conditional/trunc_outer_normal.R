

# CDF of N(mu, 1) truncated outside (-b, b)
cdf_tn_outer <- function(y, mu, b, sigma = 1, a = -b) {
  # More stable calculation of 1 - W directly
  p_below_a <- pnorm((a - mu) / sigma)
  p_above_b <- pnorm((b - mu) / sigma, lower.tail = FALSE)
  one_minus_W <- p_below_a + p_above_b

  if (y <= a) {
    return(pnorm((y - mu) / sigma) / one_minus_W)
  }
  if (y > a & y <= b) {
    return(p_below_a / one_minus_W)
  }
  if (y > b) {
    # Reformulated: 1 - P(Y > y) / (1 - W)
    # Avoids subtracting two numbers close to 1
    p_above_y <- pnorm((y - mu) / sigma, lower.tail = FALSE)
    return(1 - p_above_y / one_minus_W)
  }
}

# Quantile for N(mu,1) truncated outside (-b, b)
quantile_tn_outer <- function(p, mu, b) {
  W <- pnorm((b - mu)) - pnorm((-b - mu))
  qc <- 1 - W
  if (p == 1) {
    return(Inf)
  }
  if (p > cdf_tn_outer(b, mu, b)) {
    return(qnorm(p * qc + 1 - qc) + mu)
  } else {
    return(qnorm(p * qc) + mu)
  }
}

# Density of N(mu,1) truncated outside (-b, b)
density_tn_outer <- function(x, mu, b) {
  W <- 1 - (pnorm(b - mu) - pnorm(-b - mu))
  dens <- dnorm(x - mu) / W
  dens[abs(x) < b] <- 0
  return(dens)
}
