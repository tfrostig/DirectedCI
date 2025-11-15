# CDF of N(mu, 1) truncated outside (-b, b)
cdf_tn_outer <- function(y, mu, b, sigma = 1, a = -b) {
  W <- pnorm((b - mu) / sigma) - pnorm((a - mu) / sigma)
  if (y <= a) {
    return(pnorm((y - mu) / sigma) / (1 - W))
  }
  if (y > a & y <= b) {
    return(pnorm((a - mu) / sigma) / (1 - W))
  }
  if (y > b) {
    return((pnorm((y - mu) / sigma) - W) / (1 - W))
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
