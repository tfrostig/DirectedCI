shortest_marginal_ar <- function(theta, alpha) {
  return(c(theta + qnorm(alpha / 2), theta - qnorm(alpha / 2)))
}

