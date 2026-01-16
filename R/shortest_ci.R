shortest_ci <- function(y, alpha) {
  c(y - qnorm(1 - alpha / 2)
    y + qnorm(1 - alpha / 2))
}
