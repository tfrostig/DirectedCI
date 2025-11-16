


test_that("shortest_ar returns correct structure", {
  crit   <- qnorm(0.975)
  result <- nsnc_mcp_ci(
    y = -crit - 10^-5,
    r_l = 1.01,
    r_u = 1,
    ct = qnorm(0.975),
    alpha = 0.05
  )
  print(result)



})
