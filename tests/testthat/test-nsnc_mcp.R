# Test file: tests/testthat/test-nsnc_mcp.R

# Test shortest_ar function
test_that("shortest_ar returns correct structure", {
  result <- shortest_ar(theta = 0, ct = qnorm(0.975), alpha = 0.05)

  expect_type(result, "list")
  expect_named(result, c("A", "l"))
  expect_length(result$A, 4)
  expect_type(result$l, "double")
})

test_that("shortest_ar handles theta = 0", {
  result <- shortest_ar(theta = 0)

  # At theta = 0, should be symmetric
  expect_equal(result$A[1], -result$A[4], tolerance = 1e-10)
  expect_equal(result$A[2], -result$A[3], tolerance = 1e-10)

  # Length should be positive
  expect_true(result$l > 0)
})

test_that("shortest_ar handles positive theta", {
  result <- shortest_ar(theta = 1, ct = qnorm(0.975), alpha = 0.05)

  expect_type(result, "list")
  expect_true(result$l > 0)

  # Check that AR bounds are ordered correctly (when not NA)
  A_no_na <- result$A[!is.na(result$A)]
  expect_true(all(diff(A_no_na) > 0))
})

test_that("shortest_ar handles negative theta symmetrically", {
  result_pos <- shortest_ar(theta = 1.5)
  result_neg <- shortest_ar(theta = -1.5)

  # Results should be symmetric
  expect_equal(result_pos$A[1], -result_neg$A[4], tolerance = 1e-10)
  expect_equal(result_pos$A[2], -result_neg$A[3], tolerance = 1e-10)
  expect_equal(result_pos$A[3], -result_neg$A[2], tolerance = 1e-10)
  expect_equal(result_pos$A[4], -result_neg$A[1], tolerance = 1e-10)
  expect_equal(result_pos$l, result_neg$l, tolerance = 1e-10)
})

test_that("shortest_ar validates alpha range", {
  # Should work with valid alpha
  expect_no_error(shortest_ar(theta = 0, alpha = 0.05))
  expect_no_error(shortest_ar(theta = 0, alpha = 0.1))

  # Length should increase as alpha decreases (stricter test)
  l1 <- shortest_ar(theta = 1, alpha = 0.1)$l
  l2 <- shortest_ar(theta = 1, alpha = 0.05)$l
  expect_true(l2 > l1)
})


# Test truncated normal functions
test_that("cdf_tn_outer returns valid probabilities", {
  # Test at various points
  p1 <- cdf_tn_outer(-3, mu = 0, b = 1.96)
  p2 <- cdf_tn_outer(0, mu = 0, b = 1.96)
  p3 <- cdf_tn_outer(3, mu = 0, b = 1.96)

  expect_true(p1 >= 0 && p1 <= 1)
  expect_true(p2 >= 0 && p2 <= 1)
  expect_true(p3 >= 0 && p3 <= 1)

  # CDF should be monotone increasing
  expect_true(p1 < p3)

  # Inside truncation region, CDF should be constant
  p_left <- cdf_tn_outer(-1, mu = 0, b = 1.96)
  p_right <- cdf_tn_outer(1, mu = 0, b = 1.96)
  expect_equal(p_left, p_right, tolerance = 1e-10)
})

test_that("quantile_tn_outer is inverse of cdf_tn_outer", {
  mu <- 0
  b <- 1.96

  # Test at several probabilities
  probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)

  for (p in probs) {
    q <- quantile_tn_outer(p, mu, b)
    p_back <- cdf_tn_outer(q, mu, b)
    expect_equal(p, p_back, tolerance = 1e-6)
  }
})

test_that("density_tn_outer is zero inside truncation region", {
  mu <- 0
  b <- 1.96

  # Points inside (-b, b) should have zero density
  x_inside <- c(-1, -0.5, 0, 0.5, 1)
  dens_inside <- density_tn_outer(x_inside, mu, b)
  expect_equal(dens_inside, rep(0, 5))

  # Points outside should have positive density
  x_outside <- c(-3, -2, 2, 3)
  dens_outside <- density_tn_outer(x_outside, mu, b)
  expect_true(all(dens_outside > 0))
})

test_that("density_tn_outer integrates to approximately 1", {
  mu <- 0
  b <- 1.96

    # Numerical integration of the density
  f <- function(x) density_tn_outer(x, mu, b)

  # Left tail
  int_left <- integrate(f, -Inf, -b)$value
  # Right tail
  int_right <- integrate(f, b, Inf)$value

  total <- int_left + int_right
  expect_equal(total, 1, tolerance = 1e-3)
})


# Test nsnc_mcp_ar function
test_that("nsnc_mcp_ar returns valid interval", {
  ar <- nsnc_mcp_ar(theta = 0.5, r_l = 1.2, r_u = 1.2)

  expect_length(ar, 2)
  expect_true(ar[1] < ar[2])
  expect_type(ar, "double")
})

test_that("nsnc_mcp_ar with r=1 equals shortest_ar", {
  theta <- 1.5
  ct <- qnorm(0.975)
  alpha <- 0.05

  ar_nsnc <- nsnc_mcp_ar(theta, r_l = 1, r_u = 1, ct = ct, alpha = alpha)
  ar_short <- range(shortest_ar(theta, ct, alpha)$A, na.rm = TRUE)

  expect_equal(ar_nsnc, ar_short, tolerance = 1e-10)
})

test_that("nsnc_mcp_ar increases with extension ratios", {
  theta <- 1

  ar1 <- nsnc_mcp_ar(theta, r_l = 1.0, r_u = 1.0)
  ar2 <- nsnc_mcp_ar(theta, r_l = 1.5, r_u = 1.5)

  # Larger extension ratio should give wider interval
  width1 <- ar1[2] - ar1[1]
  width2 <- ar2[2] - ar2[1]
  expect_true(width2 > width1)
})

test_that("nsnc_mcp_ar handles asymmetric extensions", {
  theta <- 0.5

  ar_sym <- nsnc_mcp_ar(theta, r_l = 1.3, r_u = 1.3)
  ar_asym <- nsnc_mcp_ar(theta, r_l = 1.5, r_u = 1.1)

  # Both should be valid intervals
  expect_true(ar_sym[1] < ar_sym[2])
  expect_true(ar_asym[1] < ar_asym[2])

  # They should differ
  expect_false(isTRUE(all.equal(ar_sym, ar_asym)))
})


# Test nsnc_mcp_ci function
test_that("nsnc_mcp_ci returns valid confidence interval", {
  ci <- nsnc_mcp_ci(y = 2, r_l = 1.2, r_u = 1.2)

  expect_length(ci, 2)
  expect_true(ci[1] < ci[2])
  expect_type(ci, "double")
})

test_that("nsnc_mcp_ci coverage: theta in CI iff y in AR", {
  # If theta is the true parameter, and we observe y in its AR,
  # then theta should be in the CI for y

  theta <- 1.5
  r_l <- 1.2
  r_u <- 1.2
  ct <- qnorm(0.975)
  alpha <- 0.05

  # Get AR for theta
  ar <- nsnc_mcp_ar(theta, r_l, r_u, ct = ct, alpha = alpha)

  # Pick a y in the AR
  y <- mean(ar)

  # Get CI for y
  ci <- nsnc_mcp_ci(y, r_l, r_u, ct = ct, alpha = alpha)

  # theta should be in the CI
  expect_true(theta >= ci[1] && theta <= ci[2])
})

test_that("nsnc_mcp_ci handles extreme observed values", {
  # Very negative observation
  ci_neg <- nsnc_mcp_ci(y = -5, r_l = 1.2, r_u = 1.2)
  expect_true(ci_neg[1] < ci_neg[2])
  expect_true(ci_neg[2] < 0)  # CI should be in negative region

  # Very positive observation
  ci_pos <- nsnc_mcp_ci(y = 5, r_l = 1.2, r_u = 1.2)
  expect_true(ci_pos[1] < ci_pos[2])
  expect_true(ci_pos[1] > 0)  # CI should be in positive region
})

test_that("nsnc_mcp_ci width increases with extension ratios", {
  y <- 4

  ci1 <- nsnc_mcp_ci(y, r_l = 1.0, r_u = 1.0)
  ci2 <- nsnc_mcp_ci(y, r_l = 1, r_u = 1.3)

  width1 <- ci1[2] - ci1[1]
  width2 <- ci2[2] - ci2[1]

  # Larger extension should give wider CI
  expect_true(width2 > width1)
})

test_that("nsnc_mcp_ci respects significance level", {
  y <- 2

  # Different alpha should give different CI widths
  ci1 <- nsnc_mcp_ci(y, r_l = 1.2, r_u = 1.2, alpha = 0.10)
  ci2 <- nsnc_mcp_ci(y, r_l = 1.2, r_u = 1.2, alpha = 0.05)

  width1 <- ci1[2] - ci1[1]
  width2 <- ci2[2] - ci2[1]

  # Smaller alpha (more confidence) should give wider CI
  expect_true(width2 > width1)
})


# Test helper functions
test_that("theta_1_finder returns reasonable value", {
  ct <- qnorm(0.975)
  alpha <- 0.05

  theta1 <- theta_1_finder(ct, alpha)

  expect_type(theta1, "double")
  expect_true(theta1 > 0)
  expect_true(theta1 < ct + qnorm(1 - alpha))
})

test_that("find_tilde_1 returns theta1 when r=1", {
  ct <- qnorm(0.975)
  alpha <- 0.05

  theta1 <- theta_1_finder(ct, alpha)
  tilde1 <- find_tilde_1(ct, r = 1, alpha)

  expect_equal(theta1, tilde1, tolerance = 1e-10)
})

test_that("lb_finder and ub_finder are consistent", {
  theta <- 1
  ct    <- qnorm(0.975)
  alpha <- 0.05

  # Start with an upper bound
  ub <- 2.5

  # Find corresponding lower bound
  lb <- lb_finder(ub, theta, ct, alpha)

  # Find upper bound from that lower bound
  ub_back <- ub_finder(lb, theta, ct, alpha)

  # Should get back something close to original ub
  # (may not be exact due to the alpha constraint)
  expect_equal(ub, ub_back, tolerance = 0.1)
})


