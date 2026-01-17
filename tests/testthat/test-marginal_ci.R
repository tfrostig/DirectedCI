# =====================================================================
# Test suite for Marginal Confidence Interval functions
# =====================================================================

library(testthat)

# Helper small tolerance
tol <- 1e-6

# =====================================================================
# 1. Input validation tests
# =====================================================================

test_that("direction_preferring_marginal_ci rejects invalid inputs", {

  # y must be finite

expect_error(direction_preferring_marginal_ci(
    y = NA, sd = 1, r = 1.3, alpha = 0.05, direction = "positive"
  ))

  # sd must be positive
  expect_error(direction_preferring_marginal_ci(
    y = 2, sd = -1, r = 1.3, alpha = 0.05, direction = "positive"
  ))

  # r must be >= 1
  expect_error(direction_preferring_marginal_ci(
    y = 2, sd = 1, r = 0.5, alpha = 0.05, direction = "positive"
  ))

  # alpha must be in (0, 1)
  expect_error(direction_preferring_marginal_ci(
    y = 2, sd = 1, r = 1.3, alpha = 2, direction = "positive"
  ))

  # direction must be valid
  expect_error(direction_preferring_marginal_ci(
    y = 2, sd = 1, r = 1.3, alpha = 0.05, direction = "wrong"
  ))
})

test_that("modified_pratt_marginal_ci rejects invalid inputs", {

  expect_error(modified_pratt_marginal_ci(
    y = NA, sd = 1, r = 1.3, alpha = 0.05, direction = "positive"
  ))

  expect_error(modified_pratt_marginal_ci(
    y = 2, sd = -1, r = 1.3, alpha = 0.05, direction = "positive"
  ))

  expect_error(modified_pratt_marginal_ci(
    y = 2, sd = 1, r = 0.5, alpha = 0.05, direction = "positive"
  ))
})

test_that("shortest_marginal_ci_wrapper rejects invalid inputs", {

  expect_error(shortest_marginal_ci_wrapper(y = NA, sd = 1, alpha = 0.05))
  expect_error(shortest_marginal_ci_wrapper(y = 2, sd = -1, alpha = 0.05))
  expect_error(shortest_marginal_ci_wrapper(y = 2, sd = 1, alpha = 1.5))
})

# =====================================================================
# 2. Basic correctness tests - positive vs negative symmetry
# =====================================================================

test_that("dp vs dn marginal CI symmetry holds", {

  y <- 2.5
  alpha <- 0.05
  r <- 1.5

  ci_pos <- dp_marginal_ci(y, r_l = r, alpha = alpha)
  ci_neg <- dn_marginal_ci(y, r_l = r, alpha = alpha)

  # dn_marginal_ci(y) should equal reflected dp_marginal_ci(-y)
  ci_pos_reflected <- dp_marginal_ci(-y, r_l = r, alpha = alpha)
  ci_refl <- c(NA_real_, NA_real_, -ci_pos_reflected[4], -ci_pos_reflected[3])

  expect_equal(ci_neg, ci_refl, tolerance = tol)

  # bounds should be ordered (positions 3 and 4)
  expect_lte(ci_neg[3], ci_neg[4])
})

# =====================================================================
# 3. Wrapper standardization-rescaling correctness
# =====================================================================

test_that("direction_preferring_marginal_ci standardization + rescaling is correct", {

  y <- 4.6
  sd <- 2
  r <- 1.5
  alpha <- 0.05

  # Manual standardization
  z <- y / sd
  ci_z <- dp_marginal_ci(z, r_l = r, alpha = alpha)
  ci_manual <- ci_z
  ci_manual[3:4] <- sd * ci_z[3:4]

  ci_wrapper <- direction_preferring_marginal_ci(
    y = y,
    sd = sd,
    r = r,
    alpha = alpha,
    direction = "positive"
  )

  expect_equal(ci_wrapper[3:4], ci_manual[3:4], tolerance = tol)
})

test_that("modified_pratt_marginal_ci standardization + rescaling is correct", {

  y <- 3.0
  sd <- 1.5
  r <- 1.3
  alpha <- 0.05

  # Manual standardization
  z <- y / sd
  ci_z <- mp_marginal_ci(z, r = r, alpha = alpha)
  ci_manual <- ci_z
  ci_manual[3:4] <- sd * ci_z[3:4]

  ci_wrapper <- modified_pratt_marginal_ci(
    y = y,
    sd = sd,
    r = r,
    alpha = alpha,
    direction = "positive"
  )

  expect_equal(ci_wrapper[3:4], ci_manual[3:4], tolerance = tol)
})

test_that("shortest_marginal_ci_wrapper standardization + rescaling is correct", {

  y <- 2.5
  sd <- 2.0
  alpha <- 0.05

  # Manual standardization
  z <- y / sd
  ci_z <- shortest_marginal_ci(z, alpha = alpha)
  ci_manual <- ci_z
  ci_manual[3:4] <- sd * ci_z[3:4]

  ci_wrapper <- shortest_marginal_ci_wrapper(
    y = y,
    sd = sd,
    alpha = alpha
  )

  expect_equal(ci_wrapper[3:4], ci_manual[3:4], tolerance = tol)
})

# =====================================================================
# 4. Negative direction wrapper test
# =====================================================================

test_that("negative direction wrapper matches reflected positive CI", {

  y <- -3
  sd <- 1
  r <- 1.5
  alpha <- 0.05

  ci_neg <- direction_preferring_marginal_ci(
    y = y, sd = sd, r = r, alpha = alpha, direction = "negative"
  )

  # Manual reflection using dn_marginal_ci
  z <- y / sd
  ci_dn <- dn_marginal_ci(z, r_l = r, alpha = alpha)
  ci_manual <- ci_dn
  ci_manual[3:4] <- sd * ci_dn[3:4]

  expect_equal(ci_neg[3:4], ci_manual[3:4], tolerance = tol)
})

# =====================================================================
# 5. Ordering and finite output tests
# =====================================================================

test_that("Marginal CI bounds are ordered and finite", {

  y_vals <- c(-3, -1, 0, 1, 3)
  sd <- 1
  r <- 1.5
  alpha <- 0.05

  for (y in y_vals) {
    # Direction preferring
    ci_dp <- direction_preferring_marginal_ci(
      y, sd, r, alpha, direction = "positive"
    )
    expect_true(is.finite(ci_dp[3]))
    expect_true(is.finite(ci_dp[4]))
    expect_lte(ci_dp[3], ci_dp[4])

    # Modified Pratt
    ci_mp <- modified_pratt_marginal_ci(
      y, sd, r = 1.3, alpha, direction = "positive"
    )
    expect_true(is.finite(ci_mp[3]))
    expect_true(is.finite(ci_mp[4]))
    expect_lte(ci_mp[3], ci_mp[4])

    # Shortest
    ci_short <- shortest_marginal_ci_wrapper(y, sd, alpha)
    expect_true(is.finite(ci_short[3]))
    expect_true(is.finite(ci_short[4]))
    expect_lte(ci_short[3], ci_short[4])
  }
})

# =====================================================================
# 6. Shortest CI is symmetric around y
# =====================================================================

test_that("Shortest marginal CI is symmetric around y", {

  y <- 2.5
  sd <- 1.0
  alpha <- 0.05

  ci <- shortest_marginal_ci_wrapper(y, sd, alpha)

  # CI should be y +/- z_{1-alpha/2} * sd
  z_alpha <- qnorm(1 - alpha/2)
  expected_lower <- y - z_alpha * sd
  expected_upper <- y + z_alpha * sd

  expect_equal(ci[3], expected_lower, tolerance = tol)
  expect_equal(ci[4], expected_upper, tolerance = tol)
})

# =====================================================================
# 7. Coverage property check (basic)
# =====================================================================

test_that("True parameter is contained in CI for various y values", {

  # For a 95% CI, the true theta should be in the CI 95% of the time
  # Here we just check that the CI construction is valid

  theta <- 1.0  # True parameter
  sd <- 1.0
  alpha <- 0.05
  r <- 1.5

  # Generate y from N(theta, sd^2)
  set.seed(123)
  y_samples <- rnorm(100, mean = theta, sd = sd)

  # Check that a reasonable proportion contain theta
  coverage_dp <- mean(sapply(y_samples, function(y) {
    ci <- direction_preferring_marginal_ci(y, sd, r, alpha, "positive")
    theta >= ci[3] && theta <= ci[4]
  }))

  coverage_short <- mean(sapply(y_samples, function(y) {
    ci <- shortest_marginal_ci_wrapper(y, sd, alpha)
    theta >= ci[3] && theta <= ci[4]
  }))

  # Coverage should be close to nominal (with some Monte Carlo error)
  expect_gt(coverage_dp, 0.85)
  expect_gt(coverage_short, 0.85)
})

# =====================================================================
# 8. Extreme values should not crash
# =====================================================================

test_that("Extremely large y values produce finite output", {

  y <- 20
  sd <- 1
  r <- 1.5
  alpha <- 0.05

  ci_dp <- direction_preferring_marginal_ci(y, sd, r, alpha, "positive")
  expect_true(is.finite(ci_dp[3]))
  expect_true(is.finite(ci_dp[4]))
  expect_lte(ci_dp[3], ci_dp[4])

  ci_mp <- modified_pratt_marginal_ci(y, sd, r = 1.3, alpha, "positive")
  expect_true(is.finite(ci_mp[3]))
  expect_true(is.finite(ci_mp[4]))
  expect_lte(ci_mp[3], ci_mp[4])

  ci_short <- shortest_marginal_ci_wrapper(y, sd, alpha)
  expect_true(is.finite(ci_short[3]))
  expect_true(is.finite(ci_short[4]))
  expect_lte(ci_short[3], ci_short[4])
})

test_that("Extremely negative y values produce finite output", {

  y <- -20
  sd <- 1
  r <- 1.5
  alpha <- 0.05

  ci_dp <- direction_preferring_marginal_ci(y, sd, r, alpha, "positive")
  expect_true(is.finite(ci_dp[3]))
  expect_true(is.finite(ci_dp[4]))
  expect_lte(ci_dp[3], ci_dp[4])

  ci_mp <- modified_pratt_marginal_ci(y, sd, r = 1.3, alpha, "positive")
  expect_true(is.finite(ci_mp[3]))
  expect_true(is.finite(ci_mp[4]))
  expect_lte(ci_mp[3], ci_mp[4])

  ci_short <- shortest_marginal_ci_wrapper(y, sd, alpha)
  expect_true(is.finite(ci_short[3]))
  expect_true(is.finite(ci_short[4]))
  expect_lte(ci_short[3], ci_short[4])
})

# =====================================================================
# 9. r = 1 should give symmetric (shortest) CI for DP
# =====================================================================
test_that("DP marginal CI with r=1 equals shortest CI", {

  y <- 2.0
  alpha <- 0.05

  ci_dp <- dp_marginal_ci(y, r_l = 1, alpha = alpha)
  ci_short <- shortest_marginal_ci(y, alpha = alpha)

  expect_equal(ci_dp[3:4], ci_short[3:4], tolerance = tol)
})
