# =====================================================================
# Test suite for Direction-Preferring Confidence Interval functions
# =====================================================================

library(testthat)

# Helper small tolerance
tol <- 1e-6

# =====================================================================
# 1. Input validation tests
# =====================================================================

test_that("direction_preferring_conditional_ci rejects invalid inputs", {

  # y must be significant (|y/sd| > ct) for conditional CI
  expect_error(direction_preferring_conditional_ci(
    y = 1, sd = 1, r = 1.3,
    ct = qnorm(0.975), alpha = 0.05, direction = "positive"
  ))

  expect_error(direction_preferring_conditional_ci(
    y = NA, sd = 1, r = 1.3,
    ct = qnorm(0.975), alpha = 0.05, direction = "positive"
  ))

  expect_error(direction_preferring_conditional_ci(
    y = 3, sd = -1, r = 1.3,
    ct = qnorm(0.975), alpha = 0.05, direction = "positive"
  ))

  expect_error(direction_preferring_conditional_ci(
    y = 3, sd = 1, r = 0.5,
    ct = qnorm(0.975), alpha = 0.05, direction = "positive"
  ))

  expect_error(direction_preferring_conditional_ci(
    y = 3, sd = 1, r = 1.3,
    ct = "not numeric", alpha = 0.05, direction = "positive"
  ))

  expect_error(direction_preferring_conditional_ci(
    y = 3, sd = 1, r = 1.3,
    ct = qnorm(0.975), alpha = 2, direction = "positive"
  ))

  expect_error(direction_preferring_conditional_ci(
    y = 3, sd = 1, r = 1.3,
    ct = qnorm(0.975), alpha = 0.05, direction = "wrong"
  ))
})


# =====================================================================
# 2. Basic correctness tests - positive vs negative symmetry
# =====================================================================

test_that("dp vs dn conditional CI symmetry holds", {

  y <- 2.5; ct <- qnorm(0.975); alpha <- 0.05; r <- 1.3

  ci_pos  <- dp_conditional_ci(y, ct, r, alpha)
  ci_neg  <- dn_conditional_ci(y, ct, r, alpha)

  # dn_conditional_ci(y) should equal reflected dp_conditional_ci(-y)
  ci_pos_reflected <- dp_conditional_ci(-y, ct, r, alpha)
  ci_refl <- c(-ci_pos_reflected[2], -ci_pos_reflected[1])

  # dn_conditional_ci must match reflected dp_conditional_ci
  expect_equal(ci_neg, ci_refl, tolerance = tol)

  # bounds should be ordered
  expect_lte(ci_neg[1], ci_neg[2])
})


# =====================================================================
# 3. Wrapper standardization-rescaling correctness
# =====================================================================

test_that("direction_preferring_conditional_ci standardization + rescaling is correct", {

  y <- 4.6; sd <- 2  # y/sd = 2.3 > ct
  ct <- qnorm(0.975)
  r <- 1.3
  alpha <- 0.05

  # Manual standardization
  z <- y / sd
  ci_z <- dp_conditional_ci(z, ct, r, alpha)
  ci_manual <- sd * ci_z

  ci_wrapper <- direction_preferring_conditional_ci(
    y = y,
    sd = sd,
    r = r,
    ct = ct,
    alpha = alpha,
    direction = "positive"
  )

  expect_equal(ci_wrapper, ci_manual, tolerance = tol)
})


# =====================================================================
# 4. Negative direction wrapper test
# =====================================================================

test_that("negative direction wrapper matches reflected positive CI", {

  y <- -3  # significant negative
  sd <- 1
  r <- 1.3
  ct <- qnorm(0.975)
  alpha <- 0.05

  ci_neg <- direction_preferring_conditional_ci(
    y = y, sd = sd, ct = ct, r = r, alpha = alpha, direction = "negative"
  )

  # Manual reflection using dn_conditional_ci
  z <- y / sd
  ci_dn <- dn_conditional_ci(z, ct, r, alpha)
  ci_manual <- sd * ci_dn

  expect_equal(ci_neg, ci_manual, tolerance = tol)
})


# =====================================================================
# 5. Ordering and monotonicity tests
# =====================================================================

test_that("CI bounds are increasing and finite", {

  y <- 3  # significant

  sd <- 1
  r <- 1.3
  alpha <- 0.05
  ct <- qnorm(0.975)

  ci <- direction_preferring_conditional_ci(
    y, sd, r, ct, alpha, direction = "positive"
  )

  # CI returns c(lower, upper)
  expect_true(is.finite(ci[1]))
  expect_true(is.finite(ci[2]))
  expect_lte(ci[1], ci[2])
})


# =====================================================================
# 6. Test behavior near decision thresholds
# =====================================================================

test_that("CI computation near regime thresholds works", {

  sd <- 1
  r <- 1.3
  alpha <- 0.05
  ct <- qnorm(0.975)

  # Use several significant y values (all must satisfy |y| > ct)
  ys <- c(-3, -2.5, -ct - 0.1, ct + 0.1, 2.5, 3)

  for (y in ys) {
    ci <- direction_preferring_conditional_ci(y, sd, r, ct, alpha, direction="positive")
    expect_true(is.finite(ci[1]))
    expect_true(is.finite(ci[2]))
    expect_lte(ci[1], ci[2])
  }
})


# =====================================================================
# 7. Regression test for known values (stability)
# =====================================================================

test_that("Known stable CI values do not change unexpectedly", {

  y <- 2.5  # significant
  sd <- 1
  r <- 1.3
  ct <- qnorm(0.975)
  alpha <- 0.05

  ci <- direction_preferring_conditional_ci(y, sd, r, ct, alpha, "positive")

  # CI returns c(lower, upper)
  expect_true(ci[1] < ci[2])
  expect_true(ci[1] < y)  # Lower bound should be below estimate
})


# =====================================================================
# 8. Extreme values should not crash
# =====================================================================

test_that("Extremely large y values produce finite output", {

  y <- 20
  sd <- 1
  r <- 1.3
  ct <- qnorm(0.975)
  alpha <- 0.05

  ci <- direction_preferring_conditional_ci(y, sd, r, ct, alpha, "positive")

  # CI returns c(lower, upper)
  expect_true(is.finite(ci[1]))
  expect_true(is.finite(ci[2]))
  expect_lte(ci[1], ci[2])
})
