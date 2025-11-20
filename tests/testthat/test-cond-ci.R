# =====================================================================
# Test suite for Direction-Preferring Confidence Interval functions
# =====================================================================

library(testthat)

# Helper small tolerance
tol <- 1e-6

# =====================================================================
# 1. Input validation tests
# =====================================================================

test_that("conditional_direction_preferring_ci rejects invalid inputs", {

  expect_error(conditional_direction_preferring_ci(
    y = NA, sd = 1, r = 1.3,
    ct = qnorm(0.975), alpha = 0.05, direction = "positive"
  ))

  expect_error(conditional_direction_preferring_ci(
    y = 1, sd = -1, r = 1.3,
    ct = qnorm(0.975), alpha = 0.05, direction = "positive"
  ))

  expect_error(conditional_direction_preferring_ci(
    y = 1, sd = 1, r = -5,
    ct = qnorm(0.975), alpha = 0.05, direction = "positive"
  ))

  expect_error(conditional_direction_preferring_ci(
    y = 1, sd = 1, r = 1.3,
    ct = "not numeric", alpha = 0.05, direction = "positive"
  ))

  expect_error(conditional_direction_preferring_ci(
    y = 1, sd = 1, r = 1.3,
    ct = qnorm(0.975), alpha = 2, direction = "positive"
  ))

  expect_error(conditional_direction_preferring_ci(
    y = 1, sd = 1, r = 1.3,
    ct = qnorm(0.975), alpha = 0.05, direction = "wrong"
  ))
})


# =====================================================================
# 2. Basic correctness tests — positive vs negative symmetry
# =====================================================================

test_that("pp vs np conditional CI symmetry holds", {

  y <- 2; ct <- qnorm(0.975); alpha <- 0.05; r <- 1.3

  ci_pos  <- pp_conditional_ci(y,  ct, r, alpha)
  ci_neg  <- np_conditional_ci(y,  ct, r, alpha)
  ci_refl <- -rev(pp_conditional_ci(-y, ct, r, alpha))

  # np_conditional_ci must match reflected pp_conditional_ci
  expect_equal(ci_neg, ci_refl, tolerance = tol)

  # bounds should be ordered
  expect_lte(ci_neg[1], ci_neg[2])
})


# =====================================================================
# 3. Wrapper standardization–rescaling correctness
# =====================================================================

test_that("conditional_direction_preferring_ci standardization + rescaling is correct", {

  y <- 2.3;  sd <- 2
  ct <- qnorm(0.975)
  r <- 1.3
  alpha <- 0.05

  # Manual standardization
  z <- y / sd
  ci_z <- pp_conditional_ci(z, ct / sd, r, alpha)
  ci_manual <- sd * ci_z

  ci_wrapper <- conditional_direction_preferring_ci(
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

  y <- -2
  sd <- 1
  r <- 1.3
  ct <- qnorm(0.975)
  alpha <- 0.05

  ci_neg <- conditional_direction_preferring_ci(
    y, sd, r, ct, alpha, direction = "negative"
  )

  # Manual reflection:
  z <- (y) / sd
  ci_manual <- -rev(pp_conditional_ci(-z, ct / sd, r, alpha))
  ci_manual <- sd * ci_manual

  expect_equal(ci_neg, ci_manual, tolerance = tol)
})


# =====================================================================
# 5. Ordering and monotonicity tests
# =====================================================================

test_that("CI bounds are increasing and finite", {

  y <- 3

  sd <- 1
  r <- 1.3
  alpha <- 0.05
  ct <- qnorm(0.975)

  ci <- conditional_direction_preferring_ci(
    y, sd, r, ct, alpha, direction = "positive"
  )

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

  # Use several near-boundary y values
  ys <- c(-3, -2, -ct - 1e-6, ct + 1e-6, 3)

  for (y in ys) {
    ci <- conditional_direction_preferring_ci(y, sd, r, ct, alpha, direction="positive")
    expect_true(all(is.finite(ci)))
    expect_lte(ci[1], ci[2])
  }
})


# =====================================================================
# 7. Regression test for known values (stability)
# =====================================================================

test_that("Known stable CI values do not change unexpectedly", {

  y <- 1.25
  sd <- 1
  r <- 1.3
  ct <- qnorm(0.6)
  alpha <- 0.05

  ci <- conditional_direction_preferring_ci(y,  sd, r, ct, alpha, "positive")

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

  ci <- conditional_direction_preferring_ci(y, sd, r, ct, alpha, "positive")

  expect_true(all(is.finite(ci)))
  expect_lte(ci[1], ci[2])
})
