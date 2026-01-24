#' ============================================
#' MP Conditional CI Bug Investigation
#' ============================================
#'
#' Issue: mp_conditional_ci() returns CIs that do not contain the point
#' estimate for z-values near the selection threshold (|z| ~ 1.96-2.02).
#'
#' The equivalent DP conditional CI works correctly for the same values.
#'
#' Date: 2025-01-24
#' ============================================

library(DirectedCI)

# Parameters
alpha <- 0.05
ct <- qnorm(1 - alpha / 2)  # 1.96
r <- 1.3

cat("============================================\n")
cat("MP Conditional CI Bug - Reproducible Example\n")
cat("============================================\n\n")

cat("Parameters:\n")
cat("  alpha =", alpha, "\n")
cat("  ct =", round(ct, 4), "\n")
cat("  r =", r, "\n\n")

# ============================================
# CASE 1: Positive z just above threshold
# ============================================
cat("=== CASE 1: z = 1.97 (just above ct = 1.96) ===\n\n")

z1 <- 1.97

# MP Conditional CI - BUG
mp_cond1 <- mp_conditional_ci(y = z1, ct = ct, r = r, alpha = alpha)
cat("MP Conditional CI:\n")
cat("  z =", z1, "\n")
cat("  CI = [", round(mp_cond1[1], 4), ",", round(mp_cond1[2], 4), "]\n")
cat("  Contains z? ", ifelse(z1 >= mp_cond1[1] & z1 <= mp_cond1[2], "YES", "NO"), "\n")
cat("  >>> BUG: Upper bound (", round(mp_cond1[2], 4), ") < z (", z1, ") <<<\n\n")

# DP Conditional CI - WORKS
dp_cond1 <- dp_conditional_ci(y = z1, ct = ct, r = r, alpha = alpha)
cat("DP+ Conditional CI (for comparison):\n")
cat("  z =", z1, "\n")
cat("  CI = [", round(dp_cond1[1], 4), ",", round(dp_cond1[2], 4), "]\n")
cat("  Contains z? ", ifelse(z1 >= dp_cond1[1] & z1 <= dp_cond1[2], "YES", "NO"), "\n\n")

# ============================================
# CASE 2: Negative z at threshold
# ============================================
cat("=== CASE 2: z = -2.00 (at negative threshold) ===\n\n")

z2 <- -2.00

# MP Conditional CI - BUG
mp_cond2 <- mp_conditional_ci(y = z2, ct = ct, r = r, alpha = alpha)
cat("MP Conditional CI:\n")
cat("  z =", z2, "\n")
cat("  CI = [", round(mp_cond2[1], 4), ",", round(mp_cond2[2], 4), "]\n")
cat("  Contains z? ", ifelse(z2 >= mp_cond2[1] & z2 <= mp_cond2[2], "YES", "NO"), "\n")
cat("  >>> BUG: Lower bound (", round(mp_cond2[1], 4), ") > z (", z2, ") <<<\n\n")

# DP- Conditional CI - WORKS
dp_cond2 <- dn_conditional_ci(y = z2, ct = ct, r = r, alpha = alpha)
cat("DP- Conditional CI (for comparison):\n")
cat("  z =", z2, "\n")
cat("  CI = [", round(dp_cond2[1], 4), ",", round(dp_cond2[2], 4), "]\n")
cat("  Contains z? ", ifelse(z2 >= dp_cond2[1] & z2 <= dp_cond2[2], "YES", "NO"), "\n\n")

# ============================================
# TRANSITION ZONE ANALYSIS
# ============================================
cat("=== TRANSITION ZONE: z = 1.96 to 2.10 ===\n\n")

cat("Comparing MP vs DP+ Conditional CIs near threshold:\n\n")
cat("z      | MP Cond CI            | MP OK? | DP+ Cond CI           | DP OK?\n")
cat("-------|----------------------|--------|----------------------|-------\n")

z_values <- seq(1.96, 2.10, by = 0.02)

for (z in z_values) {
  # MP
 mp_ci <- tryCatch(
    mp_conditional_ci(y = z, ct = ct, r = r, alpha = alpha),
    error = function(e) c(NA, NA)
  )
  mp_ok <- !is.na(mp_ci[1]) && z >= mp_ci[1] && z <= mp_ci[2]

  # DP
  dp_ci <- tryCatch(
    dp_conditional_ci(y = z, ct = ct, r = r, alpha = alpha),
    error = function(e) c(NA, NA)
  )
  dp_ok <- !is.na(dp_ci[1]) && z >= dp_ci[1] && z <= dp_ci[2]

  cat(sprintf("%.2f   | [%7.3f, %7.3f]   | %s    | [%7.3f, %7.3f]   | %s\n",
              z,
              mp_ci[1], mp_ci[2], ifelse(mp_ok, "YES", "NO "),
              dp_ci[1], dp_ci[2], ifelse(dp_ok, "YES", "NO ")))
}

cat("\n")
cat("=== NEGATIVE TRANSITION ZONE: z = -2.10 to -1.96 ===\n\n")

cat("z      | MP Cond CI            | MP OK? | DP- Cond CI           | DP OK?\n")
cat("-------|----------------------|--------|----------------------|-------\n")

z_values <- seq(-2.10, -1.96, by = 0.02)

for (z in z_values) {
  # MP
  mp_ci <- tryCatch(
    mp_conditional_ci(y = z, ct = ct, r = r, alpha = alpha),
    error = function(e) c(NA, NA)
  )
  mp_ok <- !is.na(mp_ci[1]) && z >= mp_ci[1] && z <= mp_ci[2]

  # DP-
  dp_ci <- tryCatch(
    dn_conditional_ci(y = z, ct = ct, r = r, alpha = alpha),
    error = function(e) c(NA, NA)
  )
  dp_ok <- !is.na(dp_ci[1]) && z >= dp_ci[1] && z <= dp_ci[2]

  cat(sprintf("%.2f  | [%7.3f, %7.3f]   | %s    | [%7.3f, %7.3f]   | %s\n",
              z,
              mp_ci[1], mp_ci[2], ifelse(mp_ok, "YES", "NO "),
              dp_ci[1], dp_ci[2], ifelse(dp_ok, "YES", "NO ")))
}

# ============================================
# SUMMARY
# ============================================
cat("\n")
cat("============================================\n")
cat("SUMMARY\n")
cat("============================================\n\n")

cat("BUG: mp_conditional_ci() returns CIs that do not contain the\n")
cat("     point estimate for z-values in the range ~[1.96, 2.02]\n")
cat("     and ~[-2.02, -1.96].\n\n")

cat("The issue appears to be in the calculation of the CI bounds\n")
cat("near the selection threshold, where the acceptance region\n")
cat("transitions.\n\n")

cat("AFFECTED: mp_conditional_ci()\n")
cat("NOT AFFECTED: dp_conditional_ci(), dn_conditional_ci()\n")
cat("NOT AFFECTED: All marginal CI functions\n")
