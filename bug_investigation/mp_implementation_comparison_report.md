# Modified Pratt CI Implementation Comparison Report

## Overview

This report documents the differences between the **standalone (gold standard)** Modified Pratt CI implementations and the **DirectedCI package** implementations. The standalone implementations are based directly on the paper "Selection Adjusted Confidence Intervals with More Power to Determine the Sign" by Weinstein, Fithian, and Benjamini (JASA 2012).

## Summary of Findings

| Component | Status | Severity |
|-----------|--------|----------|
| Conditional MP CI | **Bug** | High - CI doesn't contain point estimate near threshold |
| Marginal MP CI | **Bug** | High - Lower bound not truncated to 0 when appropriate |

---

## 1. Conditional Modified Pratt CI

### Files Compared
- **Standalone (gold standard):** `standalone_modified_pratt_ci.R` → `MP.CI()`
- **Package:** `R/mp_conditional.R` → `mp_conditional_ci()`

### Observed Differences

#### Test Case: y = 4.0 (large positive observation)
| Implementation | Lower | Upper |
|----------------|-------|-------|
| Standalone     | 0.1547 | 5.6466 |
| Package        | -0.0000 | 5.6466 |

#### Test Case: y = -4.0 (large negative observation)
| Implementation | Lower | Upper |
|----------------|-------|-------|
| Standalone     | -5.6466 | -0.1547 |
| Package        | -5.6466 | 0.0000 |

### Root Cause Analysis

#### 1. Different Regime Thresholds

**Standalone Implementation** uses these key thresholds:
```r
# In MP.CI():
thetatilde1  # Critical theta value from inflated AR equation
xtilde1      # Lower critical x value (CT < x < xtilde1)
xtilde2      # Upper critical x value (xtilde1 < x < xtilde2)
zalphahalf   # qnorm(1 - 0.5 * alpha * Q(0))
ltilde1      # Shortest AR length at thetatilde1
```

**Package Implementation** uses different thresholds:
```r
# In mp_conditional_ci():
theta_1_minus_r   # -find_tilde_1(ct, r, alpha)
neg_threshold_1   # min(ar_func(theta_1_minus_r))
neg_threshold_2   # min(ar_func(0))
neg_threshold_3   # min(ar_func(EPS))
```

The package's threshold calculation doesn't match the mathematical construction in the paper.

#### 2. Upper Bound Computation Differs

**Standalone** (lines 195-200): Always computes upper bound via root finding:
```r
f <- function(theta) {
  pnorm(x + r * Shortest.AR(theta, ct, alpha)$l - theta) -
    pnorm(x - theta) - (1 - alpha) * Q(theta)
}
m <- optimize(f, c(x, x + r * 2 * qnorm(1 - alpha/2)), maximum = TRUE)$maximum
upper <- uniroot(f, c(0, m))$root
```

**Package** (lines 131-135): Sets upper bound to 0 in certain regimes:
```r
} else if (y >= neg_threshold_1 && y < neg_threshold_2) {
  lower_bound <- uniroot(lb_function, lower_interval)$root
  upper_bound <- 0  # <-- Hard-coded to 0, should use root finding
}
```

#### 3. Regime Classification Error

The package has **5 regimes** for negative y:
1. `y < neg_threshold_1`
2. `neg_threshold_1 <= y < neg_threshold_2`
3. `neg_threshold_2 <= y < neg_threshold_3`
4. `y >= neg_threshold_3`

The standalone has **5 regimes** for positive x (then applies symmetry):
1. `ct < x < xtilde1`
2. `xtilde1 <= x < zalphahalf`
3. `zalphahalf <= x < xtilde2`
4. `xtilde2 <= x < ct + r * ltilde1`
5. `x >= ct + r * ltilde1`

The package's regime boundaries don't correctly implement the paper's construction, particularly for large |y| values.

### How to Fix

1. **Replace threshold computation** in `mp_conditional_ci()`:
   - Compute `thetatilde1` using the equation for inflated AR
   - Compute `xtilde1`, `xtilde2` using the paper's equations
   - Use `Shortest.AR()` function properly

2. **Fix upper bound calculation**:
   - For the regime `x >= ct + r * ltilde1`, use root finding for both bounds
   - Never hard-code upper bound to 0 for large |y|

3. **Match the regime structure** from the standalone implementation

---

## 2. Marginal Modified Pratt CI

### Files Compared
- **Standalone (gold standard):** `standalone_non_conditional_modified_pratt_ci.R` → `ci_sym_nc()`
- **Package:** `R/mp_marginal.R` → `mp_marginal_ci()`

### Observed Differences

#### Test Cases
| y | Standalone Lower | Standalone Upper | Package Lower | Package Upper |
|---|------------------|------------------|---------------|---------------|
| 2.0 | **0.0000** | 3.6476 | -1.4483 | 3.6476 |
| 3.0 | **0.0000** | 4.6476 | -0.4483 | 4.6476 |
| -2.0 | -3.6476 | **0.0000** | -3.6476 | 1.4483 |
| -3.0 | -4.6476 | **0.0000** | -4.6476 | 0.4483 |

### Root Cause Analysis

#### 1. Missing Truncation at Zero

**Standalone** (lines 57-68) has 4 regimes with truncation:
```r
if (y >= 0 && y < z_short) {
  ci <- c(y - z_short, y + z_short)  # Full interval
}
if (y >= z_short && y < qnorm(1 - alpha / 2)) {
  ci <- c(0, y + z_short)  # Lower bound TRUNCATED to 0 (closed)
}
if (y >= qnorm(1 - alpha / 2) && y < z_long) {
  ci <- c(0 + epsilon, y + z_short)  # Lower bound TRUNCATED to 0 (open)
}
if (y >= z_long) {
  ci <- c(y - z_long, y + z_short)  # Full interval
}
```

**Package** (lines 93-103) doesn't truncate properly:
```r
if (y < neg_threshold_1) {
  lower <- y - qnorm(1 - (alpha - beta_star))
  upper <- y - qnorm(beta_star)
} else if (y >= neg_threshold_1 && y < neg_threshold_2) {
  lower <- y - qnorm(1 - (alpha - beta_star))
  upper <- 0  # open interval
} else {
  # y >= neg_threshold_2 && y <= 0
  lower <- y + qnorm(alpha - beta_star)  # <-- Should be truncated to 0
  upper <- y - qnorm(alpha - beta_star)
}
```

#### 2. Key Thresholds Differ

**Standalone** uses:
- `z_short = qnorm(1 - beta)` where beta is found from length equation
- `z_long = qnorm(1 - (alpha - beta))`

**Package** uses:
- `neg_threshold_1 = -qnorm(1 - alpha/2)`
- `neg_threshold_2 = -qnorm(1 - (alpha - beta_star))`

These are different quantities and lead to incorrect regime classification.

#### 3. The Core Issue

The Modified Pratt marginal CI is designed to **truncate at 0** for certain y values to give more power to determine the sign. Specifically:
- For `z_short <= y < z_long`, the lower bound should be 0 (not `y - z_long`)
- For `-z_long < y <= -z_short`, the upper bound should be 0 (not `y + z_long`)

The package implementation misses this truncation entirely for the regime where |y| is between `z_short` and `z_long`.

### How to Fix

1. **Compute correct thresholds**:
   ```r
   # Find beta from length equation
   len <- function(beta) qnorm(1 - (alpha - beta)) - qnorm(beta)
   find_beta <- function(beta) len(beta) - r * 2 * qnorm(1 - alpha / 2)
   beta <- uniroot(find_beta, c(alpha/2, alpha))$root

   z_short <- qnorm(1 - beta)
   z_long <- qnorm(1 - (alpha - beta))
   ```

2. **Implement 4 regimes for positive y**:
   ```r
   if (y >= 0 && y < z_short) {
     ci <- c(y - z_short, y + z_short)
   } else if (y >= z_short && y < qnorm(1 - alpha/2)) {
     ci <- c(0, y + z_short)  # Truncate lower to 0 (closed)
   } else if (y >= qnorm(1 - alpha/2) && y < z_long) {
     ci <- c(0 + epsilon, y + z_short)  # Truncate lower to 0 (open)
   } else {  # y >= z_long
     ci <- c(y - z_long, y + z_short)
   }
   ```

3. **Apply symmetry for negative y** (as the standalone does)

---

## 3. Known Bug: Conditional CI Near Threshold

As documented in `mp_conditional_ci_bug.R`, there is an additional bug where `mp_conditional_ci()` returns CIs that don't contain the point estimate for y values in the range ~[1.96, 2.02] and ~[-2.02, -1.96].

This is related to the regime classification issues described above but manifests more severely near the selection threshold.

---

## Recommendations

### Priority 1: Fix Marginal MP CI
The marginal implementation has a fundamental logic error - it doesn't implement the truncation-at-zero behavior that defines the Modified Pratt method. Replace `mp_marginal_ci()` with logic matching `ci_sym_nc()`.

### Priority 2: Fix Conditional MP CI
The conditional implementation uses incorrect regime boundaries and sometimes hard-codes bounds to 0. Rewrite to match `MP.CI()` from the standalone.

### Priority 3: Add Unit Tests
Create tests comparing package output to standalone output for a grid of y values, ensuring:
- CI always contains y (for both marginal and conditional)
- Bounds match standalone to within numerical tolerance
- Truncation at 0 occurs when appropriate

---

## Appendix: Verification Code

```r
# Compare implementations
library(DirectedCI)
source("standalone_modified_pratt_ci.R")
source("standalone_non_conditional_modified_pratt_ci.R")

# Test conditional
y_vals <- c(2.0, 2.5, 3.0, 4.0, -2.0, -2.5, -3.0, -4.0)
for (y in y_vals) {
  standalone <- MP.CI(y, r = 1.3, ct = qnorm(0.975), alpha = 0.05)
  package <- mp_conditional_ci(y, r = 1.3, ct = qnorm(0.975), alpha = 0.05)
  cat(sprintf("y=%4.1f: Standalone [%7.4f, %7.4f], Package [%7.4f, %7.4f]\n",
              y, standalone[1], standalone[2], package[1], package[2]))
}

# Test marginal
y_vals <- c(0, 1.0, 2.0, 3.0, -1.0, -2.0, -3.0)
for (y in y_vals) {
  standalone <- ci_sym_nc(y, r = 1.3, alpha = 0.05)
  package <- mp_marginal_ci(y, r = 1.3, alpha = 0.05)
  cat(sprintf("y=%4.1f: Standalone [%7.4f, %7.4f], Package [%7.4f, %7.4f]\n",
              y, standalone[1], standalone[2], package[1], package[2]))
}
```

---

**Report Generated:** 2025-01-24
**DirectedCI Package Version:** Current development branch
