# DirectedCI

Direction-Preferring Confidence Intervals for Normal Means

## Overview

`DirectedCI` implements direction-preferring confidence intervals for location parameters under normal models. The package provides both **marginal (unconditional)** and **conditional** confidence intervals using three methods:

- **Direction Preferring (DP)**: Asymmetric CIs that allocate more coverage probability to the preferred direction
- **Modified Pratt (MP)**: Alternative asymmetric CIs based on Pratt's approach
- **Shortest**: Standard symmetric CIs (for comparison)

These methods are based on acceptance region inversion and allow preferential inference in a specified direction while maintaining nominal coverage.

## Installation

The package is currently available via GitHub. Install using `devtools` or `remotes`:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install DirectedCI from GitHub
devtools::install_github("tfrostig/DirectedCI")
```

Or using `remotes`:

```r
install.packages("remotes")
remotes::install_github("tfrostig/DirectedCI")
```

## Quick Start

```r
library(DirectedCI)

# Example: Observed estimate with standard error
y <- 3.0    # Point estimate
sd <- 1.0   # Standard error
alpha <- 0.05

# ----- Marginal (Unconditional) CIs -----

# Shortest (symmetric) CI
shortest_marginal_ci_wrapper(y = y, sd = sd, alpha = alpha)

# Direction-Preferring CI (prefer positive direction)
direction_preferring_marginal_ci(
  y = y, sd = sd, r = 1.5, alpha = alpha, direction = "positive"
)

# Modified Pratt CI
modified_pratt_marginal_ci(
  y = y, sd = sd, r = 1.3, alpha = alpha, direction = "positive"
)

# ----- Conditional CIs (requires |y/sd| > ct) -----

ct <- qnorm(1 - alpha/2)  # Critical value (~1.96)

# Only valid when observation is significant: |y/sd| > ct
# In this example: |3.0/1.0| = 3.0 > 1.96 ✓

# Shortest Conditional CI
shortest_conditional_ci_wrapper(y = y, sd = sd, ct = ct, alpha = alpha)

# Direction-Preferring Conditional CI
direction_preferring_conditional_ci(
  y = y, sd = sd, ct = ct, r = 1.3, alpha = alpha, direction = "positive"
)

# Modified Pratt Conditional CI
modified_pratt_conditional_ci(
  y = y, sd = sd, ct = ct, r = 1.3, alpha = alpha, direction = "positive"
)
```

## Main Functions

### User-Friendly Wrappers

These functions handle standardization and return CIs on the original scale:

| Function | Description |
|----------|-------------|
| `direction_preferring_marginal_ci()` | Marginal DP CI |
| `direction_preferring_conditional_ci()` | Conditional DP CI |
| `modified_pratt_marginal_ci()` | Marginal MP CI |
| `modified_pratt_conditional_ci()` | Conditional MP CI |
| `shortest_marginal_ci_wrapper()` | Marginal shortest CI |
| `shortest_conditional_ci_wrapper()` | Conditional shortest CI |

### Internal Z-Scale Functions

For advanced users working on the standardized scale (Y ~ N(θ, 1)):

| Function | Description |
|----------|-------------|
| `dp_marginal_ci()` / `dn_marginal_ci()` | DP marginal CI (positive/negative) |
| `dp_conditional_ci()` / `dn_conditional_ci()` | DP conditional CI (positive/negative) |
| `mp_marginal_ci()` | MP marginal CI |
| `mp_conditional_ci()` | MP conditional CI |
| `shortest_marginal_ci()` | Shortest marginal CI |
| `shortest_conditional_ci()` | Shortest conditional CI |

### Plotting Functions

| Function | Description |
|----------|-------------|
| `plot_ci_comparison()` | Compare multiple CI methods in a forest plot |
| `plot_conditional_ar()` | Visualize conditional acceptance regions |

## Parameters

- **y**: Observed estimate
- **sd**: Standard error of the estimate
- **alpha**: Significance level (default: 0.05)
- **r**: Inflation factor for asymmetry (r >= 1; larger = more asymmetry)
- **ct**: Critical value for conditional CIs (typically `qnorm(1 - alpha/2)`)
- **direction**: Preferred direction (`"positive"` or `"negative"`)

## Return Value

All CI functions return a numeric vector of length 4: `c(ll, ul, lr, ur)`
- Positions 3 and 4 (`lr`, `ur`) contain the lower and upper CI bounds
- Positions 1 and 2 are used for split intervals (NA for single intervals)

```r
ci <- direction_preferring_marginal_ci(y = 3, sd = 1, r = 1.5, direction = "positive")
lower_bound <- ci[3]
upper_bound <- ci[4]
```

## Marginal vs Conditional CIs

| Aspect | Marginal | Conditional |
|--------|----------|-------------|
| **Use case** | General inference | Post-selection inference |
| **Requirement** | None | Observation must be significant: \|y/sd\| > ct |
| **Selection bias** | Not corrected | Corrected |
| **Multiple testing** | Use with BY adjustment for FCR control | Valid without adjustment |

## Example: Meta-Analysis with Selection

```r
library(DirectedCI)
library(metafor)

# Load example data
data(dat.konstantopoulos2011)
dat <- dat.konstantopoulos2011

# Select significant studies
z_scores <- dat$yi / sqrt(dat$vi)
alpha <- 0.05
ct <- qnorm(1 - alpha/2)
significant <- which(abs(z_scores) > ct)

# Compute conditional CIs for significant studies
for (i in significant) {
  ci <- direction_preferring_conditional_ci(
    y = dat$yi[i],
    sd = sqrt(dat$vi[i]),
    ct = ct,
    r = 1.3,
    alpha = alpha,
    direction = "positive"
  )
  cat(sprintf("Study %d: [%.3f, %.3f]\n", i, ci[3], ci[4]))
}
```

## Vignette

For a detailed example using the Konstantopoulos (2011) dataset, see the package vignette:

```r
vignette("example", package = "DirectedCI")
```

## License

MIT License

## Issues and Contributions

- Report bugs: [GitHub Issues](https://github.com/tfrostig/DirectedCI/issues)
- Repository: [https://github.com/tfrostig/DirectedCI](https://github.com/tfrostig/DirectedCI)
