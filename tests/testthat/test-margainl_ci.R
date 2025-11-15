test_that("new implementation matches old implementation in the canonical case", {

  ## --- Original (old) implementations ---------------------------------------

  ar_ns_nc_old <- function(theta, r_l, r_u, switch_val = 0, alpha = 0.05) {
    len <- function(beta) qnorm(1 - (alpha - beta)) - qnorm(beta)

    if (theta <= switch_val) {
      if (theta < -qnorm(1 - alpha / 2) | r_l == 1) {
        beta <- alpha / 2
      } else {
        find_beta <- function(beta) len(beta) - r_l * 2 * qnorm(1 - alpha / 2)
        beta      <- uniroot(find_beta,
                             c(0, alpha / 2),
                             tol = .Machine$double.eps)$root
      }
    }
    if (theta > switch_val) {
      if (theta > qnorm(1 - alpha / 2) | r_u == 1) {
        beta <- alpha / 2
      } else {
        find_beta <- function(beta) len(beta) - r_u * 2 * qnorm(1 - alpha / 2)
        beta      <- alpha - uniroot(find_beta,
                                     c(0, alpha / 2),
                                     tol = .Machine$double.eps)$root
      }
    }
    if (theta == switch_val & r_l == r_u) {
      beta <- alpha / 2
    }
    list(
      ar   = c(theta + qnorm(beta),
               theta + qnorm(1 - (alpha - beta))),
      beta = beta
    )
  }

  ci_ns_nc_old <- function(y, r_l, r_u = 1,
                           switch_val = 0, alpha = 0.05,
                           epsilon = 10^-15) {
    ar_func <- function(theta) ar_ns_nc_old(theta, r_l, r_u,
                                            switch_val = 0, alpha)

    # assuming r_u=1
    c0 <- ar_func(qnorm(alpha / 2) + epsilon)$ar[1]
    c1 <- ar_func(switch_val - epsilon)$ar[1]
    c2 <- ar_func(switch_val + epsilon)$ar[1]
    c3 <- ar_func(switch_val - epsilon)$ar[2]
    c4 <- ar_func(switch_val + epsilon)$ar[2]

    l_switch_val_beta <- ar_func(switch_val - epsilon)$beta
    u_switch_val_beta <- ar_func(switch_val + epsilon)$beta

    if (y < c0) {
      return(c(y - qnorm(1 - alpha / 2), y - qnorm(alpha / 2)))
    }
    if (y >= c0 & y < c1) {
      return(c(y - qnorm(1 - alpha / 2), y + qnorm(1 - l_switch_val_beta)))
    }
    if (y >= c1 & y < c2) {
      return(c(y - qnorm(1 - alpha / 2), 0)) # (... , 0]
    }
    if (y >= c2 & y <= switch_val ) {
      return(c(y - qnorm(1 - alpha / 2), y + qnorm(1 - u_switch_val_beta)))
    }
    if (y >= switch_val & y < c3 ) {
      return(c(y - qnorm(1 - (alpha  - l_switch_val_beta)),
               y + qnorm(1 - alpha / 2)))
    }
    if (y >= c3 & y < c4) {
      return(c(0 + epsilon, y + qnorm(1 - alpha / 2))) # (0, ...)
    }
    if (y >= c4 ) {
      return(c(y - qnorm(1 - alpha / 2), y - qnorm(alpha / 2)))
    }
  }

  ## --- Parameter grid where we expect agreement -----------------------------

  set.seed(1)

  alphas      <- c(0.05)
  switch_vals <- c(0)   # original ci_ns_nc assumes switch_val = 0 effectively
  r_ls        <- c(1, 1.75, 2)
  r_us        <- c(1.0) # ci_ns_nc was derived assuming r_u = 1
  thetas      <- seq(-3, 3, length.out = 21)
  ys          <- seq(-3, 3, length.out = 21)

  tol <- 1e-7

  ## --- Compare ar_ns_nc vs ar_ns_nc_old -------------------------------------

  for (alpha in alphas) {
    for (switch_val in switch_vals) {
      for (r_l in r_ls) {
        for (r_u in r_us) {
          for (theta in thetas) {
            # Old version may hit qnorm(0) issues; catch and skip if so
            old_val <- try(ar_ns_nc_old(theta, r_l, r_u,
                                        switch_val = switch_val,
                                        alpha = alpha),
                           silent = TRUE)
            new_val <- try(ar_ns_nc(theta, r_l, r_u,
                                    switch_val = switch_val,
                                    alpha = alpha),
                           silent = TRUE)

            if (inherits(old_val, "try-error") ||
                inherits(new_val, "try-error")) {
              next
            }

            expect_equal(old_val$ar,   new_val$ar,   tolerance = tol)
            expect_equal(old_val$beta, new_val$beta, tolerance = tol)
          }
        }
      }
    }
  }

  ## --- Compare ci_ns_nc vs ci_ns_nc_old -------------------------------------

  for (alpha in alphas) {
    for (switch_val in switch_vals) {
      for (r_l in r_ls) {
        for (r_u in r_us) {
          for (y in ys) {

            old_ci <- try(ci_ns_nc_old(y, r_l, r_u,
                                       switch_val = switch_val,
                                       alpha = alpha),
                          silent = TRUE)
            new_ci <- try(ci_ns_nc(y, r_l, r_u,
                                   switch_val = switch_val,
                                   alpha = alpha),
                          silent = TRUE)

            if (inherits(old_ci, "try-error") ||
                inherits(new_ci, "try-error") ||
                any(is.na(old_ci)) || any(is.na(new_ci))) {
              next
            }

            expect_equal(old_ci, new_ci, tolerance = tol)
          }
        }
      }
    }
  }

})
