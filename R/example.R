## example of shortest conditional
epsilon <- 0.001
alpha   <- 0.05
ct      <- 1.96
lb_quantile <- quantile_tn_outer(epsilon, mu = 0, b = ct)
ub_quantile <- quantile_tn_outer(1 - (alpha - epsilon), mu = 0, b = ct)

len_quantile <- ub_quantile - lb_quantile - 2 * ct
len_shortets <- Shortest.AR(0, ct, alpha)$l

r <- len_quantile / len_shortets
pp_conditional_ar(theta = 0, ct = ct, r = r, alpha = 0.05)

##
theta <- -0.62

lb_quantile <- quantile_tn_outer(epsilon, mu = theta, b = ct)
ub_quantile <- quantile_tn_outer(1 - (alpha - epsilon), mu = theta, b = ct)

pp_conditional_ar(theta = theta, ct = ct, r = r, alpha = 0.05)


len_quantile   <- ub_quantile - lb_quantile - 2 * ct
len_r_times_shortest <-  calculate_ar_length(pp_conditional_ar(theta = theta, ct = ct, r = r, alpha = 0.05))

len_quantile / len_r_times_shortest

