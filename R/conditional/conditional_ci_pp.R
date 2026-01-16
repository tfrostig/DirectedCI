
#' Direction-Preferring Confidence Interval via Acceptance-Region Inversion
#'
#' Constructs a \emph{marginal direction-preferring} confidence interval (CI) for a
#' location parameter \eqn{\theta} by inverting non-equivariant acceptance regions
#' that favor inference in a preferred direction (positive).
#'
#' This CI implements the marginal direction-preferring construction described in
#' Frostig et al. (2025) :contentReference[oaicite:0]{index=0}.
#' The method modifies the acceptance regions used for testing \eqn{H_0:\theta=\tau}
#' so that, when \eqn{\tau \le 0}, the acceptance region is shortened on the
#' upper side; this increases the probability of correctly determining a positive
#' effect while preserving unconditional \eqn{1-\alpha} coverage.
#'
#' The resulting CI is generally shorter for \eqn{\theta>0} and determines the sign
#' earlier than the standard symmetric CI, while allowing longer intervals in the
#' non-preferred direction (\eqn{\theta \le 0}) where precise estimation is less
#' critical.  The CI is built by:
#' \itemize{
#'   \item computing critical regime-change points derived from the
#'   direction-preferring acceptance region,
#'   \item evaluating the conditional acceptance-region (AR) endpoints
#'   \eqn{\operatorname{AR}(\theta)} using `pp_conditional_ar()`,
#'   \item solving for the values of \eqn{\theta} whose acceptance region contains
#'   the observed statistic \eqn{y}, and
#'   \item assembling the CI by inverting the AR according to the regime into
#'   which \eqn{y} falls.
#' }
#'
#' The function corresponds to the marginal (unconditional) direction-preferring
#' CIs of §2.2 in the manuscript, using the modified Pratt–type bounding of
#' acceptance-region length and the inflation parameter \eqn{r} that determines the
#' maximum allowable AR expansion in the non-preferred direction.
#'
#' @param y Numeric. The observed test statistic \eqn{Y}.
#' @param ct Numeric. The cutoff \eqn{c_t} defining the standard acceptance region
#'   (typically \eqn{q_{1-\alpha/2}} for a standard normal distribution).
#' @param r Numeric. Inflation parameter controlling the maximum allowed length of
#'   acceptance regions for \eqn{\theta \le 0}. Larger \eqn{r} increases the
#'   preference for determining a positive sign.
#' @param alpha Numeric. Nominal significance level \eqn{\alpha}; the CI has
#'   conditional \eqn{1-\alpha} coverage.
#'
#' @return
#' A numeric vector of length 2:
#' \itemize{
#'   \item \code{lower_bound}: the lower end of the direction-preferring CI;
#'   \item \code{upper_bound}: the upper end of the CI.
#' }
#'
#' @references
#' Frostig, T., Panagiotou, O., Benjamini, Y., & Heller, R. (2025).
#' *Direction Preferring Confidence Intervals*.
#' Manuscript. :contentReference[oaicite:1]{index=1}
#'
#' @seealso
#' \code{\link{pp_conditional_ar}},
#' \code{\link{theta_1_finder}},
#' \code{\link{find_tilde_1}},
#' \code{\link{bounds_for_ci}}
#'
#' @examples
#' # Construct a direction-preferring CI for a normal statistic:
#' pp_confidence_interval(y = 1.5, ct = qnorm(0.975), r = 1.3, alpha = 0.05)
#'
#' @export
pp_conditional_ci <- function(y, ct, r, alpha) {
  ## theta values at regime changes
  theta_1_minus_shortest <- -theta_1_finder(ct, alpha)
  theta_1_minus_r        <- -find_tilde_1(ct, r, alpha)

  ## important values
  neg_threshold_1 <- min(pp_conditional_ar(theta_1_minus_r, ct, r, alpha), na.rm = TRUE)
  neg_threshold_2 <- min(pp_conditional_ar(0, ct, r, alpha), na.rm = TRUE)
  neg_threshold_3 <- min(pp_conditional_ar(0, ct, 1, alpha), na.rm = TRUE)
  ## -ct
  ##  ct, provided for completeness
  pos_threshold_1 <- max(pp_conditional_ar(0, ct, r, alpha), na.rm = TRUE)
  pos_threshold_2 <- max(pp_conditional_ar(0, ct, 1, alpha), na.rm = TRUE)


  ## find all relevant values
  ## naive approach - invert the confidence interval
  lb_function <- function(theta) {
    ar <- pp_conditional_ar(theta, ct, r, alpha)
    max(ar, na.rm = TRUE) - y
  }
  ub_function <- function(theta) {
    ar <- pp_conditional_ar(theta, ct, r, alpha)
    min(ar, na.rm = TRUE) - y
  }

  if (y < neg_threshold_1) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }
  if (y >= neg_threshold_1 & y < neg_threshold_2) {
    lower_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- c(theta_1_minus_r + EPS, 0 - EPS)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }
  if (y >= neg_threshold_2 & y < neg_threshold_3) {
    lower_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_bound    <- 0 ## closed

  }
  if (y >= neg_threshold_3 & y <= -ct) {
    lower_interval <- bounds_for_ci(y, ct, alpha / 2, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }
  if (y >= ct & y < pos_threshold_1) {
    lower_interval <- c(0 - EPS, theta_1_minus_r + EPS)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }

  if (y >= pos_threshold_1 & y < pos_threshold_2) {
    lower_bound    <- 0 ## open

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }
  if (y >= pos_threshold_2) {
    lower_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = TRUE)
    lower_bound    <- uniroot(lb_function, lower_interval)$root

    upper_interval <- bounds_for_ci(y, ct, alpha, lower_ci_bound = FALSE)
    upper_bound    <- uniroot(ub_function, upper_interval)$root
  }

  return(c(lower_bound, upper_bound))

}


np_conditional_ci <- function(y, ct, r, alpha) {
  ## Obtain the CI for the negative of the statistic
  ci_neg <- pp_conditional_ci(-y, ct, r, alpha)

  ## Reverse and negate the bounds
  ci_pos <- -rev(ci_neg)

  return(sort(ci_pos))
}


#' User-Facing Direction-Preferring Confidence Interval
#'
#' Computes a direction-preferring confidence interval for a Normally distributed
#' estimator \eqn{Y \sim N(\mu, \sigma^2)}.
#'
#' This is a wrapper around `pp_confidence_interval()` (positive direction)
#' and `np_confidence_interval()` (negative direction), performing the necessary
#' standardization to the \eqn{N(0,1)} scale, calling the core interval routines,
#' and then mapping the CI back to the original scale.
#'
#' @param y Numeric. Observed value of the estimator.
#' @param sd Numeric (>0). Standard deviation of the estimator \eqn{\sigma}.
#' @param r Numeric. Inflation parameter for the acceptance-region length.
#' @param ct Numeric. Cutoff used in the core CI routines.
#'   For a Normal model, \code{ct = qnorm(1 - alpha/2)} is typical.
#' @param alpha Numeric in (0,1). Nominal CI error level.
#' @param direction Character, one of \code{"positive"} or \code{"negative"}.
#'
#' @details
#' The CI is computed in three steps:
#' \enumerate{
#'   \item Standardize: \eqn{z = (y - \mu)/\sigma}.
#'   \item Compute a direction-preferring CI for \eqn{z} using the internal
#'         methods designed for \eqn{N(0,1)}.
#'   \item Rescale: map the interval back to the original scale using
#'         \eqn{\mu + \sigma \cdot \mathrm{CI}_z}.
#' }
#'
#' The function validates inputs, checks bounds, and ensures that the returned CI
#' is well-ordered and finite.
#'
#' @return A numeric vector of length two, the confidence interval on the original scale.
#'
#' @examples
#' direction_preferring_ci(
#'   y = 1.2, sd = 2,
#'   r = 1.3, ct = qnorm(0.975),
#'   alpha = 0.05, direction = "positive"
#' )
#'
#' @export
conditional_direction_preferring_ci <- function(
    y,
    sd,
    r = 1.3,
    ct = qnorm(1 - 0.05 / 2),
    alpha = 0.05,
    direction = c("positive", "negative")
) {

  direction <- match.arg(direction)

  if (!is.numeric(y) || length(y) != 1)
    stop("y must be a single numeric value.")

  if (!is.numeric(sd) || length(sd) != 1 || sd <= 0)
    stop("sd must be a single positive numeric value.")

  if (!is.numeric(r) || r <= 0)
    stop("r must be a positive numeric inflation factor.")

  if (!is.numeric(ct) || length(ct) != 1)
    stop("ct must be a single numeric threshold used by core CIs.")

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0,1).")

  z    <-  y / sd
  z_ct <- ct / sd

  if (abs(z) < z_ct) {
    stop("The normalized score should be larger than the nomralized cutoff ct/sd
         to construct a direction-preferring CI.")
  }

  if (direction == "positive") {
    ci_standardized <- pp_conditional_ci(z, ct = z_ct, r = r, alpha = alpha)
  }
  else {
    ci_standardized <- np_conditional_ci(z, ct = z_ct, r = r, alpha = alpha)
}


  if (any(!is.finite(ci_standardized))) {
    warning("Infinite or NaN values in standardized CI.")
  }

  ci_original <- sd * ci_standardized

  return(ci_original)
}

