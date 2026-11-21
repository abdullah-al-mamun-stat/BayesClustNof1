#' Bayesian hierarchical model for series of N-of-1 trials
#'
#' Fits the hierarchical model, in which the patient means and treatment
#' effects are drawn from normal distributions,
#' \eqn{\mu_i \sim N(\tilde\mu, \sigma^2_\mu)} and
#' \eqn{\beta_i \sim N(\tilde\beta, \sigma^2_\beta)}, with flat priors on the
#' centres and half-Cauchy priors on \eqn{\sigma_\mu} and \eqn{\sigma_\beta}
#' (or fixed variances). Patient effects are drawn exactly by Gibbs sampling,
#' spike-and-slab indicators by collapsed Gibbs sampling, the error variance
#' has prior \eqn{1/\sigma^2}, and AR(1) errors within treatment periods are
#' available (\code{ar1 = TRUE}).
#'
#' @param iter Total number of MCMC iterations; the first half is burn-in.
#' @param ID Subject IDs (one per observation). Data must be sorted by subject and time.
#' @param Y Outcome.
#' @param A Treatment indicator (lowest level = reference).
#' @param t Time within period.
#' @param B Block indicator (lowest level = reference).
#' @param P Period position within a block, e.g. 1 and 2 (lowest level = reference).
#' @param p_zeta_alpha,p_zeta_gamma,p_zeta_delta Prior inclusion probabilities.
#' @param sigma_sq_alpha,sigma_sq_gamma,sigma_sq_delta Slab variances.
#' @param sigma_sq_mu,sigma_sq_beta \code{NULL} (default) estimates the
#'   between-subject variance with a half-Cauchy prior; a number fixes it
#'   (e.g. 25 reproduces the original comparator).
#' @param A_mu,A_beta Half-Cauchy scales for \eqn{\sigma_\mu} and \eqn{\sigma_\beta}, in outcome units.
#' @param ar1 Logical; AR(1) errors within each treatment period.
#' @param rho_jump Random-walk standard deviation for \code{atanh(rho)}.
#' @return A list with the same components as before (\code{subject_summary},
#'   \code{other_summary}, \code{bar_summary}, \code{alpha_summary},
#'   \code{gamma_summary}, \code{delta_summary}, \code{posterior_samples}, ...),
#'   plus \code{sigma_sq_mu}, \code{sigma_sq_beta} and \code{rho} in
#'   \code{other_summary}, Rao-Blackwellised inclusion probabilities, and the
#'   raw sampler output in \code{raw}.
#' @examples
#' \dontrun{
#' fit <- hierarchical_sampler(iter = 20000, ID = dat$ID, Y = dat$Y, A = dat$A,
#'                             t = dat$t, B = dat$B, P = dat$P)
#' fit$subject_summary
#' }
#' @export
hierarchical_sampler <- function(iter, ID, Y, A, t, B, P,
                                 p_zeta_alpha = 0.1, p_zeta_gamma = 0.1, p_zeta_delta = 0.1,
                                 sigma_sq_alpha = 1, sigma_sq_gamma = 1, sigma_sq_delta = 1,
                                 sigma_sq_mu = NULL, sigma_sq_beta = NULL,
                                 A_mu = 5, A_beta = 1, ar1 = FALSE, rho_jump = 0.05) {
  results <- hierarchical_sampler_cpp_corrected(
    iter, ID, Y, A, t, B, P, p_zeta_alpha, p_zeta_gamma, p_zeta_delta,
    sigma_sq_alpha, sigma_sq_gamma, sigma_sq_delta,
    if (is.null(sigma_sq_mu)) -1 else sigma_sq_mu,         # -1 tells C++ to estimate it
    if (is.null(sigma_sq_beta)) -1 else sigma_sq_beta,
    A_mu, A_beta, ar1, rho_jump)
  out <- summarise_common(results, iter, ids = sort(unique(ID)),
                          extra_other = list(mu_tilde = results$mu_tilde, beta_tilde = results$beta_tilde,
                                             sigma_sq_mu = results$sigma_sq_mu,
                                             sigma_sq_beta = results$sigma_sq_beta, rho = results$rho))
  out$hyperparameters <- list(p_zeta_alpha = p_zeta_alpha, p_zeta_gamma = p_zeta_gamma,
                              p_zeta_delta = p_zeta_delta, sigma_sq_alpha = sigma_sq_alpha,
                              sigma_sq_gamma = sigma_sq_gamma, sigma_sq_delta = sigma_sq_delta,
                              sigma_sq_mu = sigma_sq_mu, sigma_sq_beta = sigma_sq_beta,
                              A_mu = A_mu, A_beta = A_beta, ar1 = ar1)
  out$raw <- results
  out
}
