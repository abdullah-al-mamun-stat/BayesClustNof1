#' Bayesian home-host clustering model for series of N-of-1 trials
#'
#' Fits the home-host clustering model (indicators \code{zeta_mu},
#' \code{zeta_beta}; prior \eqn{p(\zeta) \propto r^{N-K}},
#' \eqn{r = p_c/(1-p_c)}). Home values have a
#' proper base distribution \eqn{N(m, \tau^2)} (m flat, \eqn{\tau} half-Cauchy),
#' so the hierarchical model is the limit \eqn{p_c \to 0}. The merge, switch
#' and return-home moves use acceptance probabilities derived from detailed
#' balance, and home values are drawn exactly. Design effects, error variance and AR(1) as in
#' \code{\link{hierarchical_sampler}}.
#'
#' @inheritParams hierarchical_sampler
#' @param clust_prob Clustering probability \eqn{p_c} (strictly between 0 and 1).
#' @param A_tau_mu,A_tau_beta Half-Cauchy scales for the base standard deviations.
#' @param tau_sq_mu_fixed,tau_sq_beta_fixed \code{NULL} = estimate the base
#'   variance (default, as in the paper); a number fixes it.
#' @param q_sd_mu,q_sd_beta Standard deviation of the return-home proposal.
#' @return As \code{\link{hierarchical_sampler}}, plus \code{coclustering_mu},
#'   \code{coclustering_beta}, \code{n_clusters_summary} and
#'   \code{clustering_summary} (posterior probability that each subject is a home).
#' @export
clustering_sampler <- function(iter, ID, Y, A, t, B, P, clust_prob = 0.3,
                               p_zeta_alpha = 0.1, p_zeta_gamma = 0.1, p_zeta_delta = 0.1,
                               sigma_sq_alpha = 1, sigma_sq_gamma = 1, sigma_sq_delta = 1,
                               A_tau_mu = 5, A_tau_beta = 1,
                               tau_sq_mu_fixed = NULL, tau_sq_beta_fixed = NULL,
                               q_sd_mu = 1, q_sd_beta = 1, ar1 = FALSE, rho_jump = 0.05) {
  results <- clustering_sampler_cpp_corrected(
    iter, ID, Y, A, t, B, P, clust_prob, p_zeta_alpha, p_zeta_gamma, p_zeta_delta,
    sigma_sq_alpha, sigma_sq_gamma, sigma_sq_delta, A_tau_mu, A_tau_beta,
    if (is.null(tau_sq_mu_fixed)) -1 else tau_sq_mu_fixed,
    if (is.null(tau_sq_beta_fixed)) -1 else tau_sq_beta_fixed,
    q_sd_mu, q_sd_beta, ar1, rho_jump)
  out <- summarise_common(results, iter, ids = sort(unique(ID)),
                          extra_other = list(m_mu = results$m_mu, tau_mu = results$tau_mu,
                                             m_beta = results$m_beta, tau_beta = results$tau_beta,
                                             rho = results$rho))
  out <- add_cluster_summaries(out, results, "zeta_mu", "zeta_beta", home_labels = TRUE)
  out$hyperparameters <- list(clust_prob = clust_prob, p_zeta_alpha = p_zeta_alpha,
                              p_zeta_gamma = p_zeta_gamma, p_zeta_delta = p_zeta_delta,
                              sigma_sq_alpha = sigma_sq_alpha, sigma_sq_gamma = sigma_sq_gamma,
                              sigma_sq_delta = sigma_sq_delta, A_tau_mu = A_tau_mu, A_tau_beta = A_tau_beta,
                              tau_sq_mu_fixed = tau_sq_mu_fixed, tau_sq_beta_fixed = tau_sq_beta_fixed, ar1 = ar1)
  out$raw <- results
  out
}
