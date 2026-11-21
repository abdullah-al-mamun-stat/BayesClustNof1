// =====================================================================
// CLUSTERING_SAMPLER_CORRECTED.cpp
//
// C++ version of clustering_sampler_corrected() in
// R CODE/CLUSTERING_SAMPLER_CORRECTED.R (your home-host clustering model
// with corrected acceptance ratios). Same seed -> same draws as R.
//
//   Rcpp::sourceCpp("CLUSTERING_SAMPLER_CORRECTED.cpp")   # NOF1_COMMON_CORRECTED.h in the same folder
//   fit <- clustering_sampler_cpp_corrected(iter, ID, Y, A, t, B, P, clust_prob = 0.3)
//
// tau_sq_mu_fixed / tau_sq_beta_fixed: negative (default) = estimate the
// base variance (half-Cauchy); positive = fix it. The corrected MIXED
// sampler is this function with tau_sq_*_fixed = 25 (see mixed_sampler_cpp_corrected).
// zeta_mu / zeta_beta are returned 1-based, as in R.
// =====================================================================
#include "NOF1_COMMON_CORRECTED.h"   // shared helpers (in this src folder)

// [[Rcpp::export]]
List clustering_sampler_cpp_corrected(const int iter,
                                      SEXP ID, SEXP Y, SEXP A, SEXP t, SEXP B, SEXP P,
                                      const double clust_prob     = 0.3,
                                      const double p_zeta_alpha   = 0.1,
                                      const double p_zeta_gamma   = 0.1,
                                      const double p_zeta_delta   = 0.1,
                                      const double sigma_sq_alpha = 1.0,
                                      const double sigma_sq_gamma = 1.0,
                                      const double sigma_sq_delta = 1.0,
                                      const double A_tau_mu   = 5.0,
                                      const double A_tau_beta = 1.0,
                                      const double tau_sq_mu_fixed   = -1.0,   // < 0: estimate
                                      const double tau_sq_beta_fixed = -1.0,   // < 0: estimate
                                      const double q_sd_mu   = 1.0,
                                      const double q_sd_beta = 1.0,
                                      const bool   ar1 = false,
                                      const double rho_jump = 0.05) {

  if (clust_prob <= 0.0 || clust_prob >= 1.0) stop("clust_prob must be strictly between 0 and 1.");
  const double r_odds = clust_prob / (1.0 - clust_prob);

  const Design d = make_design_corrected(ID, Y, A, t, B, P);
  const int N_id = d.N_id, N_tot = d.N_tot, K1 = d.K_minus_1;
  const bool fix_mu = (tau_sq_mu_fixed > 0.0), fix_beta = (tau_sq_beta_fixed > 0.0);

  // ---------------- storage ----------------
  NumericMatrix mu_store(iter, N_id), beta_store(iter, N_id);
  IntegerMatrix zeta_mu_store(iter, N_id), zeta_beta_store(iter, N_id);
  NumericVector alpha_store(iter), prob_alpha_store(iter), delta_store(iter), prob_delta_store(iter);
  IntegerVector zeta_alpha_store(iter), zeta_delta_store(iter), accept_mu_store(iter), accept_beta_store(iter);
  NumericMatrix gamma_store(iter, K1), prob_gamma_store(iter, K1);
  IntegerMatrix zeta_gamma_store(iter, K1);
  NumericVector sigma_sq_store(iter), rho_store(iter), m_mu_store(iter), tau_mu_store(iter),
                m_beta_store(iter), tau_beta_store(iter);

  // ---------------- initial values: everybody is a home ----------------
  dvec mu(N_id), beta(N_id, 0.0);
  for (int i = 0; i < N_id; ++i) {
    long double s = 0.0;
    for (std::size_t h = 0; h < d.obs_by_subj[i].size(); ++h) s += d.Y[d.obs_by_subj[i][h]];
    mu[i] = (double) s / (double) d.obs_by_subj[i].size();
  }
  ivec zeta_mu(N_id), zeta_beta(N_id);
  for (int i = 0; i < N_id; ++i) { zeta_mu[i] = i; zeta_beta[i] = i; }
  Nuisance nu;
  nu.gamma.assign(K1, 0.0); nu.zeta_gamma.assign(K1, 0); nu.prob_gamma.assign(K1, NA_REAL);
  double rho = 0.0;
  dvec res0(N_tot);
  for (int j = 0; j < N_tot; ++j) res0[j] = d.Y[j] - mu[d.ID0[j]];
  double sig2 = sum_prod_ld(res0, res0) / (double) N_tot;
  double m_mu = sum_ld(mu) / (double) N_id;
  double tau_mu = fix_mu ? std::sqrt(tau_sq_mu_fixed) : A_tau_mu, a_mu = 1.0;
  double m_beta = 0.0;
  double tau_beta = fix_beta ? std::sqrt(tau_sq_beta_fixed) : A_tau_beta, a_beta = 1.0;

  dvec fit_mb(N_tot), r(N_tot), Tv, Zv;

  // ---------------- MCMC ----------------
  for (int jj = 0; jj < iter; ++jj) {
    if (jj % 500 == 0) Rcpp::checkUserInterrupt();

    // (1) transformed data
    const Transformed tr = transform_design(d, rho);
    for (int j = 0; j < N_tot; ++j) fit_mb[j] = mu[d.ID0[j]] * tr.X0[j] + beta[d.ID0[j]] * tr.As[j];

    // (2) spike-and-slab design effects
    update_nuisance(nu, tr, fit_mb, sig2, N_tot, K1, p_zeta_alpha, p_zeta_gamma, p_zeta_delta,
                    sigma_sq_alpha, sigma_sq_gamma, sigma_sq_delta);
    const dvec nuis = nuisance_term(nu, tr, N_tot);

    // (3) error variance
    sig2 = update_sigma_sq(tr, fit_mb, nuis, N_tot);

    // (4) clustering of mu
    for (int j = 0; j < N_tot; ++j) r[j] = tr.Ys[j] - beta[d.ID0[j]] * tr.As[j] - nuis[j];
    subject_TZ(tr.X0, r, d.obs_by_subj, Tv, Zv);
    accept_mu_store[jj] = update_clusters_corrected(zeta_mu, mu, Tv, Zv, sig2, m_mu, tau_mu, r_odds, q_sd_mu);
    { dvec hv; for (int i = 0; i < N_id; ++i) if (zeta_mu[i] == i) hv.push_back(mu[i]);
      update_base_corrected(hv, m_mu, tau_mu, a_mu, A_tau_mu, fix_mu); }

    // (5) clustering of beta
    for (int j = 0; j < N_tot; ++j) r[j] = tr.Ys[j] - mu[d.ID0[j]] * tr.X0[j] - nuis[j];
    subject_TZ(tr.As, r, d.obs_by_subj, Tv, Zv);
    accept_beta_store[jj] = update_clusters_corrected(zeta_beta, beta, Tv, Zv, sig2, m_beta, tau_beta, r_odds, q_sd_beta);
    { dvec hv; for (int i = 0; i < N_id; ++i) if (zeta_beta[i] == i) hv.push_back(beta[i]);
      update_base_corrected(hv, m_beta, tau_beta, a_beta, A_tau_beta, fix_beta); }

    // (6) AR(1)
    if (ar1) rho = update_rho(rho, rho_jump, raw_residuals(d, mu, beta, nu), d, sig2);

    // (7) store (zeta 1-based)
    store_row(mu_store, jj, mu); store_row_int(zeta_mu_store, jj, zeta_mu, 1);
    store_row(beta_store, jj, beta); store_row_int(zeta_beta_store, jj, zeta_beta, 1);
    alpha_store[jj] = nu.alpha; zeta_alpha_store[jj] = nu.zeta_alpha; prob_alpha_store[jj] = nu.prob_alpha;
    if (K1 > 0) { store_row(gamma_store, jj, nu.gamma); store_row_int(zeta_gamma_store, jj, nu.zeta_gamma);
                  store_row(prob_gamma_store, jj, nu.prob_gamma); }
    delta_store[jj] = nu.delta; zeta_delta_store[jj] = nu.zeta_delta; prob_delta_store[jj] = nu.prob_delta;
    sigma_sq_store[jj] = sig2; rho_store[jj] = rho;
    m_mu_store[jj] = m_mu; tau_mu_store[jj] = tau_mu; m_beta_store[jj] = m_beta; tau_beta_store[jj] = tau_beta;
  }

  // (List::create takes at most 20 elements, so the list is filled by name)
  List lev = design_levels(d);
  List out;
  out["mu"] = mu_store;             out["zeta_mu"] = zeta_mu_store;
  out["beta"] = beta_store;         out["zeta_beta"] = zeta_beta_store;
  out["alpha"] = alpha_store;       out["zeta_alpha"] = zeta_alpha_store; out["prob_alpha"] = prob_alpha_store;
  out["gamma"] = gamma_store;       out["zeta_gamma"] = zeta_gamma_store; out["prob_gamma"] = prob_gamma_store;
  out["delta"] = delta_store;       out["zeta_delta"] = zeta_delta_store; out["prob_delta"] = prob_delta_store;
  out["sigma_sq"] = sigma_sq_store; out["rho"] = rho_store;
  out["m_mu"] = m_mu_store;         out["tau_mu"] = tau_mu_store;
  out["m_beta"] = m_beta_store;     out["tau_beta"] = tau_beta_store;
  out["accept_mu"] = accept_mu_store; out["accept_beta"] = accept_beta_store;
  out["trt_levels"] = lev["trt_levels"]; out["block_levels"] = lev["block_levels"];
  out["period_levels"] = lev["period_levels"];
  return out;
}

// Corrected MIXED sampler: clustering model with base variances fixed (25)
// [[Rcpp::export]]
List mixed_sampler_cpp_corrected(const int iter,
                                 SEXP ID, SEXP Y, SEXP A, SEXP t, SEXP B, SEXP P,
                                 const double clust_prob     = 0.3,
                                 const double p_zeta_alpha   = 0.1,
                                 const double p_zeta_gamma   = 0.1,
                                 const double p_zeta_delta   = 0.1,
                                 const double sigma_sq_alpha = 1.0,
                                 const double sigma_sq_gamma = 1.0,
                                 const double sigma_sq_delta = 1.0,
                                 const double sigma_sq_mu    = 25.0,
                                 const double sigma_sq_beta  = 25.0,
                                 const double q_sd_mu   = 1.0,
                                 const double q_sd_beta = 1.0,
                                 const bool   ar1 = false,
                                 const double rho_jump = 0.05) {
  return clustering_sampler_cpp_corrected(iter, ID, Y, A, t, B, P, clust_prob,
                                          p_zeta_alpha, p_zeta_gamma, p_zeta_delta,
                                          sigma_sq_alpha, sigma_sq_gamma, sigma_sq_delta,
                                          5.0, 1.0, sigma_sq_mu, sigma_sq_beta,
                                          q_sd_mu, q_sd_beta, ar1, rho_jump);
}
