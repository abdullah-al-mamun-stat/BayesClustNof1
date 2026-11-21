// =====================================================================
// HIERARCHICAL_SAMPLER_CORRECTED.cpp
//
// C++ version of hierarchical_sampler_corrected() in
// R CODE/HIERARCHICAL_SAMPLER_CORRECTED.R. With the same set.seed() it
// returns the same draws as the R function (checked in
// TESTS/test_R_vs_Cpp_identical.R).
//
// Use from R:
//   Rcpp::sourceCpp("HIERARCHICAL_SAMPLER_CORRECTED.cpp")   # needs NOF1_COMMON_CORRECTED.h in the same folder
//   fit <- hierarchical_sampler_cpp_corrected(iter, ID, Y, A, t, B, P)
//
// sigma_sq_mu / sigma_sq_beta: a negative value (default -1) means
// "estimate with a half-Cauchy prior"; a positive value fixes it (e.g. 25).
// =====================================================================
#include "NOF1_COMMON_CORRECTED.h"   // shared helpers (in this src folder)

// [[Rcpp::export]]
List hierarchical_sampler_cpp_corrected(const int iter,
                                        SEXP ID, SEXP Y, SEXP A, SEXP t, SEXP B, SEXP P,
                                        const double p_zeta_alpha   = 0.1,
                                        const double p_zeta_gamma   = 0.1,
                                        const double p_zeta_delta   = 0.1,
                                        const double sigma_sq_alpha = 1.0,
                                        const double sigma_sq_gamma = 1.0,
                                        const double sigma_sq_delta = 1.0,
                                        const double sigma_sq_mu    = -1.0,  // < 0: estimate
                                        const double sigma_sq_beta  = -1.0,  // < 0: estimate
                                        const double A_mu   = 5.0,
                                        const double A_beta = 1.0,
                                        const bool   ar1 = false,
                                        const double rho_jump = 0.05) {

  const Design d = make_design_corrected(ID, Y, A, t, B, P);
  const int N_id = d.N_id, N_tot = d.N_tot, K1 = d.K_minus_1;
  const bool estimate_s2_mu = (sigma_sq_mu < 0.0), estimate_s2_beta = (sigma_sq_beta < 0.0);

  // ---------------- storage (same names as the R function) ----------------
  NumericMatrix mu_store(iter, N_id), beta_store(iter, N_id);
  NumericVector alpha_store(iter), prob_alpha_store(iter), delta_store(iter), prob_delta_store(iter);
  IntegerVector zeta_alpha_store(iter), zeta_delta_store(iter);
  NumericMatrix gamma_store(iter, K1), prob_gamma_store(iter, K1);
  IntegerMatrix zeta_gamma_store(iter, K1);
  NumericVector sigma_sq_store(iter), rho_store(iter), mu_tilde_store(iter), beta_tilde_store(iter),
                sigma_sq_mu_store(iter), sigma_sq_beta_store(iter);

  // ---------------- initial values (deterministic) ----------------
  dvec mu(N_id), beta(N_id, 0.0);
  for (int i = 0; i < N_id; ++i) {
    long double s = 0.0;
    for (std::size_t h = 0; h < d.obs_by_subj[i].size(); ++h) s += d.Y[d.obs_by_subj[i][h]];
    mu[i] = (double) s / (double) d.obs_by_subj[i].size();
  }
  Nuisance nu;
  nu.gamma.assign(K1, 0.0); nu.zeta_gamma.assign(K1, 0); nu.prob_gamma.assign(K1, NA_REAL);
  double rho = 0.0;
  dvec res0(N_tot);
  for (int j = 0; j < N_tot; ++j) res0[j] = d.Y[j] - mu[d.ID0[j]];
  double sig2 = sum_prod_ld(res0, res0) / (double) N_tot;
  double mu_tilde = sum_ld(mu) / (double) N_id, beta_tilde = 0.0;
  double s2_mu   = estimate_s2_mu   ? A_mu * A_mu     : sigma_sq_mu;
  double s2_beta = estimate_s2_beta ? A_beta * A_beta : sigma_sq_beta;
  double a_mu = 1.0, a_beta = 1.0;

  dvec fit_mb(N_tot), r(N_tot), Tv, Zv;

  // ---------------- MCMC ----------------
  for (int jj = 0; jj < iter; ++jj) {
    if (jj % 500 == 0) Rcpp::checkUserInterrupt();

    // (1) transformed data for the current rho
    const Transformed tr = transform_design(d, rho);
    for (int j = 0; j < N_tot; ++j) fit_mb[j] = mu[d.ID0[j]] * tr.X0[j] + beta[d.ID0[j]] * tr.As[j];

    // (2) spike-and-slab design effects (collapsed Gibbs)
    update_nuisance(nu, tr, fit_mb, sig2, N_tot, K1, p_zeta_alpha, p_zeta_gamma, p_zeta_delta,
                    sigma_sq_alpha, sigma_sq_gamma, sigma_sq_delta);
    const dvec nuis = nuisance_term(nu, tr, N_tot);

    // (3) error variance, prior 1/sigma^2
    sig2 = update_sigma_sq(tr, fit_mb, nuis, N_tot);

    // (4) mu_i | rest (exact normal draws), mu_tilde, sigma_sq_mu
    for (int j = 0; j < N_tot; ++j) r[j] = tr.Ys[j] - beta[d.ID0[j]] * tr.As[j] - nuis[j];
    subject_TZ(tr.X0, r, d.obs_by_subj, Tv, Zv);
    for (int i = 0; i < N_id; ++i) {
      const double v = 1.0 / (Tv[i] / sig2 + 1.0 / s2_mu);
      const double m = v * (Zv[i] / sig2 + mu_tilde / s2_mu);
      mu[i] = R::rnorm(m, std::sqrt(v));
    }
    mu_tilde = R::rnorm(sum_ld(mu) / (double) N_id, std::sqrt(s2_mu / (double) N_id));
    if (estimate_s2_mu) hc_update(mu, mu_tilde, a_mu, A_mu, s2_mu);

    // (5) beta_i | rest, beta_tilde, sigma_sq_beta
    for (int j = 0; j < N_tot; ++j) r[j] = tr.Ys[j] - mu[d.ID0[j]] * tr.X0[j] - nuis[j];
    subject_TZ(tr.As, r, d.obs_by_subj, Tv, Zv);
    for (int i = 0; i < N_id; ++i) {
      const double v = 1.0 / (Tv[i] / sig2 + 1.0 / s2_beta);
      const double m = v * (Zv[i] / sig2 + beta_tilde / s2_beta);
      beta[i] = R::rnorm(m, std::sqrt(v));
    }
    beta_tilde = R::rnorm(sum_ld(beta) / (double) N_id, std::sqrt(s2_beta / (double) N_id));
    if (estimate_s2_beta) hc_update(beta, beta_tilde, a_beta, A_beta, s2_beta);

    // (6) AR(1) coefficient
    if (ar1) rho = update_rho(rho, rho_jump, raw_residuals(d, mu, beta, nu), d, sig2);

    // (7) store
    store_row(mu_store, jj, mu); store_row(beta_store, jj, beta);
    alpha_store[jj] = nu.alpha; zeta_alpha_store[jj] = nu.zeta_alpha; prob_alpha_store[jj] = nu.prob_alpha;
    if (K1 > 0) { store_row(gamma_store, jj, nu.gamma); store_row_int(zeta_gamma_store, jj, nu.zeta_gamma);
                  store_row(prob_gamma_store, jj, nu.prob_gamma); }
    delta_store[jj] = nu.delta; zeta_delta_store[jj] = nu.zeta_delta; prob_delta_store[jj] = nu.prob_delta;
    sigma_sq_store[jj] = sig2; rho_store[jj] = rho;
    mu_tilde_store[jj] = mu_tilde; beta_tilde_store[jj] = beta_tilde;
    sigma_sq_mu_store[jj] = s2_mu; sigma_sq_beta_store[jj] = s2_beta;
  }

  List lev = design_levels(d);
  return List::create(
    _["mu"] = mu_store, _["beta"] = beta_store,
    _["alpha"] = alpha_store, _["zeta_alpha"] = zeta_alpha_store, _["prob_alpha"] = prob_alpha_store,
    _["gamma"] = gamma_store, _["zeta_gamma"] = zeta_gamma_store, _["prob_gamma"] = prob_gamma_store,
    _["delta"] = delta_store, _["zeta_delta"] = zeta_delta_store, _["prob_delta"] = prob_delta_store,
    _["sigma_sq"] = sigma_sq_store, _["rho"] = rho_store,
    _["mu_tilde"] = mu_tilde_store, _["beta_tilde"] = beta_tilde_store,
    _["sigma_sq_mu"] = sigma_sq_mu_store, _["sigma_sq_beta"] = sigma_sq_beta_store,
    _["trt_levels"] = lev["trt_levels"], _["block_levels"] = lev["block_levels"],
    _["period_levels"] = lev["period_levels"]);
}
