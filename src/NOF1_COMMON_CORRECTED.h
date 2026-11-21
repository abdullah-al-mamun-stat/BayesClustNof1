// =====================================================================
// NOF1_COMMON_CORRECTED.h
//
// Helper functions shared by the corrected C++ samplers:
//   HIERARCHICAL_SAMPLER_CORRECTED.cpp
//   CLUSTERING_SAMPLER_CORRECTED.cpp   (also the corrected MIXED sampler)
//   DP_SAMPLER.cpp                     (Dirichlet-process comparison model)
//
// Each function is a line-by-line translation of the R function with the
// same name in NOF1_COMMON_FUNCTIONS_CORRECTED.R, so that R and C++ use the
// random numbers in the same order and give the same draws for the same seed:
//   * random numbers only through R::runif, R::rnorm, R::rgamma, R::rbeta
//     (the same generators R uses);
//   * a random index is floor(runif * n), never sample();
//   * sums are accumulated in long double, in the same order as R's sum().
//
// Indices: C++ uses 0-based indices internally; zeta is returned 1-based.
// =====================================================================
#ifndef NOF1_COMMON_CORRECTED_H
#define NOF1_COMMON_CORRECTED_H

#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

typedef std::vector<double> dvec;
typedef std::vector<int>    ivec;

// ---------------------------------------------------------------------
// Sums in long double, exactly like R's sum()
// ---------------------------------------------------------------------
inline double sum_ld(const dvec& x) {
  long double s = 0.0;
  for (std::size_t j = 0; j < x.size(); ++j) s += x[j];
  return (double) s;
}
inline double sum_prod_ld(const dvec& x, const dvec& y) {        // sum(x * y)
  long double s = 0.0;
  for (std::size_t j = 0; j < x.size(); ++j) s += x[j] * y[j];
  return (double) s;
}

// ---------------------------------------------------------------------
// 1) Design (same coding as your original samplers)
// ---------------------------------------------------------------------
struct Design {
  int N_id, N_tot, K_minus_1, n_periods;
  ivec ID0;                         // subject index 0..N_id-1 for every row
  dvec Y, t, A_dummy, P_dummy;
  std::vector<dvec> B_cols;         // block dummy columns
  std::vector<bool> first_in_period;
  std::vector<ivec> obs_by_subj;    // row numbers of each subject
  double trt_ref, trt_other, block_ref, per_ref, per_other;
  NumericVector block_nonref;
  IntegerVector ids;
};

inline NumericVector coerce_numeric_like_r(SEXP x) {
  Function as_numeric("as.numeric");
  Function as_character("as.character");
  if (Rf_isFactor(x)) return as<NumericVector>(as_numeric(as_character(x)));
  return as<NumericVector>(as_numeric(x));
}
inline NumericVector sort_unique_num(SEXP x) {
  Function sort_f("sort"); Function unique_f("unique");
  return as<NumericVector>(sort_f(unique_f(x)));
}

inline Design make_design_corrected(SEXP ID, SEXP Y, SEXP A, SEXP t, SEXP B, SEXP P) {
  Design d;
  Function sort_f("sort"); Function unique_f("unique"); Function match_f("match");
  d.ids = as<IntegerVector>(sort_f(unique_f(ID)));
  IntegerVector mapped = as<IntegerVector>(match_f(ID, d.ids));
  NumericVector Yv = coerce_numeric_like_r(Y), tv = coerce_numeric_like_r(t);
  d.N_tot = Yv.size(); d.N_id = d.ids.size();
  if (mapped.size() != d.N_tot || tv.size() != d.N_tot) stop("ID, Y and t must have the same length.");
  d.ID0.resize(d.N_tot); d.Y.resize(d.N_tot); d.t.resize(d.N_tot);
  for (int j = 0; j < d.N_tot; ++j) { d.ID0[j] = mapped[j] - 1; d.Y[j] = Yv[j]; d.t[j] = tv[j]; }
  for (int j = 1; j < d.N_tot; ++j)
    if (d.ID0[j] < d.ID0[j - 1]) stop("Sort the data by subject (and by time within subject) first.");

  // treatment dummy (ref = min level)
  NumericVector A_num = coerce_numeric_like_r(A);
  NumericVector trt_levels = sort_unique_num(A_num);
  d.trt_ref = trt_levels[0];
  d.trt_other = (trt_levels.size() > 1) ? trt_levels[1] : NA_REAL;
  d.A_dummy.assign(d.N_tot, 0.0);
  if (trt_levels.size() > 1) for (int j = 0; j < d.N_tot; ++j) d.A_dummy[j] = (A_num[j] == d.trt_other) ? 1.0 : 0.0;

  // block dummies (ref = min level)
  NumericVector B_num = coerce_numeric_like_r(B);
  NumericVector block_levels = sort_unique_num(B_num);
  d.block_ref = block_levels[0];
  d.K_minus_1 = block_levels.size() - 1;
  d.block_nonref = NumericVector(d.K_minus_1);
  d.B_cols.assign(d.K_minus_1, dvec(d.N_tot, 0.0));
  for (int k = 0; k < d.K_minus_1; ++k) {
    d.block_nonref[k] = block_levels[k + 1];
    for (int j = 0; j < d.N_tot; ++j) d.B_cols[k][j] = (B_num[j] == block_levels[k + 1]) ? 1.0 : 0.0;
  }

  // period dummy (ref = min level)
  NumericVector P_num = coerce_numeric_like_r(P);
  NumericVector per_levels = sort_unique_num(P_num);
  d.per_ref = per_levels[0];
  d.per_other = (per_levels.size() > 1) ? per_levels[1] : NA_REAL;
  d.P_dummy.assign(d.N_tot, 0.0);
  if (per_levels.size() > 1) for (int j = 0; j < d.N_tot; ++j) d.P_dummy[j] = (P_num[j] == d.per_other) ? 1.0 : 0.0;
  if (per_levels.size() > 2)
    Rf_warning("P has more than 2 levels; as in the original code only the second-lowest level gets the period dummy. P should be the period position within a block.");

  // first day of each treatment period (AR(1) restarts there)
  d.first_in_period.assign(d.N_tot, true);
  d.n_periods = 1;
  for (int j = 1; j < d.N_tot; ++j) {
    d.first_in_period[j] = (d.ID0[j] != d.ID0[j - 1]) || (B_num[j] != B_num[j - 1]) || (P_num[j] != P_num[j - 1]);
    if (d.first_in_period[j]) d.n_periods++;
  }

  d.obs_by_subj.assign(d.N_id, ivec());
  for (int j = 0; j < d.N_tot; ++j) d.obs_by_subj[d.ID0[j]].push_back(j);
  return d;
}

// ---------------------------------------------------------------------
// 2) AR(1): Prais-Winsten transform within each treatment period
// ---------------------------------------------------------------------
inline dvec pw_transform(const dvec& x, const std::vector<bool>& first, const double rho) {
  if (rho == 0.0) return x;
  dvec out(x.size());
  const double c1 = std::sqrt(1.0 - rho * rho);
  for (std::size_t j = 0; j < x.size(); ++j) {
    if (first[j]) out[j] = c1 * x[j];
    else          out[j] = x[j] - rho * x[j - 1];
  }
  return out;
}

struct Transformed {                // transformed outcome and design columns
  dvec Ys, X0, As, Ts, Ps;
  std::vector<dvec> Bs;
};

inline Transformed transform_design(const Design& d, const double rho) {
  Transformed tr;
  tr.Ys = pw_transform(d.Y, d.first_in_period, rho);
  tr.X0 = pw_transform(dvec(d.N_tot, 1.0), d.first_in_period, rho);
  tr.As = pw_transform(d.A_dummy, d.first_in_period, rho);
  tr.Ts = pw_transform(d.t, d.first_in_period, rho);
  tr.Ps = pw_transform(d.P_dummy, d.first_in_period, rho);
  tr.Bs.resize(d.K_minus_1);
  for (int k = 0; k < d.K_minus_1; ++k) tr.Bs[k] = pw_transform(d.B_cols[k], d.first_in_period, rho);
  return tr;
}

// block term sum_k Bs[, k] * gamma[k], added column by column (same as R)
inline dvec block_term(const std::vector<dvec>& Bs, const dvec& gamma_vec, const int n) {
  dvec out(n, 0.0);
  for (std::size_t k = 0; k < gamma_vec.size(); ++k)
    for (int j = 0; j < n; ++j) out[j] = out[j] + Bs[k][j] * gamma_vec[k];
  return out;
}

// ---------------------------------------------------------------------
// 3) Random index 0..n-1 with runif (R: floor(runif(1) * n) + 1)
// ---------------------------------------------------------------------
inline int draw_index(const int n) {
  return (int) std::floor(R::runif(0.0, 1.0) * (double) n);
}

// ---------------------------------------------------------------------
// 4) Collapsed Gibbs update for ONE spike-and-slab coefficient
// ---------------------------------------------------------------------
struct SSResult { int z; double coef; double prob; };

inline SSResult ss_gibbs_update(const dvec& r, const dvec& x, const double sig2,
                                const double sig2_slab, const double p_incl) {
  const double Q   = sum_prod_ld(x, x) / sig2 + 1.0 / sig2_slab;
  const double b   = sum_prod_ld(x, r) / sig2;
  const double lBF = -0.5 * std::log(sig2_slab * Q) + (b * b) / (2.0 * Q);
  const double lo  = std::log(p_incl) - std::log(1.0 - p_incl) + lBF;
  SSResult out;
  out.prob = 1.0 / (1.0 + std::exp(-lo));
  out.z    = (R::runif(0.0, 1.0) < out.prob) ? 1 : 0;
  out.coef = (out.z == 1) ? R::rnorm(b / Q, std::sqrt(1.0 / Q)) : 0.0;
  return out;
}

// ---------------------------------------------------------------------
// 5) Half-Cauchy(0, A) scale update (inverse-gamma mixture)
// ---------------------------------------------------------------------
inline void hc_update(const dvec& values, const double centre, double& a, const double A, double& s2) {
  long double S_ld = 0.0;
  for (std::size_t j = 0; j < values.size(); ++j) S_ld += (values[j] - centre) * (values[j] - centre);
  const double S = (double) S_ld;
  const double n = (double) values.size();
  s2 = 1.0 / R::rgamma((n + 1.0) / 2.0, 1.0 / (1.0 / a + S / 2.0));   // R::rgamma uses scale
  a  = 1.0 / R::rgamma(1.0, 1.0 / (1.0 / (A * A) + 1.0 / s2));
}

// ---------------------------------------------------------------------
// 6) Per-subject summaries T_i = sum x^2, Z_i = sum x r
// ---------------------------------------------------------------------
inline void subject_TZ(const dvec& x, const dvec& r, const std::vector<ivec>& obs, dvec& Tv, dvec& Zv) {
  const int N = obs.size();
  Tv.assign(N, 0.0); Zv.assign(N, 0.0);
  for (int i = 0; i < N; ++i) {
    long double sT = 0.0, sZ = 0.0;
    for (std::size_t h = 0; h < obs[i].size(); ++h) {
      const int j = obs[i][h];
      sT += x[j] * x[j];
      sZ += x[j] * r[j];
    }
    Tv[i] = (double) sT; Zv[i] = (double) sZ;
  }
}

// ---------------------------------------------------------------------
// 7) Exact AR(1) log-likelihood (+ Jacobian of rho = tanh(theta))
// ---------------------------------------------------------------------
inline double ar1_log_target(const double rho, const dvec& resid_raw, const std::vector<bool>& first,
                             const int n_periods, const double sig2) {
  const dvec es = pw_transform(resid_raw, first, rho);
  const double ll = 0.5 * (double) n_periods * std::log(1.0 - rho * rho) - 0.5 * sum_prod_ld(es, es) / sig2;
  return ll + std::log(1.0 - rho * rho);
}

// random-walk update of rho on the atanh scale; returns the new rho
inline double update_rho(const double rho, const double rho_jump, const dvec& resid_raw,
                         const Design& d, const double sig2) {
  const double rho_prop = std::tanh(std::atanh(rho) + R::rnorm(0.0, rho_jump));
  const double diff = ar1_log_target(rho_prop, resid_raw, d.first_in_period, d.n_periods, sig2) -
                      ar1_log_target(rho,      resid_raw, d.first_in_period, d.n_periods, sig2);
  return (std::log(R::runif(0.0, 1.0)) < diff) ? rho_prop : rho;
}

// ---------------------------------------------------------------------
// 8) YOUR home-host clustering moves, corrected (N moves per call)
//    zeta: 0-based here (zeta[i] == i means home)
// ---------------------------------------------------------------------
inline double ll_i(const int i, const double v, const dvec& Tv, const dvec& Zv, const double sig2) {
  return (v * Zv[i] - 0.5 * v * v * Tv[i]) / sig2;
}

inline int update_clusters_corrected(ivec& zeta, dvec& val, const dvec& Tv, const dvec& Zv,
                                     const double sig2, const double m, const double tau,
                                     const double r, const double q_sd) {
  const int N = zeta.size();
  int accepted = 0;
  for (int mv = 0; mv < N; ++mv) {
    const int i = draw_index(N);                        // patient chosen uniformly
    ivec H;                                             // current homes, ascending
    for (int h = 0; h < N; ++h) if (zeta[h] == h) H.push_back(h);
    const int K = H.size();

    if (zeta[i] == i) {
      // ---- i is a home ----
      int n_members = 0;
      for (int h = 0; h < N; ++h) if (zeta[h] == i) n_members++;
      if (n_members > 1 || K == 1) continue;            // home with guests / only cluster: no move
      ivec others;
      for (int h = 0; h < K; ++h) if (H[h] != i) others.push_back(H[h]);
      const int j = others[draw_index(K - 1)];          // merge target
      const double logR = ll_i(i, val[j], Tv, Zv, sig2) - ll_i(i, val[i], Tv, Zv, sig2) +
        std::log(r) +
        std::log((double)(K - 1)) +
        R::dnorm4(val[i], val[j], q_sd, 1) -
        R::dnorm4(val[i], m, tau, 1) -
        std::log(2.0);
      if (std::log(R::runif(0.0, 1.0)) < logR) { zeta[i] = j; val[i] = val[j]; accepted++; }
    } else {
      // ---- i is a guest of home h ----
      const int h = zeta[i];
      if (R::runif(0.0, 1.0) < 0.5) {
        // switch host
        if (K == 1) continue;
        ivec others;
        for (int q = 0; q < K; ++q) if (H[q] != h) others.push_back(H[q]);
        const int k = others[draw_index(K - 1)];
        const double logR = ll_i(i, val[k], Tv, Zv, sig2) - ll_i(i, val[h], Tv, Zv, sig2);
        if (std::log(R::runif(0.0, 1.0)) < logR) { zeta[i] = k; val[i] = val[k]; accepted++; }
      } else {
        // return home with a new value u ~ N(val[h], q_sd^2)
        const double u = R::rnorm(val[h], q_sd);
        const double logR = ll_i(i, u, Tv, Zv, sig2) - ll_i(i, val[h], Tv, Zv, sig2) -
          std::log(r) +
          R::dnorm4(u, m, tau, 1) -
          std::log((double) K) -
          R::dnorm4(u, val[h], q_sd, 1) +
          std::log(2.0);
        if (std::log(R::runif(0.0, 1.0)) < logR) { zeta[i] = i; val[i] = u; accepted++; }
      }
    }
  }

  // exact (conjugate) draw of each home value, shared by its guests
  for (int h = 0; h < N; ++h) {
    if (zeta[h] != h) continue;
    long double sT = 0.0, sZ = 0.0;
    for (int q = 0; q < N; ++q) if (zeta[q] == h) { sT += Tv[q]; sZ += Zv[q]; }
    const double prec = (double) sT / sig2 + 1.0 / (tau * tau);
    const double mn   = ((double) sZ / sig2 + m / (tau * tau)) / prec;
    const double draw = R::rnorm(mn, std::sqrt(1.0 / prec));
    for (int q = 0; q < N; ++q) if (zeta[q] == h) val[q] = draw;
  }
  return accepted;
}

// ---------------------------------------------------------------------
// 9) Base distribution update (m flat; tau half-Cauchy unless fixed)
// ---------------------------------------------------------------------
inline void update_base_corrected(const dvec& home_values, double& m, double& tau, double& a,
                                  const double A_tau, const bool tau_fixed) {
  const double K = (double) home_values.size();
  m = R::rnorm(sum_ld(home_values) / K, tau / std::sqrt(K));
  if (!tau_fixed) {
    double s2;
    hc_update(home_values, m, a, A_tau, s2);
    tau = std::sqrt(s2);
  }
}

// ---------------------------------------------------------------------
// Spike-and-slab block shared by all samplers (alpha, gamma_k, delta)
// fit_mb = subject part of the mean on the transformed scale
// ---------------------------------------------------------------------
struct Nuisance {
  double alpha = 0.0, delta = 0.0, prob_alpha = NA_REAL, prob_delta = NA_REAL;
  int zeta_alpha = 0, zeta_delta = 0;
  dvec gamma, prob_gamma; ivec zeta_gamma;
};

inline void update_nuisance(Nuisance& nu, const Transformed& tr, const dvec& fit_mb, const double sig2,
                            const int n, const int K1,
                            const double p_zeta_alpha, const double p_zeta_gamma, const double p_zeta_delta,
                            const double sigma_sq_alpha, const double sigma_sq_gamma, const double sigma_sq_delta) {
  dvec r(n), bt;
  // alpha
  bt = block_term(tr.Bs, nu.gamma, n);
  for (int j = 0; j < n; ++j) r[j] = tr.Ys[j] - fit_mb[j] - bt[j] - nu.delta * tr.Ps[j];
  SSResult u = ss_gibbs_update(r, tr.Ts, sig2, sigma_sq_alpha, p_zeta_alpha);
  nu.alpha = u.coef; nu.zeta_alpha = u.z; nu.prob_alpha = u.prob;
  // gamma_k
  for (int k = 0; k < K1; ++k) {
    dvec g0 = nu.gamma; g0[k] = 0.0;
    bt = block_term(tr.Bs, g0, n);
    for (int j = 0; j < n; ++j) r[j] = tr.Ys[j] - fit_mb[j] - nu.alpha * tr.Ts[j] - bt[j] - nu.delta * tr.Ps[j];
    u = ss_gibbs_update(r, tr.Bs[k], sig2, sigma_sq_gamma, p_zeta_gamma);
    nu.gamma[k] = u.coef; nu.zeta_gamma[k] = u.z; nu.prob_gamma[k] = u.prob;
  }
  // delta
  bt = block_term(tr.Bs, nu.gamma, n);
  for (int j = 0; j < n; ++j) r[j] = tr.Ys[j] - fit_mb[j] - nu.alpha * tr.Ts[j] - bt[j];
  u = ss_gibbs_update(r, tr.Ps, sig2, sigma_sq_delta, p_zeta_delta);
  nu.delta = u.coef; nu.zeta_delta = u.z; nu.prob_delta = u.prob;
}

// nuisance part of the mean on the transformed scale
inline dvec nuisance_term(const Nuisance& nu, const Transformed& tr, const int n) {
  dvec bt = block_term(tr.Bs, nu.gamma, n), out(n);
  for (int j = 0; j < n; ++j) out[j] = nu.alpha * tr.Ts[j] + bt[j] + nu.delta * tr.Ps[j];
  return out;
}

// untransformed residuals for the rho update
inline dvec raw_residuals(const Design& d, const dvec& mu, const dvec& beta, const Nuisance& nu) {
  dvec bt = block_term(d.B_cols, nu.gamma, d.N_tot), out(d.N_tot);
  for (int j = 0; j < d.N_tot; ++j) {
    const int s = d.ID0[j];
    out[j] = d.Y[j] - (mu[s] + beta[s] * d.A_dummy[j]) - (nu.alpha * d.t[j] + bt[j] + nu.delta * d.P_dummy[j]);
  }
  return out;
}

// sigma^2 | rest with prior 1/sigma^2
inline double update_sigma_sq(const Transformed& tr, const dvec& fit_mb, const dvec& nuis, const int n) {
  dvec res(n);
  for (int j = 0; j < n; ++j) res[j] = tr.Ys[j] - fit_mb[j] - nuis[j];
  const double rate = sum_prod_ld(res, res) / 2.0;
  return 1.0 / R::rgamma((double) n / 2.0, 1.0 / rate);
}

// helpers to copy into R matrices
inline void store_row(NumericMatrix& M, const int row, const dvec& x) {
  for (std::size_t c = 0; c < x.size(); ++c) M(row, c) = x[c];
}
inline void store_row_int(IntegerMatrix& M, const int row, const ivec& x, const int add = 0) {
  for (std::size_t c = 0; c < x.size(); ++c) M(row, c) = x[c] + add;
}
inline List design_levels(const Design& d) {
  return List::create(
    _["trt_levels"]    = List::create(_["ref"] = d.trt_ref, _["other"] = d.trt_other),
    _["block_levels"]  = List::create(_["ref"] = d.block_ref, _["nonref"] = d.block_nonref),
    _["period_levels"] = List::create(_["ref"] = d.per_ref, _["other"] = d.per_other));
}

#endif
