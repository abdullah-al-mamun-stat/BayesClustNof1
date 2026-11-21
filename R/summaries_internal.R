# =====================================================================
# Internal helpers used by the sampler wrappers (not exported).
# Same summaries as the original package wrappers; the first half of the
# iterations is discarded as burn-in.
# =====================================================================

ci95   <- function(x) unname(stats::quantile(x, c(0.025, 0.975), na.rm = TRUE))
fmt_ci <- function(ci) paste0("(", round(ci[1], 3), ", ", round(ci[2], 3), ")")

# posterior co-clustering matrix from a matrix of cluster labels (draws x subjects)
coclustering_matrix <- function(label_draws) {
  N <- ncol(label_draws); P <- matrix(0, N, N)
  for (s in seq_len(nrow(label_draws))) P <- P + outer(label_draws[s, ], label_draws[s, ], "==")
  P / nrow(label_draws)
}

summarise_common <- function(results, iter, ids, extra_other = list()) {
  half_iter <- ceiling(iter / 2)
  keep <- (half_iter + 1):iter

  mu_last    <- results$mu[keep, , drop = FALSE]
  beta_last  <- results$beta[keep, , drop = FALSE]
  gamma_last <- results$gamma[keep, , drop = FALSE]
  zeta_gamma_l <- results$zeta_gamma[keep, , drop = FALSE]
  if (ncol(gamma_last) > 0) {
    gnames <- paste0("gamma_", results$block_levels$nonref)
    colnames(gamma_last) <- gnames; colnames(zeta_gamma_l) <- gnames
  }

  subject_df <- data.frame(
    ID = ids, mu = colMeans(mu_last), beta = colMeans(beta_last),
    mu_95CI = apply(apply(mu_last, 2, ci95), 2, fmt_ci),
    beta_95CI = apply(apply(beta_last, 2, ci95), 2, fmt_ci))

  oth <- c(list(sigma_sq_epsilon = results$sigma_sq[keep]), lapply(extra_other, function(x) x[keep]))
  others_df <- data.frame(parameter = names(oth), mean = sapply(oth, mean),
                          lower_95 = sapply(oth, function(x) ci95(x)[1]),
                          upper_95 = sapply(oth, function(x) ci95(x)[2]), row.names = NULL)

  mu_bar <- rowMeans(mu_last); beta_bar <- rowMeans(beta_last)
  bar_df <- data.frame(parameter = c("mu_bar", "beta_bar"), mean = c(mean(mu_bar), mean(beta_bar)),
                       lower_95 = c(ci95(mu_bar)[1], ci95(beta_bar)[1]),
                       upper_95 = c(ci95(mu_bar)[2], ci95(beta_bar)[2]))

  alpha_df <- data.frame(parameter = c("alpha", "Pr(zeta_alpha=1)"),
                         mean = c(mean(results$alpha[keep]), mean(results$zeta_alpha[keep] == 1L)),
                         lower_95 = c(ci95(results$alpha[keep])[1], NA),
                         upper_95 = c(ci95(results$alpha[keep])[2], NA),
                         RaoBlackwell_Pr = c(NA, mean(results$prob_alpha[keep])), row.names = NULL)
  if (ncol(gamma_last) > 0) {
    gamma_df <- data.frame(parameter = colnames(gamma_last), mean = colMeans(gamma_last),
                           lower_95 = apply(gamma_last, 2, function(x) ci95(x)[1]),
                           upper_95 = apply(gamma_last, 2, function(x) ci95(x)[2]),
                           Pr_included = colMeans(zeta_gamma_l == 1L),
                           RaoBlackwell_Pr = colMeans(results$prob_gamma[keep, , drop = FALSE]), row.names = NULL)
  } else {
    gamma_df <- data.frame(parameter = character(0), mean = numeric(0), lower_95 = numeric(0),
                           upper_95 = numeric(0), Pr_included = numeric(0), RaoBlackwell_Pr = numeric(0))
  }
  delta_df <- data.frame(parameter = c("delta", "Pr(zeta_delta=1)"),
                         mean = c(mean(results$delta[keep]), mean(results$zeta_delta[keep] == 1L)),
                         lower_95 = c(ci95(results$delta[keep])[1], NA),
                         upper_95 = c(ci95(results$delta[keep])[2], NA),
                         RaoBlackwell_Pr = c(NA, mean(results$prob_delta[keep])), row.names = NULL)

  posterior_samples <- list(mu = mu_last, beta = beta_last,
                            alpha = results$alpha[keep], zeta_alpha = results$zeta_alpha[keep],
                            gamma = gamma_last, zeta_gamma = zeta_gamma_l,
                            delta = results$delta[keep], zeta_delta = results$zeta_delta[keep],
                            sigma_sq = results$sigma_sq[keep], rho = results$rho[keep])
  for (nm in names(extra_other)) posterior_samples[[nm]] <- extra_other[[nm]][keep]

  list(subject_summary = subject_df, other_summary = others_df, bar_summary = bar_df,
       alpha_summary = alpha_df, gamma_summary = gamma_df, delta_summary = delta_df,
       posterior_samples = posterior_samples, mu_bar = mu_bar, beta_bar = beta_bar,
       block_levels = results$block_levels, period_levels = results$period_levels, keep = keep)
}

# add clustering-specific summaries (co-clustering, number of clusters, Pr(home))
add_cluster_summaries <- function(out, results, lab_mu, lab_beta, home_labels = TRUE) {
  keep <- out$keep
  Lmu <- results[[lab_mu]][keep, , drop = FALSE]; Lbeta <- results[[lab_beta]][keep, , drop = FALSE]
  out$posterior_samples[[lab_mu]] <- Lmu; out$posterior_samples[[lab_beta]] <- Lbeta
  out$coclustering_mu <- coclustering_matrix(Lmu)
  out$coclustering_beta <- coclustering_matrix(Lbeta)
  N <- ncol(Lbeta)
  n_clusters <- function(L) apply(L, 1, function(z) length(unique(z)))
  out$n_clusters_summary <- data.frame(
    parameter = c("mu", "beta"),
    mean = c(mean(n_clusters(Lmu)), mean(n_clusters(Lbeta))),
    median = c(stats::median(n_clusters(Lmu)), stats::median(n_clusters(Lbeta))))
  if (home_labels) {
    idx <- matrix(seq_len(N), nrow(Lbeta), N, byrow = TRUE)
    out$clustering_summary <- data.frame(ID = out$subject_summary$ID,
                                         Pr_home_mu = colMeans(Lmu == idx),
                                         Pr_home_beta = colMeans(Lbeta == idx))
  }
  out
}
