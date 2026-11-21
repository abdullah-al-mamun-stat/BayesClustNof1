# BayesClustNof1 0.1.0

* First release. This version contains the code used for all results in
  Al-Mamun and Chapple (in preparation): `clustering_sampler()`,
  `hierarchical_sampler()` (including the no-pooling reference with fixed
  variances), `estimate_partition()`, and the data `nof1_example` and
  `raynaud`.
* Compared with the samplers in the dissertation (development package
  `mynof1pkg` 0.0.0.9000):
  * the home values of the clustering model have a proper base distribution
    N(m, tau^2) with a half-Cauchy prior on tau, so the hierarchical model is
    the limit p_c -> 0;
  * the acceptance probabilities of the merge, switch and return-home moves are
    derived from detailed balance (checked against exact posteriors);
  * the hierarchical model estimates the between-patient variances with
    half-Cauchy priors (fixed variances are still available);
  * spike-and-slab indicators are drawn by collapsed Gibbs sampling, the error
    variance has prior 1/sigma^2, and AR(1) errors are available (`ar1 = TRUE`).
