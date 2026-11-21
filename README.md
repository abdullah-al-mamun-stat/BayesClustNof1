# BayesClustNof1

Bayesian clustering and hierarchical models for series of N-of-1 trials.

When several patients complete N-of-1 trials with a common protocol, the
standard hierarchical model pulls every patient's treatment effect toward the
same population mean. **BayesClustNof1** fits the *home–host clustering
model*, which groups patients with similar treatment effects and borrows
information mainly within those groups, together with the hierarchical model
for comparison. Both models include spike-and-slab priors for time, block and
period effects and optional AR(1) errors within treatment periods. The
samplers are written in C++ (Rcpp).

## Installation

```r
# install.packages("remotes")
remotes::install_github("abdullah-al-mamun-stat/BayesClustNof1", build_vignettes = TRUE)
```

A C++ compiler is needed (Rtools on Windows, Xcode command line tools on macOS).

## Example

```r
library(BayesClustNof1)
d <- nof1_example                      # 10 simulated patients; patients 3 and 7 respond

set.seed(1)
fit <- clustering_sampler(iter = 20000, ID = d$ID, Y = d$Y, A = d$A, t = d$t,
                          B = d$B, P = d$P, clust_prob = 0.3, ar1 = TRUE)
fit$subject_summary                    # patient-specific effects with 95% intervals
colMeans(fit$posterior_samples$beta > 0)   # posterior probability of benefit
fit$coclustering_beta                  # posterior probability that two patients share an effect

set.seed(1)
fit_h <- hierarchical_sampler(iter = 20000, ID = d$ID, Y = d$Y, A = d$A, t = d$t,
                              B = d$B, P = d$P, ar1 = TRUE)
```

The user guide, `vignette("BayesClustNof1")`, covers the data format, the
output, the choice of the clustering probability, convergence checks and the
reanalysis of the Raynaud data (`raynaud`) from the paper.

## Main functions

| Function | Model |
|---|---|
| `clustering_sampler()` | home–host clustering model |
| `hierarchical_sampler()` | hierarchical model (half-Cauchy priors on the between-patient SDs); with `sigma_sq_mu = 25, sigma_sq_beta = 25`, the no-pooling reference |
| `estimate_partition()` | partition minimising the posterior expected variation of information |

## Reproducibility

Version 0.1.0 contains the code that produced all results in the paper. With
the same seeds, the package gives posterior draws identical to those reported.

## Citation

Al-Mamun, M. A. and Chapple, A. G. (2026). Local versus global borrowing of
patient-specific treatment effects in series of N-of-1 trials: a Bayesian
home–host clustering model. Manuscript in preparation.

See `citation("BayesClustNof1")`.

## License

MIT
