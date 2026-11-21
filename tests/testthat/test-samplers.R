# small simulated N-of-1 dataset: 6 subjects, 3 blocks, 2 periods of 5 days
make_data <- function(seed = 1) {
  set.seed(seed)
  rows <- list()
  for (i in 1:6) for (b in 1:3) { ord <- sample(0:1)
    for (p in 1:2) rows[[length(rows) + 1]] <- data.frame(ID = i, B = b, P = p, t = 0:4, A = ord[p]) }
  d <- do.call(rbind, rows)
  beta_true <- c(0, 0, 0, 1.5, 1.5, 1.5)
  d$Y <- c(2, 4, 6, 3, 5, 7)[d$ID] + beta_true[d$ID] * d$A + rnorm(nrow(d), 0, 0.5)
  d
}

test_that("samplers run and return the documented components", {
  d <- make_data()
  for (f in list(hierarchical_sampler, clustering_sampler)) {
    fit <- f(400, d$ID, d$Y, d$A, d$t, d$B, d$P)
    expect_equal(nrow(fit$subject_summary), 6)
    expect_true(all(c("subject_summary", "other_summary", "alpha_summary", "gamma_summary",
                      "delta_summary", "posterior_samples") %in% names(fit)))
    expect_equal(nrow(fit$posterior_samples$beta), 200)
  }
})

test_that("same seed gives identical draws", {
  d <- make_data()
  set.seed(5); a <- clustering_sampler(300, d$ID, d$Y, d$A, d$t, d$B, d$P, clust_prob = 0.5, ar1 = TRUE)
  set.seed(5); b <- clustering_sampler(300, d$ID, d$Y, d$A, d$t, d$B, d$P, clust_prob = 0.5, ar1 = TRUE)
  expect_identical(a$raw$beta, b$raw$beta)
  expect_identical(a$raw$zeta_beta, b$raw$zeta_beta)
})

test_that("treatment effects are recovered and responders cluster together", {
  d <- make_data(2)
  set.seed(3); fit <- clustering_sampler(4000, d$ID, d$Y, d$A, d$t, d$B, d$P, clust_prob = 0.5)
  expect_lt(max(abs(fit$subject_summary$beta - c(0, 0, 0, 1.5, 1.5, 1.5))), 0.5)
  expect_gt(fit$coclustering_beta[4, 5], 0.5)
  expect_lt(fit$coclustering_beta[1, 4], 0.2)
})

test_that("partition estimate behaves", {
  L <- rbind(c(1, 1, 2, 2), c(1, 1, 2, 2), c(1, 2, 2, 2))
  expect_equal(estimate_partition(L, "Binder")$labels, c(1, 1, 2, 2))
  expect_equal(estimate_partition(L, "VI")$labels, c(1, 1, 2, 2))
})

test_that("included data sets have the documented structure", {
  expect_equal(dim(nof1_example), c(600L, 7L))
  expect_equal(length(unique(raynaud$ID)), 38L)
  expect_equal(nrow(raynaud), 1537L)
  expect_true(all(raynaud$A %in% 0:1) && all(raynaud$P %in% 1:2))
  expect_false(is.unsorted(raynaud$ID))
})

test_that("no-pooling reference (fixed variances) runs", {
  d <- make_data()
  fit <- hierarchical_sampler(400, d$ID, d$Y, d$A, d$t, d$B, d$P, sigma_sq_mu = 25, sigma_sq_beta = 25)
  expect_equal(nrow(fit$subject_summary), 6)
})
