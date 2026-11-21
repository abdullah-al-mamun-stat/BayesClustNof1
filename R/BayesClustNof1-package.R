#' BayesClustNof1: Bayesian clustering and hierarchical models for series of N-of-1 trials
#'
#' The models of the accompanying paper:
#' \code{\link{clustering_sampler}} (home-host clustering model) and
#' \code{\link{hierarchical_sampler}} (hierarchical model; with fixed
#' variances it gives the no-pooling reference). The partition that minimises
#' the posterior expected variation of information is given by
#' \code{\link{estimate_partition}}.
#' Data: \code{\link{nof1_example}} (simulated) and \code{\link{raynaud}}.
#' See \code{vignette("BayesClustNof1")} for a worked example.
#' @keywords internal
#' @useDynLib BayesClustNof1, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats quantile median
"_PACKAGE"
