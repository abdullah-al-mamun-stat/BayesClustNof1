# Point estimate of the partition (used in the paper: the partition that
# minimises the posterior expected variation of information).

#' Point estimate of the partition from posterior cluster labels
#'
#' Chooses, among the sampled partitions, the one with the smallest posterior
#' expected loss: Binder's loss (Dahl's least-squares method) or the variation
#' of information (Wade and Ghahramani, 2018).
#'
#' @param label_draws Matrix of cluster labels, draws in rows and patients in
#'   columns (e.g. \code{fit$posterior_samples$zeta_beta} from \code{\link{clustering_sampler}}).
#' @param loss \code{"VI"} or \code{"Binder"}.
#' @param max_draws Maximum number of draws used (thinned evenly) to keep the
#'   VI computation fast.
#' @return A list with \code{labels} (1..K), \code{groups} and \code{expected_loss}.
#' @export
estimate_partition <- function(label_draws, loss = c("VI", "Binder"), max_draws = 500) {
  loss <- match.arg(loss)
  if (nrow(label_draws) > max_draws)
    label_draws <- label_draws[round(seq(1, nrow(label_draws), length.out = max_draws)), , drop = FALSE]
  relabel <- function(z) match(z, unique(z))
  L <- t(apply(label_draws, 1, relabel))
  cand <- L[!duplicated(apply(L, 1, paste, collapse = ",")), , drop = FALSE]
  if (loss == "Binder") {
    P <- coclustering_matrix(L); ut <- upper.tri(P)
    el <- apply(cand, 1, function(z) sum(abs(outer(z, z, "==")[ut] - P[ut])))
  } else {
    # VI(a, b) = 2 H(a, b) - H(a) - H(b); entropies from fast integer counts
    N  <- ncol(L)
    H  <- function(counts) { p <- counts[counts > 0] / N; -sum(p * log(p)) }
    HL <- apply(L, 1, function(w) H(tabulate(w, N)))
    el <- apply(cand, 1, function(z) {
      Hz <- H(tabulate(z, N))
      mean(sapply(seq_len(nrow(L)), function(s) 2 * H(tabulate((z - 1L) * N + L[s, ], N * N)) - Hz - HL[s]))
    })
  }
  best <- cand[which.min(el), ]
  list(labels = best, groups = unname(split(seq_along(best), best)), expected_loss = min(el), loss = loss)
}
