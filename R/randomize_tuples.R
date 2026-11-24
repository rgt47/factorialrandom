#' Randomize treatment assignments within matched tuples
#'
#' Within each tuple of size k, the k treatment combinations
#' are assigned uniformly at random to the k subjects.
#'
#' @param tuples Matching result from match_tuples.
#' @param treatment_labels Character or integer vector of
#'   length k naming the treatment combinations. Default is
#'   0:(k-1).
#' @return An integer vector of treatment assignments
#'   (length = total subjects, NA for unmatched).
#' @export
randomize_tuples <- function(tuples,
                             treatment_labels = NULL) {
  k <- tuples$k
  if (is.null(treatment_labels)) {
    treatment_labels <- seq(0L, k - 1L)
  }
  stopifnot(length(treatment_labels) == k)

  n <- length(tuples$tuple_id)
  trt <- rep(NA_integer_, n)

  for (m in seq_len(tuples$n_tuples)) {
    members <- which(tuples$tuple_id == m)
    if (length(members) != k) next
    perm <- sample(treatment_labels)
    trt[members] <- perm
  }

  trt
}

#' Generate factorial treatment labels for a 2^F design
#'
#' @param n_factors Integer. Number of factors (F).
#' @return A data.frame with k = 2^F rows and F columns of
#'   0/1 indicators, plus a 'label' column with the
#'   treatment combination index (0 to k-1).
#' @export
factorial_labels <- function(n_factors) {
  k <- 2^n_factors
  grid <- expand.grid(
    replicate(n_factors, c(0L, 1L), simplify = FALSE)
  )
  names(grid) <- paste0("F", seq_len(n_factors))
  grid$label <- seq(0L, k - 1L)
  grid
}

#' Simple randomization for factorial design
#'
#' Assigns each of N subjects to one of k treatment
#' combinations with equal probability.
#'
#' @param n Integer. Number of subjects.
#' @param k Integer. Number of treatment combinations.
#' @return Integer vector of assignments (0 to k-1).
#' @export
randomize_simple_factorial <- function(n, k) {
  sample(0L:(k - 1L), n, replace = TRUE)
}

#' Stratified block randomization for factorial design
#'
#' Permuted block randomization within strata, with blocks
#' of size k.
#'
#' @param data data.frame of subject data.
#' @param strat_vars Character vector of stratification
#'   variable names.
#' @param k Integer. Number of treatment combinations.
#' @return Integer vector of assignments (0 to k-1).
#' @export
randomize_stratified_factorial <- function(data,
                                           strat_vars,
                                           k) {
  n <- nrow(data)
  trt <- integer(n)

  if (length(strat_vars) == 0) {
    return(randomize_block_factorial(n, k))
  }

  strata <- interaction(data[strat_vars], drop = TRUE)
  for (stratum in levels(strata)) {
    idx <- which(strata == stratum)
    trt[idx] <- randomize_block_factorial(
      length(idx), k
    )
  }
  trt
}

#' Permuted block randomization for a single stratum
#'
#' @param n Integer. Number of subjects.
#' @param k Integer. Number of treatment arms.
#' @return Integer vector of assignments (0 to k-1).
#' @keywords internal
randomize_block_factorial <- function(n, k) {
  n_full <- n %/% k
  remainder <- n %% k

  trt <- integer(0)
  for (b in seq_len(n_full)) {
    trt <- c(trt, sample(0L:(k - 1L)))
  }

  if (remainder > 0) {
    partial <- sample(0L:(k - 1L))
    trt <- c(trt, partial[seq_len(remainder)])
  }
  trt
}

#' Rerandomization for factorial design
#'
#' Repeatedly draws complete (equal-size) randomizations
#' and accepts the first whose k-arm Mahalanobis balance
#' statistic (see `mahalanobis_balance()`) falls below a
#' threshold. Under complete randomization the statistic is
#' asymptotically chi-squared with `p * (k - 1)` degrees of
#' freedom, so the default threshold
#' `qchisq(accept_prob, p * (k - 1))` gives an acceptance
#' probability of about `accept_prob` (Morgan and Rubin,
#' 2012; Branson, Dasgupta and Rubin, 2016).
#'
#' @param data data.frame.
#' @param covariates Character vector of covariate names.
#' @param k Integer. Number of treatment arms.
#' @param threshold Numeric. Acceptance threshold. If
#'   `NULL` (default), calibrated from `accept_prob`.
#' @param accept_prob Numeric in (0, 1). Target acceptance
#'   probability used when `threshold` is `NULL`.
#' @param max_iter Integer. Maximum randomization attempts.
#'   If none is accepted, the best draw is returned.
#' @return Integer vector of assignments (0 to k-1).
#' @export
randomize_rerand_factorial <- function(data,
                                       covariates,
                                       k,
                                       threshold = NULL,
                                       accept_prob = 0.001,
                                       max_iter = 10000L) {
  n <- nrow(data)
  p <- length(covariates)
  X <- as.matrix(data[, covariates, drop = FALSE])

  if (is.null(threshold)) {
    threshold <- stats::qchisq(accept_prob, df = p * (k - 1))
  }

  S <- stats::cov(X)
  S_inv <- tryCatch(
    solve(S),
    error = function(e) MASS::ginv(S)
  )
  base <- rep(0L:(k - 1L), length.out = n)
  # Whitened, centered covariates: the statistic is then
  # sum_a |colSums(Z[arm a, ])|^2 / n_a, computed per draw
  # without re-deriving S or the grand mean.
  Z <- scale(X, center = TRUE, scale = FALSE) %*%
    t(chol(S_inv))
  n_a <- tabulate(base + 1L, nbins = k)
  slots <- lapply(seq_len(k), function(a) which(base == a - 1L))

  # Draws are scored in batches for speed; the first accepted
  # draw in sequence order is returned, so the accepted
  # allocation has the same distribution as one-at-a-time
  # rejection sampling.
  batch <- 250L
  best_perm <- NULL
  best_M <- Inf
  n_done <- 0L

  while (n_done < max_iter) {
    b <- min(batch, max_iter - n_done)
    perms <- replicate(b, sample.int(n))
    stat <- numeric(b)
    for (a in seq_len(k)) {
      idx <- perms[slots[[a]], , drop = FALSE]
      for (j in seq_len(ncol(Z))) {
        s_aj <- colSums(matrix(Z[idx, j], nrow = n_a[a]))
        stat <- stat + s_aj^2 / n_a[a]
      }
    }
    hit <- which(stat < threshold)
    if (length(hit) > 0) {
      best_perm <- perms[, hit[1]]
      best_M <- stat[hit[1]]
      break
    }
    if (min(stat) < best_M) {
      best_M <- min(stat)
      best_perm <- perms[, which.min(stat)]
    }
    n_done <- n_done + b
  }

  trt <- integer(n)
  trt[best_perm] <- base
  trt
}

#' Compute the k-arm Mahalanobis balance criterion
#'
#' Computes `sum_a n_a (xbar_a - xbar)' S^{-1}
#' (xbar_a - xbar)`, where `S` is the covariance of the
#' individual covariates. Because `Var(xbar_a)` is about
#' `S / n_a`, this is the Morgan-Rubin statistic generalized
#' to k arms: under complete randomization it is
#' asymptotically chi-squared on `p * (k - 1)` df. For
#' `k = 2` it equals the Morgan-Rubin `M` statistic up to a
#' finite-population factor.
#'
#' @param X Numeric matrix (n x p).
#' @param trt Integer vector of treatment assignments.
#' @param k Integer. Number of arms.
#' @param S_inv Optional precomputed inverse covariance.
#' @return Numeric scalar.
#' @keywords internal
mahalanobis_balance <- function(X, trt, k, S_inv = NULL) {
  if (is.null(S_inv)) {
    S <- stats::cov(X)
    S_inv <- tryCatch(
      solve(S),
      error = function(e) MASS::ginv(S)
    )
  }
  grand <- colMeans(X)
  sums <- rowsum(X, trt, reorder = TRUE)
  n_a <- tabulate(match(trt, sort(unique(trt))))
  dev <- sums / n_a - matrix(grand, nrow = nrow(sums),
                             ncol = ncol(X), byrow = TRUE)
  sum(n_a * rowSums((dev %*% S_inv) * dev))
}
