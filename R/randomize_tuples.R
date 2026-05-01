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
#' Repeatedly generate simple randomizations and accept
#' only those where the Mahalanobis distance between arm
#' covariate means is below a threshold.
#'
#' @param data data.frame.
#' @param covariates Character vector of covariate names.
#' @param k Integer. Number of treatment arms.
#' @param threshold Numeric. Acceptance threshold for the
#'   Mahalanobis criterion.
#' @param max_iter Integer. Maximum randomization attempts.
#' @return Integer vector of assignments (0 to k-1).
#' @export
randomize_rerand_factorial <- function(data,
                                       covariates,
                                       k,
                                       threshold = NULL,
                                       max_iter = 10000L) {
  n <- nrow(data)
  p <- length(covariates)
  X <- as.matrix(data[, covariates, drop = FALSE])

  if (is.null(threshold)) {
    threshold <- stats::qchisq(0.001, df = p * (k - 1))
  }

  best_trt <- NULL
  best_M <- Inf

  for (iter in seq_len(max_iter)) {
    trt <- sample(rep(0L:(k - 1L), length.out = n))
    M_stat <- mahalanobis_balance(X, trt, k)

    if (M_stat < threshold) {
      return(trt)
    }

    if (M_stat < best_M) {
      best_M <- M_stat
      best_trt <- trt
    }
  }

  best_trt
}

#' Compute Mahalanobis balance criterion for k arms
#'
#' @param X Numeric matrix (n x p).
#' @param trt Integer vector of treatment assignments.
#' @param k Integer. Number of arms.
#' @return Numeric. Sum of pairwise Mahalanobis distances
#'   between arm means.
#' @keywords internal
mahalanobis_balance <- function(X, trt, k) {
  S <- stats::cov(X)
  S_inv <- tryCatch(
    solve(S),
    error = function(e) MASS::ginv(S)
  )

  arm_means <- lapply(0:(k - 1), function(a) {
    colMeans(X[trt == a, , drop = FALSE])
  })

  total <- 0
  for (a in seq_len(k - 1)) {
    for (b in seq(a + 1, k)) {
      diff <- arm_means[[a]] - arm_means[[b]]
      total <- total +
        as.numeric(t(diff) %*% S_inv %*% diff)
    }
  }
  total
}
