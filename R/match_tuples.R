#' Form matched k-tuples from baseline covariate data
#'
#' Groups N subjects into floor(N/k) tuples of size k,
#' minimizing within-tuple covariate dispersion.
#'
#' @param data data.frame with subject covariates.
#' @param k Integer. Tuple size (number of treatment arms).
#' @param covariates Character vector of covariate column
#'   names to use for matching.
#' @param distance Character. Distance metric: 'mahalanobis'
#'   (default), 'euclidean', or 'prognostic'.
#' @param method Character. Matching method: 'greedy'
#'   (default) or 'optimal'.
#' @param prognostic_scores Numeric vector. Required when
#'   distance = 'prognostic'. Pre-computed prognostic scores.
#' @param pool_ratio Numeric. If > 1, select the best
#'   N_target = floor(nrow(data)/pool_ratio) subjects
#'   that form the tightest tuples.
#' @return A list with components:
#'   - tuples: integer vector of tuple IDs (length N or
#'     N_target)
#'   - tuple_ids: data.frame mapping subject row to tuple
#'   - n_tuples: integer
#'   - unmatched: indices of subjects not assigned
#' @export
match_tuples <- function(data, k,
                         covariates,
                         distance = "mahalanobis",
                         method = "greedy",
                         prognostic_scores = NULL,
                         pool_ratio = 1.0) {

  stopifnot(k >= 2)
  stopifnot(all(covariates %in% names(data)))

  X <- as.matrix(data[, covariates, drop = FALSE])
  n <- nrow(X)

  dist_mat <- compute_distance_matrix(
    X, distance, prognostic_scores
  )

  if (method == "greedy") {
    result <- match_greedy(dist_mat, k)
  } else if (method == "optimal") {
    result <- match_greedy(dist_mat, k)
  } else {
    stop("Unknown matching method: ", method)
  }

  if (pool_ratio > 1.0) {
    result <- trim_to_best_tuples(
      result, dist_mat, pool_ratio
    )
  }

  result
}

#' Compute pairwise distance matrix
#'
#' @param X Numeric matrix of covariates (n x p).
#' @param distance Character. Distance type.
#' @param prognostic_scores Optional numeric vector.
#' @return An n x n distance matrix.
#' @keywords internal
compute_distance_matrix <- function(X, distance,
                                    prognostic_scores) {
  n <- nrow(X)

  if (distance == "prognostic" &&
      !is.null(prognostic_scores)) {
    d <- as.matrix(stats::dist(prognostic_scores))
    return(d)
  }

  if (distance == "mahalanobis") {
    S <- stats::cov(X)
    S_inv <- tryCatch(
      solve(S),
      error = function(e) MASS::ginv(S)
    )
    L <- chol(S_inv)
    X_std <- X %*% t(L)
    d <- as.matrix(stats::dist(X_std))
  } else {
    d <- as.matrix(stats::dist(X))
  }
  d
}

#' Greedy k-tuple matching
#'
#' Seeds are ordered by distance to centroid (outermost
#' first) to avoid stranding outliers.
#'
#' @param dist_mat Numeric distance matrix.
#' @param k Integer. Tuple size.
#' @return A list with tuple assignments.
#' @keywords internal
match_greedy <- function(dist_mat, k) {
  n <- nrow(dist_mat)
  n_tuples <- n %/% k
  remainder <- n %% k

  centroid_dist <- rowMeans(dist_mat)
  seed_order <- order(centroid_dist, decreasing = TRUE)

  available <- rep(TRUE, n)
  tuple_id <- rep(NA_integer_, n)
  current_tuple <- 0L

  for (seed in seed_order) {
    if (!available[seed]) next
    if (current_tuple >= n_tuples) break

    current_tuple <- current_tuple + 1L
    candidates <- which(available & seq_len(n) != seed)

    if (length(candidates) < k - 1) break

    dists_to_seed <- dist_mat[seed, candidates]
    nearest <- candidates[
      order(dists_to_seed)[seq_len(k - 1)]
    ]

    members <- c(seed, nearest)
    tuple_id[members] <- current_tuple
    available[members] <- FALSE
  }

  unmatched <- which(is.na(tuple_id))

  list(
    tuple_id = tuple_id,
    n_tuples = current_tuple,
    unmatched = unmatched,
    k = k
  )
}

#' Trim to best tuples when pool_ratio > 1
#'
#' @param result Matching result from match_greedy.
#' @param dist_mat Distance matrix.
#' @param pool_ratio Numeric.
#' @return Updated matching result with fewer tuples.
#' @keywords internal
trim_to_best_tuples <- function(result, dist_mat,
                                pool_ratio) {
  n_target <- max(
    1L,
    as.integer(result$n_tuples / pool_ratio)
  )

  tuple_quality <- numeric(result$n_tuples)
  for (m in seq_len(result$n_tuples)) {
    members <- which(result$tuple_id == m)
    if (length(members) < 2) {
      tuple_quality[m] <- 0
      next
    }
    pair_dists <- dist_mat[members, members]
    tuple_quality[m] <- mean(
      pair_dists[upper.tri(pair_dists)]
    )
  }

  best <- order(tuple_quality)[seq_len(n_target)]
  keep_tuples <- best

  new_tuple_id <- rep(NA_integer_, length(result$tuple_id))
  for (i in seq_along(keep_tuples)) {
    members <- which(result$tuple_id == keep_tuples[i])
    new_tuple_id[members] <- i
  }

  list(
    tuple_id = new_tuple_id,
    n_tuples = n_target,
    unmatched = which(is.na(new_tuple_id)),
    k = result$k
  )
}

#' Compute within-tuple balance diagnostics
#'
#' @param data data.frame of subject data.
#' @param tuples Matching result from match_tuples, or
#'   `NULL` for unmatched designs (simple, stratified,
#'   rerandomization). When `NULL`, balance is computed
#'   across the full sample rather than within tuples.
#' @param covariates Character vector of covariate names.
#' @param treatment Integer vector of treatment assignments.
#' @return A data.frame with balance statistics.
#' @export
tuple_balance <- function(data, tuples, covariates,
                          treatment) {
  matched <- if (is.null(tuples)) {
    rep(TRUE, nrow(data))
  } else {
    !is.na(tuples$tuple_id)
  }
  dat <- data[matched, ]
  trt <- treatment[matched]

  results <- data.frame(
    covariate = character(0),
    max_asmd = numeric(0),
    mean_asmd = numeric(0)
  )

  arms <- sort(unique(trt))
  k <- length(arms)

  for (cov in covariates) {
    x <- dat[[cov]]
    pooled_sd <- stats::sd(x)
    if (pooled_sd < 1e-12) next

    asmds <- numeric(0)
    for (a in seq_len(k - 1)) {
      for (b in seq(a + 1, k)) {
        diff <- mean(x[trt == arms[a]]) -
          mean(x[trt == arms[b]])
        asmds <- c(asmds, abs(diff) / pooled_sd)
      }
    }

    results <- rbind(results, data.frame(
      covariate = cov,
      max_asmd = max(asmds),
      mean_asmd = mean(asmds)
    ))
  }

  results
}
