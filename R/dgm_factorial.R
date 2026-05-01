#' Generate a simulated factorial trial dataset
#'
#' Creates baseline covariates, assigns treatment via
#' the specified randomization method, and generates
#' outcomes.
#'
#' @param n Integer. Total sample size.
#' @param n_factors Integer. Number of factors.
#' @param p_covariates Integer. Number of baseline
#'   covariates.
#' @param r_squared Numeric. Target prognostic R-squared.
#' @param tau_main Numeric vector of length n_factors.
#'   Main effects in SD units.
#' @param tau_interaction Numeric. Two-way interaction
#'   effect in SD units (for the first pair of factors).
#' @param sigma_y Numeric. Outcome SD (default 1).
#' @param cov_rho Numeric. AR(1) correlation among
#'   covariates (default 0.3).
#' @param randomization Character. One of 'simple',
#'   'stratified', 'rerandomization', 'matched'.
#' @param pool_ratio Numeric. Pool/sample ratio for
#'   matched designs (default 1.0).
#' @param rerand_threshold Numeric. Threshold for
#'   rerandomization (default NULL = use chi-sq 0.001).
#' @return A list with: data (data.frame), tuples
#'   (matching result or NULL), treatment (integer vector).
#' @export
generate_factorial_trial <- function(
    n,
    n_factors = 2L,
    p_covariates = 3L,
    r_squared = 0.3,
    tau_main = c(0.2, 0.2),
    tau_interaction = 0,
    sigma_y = 1,
    cov_rho = 0.3,
    randomization = "simple",
    pool_ratio = 1.0,
    rerand_threshold = NULL) {

  k <- 2^n_factors
  n_pool <- ceiling(n * pool_ratio)

  Sigma_x <- outer(
    seq_len(p_covariates),
    seq_len(p_covariates),
    function(i, j) cov_rho^abs(i - j)
  )
  X <- MASS::mvrnorm(n_pool, mu = rep(0, p_covariates),
                     Sigma = Sigma_x)
  cov_names <- paste0("X", seq_len(p_covariates))
  colnames(X) <- cov_names
  data <- as.data.frame(X)
  data$id <- seq_len(n_pool)

  labels <- factorial_labels(n_factors)
  factor_cols <- paste0("F", seq_len(n_factors))

  tuples_result <- NULL

  if (randomization == "matched") {
    tuples_result <- match_tuples(
      data, k = k,
      covariates = cov_names,
      distance = "mahalanobis",
      method = "greedy",
      pool_ratio = pool_ratio
    )
    trt <- randomize_tuples(tuples_result)
    matched <- !is.na(trt)
    data <- data[matched, ]
    X <- X[matched, , drop = FALSE]
    trt <- trt[matched]
    if (!is.null(tuples_result)) {
      tuples_result$tuple_id <-
        tuples_result$tuple_id[matched]
    }

  } else if (randomization == "stratified") {
    data <- data[seq_len(n), ]
    X <- X[seq_len(n), , drop = FALSE]
    data$X1_cat <- as.integer(
      cut(data$X1, breaks = 4)
    )
    trt <- randomize_stratified_factorial(
      data, strat_vars = "X1_cat", k = k
    )
    data$X1_cat <- NULL

  } else if (randomization == "rerandomization") {
    data <- data[seq_len(n), ]
    X <- X[seq_len(n), , drop = FALSE]
    trt <- randomize_rerand_factorial(
      data, covariates = cov_names, k = k,
      threshold = rerand_threshold
    )

  } else {
    data <- data[seq_len(n), ]
    X <- X[seq_len(n), , drop = FALSE]
    trt <- randomize_simple_factorial(n, k)
  }

  n_actual <- nrow(data)
  trt_factors <- labels[trt + 1, factor_cols,
                        drop = FALSE]

  beta <- compute_beta_for_r2(p_covariates, r_squared,
                              sigma_y)
  y <- as.numeric(X %*% beta)

  for (j in seq_len(min(n_factors, length(tau_main)))) {
    coded_j <- 2 * trt_factors[[paste0("F", j)]] - 1
    y <- y + tau_main[j] * sigma_y * coded_j / 2
  }

  if (n_factors >= 2 && tau_interaction != 0) {
    coded_1 <- 2 * trt_factors$F1 - 1
    coded_2 <- 2 * trt_factors$F2 - 1
    y <- y + tau_interaction * sigma_y *
      coded_1 * coded_2 / 4
  }

  y <- y + stats::rnorm(n_actual, 0,
                        sigma_y * sqrt(1 - r_squared))

  data$y <- y
  data$trt <- trt

  for (f in factor_cols) {
    data[[f]] <- trt_factors[[f]]
  }

  list(
    data = data,
    tuples = tuples_result,
    treatment = trt
  )
}

#' Compute beta coefficients to achieve target R-squared
#'
#' @param p Integer. Number of covariates.
#' @param r_squared Numeric. Target R-squared.
#' @param sigma_y Numeric. Outcome SD.
#' @return Numeric vector of beta coefficients.
#' @keywords internal
compute_beta_for_r2 <- function(p, r_squared, sigma_y) {
  if (r_squared <= 0) return(rep(0, p))
  raw <- rep(1 / sqrt(p), p)
  raw * sigma_y * sqrt(r_squared)
}
