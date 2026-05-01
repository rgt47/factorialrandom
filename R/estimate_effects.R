#' Estimate factorial effects under matched-tuple design
#'
#' Computes main effects and interactions for a 2^F
#' factorial design, with variance estimation appropriate
#' for the matched-tuple randomization distribution.
#'
#' @param outcome Numeric vector of outcomes.
#' @param treatment Integer vector of treatment assignments
#'   (0 to k-1).
#' @param tuples Matching result from match_tuples (NULL
#'   for unmatched designs).
#' @param n_factors Integer. Number of factors.
#' @param covariates Optional matrix of covariates for
#'   regression adjustment.
#' @param method Character. 'neyman' (default) for
#'   Neyman-type variance, 'paired' for within-tuple.
#' @return A data.frame with columns: effect, estimate,
#'   se, ci_lower, ci_upper, p_value.
#' @export
estimate_effects <- function(outcome,
                             treatment,
                             tuples = NULL,
                             n_factors = 2L,
                             covariates = NULL,
                             method = "neyman") {

  k <- 2^n_factors
  labels <- factorial_labels(n_factors)
  contrasts <- factorial_contrasts(n_factors)

  results <- data.frame(
    effect = character(0),
    estimate = numeric(0),
    se = numeric(0),
    ci_lower = numeric(0),
    ci_upper = numeric(0),
    p_value = numeric(0)
  )

  for (i in seq_len(nrow(contrasts))) {
    cname <- contrasts$name[i]
    weights <- as.numeric(contrasts[i, -1])

    arm_means <- numeric(k)
    arm_vars <- numeric(k)
    arm_ns <- integer(k)

    matched <- if (!is.null(tuples)) {
      !is.na(tuples$tuple_id)
    } else {
      rep(TRUE, length(outcome))
    }

    for (a in seq_len(k)) {
      idx <- which(treatment == (a - 1) & matched)
      arm_ns[a] <- length(idx)
      if (arm_ns[a] > 0) {
        arm_means[a] <- mean(outcome[idx])
        arm_vars[a] <- if (arm_ns[a] > 1) {
          stats::var(outcome[idx])
        } else {
          0
        }
      }
    }

    est <- sum(weights * arm_means)

    if (!is.null(tuples) && method == "paired") {
      se <- paired_se(outcome, treatment, tuples,
                      weights, k)
    } else {
      se <- sqrt(sum(weights^2 * arm_vars / arm_ns))
    }

    if (!is.null(covariates) && !is.null(tuples)) {
      adj <- regression_adjustment(
        outcome[matched], treatment[matched],
        covariates[matched, , drop = FALSE],
        weights, k
      )
      est <- est - adj$correction
      se <- adj$se
    }

    p_val <- if (se > 0) {
      2 * stats::pnorm(-abs(est / se))
    } else {
      NA_real_
    }

    results <- rbind(results, data.frame(
      effect = cname,
      estimate = est,
      se = se,
      ci_lower = est - 1.96 * se,
      ci_upper = est + 1.96 * se,
      p_value = p_val
    ))
  }

  results
}

#' Generate factorial contrast definitions
#'
#' @param n_factors Integer. Number of factors.
#' @return A data.frame with contrast name and weights
#'   for each treatment combination.
#' @export
factorial_contrasts <- function(n_factors) {
  labels <- factorial_labels(n_factors)
  k <- nrow(labels)
  factor_cols <- paste0("F", seq_len(n_factors))

  coded <- labels[, factor_cols, drop = FALSE]
  coded[coded == 0] <- -1

  contrasts <- data.frame(name = character(0))
  for (j in seq_len(n_factors)) {
    contrasts <- rbind(
      contrasts,
      data.frame(name = paste0("F", j))
    )
  }

  if (n_factors >= 2) {
    for (a in seq_len(n_factors - 1)) {
      for (b in seq(a + 1, n_factors)) {
        contrasts <- rbind(
          contrasts,
          data.frame(
            name = paste0("F", a, ":F", b)
          )
        )
      }
    }
  }

  weight_mat <- matrix(NA, nrow = nrow(contrasts),
                        ncol = k)

  row_idx <- 0
  for (j in seq_len(n_factors)) {
    row_idx <- row_idx + 1
    weight_mat[row_idx, ] <- coded[[paste0("F", j)]] / k
  }

  if (n_factors >= 2) {
    for (a in seq_len(n_factors - 1)) {
      for (b in seq(a + 1, n_factors)) {
        row_idx <- row_idx + 1
        weight_mat[row_idx, ] <-
          coded[[paste0("F", a)]] *
          coded[[paste0("F", b)]] / k
      }
    }
  }

  result <- cbind(contrasts,
                  as.data.frame(weight_mat))
  names(result)[-1] <- paste0("arm", 0:(k - 1))
  result
}

#' Paired (within-tuple) SE for factorial contrasts
#'
#' @param outcome Numeric vector.
#' @param treatment Integer vector.
#' @param tuples Matching result.
#' @param weights Numeric vector of contrast weights.
#' @param k Integer. Number of arms.
#' @return Numeric. Standard error.
#' @keywords internal
paired_se <- function(outcome, treatment, tuples,
                      weights, k) {
  M <- tuples$n_tuples
  tuple_contrasts <- numeric(M)

  for (m in seq_len(M)) {
    members <- which(tuples$tuple_id == m)
    if (length(members) != k) {
      tuple_contrasts[m] <- NA
      next
    }
    trt_m <- treatment[members]
    y_m <- outcome[members]

    arm_vals <- numeric(k)
    for (a in seq_len(k)) {
      idx <- which(trt_m == (a - 1))
      arm_vals[a] <- if (length(idx) > 0) {
        y_m[idx[1]]
      } else {
        NA
      }
    }
    tuple_contrasts[m] <- sum(weights * arm_vals)
  }

  valid <- !is.na(tuple_contrasts)
  M_valid <- sum(valid)

  if (M_valid < 2) return(NA_real_)

  tc <- tuple_contrasts[valid]
  sqrt(stats::var(tc) / M_valid)
}

#' Lin-style regression adjustment
#'
#' @param outcome Numeric vector (matched only).
#' @param treatment Integer vector (matched only).
#' @param covariates Matrix of covariates.
#' @param weights Contrast weights.
#' @param k Number of arms.
#' @return List with correction and adjusted se.
#' @keywords internal
regression_adjustment <- function(outcome, treatment,
                                  covariates, weights,
                                  k) {
  X <- as.matrix(covariates)
  X_centered <- scale(X, center = TRUE, scale = FALSE)

  trt_dummies <- matrix(0, nrow = length(treatment),
                        ncol = k)
  for (a in seq_len(k)) {
    trt_dummies[, a] <- as.numeric(
      treatment == (a - 1)
    )
  }

  interactions <- do.call(cbind, lapply(
    seq_len(k), function(a) {
      trt_dummies[, a] * X_centered
    }
  ))

  design <- cbind(trt_dummies, interactions)

  fit <- tryCatch(
    stats::lm.fit(design, outcome),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(correction = 0, se = NA_real_))
  }

  arm_coefs <- fit$coefficients[seq_len(k)]
  arm_coefs[is.na(arm_coefs)] <- 0
  est_adj <- sum(weights * arm_coefs)

  resid <- fit$residuals
  arm_resid_vars <- numeric(k)
  arm_ns <- integer(k)
  for (a in seq_len(k)) {
    idx <- which(treatment == (a - 1))
    arm_ns[a] <- length(idx)
    arm_resid_vars[a] <- if (length(idx) > 1) {
      stats::var(resid[idx])
    } else {
      0
    }
  }

  se_adj <- sqrt(sum(weights^2 * arm_resid_vars / arm_ns))

  list(
    correction = est_adj - sum(weights * arm_coefs),
    se = se_adj
  )
}
