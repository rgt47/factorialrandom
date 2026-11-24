library(tinytest)
library(factorialrandom)

## ── factorial_labels() / factorial_contrasts() ──────────────

labels <- factorial_labels(2L)
expect_equal(nrow(labels), 4L,
  info = "2^2 design has 4 treatment combinations")
expect_equal(labels$label, 0:3,
  info = "labels run 0 to k-1")

contrasts <- factorial_contrasts(2L)
expect_equal(nrow(contrasts), 3L,
  info = "2 main effects + 1 interaction for n_factors = 2")
expect_true(all(abs(rowSums(as.matrix(
  contrasts[, -1]
))) < 1e-10), info = "each contrast's weights sum to zero")

## ── match_tuples() ───────────────────────────────────────────

set.seed(1)
dat <- data.frame(
  X1 = rnorm(20), X2 = rnorm(20)
)
res <- match_tuples(dat, k = 4L, covariates = c("X1", "X2"))
expect_equal(res$n_tuples, 5L,
  info = "20 subjects / k=4 forms 5 complete tuples")
expect_equal(length(res$unmatched), 0L,
  info = "no remainder when n is a multiple of k")
expect_true(all(table(res$tuple_id) == 4L),
  info = "every tuple has exactly k members")

## unmatched remainder when n is not a multiple of k
dat_rem <- data.frame(X1 = rnorm(22), X2 = rnorm(22))
res_rem <- match_tuples(dat_rem, k = 4L, covariates = c("X1", "X2"))
expect_equal(res_rem$n_tuples, 5L,
  info = "22 subjects / k=4 forms 5 tuples with 2 unmatched")
expect_equal(length(res_rem$unmatched), 2L,
  info = "remainder subjects are reported as unmatched")

## pool_ratio > 1 trims to the tightest tuples
res_pool <- match_tuples(
  dat_rem, k = 4L, covariates = c("X1", "X2"), pool_ratio = 2.0
)
expect_true(res_pool$n_tuples <= res_rem$n_tuples,
  info = "pool_ratio > 1 keeps no more tuples than pool_ratio = 1")

## ── randomize_tuples() ───────────────────────────────────────

set.seed(2)
trt <- randomize_tuples(res)
expect_equal(length(trt), length(res$tuple_id),
  info = "treatment vector length matches tuple_id length")
tab <- table(trt[!is.na(trt)], res$tuple_id[!is.na(trt)])
expect_true(all(colSums(tab > 0) == 4L),
  info = "each tuple receives all k distinct treatment labels")

## ── tuple_balance(): matched and tuples = NULL paths ─────────

set.seed(3)
n <- 20
dat_bal <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
tup <- match_tuples(dat_bal, k = 4L, covariates = c("X1", "X2"))
trt_bal <- randomize_tuples(tup)

bal <- tuple_balance(dat_bal, tup, c("X1", "X2"), trt_bal)
expect_true(is.data.frame(bal),
  info = "tuple_balance returns a data.frame for matched designs")
expect_true(all(bal$max_asmd >= bal$mean_asmd),
  info = "max ASMD is never less than mean ASMD per covariate")

## tuples = NULL must not crash (regression test for the
## tuple_balance() NULL-tuples bug found in the 2026-08-16
## pub_review: previously threw "missing value where TRUE/FALSE
## needed" for every non-matched randomization method)
trt_simple <- sample(0:3, n, replace = TRUE)
expect_silent(
  tuple_balance(dat_bal, NULL, c("X1", "X2"), trt_simple)
)
bal_null <- tuple_balance(dat_bal, NULL, c("X1", "X2"), trt_simple)
expect_true(is.data.frame(bal_null),
  info = "tuple_balance(tuples = NULL) returns a data.frame")
expect_true(nrow(bal_null) > 0,
  info = "tuple_balance(tuples = NULL) computes balance across the full sample")

## degenerate covariance (a constant covariate) is skipped, not
## an error
dat_degenerate <- dat_bal
dat_degenerate$X1 <- 5
expect_silent(
  tuple_balance(dat_degenerate, NULL, c("X1", "X2"), trt_simple)
)
bal_degenerate <- tuple_balance(
  dat_degenerate, NULL, c("X1", "X2"), trt_simple
)
expect_false("X1" %in% bal_degenerate$covariate,
  info = "zero-variance covariate is excluded rather than producing NA/Inf")

## ── estimate_effects(): unbiasedness and regression adjustment ──

set.seed(4)
n_reps <- 300L
ests_f1 <- numeric(n_reps)
ests_int <- numeric(n_reps)
for (i in seq_len(n_reps)) {
  trial <- generate_factorial_trial(
    n = 200, n_factors = 2L, p_covariates = 3L,
    r_squared = 0.3, tau_main = c(0.2, 0.2),
    tau_interaction = 0.15, sigma_y = 1, cov_rho = 0.3,
    randomization = "simple", pool_ratio = 1.0
  )
  eff <- estimate_effects(
    trial$data$y, trial$treatment, trial$tuples, 2L,
    method = "neyman"
  )
  ests_f1[i] <- eff$estimate[eff$effect == "F1"]
  ests_int[i] <- eff$estimate[eff$effect == "F1:F2"]
}
## Regression test for the factorial_contrasts()/dgm_factorial()
## normalization mismatch found in the 2026-08-16 pub_review: the
## estimator previously recovered only ~50% of the true main
## effect and ~25% of the true interaction effect on average.
expect_true(
  abs(mean(ests_f1) - 0.2) < 0.03,
  info = "main-effect estimator is unbiased for tau_main = 0.2 (300 reps, n = 200)"
)
expect_true(
  abs(mean(ests_int) - 0.15) < 0.03,
  info = "interaction estimator is unbiased for tau_interaction = 0.15 (300 reps, n = 200)"
)

## Regression test for the regression_adjustment() zero-correction
## bug found in the 2026-08-16 pub_review: adjusted and unadjusted
## point estimates were previously identical by construction.
set.seed(5)
trial_adj <- generate_factorial_trial(
  n = 200, n_factors = 2L, p_covariates = 3L,
  r_squared = 0.5, tau_main = c(0.2, 0.2),
  tau_interaction = 0.15, sigma_y = 1, cov_rho = 0.3,
  randomization = "matched", pool_ratio = 1.0
)
cov_mat <- as.matrix(
  trial_adj$data[, paste0("X", 1:3), drop = FALSE]
)
eff_unadj <- estimate_effects(
  trial_adj$data$y, trial_adj$treatment, trial_adj$tuples, 2L,
  method = "neyman"
)
eff_adj <- estimate_effects(
  trial_adj$data$y, trial_adj$treatment, trial_adj$tuples, 2L,
  covariates = cov_mat, method = "neyman"
)
expect_false(
  identical(eff_unadj$estimate, eff_adj$estimate),
  info = "regression-adjusted point estimates differ from unadjusted (correction is not identically zero)"
)

## ── generate_factorial_trial(): tuples = NULL for non-matched ──

trial_simple <- generate_factorial_trial(
  n = 40, n_factors = 2L, p_covariates = 2L,
  randomization = "simple"
)
expect_true(is.null(trial_simple$tuples),
  info = "simple randomization returns tuples = NULL")

trial_matched <- generate_factorial_trial(
  n = 40, n_factors = 2L, p_covariates = 2L,
  randomization = "matched", pool_ratio = 1.0
)
expect_false(is.null(trial_matched$tuples),
  info = "matched randomization returns a non-NULL tuples object")
expect_equal(nrow(trial_matched$data), length(trial_matched$treatment),
  info = "data and treatment vector stay aligned after matching")

## ── 2026-11 referee remediation (M4, M5, M7, DGM R^2) ────────

## M4: the k-arm balance statistic is calibrated to chi^2 on
## p (k - 1) df, so the nominal acceptance rate is achieved. The
## earlier pairwise statistic on the individual-level covariance
## accepted nearly every draw at a nominal rate of 0.001.
set.seed(11)
acc_rate <- function(n, p, k = 4L, pa = 0.05, draws = 1500L) {
  X <- matrix(rnorm(n * p), n, p)
  base <- rep(0L:(k - 1L), length.out = n)
  thr <- qchisq(pa, p * (k - 1))
  mean(replicate(draws,
    factorialrandom:::mahalanobis_balance(X, sample(base), k) < thr))
}
for (cell in list(c(100, 3), c(500, 8))) {
  r <- acc_rate(cell[1], cell[2])
  expect_true(r > 0.025 && r < 0.08,
    info = sprintf("acceptance rate near 0.05 at n=%d p=%d (got %.3f)",
                   cell[1], cell[2], r))
}

set.seed(12)
dat_rr <- data.frame(X1 = rnorm(100), X2 = rnorm(100), X3 = rnorm(100))
trt_rr <- randomize_rerand_factorial(dat_rr, c("X1", "X2", "X3"),
                                     k = 4L, accept_prob = 0.01)
expect_equal(as.vector(table(trt_rr)), rep(25L, 4),
  info = "rerandomization returns a complete (equal-size) allocation")
stat_rr <- factorialrandom:::mahalanobis_balance(
  as.matrix(dat_rr), trt_rr, 4L)
expect_true(stat_rr < qchisq(0.01, 9),
  info = "accepted rerandomization meets the calibrated threshold")

## M5: stratified randomization uses sample quartiles of X1
set.seed(13)
tr_st <- generate_factorial_trial(n = 100, p_covariates = 3L,
                                  randomization = "stratified")
q_x1 <- cut(tr_st$data$X1, quantile(tr_st$data$X1, 0:4 / 4),
            include.lowest = TRUE)
expect_equal(as.vector(table(q_x1)), rep(25L, 4),
  info = "X1 strata are equal-count quartiles")
spread <- tapply(tr_st$treatment, q_x1,
                 function(z) diff(range(tabulate(z + 1L, 4))))
expect_true(all(spread <= 1),
  info = "permuted blocks keep arms within 1 of each other per stratum")

## DGM: beta is scaled so that Var(X'beta) = R^2 under AR(1)
Sig <- outer(1:3, 1:3, function(i, j) 0.3^abs(i - j))
b <- factorialrandom:::compute_beta_for_r2(3, 0.5, 1, Sigma_x = Sig)
expect_equal(as.numeric(t(b) %*% Sig %*% b), 0.5,
  info = "linear predictor variance equals target R^2")

## M7: paired intervals use t on M - 1 df; adjusted SE is HC2
set.seed(14)
tr_p <- generate_factorial_trial(n = 100, p_covariates = 3L,
                                 randomization = "matched")
ep <- estimate_effects(tr_p$data$y, tr_p$treatment, tr_p$tuples,
                       2L, method = "paired")
expect_equal((ep$ci_upper - ep$estimate) / ep$se,
             rep(qt(0.975, 24), 3),
  info = "paired CI half-width uses t_{M-1} with M = 25")
tr_s <- generate_factorial_trial(n = 100, p_covariates = 3L,
                                 randomization = "simple")
ea <- estimate_effects(tr_s$data$y, tr_s$treatment, NULL, 2L,
  covariates = as.matrix(tr_s$data[, c("X1", "X2", "X3")]))
eu <- estimate_effects(tr_s$data$y, tr_s$treatment, NULL, 2L)
expect_true(all(is.finite(ea$se)) && !identical(ea$estimate,
                                                 eu$estimate),
  info = "Lin adjustment also applies to unmatched designs")
