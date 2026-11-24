## ──────────────────────────────────────────────────────────
## 01_run_simulation.R
## Monte Carlo simulation comparing randomization strategies
## for 2x2 factorial designs
## ──────────────────────────────────────────────────────────

## Usage (from the package root):
##   Rscript analysis/scripts/01_run_simulation.R [n_reps] [n_cores]
## Defaults: n_reps = 2000, n_cores = 1. The package is loaded
## from source with pkgload so that the run uses the R/ code in
## this repository, not a possibly stale installed copy.

if (file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(factorialrandom)
}
library(MASS)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

# Morris, White, and Crowther (2019) sec. 5.3: the MCSE of a
# rejection rate p is sqrt(p (1 - p) / n_reps). At p = 0.5 an
# MCSE <= 1 percentage point needs n_reps >= 2500; 2000 reps give
# 1.1 points. At p = 0.05 (nominal type I error of a two-sided
# 0.05 test) 2000 reps give an MCSE of 0.49 points.
n_replicates <- if (length(args) >= 1) as.integer(args[1]) else 2000L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
base_seed <- 20260312L
n_factors <- 2L
k <- 2^n_factors
sigma_y <- 1

scenarios <- expand.grid(
  n = c(100L, 200L, 500L),
  p_covariates = c(3L, 8L),
  r_squared = c(0.1, 0.3, 0.5),
  tau_main_a = c(0, 0.2),
  tau_interaction = c(0, 0.15),
  stringsAsFactors = FALSE
)

randomization_methods <- c(
  "simple",
  "stratified",
  "rerandomization",
  "matched_1.0",
  "matched_2.0"
)

parse_method <- function(method_label) {
  if (grepl("^matched_", method_label)) {
    pool_ratio <- as.numeric(sub("matched_", "", method_label))
    list(randomization = "matched", pool_ratio = pool_ratio)
  } else {
    list(randomization = method_label, pool_ratio = 1.0)
  }
}

run_one_replicate <- function(scenario, method_label) {
  mp <- parse_method(method_label)

  tau_main <- c(scenario$tau_main_a, scenario$tau_main_a)

  trial <- generate_factorial_trial(
    n = scenario$n,
    n_factors = n_factors,
    p_covariates = scenario$p_covariates,
    r_squared = scenario$r_squared,
    tau_main = tau_main,
    tau_interaction = scenario$tau_interaction,
    sigma_y = sigma_y,
    cov_rho = 0.3,
    randomization = mp$randomization,
    pool_ratio = mp$pool_ratio
  )

  dat <- trial$data
  cov_names <- paste0("X", seq_len(scenario$p_covariates))

  bal <- tuple_balance(
    dat, trial$tuples, cov_names, trial$treatment
  )
  max_asmd <- if (nrow(bal) > 0) max(bal$max_asmd) else NA_real_

  eff <- estimate_effects(
    outcome = dat$y,
    treatment = trial$treatment,
    tuples = trial$tuples,
    n_factors = n_factors,
    method = "neyman"
  )

  eff_paired <- if (!is.null(trial$tuples)) {
    estimate_effects(
      outcome = dat$y,
      treatment = trial$treatment,
      tuples = trial$tuples,
      n_factors = n_factors,
      method = "paired"
    )
  } else {
    NULL
  }

  # Lin regression adjustment for every design, so that the
  # adjusted matched designs can be compared with Lin-adjusted
  # unmatched designs.
  eff_adj <- {
    cov_mat <- as.matrix(dat[, cov_names, drop = FALSE])
    estimate_effects(
      outcome = dat$y,
      treatment = trial$treatment,
      tuples = trial$tuples,
      n_factors = n_factors,
      covariates = cov_mat,
      method = "neyman"
    )
  }

  results <- data.frame(
    effect = eff$effect,
    estimate = eff$estimate,
    se = eff$se,
    ci_lower = eff$ci_lower,
    ci_upper = eff$ci_upper,
    p_value = eff$p_value,
    estimator = "neyman",
    max_asmd = max_asmd,
    stringsAsFactors = FALSE
  )

  if (!is.null(eff_paired)) {
    paired_rows <- data.frame(
      effect = eff_paired$effect,
      estimate = eff_paired$estimate,
      se = eff_paired$se,
      ci_lower = eff_paired$ci_lower,
      ci_upper = eff_paired$ci_upper,
      p_value = eff_paired$p_value,
      estimator = "paired",
      max_asmd = max_asmd,
      stringsAsFactors = FALSE
    )
    results <- rbind(results, paired_rows)
  }

  if (!is.null(eff_adj)) {
    adj_rows <- data.frame(
      effect = eff_adj$effect,
      estimate = eff_adj$estimate,
      se = eff_adj$se,
      ci_lower = eff_adj$ci_lower,
      ci_upper = eff_adj$ci_upper,
      p_value = eff_adj$p_value,
      estimator = "adjusted",
      max_asmd = max_asmd,
      stringsAsFactors = FALSE
    )
    results <- rbind(results, adj_rows)
  }

  results
}

run_condition <- function(scenario_idx, method_label,
                          n_reps = n_replicates) {
  scenario <- scenarios[scenario_idx, ]
  cat(sprintf(
    "  Scenario %d/%d | Method: %s | n=%d p=%d R2=%.1f tau=%.1f int=%.2f\n",
    scenario_idx, nrow(scenarios), method_label,
    scenario$n, scenario$p_covariates, scenario$r_squared,
    scenario$tau_main_a, scenario$tau_interaction
  ))

  rep_results <- vector("list", n_reps)
  n_errors <- 0L
  first_error <- NULL
  for (r in seq_len(n_reps)) {
    rep_results[[r]] <- tryCatch(
      run_one_replicate(scenario, method_label),
      error = function(e) {
        n_errors <<- n_errors + 1L
        if (is.null(first_error)) first_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  if (n_errors > 0L) {
    cat(sprintf(
      "  WARNING: %d/%d replicates failed for %s (first error: %s)\n",
      n_errors, n_reps, method_label, first_error
    ))
  }

  valid <- !vapply(rep_results, is.null, logical(1))
  if (sum(valid) == 0) {
    warning(sprintf(
      "run_condition: all %d replicates failed for scenario %d, method %s",
      n_reps, scenario_idx, method_label
    ))
    return(NULL)
  }

  all_reps <- do.call(rbind, rep_results[valid])
  all_reps$replicate <- rep(
    which(valid),
    each = nrow(rep_results[valid][[1]])
  )
  all_reps$scenario_idx <- scenario_idx
  all_reps$method <- method_label
  all_reps$n <- scenario$n
  all_reps$p_covariates <- scenario$p_covariates
  all_reps$r_squared <- scenario$r_squared
  all_reps$tau_main_a <- scenario$tau_main_a
  all_reps$tau_interaction <- scenario$tau_interaction

  all_reps
}

compute_performance <- function(results) {
  true_effects <- unique(
    results[, c("scenario_idx", "effect", "tau_main_a",
                "tau_interaction")]
  )
  true_effects$true_tau <- ifelse(
    true_effects$effect %in% c("F1", "F2"),
    true_effects$tau_main_a * sigma_y,
    true_effects$tau_interaction * sigma_y
  )

  results <- merge(results, true_effects[, c(
    "scenario_idx", "effect", "true_tau"
  )])

  perf <- results |>
    split(interaction(
      results$scenario_idx, results$method,
      results$effect, results$estimator
    )) |>
    lapply(function(df) {
      if (nrow(df) < 10) return(NULL)
      B <- nrow(df)

      bias <- mean(df$estimate) - df$true_tau[1]
      emp_se <- sd(df$estimate)
      model_se <- mean(df$se, na.rm = TRUE)
      rel_se_error <- (model_se - emp_se) / emp_se

      covers <- df$ci_lower <= df$true_tau[1] &
        df$ci_upper >= df$true_tau[1]
      coverage <- mean(covers, na.rm = TRUE)

      rejects <- df$p_value < 0.05
      rejection_rate <- mean(rejects, na.rm = TRUE)

      mean_asmd <- mean(df$max_asmd, na.rm = TRUE)

      mcse_bias <- emp_se / sqrt(B)
      mcse_coverage <- sqrt(
        coverage * (1 - coverage) / B
      )
      mcse_rejection <- sqrt(
        rejection_rate * (1 - rejection_rate) / B
      )

      data.frame(
        scenario_idx = df$scenario_idx[1],
        method = df$method[1],
        effect = df$effect[1],
        estimator = df$estimator[1],
        n = df$n[1],
        p_covariates = df$p_covariates[1],
        r_squared = df$r_squared[1],
        tau_main_a = df$tau_main_a[1],
        tau_interaction_dgm = df$tau_interaction[1],
        true_tau = df$true_tau[1],
        n_reps = B,
        bias = bias,
        emp_se = emp_se,
        model_se = model_se,
        rel_se_error = rel_se_error,
        coverage = coverage,
        rejection_rate = rejection_rate,
        mean_max_asmd = mean_asmd,
        mcse_bias = mcse_bias,
        mcse_coverage = mcse_coverage,
        mcse_rejection = mcse_rejection,
        stringsAsFactors = FALSE
      )
    })

  do.call(rbind, perf[!vapply(perf, is.null, logical(1))])
}

## ── Main execution ──────────────────────────────────────

cat("Factorial matched randomization simulation\n")
cat(sprintf("Scenarios: %d\n", nrow(scenarios)))
cat(sprintf("Methods: %d\n", length(randomization_methods)))
cat(sprintf("Replicates per condition: %d\n", n_replicates))
cat(sprintf(
  "Total conditions: %d\n",
  nrow(scenarios) * length(randomization_methods)
))

# Morris, White, and Crowther (2019) sec. 4.1: pin RNGkind and
# seed deterministically. Each scenario-method condition gets its
# own seed (base_seed + condition index), set inside the worker,
# so any condition can be rerun alone and results do not depend
# on n_cores. Per-replicate `.Random.seed` states are not stored.
RNGkind("L'Ecuyer-CMRG")

conditions <- expand.grid(
  m = randomization_methods,
  s = seq_len(nrow(scenarios)),
  stringsAsFactors = FALSE
)
conditions$seed <- base_seed + seq_len(nrow(conditions))

run_seeded <- function(i) {
  set.seed(conditions$seed[i], kind = "L'Ecuyer-CMRG")
  run_condition(conditions$s[i], conditions$m[i])
}

t_start <- Sys.time()
all_results <- parallel::mclapply(
  seq_len(nrow(conditions)), run_seeded,
  mc.cores = n_cores, mc.preschedule = FALSE
)
cat(sprintf("Elapsed: %.1f min\n",
  as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
failed <- vapply(all_results, inherits, logical(1), "try-error")
stopifnot(!any(failed))

raw_results <- do.call(rbind,
  all_results[!vapply(all_results, is.null, logical(1))]
)

performance <- compute_performance(raw_results)

out_dir <- file.path("analysis", "data", "derived_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  raw_results,
  file.path(out_dir, "sim_raw_results.rds")
)
saveRDS(
  performance,
  file.path(out_dir, "sim_performance.rds")
)

cat("\nPerformance summary (main effects, neyman estimator):\n")
sub <- performance[
  performance$effect == "F1" &
  performance$estimator == "neyman",
]
print(
  sub[, c("method", "n", "r_squared", "tau_main_a",
          "bias", "coverage", "rejection_rate",
          "mean_max_asmd")],
  digits = 3
)

cat("\nSimulation complete.\n")
cat(sprintf("Results saved to %s/\n", out_dir))
