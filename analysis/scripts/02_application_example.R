## ──────────────────────────────────────────────────────────
## 02_application_example.R
## Design-specific simulation for the AHEAD-3-45-like worked
## example in the Application section of report.Rmd.
##
## Reduced-replicate run (see n_replicates below) added during
## the 2026-08-16 pub_review remediation pass to replace the
## previously-empty Application subsections with real,
## reproducible numbers. This is NOT the full 2,000-replicate
## sweep used elsewhere in the manuscript; it exists to give
## the Application section concrete, non-fabricated numbers at
## the specific design point described in the text
## (N_pool = 2000, r = 4.0, k = 4). See
## docs/pub_review_remediation_2026-08-16.md for rationale and
## a TODO to raise n_replicates once time allows.
## ──────────────────────────────────────────────────────────

library(factorialrandom)
library(MASS)

n_replicates <- 300L
n_factors <- 2L
k <- 2^n_factors
sigma_y <- 1

application_scenario <- list(
  n = 500L,
  p_covariates = 5L,
  r_squared = 0.3,
  tau_main_a = 0.2,
  tau_interaction = 0.15
)

randomization_methods <- c(
  "simple",
  "stratified",
  "rerandomization",
  "matched_1.0",
  "matched_4.0"
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

  data.frame(
    effect = eff$effect,
    estimate = eff$estimate,
    se = eff$se,
    ci_lower = eff$ci_lower,
    ci_upper = eff$ci_upper,
    p_value = eff$p_value,
    max_asmd = max_asmd,
    stringsAsFactors = FALSE
  )
}

run_condition <- function(method_label, n_reps = n_replicates) {
  cat(sprintf("  Method: %s\n", method_label))
  rep_results <- vector("list", n_reps)
  n_errors <- 0L
  for (r in seq_len(n_reps)) {
    rep_results[[r]] <- tryCatch(
      run_one_replicate(application_scenario, method_label),
      error = function(e) {
        n_errors <<- n_errors + 1L
        NULL
      }
    )
  }
  if (n_errors > 0L) {
    cat(sprintf(
      "  WARNING: %d/%d replicates failed for %s\n",
      n_errors, n_reps, method_label
    ))
  }
  valid <- !vapply(rep_results, is.null, logical(1))
  all_reps <- do.call(rbind, rep_results[valid])
  all_reps$method <- method_label
  all_reps
}

RNGkind("L'Ecuyer-CMRG")
set.seed(20260816)

raw <- do.call(rbind, lapply(randomization_methods, run_condition))

true_tau <- c(
  F1 = application_scenario$tau_main_a * sigma_y,
  F2 = application_scenario$tau_main_a * sigma_y,
  "F1:F2" = application_scenario$tau_interaction * sigma_y
)
raw$true_tau <- true_tau[raw$effect]

performance <- do.call(rbind, lapply(
  split(raw, interaction(raw$method, raw$effect, drop = TRUE)),
  function(df) {
    if (nrow(df) < 10) return(NULL)
    B <- nrow(df)
    bias <- mean(df$estimate) - df$true_tau[1]
    emp_se <- sd(df$estimate)
    covers <- df$ci_lower <= df$true_tau[1] &
      df$ci_upper >= df$true_tau[1]
    coverage <- mean(covers, na.rm = TRUE)
    rejection_rate <- mean(df$p_value < 0.05, na.rm = TRUE)
    data.frame(
      method = df$method[1],
      effect = df$effect[1],
      n_reps = B,
      bias = bias,
      emp_se = emp_se,
      coverage = coverage,
      rejection_rate = rejection_rate,
      mean_max_asmd = mean(df$max_asmd, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))

out_dir <- file.path("analysis", "data", "derived_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  performance,
  file.path(out_dir, "application_example.rds")
)

cat("\nApplication example performance:\n")
print(performance, digits = 3)
cat(sprintf("\nSaved to %s/application_example.rds\n", out_dir))
