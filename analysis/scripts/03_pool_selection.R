## ──────────────────────────────────────────────────────────
## 03_pool_selection.R
## Who is retained under pool oversampling, and what that does
## to the estimand under an illustrative effect modifier.
##
## Usage (from the package root):
##   Rscript analysis/scripts/03_pool_selection.R [n_pools]
##
## For each pool ratio r, draw a pool of N * r subjects with
## AR(1) normal covariates, form greedy Mahalanobis tuples,
## keep the tightest N / k tuples (match_tuples(pool_ratio = r))
## and summarize the retained subjects against the full pool.
##
## The effect modifier below is ILLUSTRATIVE, chosen for this
## script only: tau(x) = 0.2 + 0.1 * (x1^2 - 1), so the pool
## (population) average effect is 0.2 SD. Under within-tuple
## randomization the design estimates the average effect among
## retained subjects, so the gap between the two averages is
## the estimand shift caused by oversampling.
## ──────────────────────────────────────────────────────────

if (file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(factorialrandom)
}

args <- commandArgs(trailingOnly = TRUE)
n_pools <- if (length(args) >= 1) as.integer(args[1]) else 50L

n_trial <- 500L
k <- 4L
p <- 5L
cov_rho <- 0.3
pool_ratios <- c(1, 2, 4)
tau_0 <- 0.2
gamma <- 0.1

Sigma_x <- outer(seq_len(p), seq_len(p),
                 function(i, j) cov_rho^abs(i - j))
cov_names <- paste0("X", seq_len(p))

one_pool <- function(r) {
  n_pool <- ceiling(n_trial * r)
  X <- MASS::mvrnorm(n_pool, rep(0, p), Sigma_x)
  colnames(X) <- cov_names
  tup <- match_tuples(as.data.frame(X), k = k,
                      covariates = cov_names,
                      pool_ratio = r)
  keep <- !is.na(tup$tuple_id)
  tau_i <- tau_0 + gamma * (X[, 1]^2 - 1)
  d2 <- stats::mahalanobis(X, rep(0, p), Sigma_x)
  outer_pool <- d2 > stats::qchisq(0.9, p)
  data.frame(
    pool_ratio = r,
    n_retained = sum(keep),
    sd_x1_retained = stats::sd(X[keep, 1]),
    frac_outer_retained = mean(keep[outer_pool]),
    frac_inner_retained = mean(keep[!outer_pool]),
    tau_pool = mean(tau_i),
    tau_retained = mean(tau_i[keep])
  )
}

RNGkind("L'Ecuyer-CMRG")
set.seed(20261119)

res <- do.call(rbind, lapply(pool_ratios, function(r) {
  do.call(rbind, replicate(n_pools, one_pool(r), simplify = FALSE))
}))
stopifnot(all(res$n_retained == n_trial))

summ <- do.call(rbind, lapply(split(res, res$pool_ratio), function(d) {
  data.frame(
    pool_ratio = d$pool_ratio[1],
    n_pools = nrow(d),
    sd_x1_retained = mean(d$sd_x1_retained),
    frac_outer_retained = mean(d$frac_outer_retained),
    frac_inner_retained = mean(d$frac_inner_retained),
    tau_pool = mean(d$tau_pool),
    tau_retained = mean(d$tau_retained),
    shift = mean(d$tau_retained - d$tau_pool),
    mcse_shift = stats::sd(d$tau_retained - d$tau_pool) /
      sqrt(nrow(d))
  )
}))

out_dir <- file.path("analysis", "data", "derived_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(summ, file.path(out_dir, "pool_selection.rds"))
print(summ, digits = 3)
