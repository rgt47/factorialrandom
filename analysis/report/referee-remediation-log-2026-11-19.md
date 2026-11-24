# Referee Remediation Log: report.Rmd
*2026-11-19 00:45 PST*

Work list: `referee-report-2026-11-18.md` (recommendation:
reject). Built on the 2026-11-18 voice pass. Nothing was
committed, staged, or rendered. Pre-edit copies of
`report.Rmd`, `references.bib`, `R/`, the scripts, the tests
and the two small `.rds` files are in
`/private/tmp/claude-501/rem-16/`.

## Fixed

### Code (all [verified] by tinytest, 34/34 pass)

- **M4 rerandomization.** `R/randomize_tuples.R`:
  `mahalanobis_balance()` is now the k-arm statistic
  `sum_a n_a (xbar_a - xbar)' S^-1 (xbar_a - xbar)`, which is
  approximately chi-squared on p(k-1) df. The threshold is
  `qchisq(accept_prob, p(k-1))`, new argument `accept_prob =
  0.001`. Draws are scored in vectorized batches; the first
  accepted draw in order is returned. New tests: empirical
  acceptance near nominal (0.05) at (n, p) = (100, 3) and
  (500, 8); accepted draw meets the threshold. Direct check
  at the nominal 0.001 with 30,000 draws: 0.00083 at
  (100, 3) and 0.00097 at (500, 8).
- **M5 stratification.** `R/dgm_factorial.R`: strata are now
  sample quartiles of X1 (was equal-width `cut(breaks = 4)`).
  Comparator table now says 'quartiles of X1'. Test added.
- **M7 inference.** `R/estimate_effects.R`: paired estimator
  uses t on M-1 df for CIs and p-values; Lin-adjusted SE is
  now HC2 sandwich (was unadjusted per-arm residual
  variances); adjustment now also applies to unmatched
  designs, so Lin-adjusted simple randomization is a
  comparator. Tests added.
- **New: R-squared calibration bug.** `compute_beta_for_r2()`
  ignored the AR(1) covariance, so realized R-squared exceeded
  the target (Var(X'beta) 1.46x target at p = 3). Fixed;
  test added.
- `estimate_effects()` no longer grows a data.frame by
  `rbind` (speed only; outputs unchanged in form).

### Simulation rerun (M1, M8) [verified]

- `01_run_simulation.R`: loads the package from source
  (`pkgload::load_all`), takes `[n_reps] [n_cores]`, seeds each
  condition as `20260312 + condition index` inside the worker,
  runs conditions with `mclapply`. `pkgload` is in `renv.lock`
  but not in `DESCRIPTION`. Determinism checked: two
  conditions rerun alone reproduced the saved estimates
  exactly. Comment errors fixed (2,500 vs 2,000 MCSE wording;
  0.025 called type I error; false '.Random.seed stored').
- Command run: `Rscript analysis/scripts/01_run_simulation.R
  60 8` (4.4 min). Output: 2,592 rows, 5 methods x 72
  scenarios, 0 failed replicates, max |bias| 0.076 (max 3.5
  MCSE). `sim_performance.rds` and `sim_raw_results.rds`
  regenerated.
- `02_application_example.R`: same loading/seeding changes,
  MCSE columns added; rerun with `300 5`.
- New `03_pool_selection.R` (M6), run with `50`.

### Manuscript

- **M1**: Provenance note removed; replaced by 'Simulation
  scale' with replicate count read from the data and a
  visible TODO. All Results numbers are inline R from the new
  `.rds`. Coverage plot no longer clipped; axes say 'Total
  sample size N'. Type I error and power reported separately
  (minor 10). New RE table and power table. [verified: purled
  code evaluated, every inline value printed and checked
  against the prose]
- **M2**: RE limit `k/(k-1) M` withdrawn (abstract,
  contributions, Methods, Discussion). Replaced with a
  labeled heuristic bound `RE <= 1/(1-R^2)` for the linear
  model, compared with empirical RE. [applied; bound
  inspected, not proved]
- **M3**: Contribution 2 narrowed. Added unbiasedness proof,
  explicit V_m(c) over the k! permutations, conservativeness
  identity for the between-tuple estimator (the text wrongly
  said 'within-tuple'). Matched-pair Eq. (3) replaced with a
  direct calculation, citing Imai (2008) as 'see also' (not
  re-read). Asymptotics deferred to
  Bai, Liu and Tabord-Meehan (2024). [applied, derivations
  checked by hand]
- **M4/M5/M7**: prose matches the corrected code.
  [verified against code]
- **M6**: estimand paragraph in Methods; new results
  subsection 'Pool oversampling and the estimand'.
  [verified: values from `pool_selection.rds`]
- **M8**: Application now described as a one-point
  head-to-head comparison with MCSE. [verified]
- Minor 1, 2, 5, 8, 9, 10, 11, 12, 13, 14, 15 [verified
  against code or data]. Minor 3 (nonbipartite), 4
  (convolution), 6 (Higgins 'at least k units') and the
  Bai et al. 'consistent estimator' sentence are
  [applied, unverified]: they rest on the referee's recall,
  my recall, or abstracts. Minor 7: the outcome equation
  was indicator-coded and did not match the code; it is now
  in +-1 coding, and the weights +-1/2^(F-1) are tied to
  equation (1). I believe this is the Dasgupta et al.
  convention, so the referee's 'twice' would not apply, but
  I did not open the paper; the attribution was removed.
- Knitr `cache` set to FALSE so a rerun is not masked by
  cached 60-replicate chunks. Two kable captions escaped
  (`%`, `^`).
  Data availability and ADEMP moved out of 'References'
  (unnumbered sections); internal paths and remediation
  narration removed.

### References

- Added (each DOI checked on Crossref): Bai, Liu,
  Tabord-Meehan 2024 (10.3982/QE2354); Higgins, Savje,
  Sekhon 2016; Branson, Dasgupta, Rubin 2016; Li, Ding,
  Rubin 2020 (10.1214/18-AOS1790); Fogarty 2018 JRSS-B;
  Bai, Romano, Shaikh 2022 JASA.
- `solomonEPAD2022`: authors corrected to Saunders et al.
  (key kept).
- `mcdonoughTRCScreening2024`: DOI resolves to an unrelated
  nanolithium review. Removed. Replaced by
  `walterAlzMatch2024` (verified) and the bullet narrowed.
- `baiCovariateAdaptive2024`: uncited, wrong DOI. Removed.
- 80% screen failure and $10,000 cost: removed; Sperling
  2011 now cited only for the biomarker definition.
- Cummings 2024: claim narrowed to 'combination therapies
  appear in the pipeline'; examples marked as ours.
- New entries [verified on Crossref]; removals and the
  Saunders author fix [verified on Crossref]; Cummings
  narrowing [applied, unverified].

## Deferred

- **Full 2,000-replicate rerun** (about 2.5 h at 8 cores on
  the load seen here): `Rscript
  analysis/scripts/01_run_simulation.R 2000 8` from the
  package root, then re-render. Manuscript TODO in place.
  After the rerun, re-read these hardcoded sentences, which
  do not update from the data: the abstract findings; the
  Discussion summary; 'no replicate failed' (Simulation
  scale); 'best balance with p = 3' and 'as well as or
  better' (Covariate balance); 'small and inconsistent'
  (oversampling); 'did not outperform rerandomization' and
  'appeared only when it was also adjusted' (Results).
- Monte Carlo SEs for the relative-efficiency table were not
  computed; the text says so.
- Title still says 'Theory and Application'; with
  Contribution 2 narrowed, the author may wish to change it.
- Bai et al. variance estimator and formal RE derivation
  (needs the full paper; read from abstracts only).
- Heterogeneous-effect scenarios inside the main grid,
  skewed outcomes, k = 8: not done.
- EPAD '2,000 participants, 39 institutions' and TRC-PAD
  '30,000' figures: not checked against sources.
- Cummings 2024 support for combination therapies: from
  recall, not re-read.
- `zotero-sync/missing-references-for-zotero-import.bib` is
  an import file, not a source of `references.bib`; the
  Zotero library itself may still hold the bad entries.
- `man/` not regenerated for the changed roxygen
  (`accept_prob`, `Sigma_x`): run `Rscript -e "devtools::document()"`.
- Render not attempted (instructions). The vignette code was
  purled and run against the changed functions: it runs.

## New issues found while fixing

- R-squared calibration bug (fixed, above).
- Under a linear model, Lin-adjusted simple randomization
  reaches the same asymptotic bound as perfect matching. The
  data agree: adjusted simple RE 1.70 at R2 = 0.5, equal to
  or above unadjusted matched tuples. The paper's case for
  matching is weaker than its framing implied.
- Pool oversampling is highly selective: at r = 2 only 5%
  of outer subjects are kept vs 55% of inner subjects.

## Abstract voice pass

- Words: 178 -> 171. YAML parses (`yaml::yaml.load`).
- Factual changes: (1) removed the RE result (M2);
  (2) removed 'derive ... variance estimators' in favor of
  unbiasedness and conservativeness (M3); (3) removed 'not
  yet reported' for the comparison and added the findings
  (M1, M8); (4) added credit to Bai, Liu and Tabord-Meehan
  (novelty); (5) added the estimand change (M6); (6) noted
  the reduced replicate count.
- Longest sentence 19 words. None over 25.

## Illustrative values chosen in this pass

- Effect modifier `tau(x) = 0.2 + 0.1 (x1^2 - 1)`, outer
  cutoff at the chi-squared(5) 90th percentile, p = 5,
  N = 500, r in {1, 2, 4}, 50 pools (`03_pool_selection.R`).
- Replicate counts 60 (main) and 300 (Application), set by
  the 10-minute rule.
- Test-only acceptance rate 0.05.
