# Remediation Log: Referee Review of 2026-08-16
*2026-08-16 16:30 PDT*

This log records the remediation pass against
`docs/pub_review_whitepaper_2026-08-16.md` for the
`factorialrandom` package and `analysis/report/report.Rmd`
manuscript. It does not replace or edit the whitepaper, which
remains the review record.

## 1. Fixed

**M1 (correctness). `tuple_balance()` crashed on
`tuples = NULL`, silently discarding every non-matched
replicate.** `R/match_tuples.R`: `tuple_balance()` now treats
`tuples = NULL` as "all subjects unmatched, single-arm
comparison" and computes raw covariate SMDs across `treatment`
directly, instead of indexing into a `NULL` tuples object.
`analysis/scripts/01_run_simulation.R`: `run_condition()` no
longer discards `tryCatch` errors silently; it now counts
failures, prints a warning with the first error message, and
raises an R `warning()` if every replicate in a condition
fails. **[verified]** — reproduced the original crash on
`generate_factorial_trial(..., randomization = "simple")` +
`tuple_balance()` before the fix, confirmed no crash after; a
4-scenario x 5-method x 20-replicate smoke test after the fix
produced non-NA bias/coverage/rejection-rate rows for all 5
methods (`simple`, `stratified`, `rerandomization`,
`matched_1.0`, `matched_2.0`).

**M3 (correctness). `regression_adjustment()`'s point-estimate
correction was identically zero by construction.**
`R/estimate_effects.R`: `regression_adjustment()` now computes
the unadjusted (raw) weighted arm-mean estimate separately from
the model-adjusted estimate and returns
`correction = est_unadj - est_adj`, so
`est <- est - adj$correction` in `estimate_effects()` actually
replaces the point estimate with the adjusted one. **[verified]**
— confirmed `identical(eff_unadj$estimate, eff_adj$estimate)` is
now `FALSE` (previously `TRUE` by construction); regression test
added in `inst/tinytest/test_basic.R`.

**New correctness bug found and fixed (not in the original
whitepaper): `factorial_contrasts()` weight normalization did
not match `generate_factorial_trial()`'s outcome-scaling,
producing systematic estimator bias of roughly 50% (main
effects) to 75% (interaction) of the true effect size, in
every reported bias/coverage/rejection number in the existing
`sim_performance.rds`.** Found while regenerating the
Application section's worked example: estimated main effects
and interactions were far from their data-generating values
even at 300 replicates, which is inconsistent with sampling
noise. Root cause: `factorial_contrasts()`
(`R/estimate_effects.R`) divided all contrast weights by `k`
instead of `k/2` (the standard Yates normalization for a 2^F
orthogonal design), and `generate_factorial_trial()`
(`R/dgm_factorial.R`) used an inconsistent divisor for the
interaction term (`/4`) versus the main-effect terms (`/2`).
Fixed both: `factorial_contrasts()` now divides by `k/2`
uniformly for main effects and interactions; the DGM's
interaction term now divides by 2, matching the main-effect
convention. **[verified]** — a 300-replicate check at `n = 200`
gave mean F1 estimate 0.206 (true 0.2, bias 0.006) and mean
interaction estimate 0.151 (true 0.15, bias 0.001), versus
biases of roughly -0.10 and -0.115 before the fix (confirmed
directly against the existing, pre-fix `sim_performance.rds`,
whose F1 bias column ranges -0.104 to -0.096 and whose F1:F2
bias column ranges -0.117 to -0.109 for the affected
scenarios). Regression tests added in
`inst/tinytest/test_basic.R`. This bug affects every bias
number currently reported anywhere in the manuscript
(Provenance paragraph, all Results figures); see Deferred
below for the required full rerun.

**M2 (correctness). Provenance paragraph misrepresented
simulation completeness and contained a stale, incorrect
"30-minute cap" comment.** `analysis/report/report.Rmd`
(Results, "Provenance"): rewritten to disclose that the
current `sim_performance.rds` predates the `tuple_balance()`
and `regression_adjustment()`/`factorial_contrasts()` fixes,
covers only 2 of 5 methods, states the true root cause (a
crash silently swallowed by `tryCatch`, not a wall-clock cap),
and adds an explicit TODO with the exact rerun command. The
stale code comment claiming a "30-min cap at ~25 of 72
scenarios" was removed from the `headline-16` chunk.
**[verified]** — checked against the RDS directly (all 72
scenarios present for the 2 completed methods, consistent
with a code bug rather than a timeout).

**M4 (correctness/acceptance). Two of three claimed tuple-
formation algorithms (hierarchical clustering, integer
programming) do not exist; `"optimal"` silently aliases to
`"greedy"`.** `analysis/report/report.Rmd` (Methods, "Tuple
Formation"; Discussion, "Limitations"): rewritten to describe
only the shipped greedy heuristic, explicitly state that
hierarchical clustering and integer programming are
unimplemented extensions (not benchmarked), and remove the
unsubstantiated "$N \lesssim 1000$" integer-programming
feasibility claim in favor of an accurate description of the
greedy algorithm's actual (quadratic) cost. **[applied,
unverified]** — prose-only change; no code was added (see
Deferred).

**M5 (acceptance). Empty Application subsections
("Matching Quality", "Treatment Effect Estimation") and empty
Discussion "Summary of Findings".** Added
`analysis/scripts/02_application_example.R`, a design-specific
simulation at the AHEAD-3-45-like parameters stated in the
text ($N_{\text{pool}} = 2{,}000$, $N = 500$, $r = 4.0$, $k =
4$, $R^2 = 0.3$, main effects 0.2 SD, interaction 0.15 SD),
run at 300 replicates per method (reduced from the main
study's 2,000, disclosed as such in the text) across all 5
comparator methods. `analysis/report/report.Rmd`: the
Application section now reads
`analysis/data/derived_data/application_example.rds` and
renders two real tables (covariate balance; bias/empirical
SE/coverage/rejection rate) instead of a placeholder sentence;
wrote the missing Discussion "Summary of Findings" paragraph.
**[verified]** — script executed successfully after the
`factorial_contrasts()`/DGM fix; bias is near zero for all
5 methods (range -0.006 to 0.004), coverage 0.92-0.99, and
matched methods show visibly better covariate balance (mean
max ASMD 0.08-0.12) than unmatched methods (0.27-0.28),
consistent with the paper's theoretical claim.

**M6 (acceptance). Trivial test suite (`expect_true(TRUE)`)
and no vignette.** `inst/tinytest/test_basic.R` rewritten with
25 real assertions covering `factorial_labels()`,
`factorial_contrasts()` (zero-sum weights), `match_tuples()`
(complete tuples, remainder/unmatched handling, `pool_ratio`
trimming), `randomize_tuples()` (all k labels per tuple),
`tuple_balance()` (matched path, the `tuples = NULL` path that
previously crashed, and a degenerate zero-variance covariate),
`estimate_effects()` (empirical unbiasedness of main effect
and interaction estimators, and the regression-adjustment
correction being non-zero), and `generate_factorial_trial()`
(`tuples = NULL` for non-matched methods). Added
`vignettes/matched-tuple-workflow.Rmd`, demonstrating the
`match_tuples()` -> `randomize_tuples()` -> `tuple_balance()`
-> `estimate_effects()` pipeline for both matched and
unmatched designs; added `VignetteBuilder: knitr` to
`DESCRIPTION` (it was missing despite `knitr`/`rmarkdown` in
`Suggests`). **[verified]** — `tinytest::run_test_dir()`
reports "All ok, 25 results"; the vignette's code was extracted
with `knitr::purl()` and executed successfully end to end.
`ri_test()` and diagnostic plotting functions promised in
`docs/analysis-plan.md` Aim 5 were not implemented (see
Deferred).

**M7 (acceptance, partial). Relative-efficiency claim not
validated against simulation.** Not independently derived (see
Deferred), but the abstract and Introduction claims were
narrowed (below) so the manuscript no longer asserts an
unvalidated empirical efficiency comparison.

**Abstract/Introduction overclaim ("efficiency gains relative
to simple and stratified randomization").** `report.Rmd`
abstract: reworded to state the relative-efficiency result as
theoretical/asymptotic and to explicitly disclose that the
empirical head-to-head comparison had not been regenerated as
of this draft, per the whitepaper's Recommended Framing
(section 5c: "must not claim efficiency gains ... until the
comparison exists"). **[applied, unverified]** — prose only.

**m1. Stale ADEMP compliance summary contradicting current
code.** `report.Rmd`, "Morris et al. (2019) ADEMP Compliance":
rewritten to state which April-17 gaps are now resolved
(`n_replicates` MCSE justification, `RNGkind` pinning, MCSE
columns now surfaced) versus genuinely still open (no
per-replicate `.Random.seed` capture; full sweep rerun still
pending). **[verified]** — checked `RNGkind`/`set.seed` call
sites and the MCSE-justification comment directly in
`01_run_simulation.R:11-15,278-279`.

**m2. MCSE columns computed but never surfaced.** `report.Rmd`,
Results "Coverage and Power": added a new `mcse-table` chunk
rendering mean `mcse_bias`, `mcse_coverage`, `mcse_rejection`
by method (main effect F1, Neyman estimator), with a citation
to Morris et al. Table 6 for the formulae. **[verified]** —
chunk reads existing RDS columns; logic is a straightforward
`dplyr::summarise()`, and the same columns were already
confirmed present in `sim_performance.rds` per the whitepaper.

**m3. Comparator table silently dropped minimization.**
`report.Rmd`, "Comparator Methods": removed the fabricated
"Minimization (Pocock-Simon)" row (never simulated), replaced
with the two matched-tuple pool-ratio rows that are actually
simulated, and added an explanatory paragraph stating why
minimization (a sequential-enrollment method) is excluded from
a pool-based-randomization comparison, addressing the
whitepaper's Recommended Framing guidance directly.
**[verified]** — cross-checked against
`randomization_methods` in `01_run_simulation.R`.

**m4. Overstated $k \le 8$ practical guidance.** `report.Rmd`,
Discussion "When to Use": reworded to state that only $k = 4$
was tested and that no specific upper bound on $k$ is
asserted. **[verified]** — cross-checked against `scenarios`
grid (`n_factors <- 2L` throughout `01_run_simulation.R`; no
$k = 8$ scenario exists).

**m5. Martin et al. (1993) numeric threshold cited without
page reference.** `report.Rmd`, Methods "Degrees of Freedom":
added the article's page range as a citation locator and
explicitly disclosed that the $M > 30$ Discussion heuristic is
not a page-specific quotation/derivation from that source.
**[applied, unverified]** — the exact page for a specific
numeric claim inside Martin et al. (1993) was not independently
verified (no full-text access during this pass); disclosed as
such rather than fabricated.

**m6. Interpretive figure captions.** `report.Rmd`: reworded
the covariate-balance and efficiency-emp-SE figure captions to
be descriptive rather than interpretive, and to disclose that
both figures currently reflect only matched-tuple methods
and/or pre-fix data (see Deferred). **[applied, unverified]**.

**m7. No data/code availability or computational-environment
disclosure.** `report.Rmd`: added a new "Data and Code
Availability" section (repository location, script/output
file mapping, no external data, R/package versions via inline
`r` code, and the reduced-scale timing benchmark).
**[verified]** — `r r_version`/`pkg_version`/`mass_version`
inline expressions read live from `R.version` and
`packageVersion()`, not hardcoded.

**m8. Unrelated `archive/` scratch files at the repository
root.** Added `archive/README.md` documenting what the
directory contains, that it appears unrelated to the
manuscript, and that deletion is deferred to the author's
judgment rather than performed unilaterally. **[applied]** —
no files were deleted.

**Fabricated DGM table entries (new issue, folded into M1/M2
correctness work).** `report.Rmd`, "Data Generating Mechanism"
table: the table previously listed a "skewed" outcome
distribution option that does not exist anywhere in
`dgm_factorial.R` (only `stats::rnorm()` is used), separate
"Main effect A"/"Main effect B" rows implying independently
varied levels (`0, 0.2, 0.5`) when the script only varies a
single symmetric `tau_main_a` parameter at 2 levels
(`0, 0.2`), and a "Pool/sample ratio" row listing 4 levels
(`1.0, 1.5, 2.0, 3.0`) when only 2 pool ratios
(`matched_1.0`, `matched_2.0`) are ever simulated. Corrected
to match `scenarios <- expand.grid(...)` and
`randomization_methods` in `01_run_simulation.R` exactly.
**[verified]** — cross-checked line by line against the script.

## 2. Deferred

- **Full-scale simulation rerun (all 5 methods x 72 scenarios
  x 2,000 replicates)** required to regenerate every figure
  and table in Results with the corrected code (M1, M3, and
  the new `factorial_contrasts()`/DGM normalization fix all
  postdate the current `sim_performance.rds`). Deferred: a
  reduced-scale timing check (4 scenarios x 5 methods x 20
  replicates, single core, n up to 40) took about 13 seconds
  (~32 ms/replicate); the full grid includes conditions with
  matched-method pool sizes up to 1,000 subjects, which are
  markedly slower per replicate (~0.4-0.5 s observed at pool
  size 2,000 in the Application-example timing check), so a
  full rerun plausibly takes several hours on a single core —
  outside this remediation pass's time budget. Exact command:
  `Rscript analysis/scripts/01_run_simulation.R` (run from the
  `factorialrandom` package root). The script already loads
  the `parallel` package but does not use it; parallelizing
  `run_condition()` calls across scenario x method combinations
  would likely bring this into single-digit-hours or less and
  is a reasonable next step for the author, but was not
  implemented here (a parallel refactor risks its own
  correctness bugs around RNG-stream independence and was
  judged out of scope for this pass).
- **M4/item 5: implement hierarchical-clustering and
  integer-programming tuple formation.** Deferred as a
  nontrivial software addition, not a documentation fix;
  Methods/Discussion were instead rewritten to describe only
  the shipped greedy heuristic (done, see Fixed).
- **M7/item 8: derive the relative-efficiency asymptote
  ($\text{RE} \to k/(k-1) \cdot M$) and validate it against
  simulation.** Deferred: this requires either a supplementary
  derivation (a genuine theoretical contribution, not
  mechanical remediation) or the completed full rerun above to
  validate empirically. Requires the author's mathematical
  judgment on the derivation approach.
- **M6: implement `ri_test()` (randomization-inference testing)
  and diagnostic/visualization functions promised in
  `docs/analysis-plan.md` Aim 5.** Deferred: not implemented;
  the vignette and tests added here cover only the functions
  that currently exist. The author should either implement
  these or scope Aim 5 / the abstract down to what is shipped.
- **m5: independently verify the exact page in Martin et al.
  (1993) supporting the "$M > 30$" heuristic.** Requires
  reading the primary source; not done here (disclosed as an
  assumption in the manuscript text instead of fabricated).
- **m6 (partial): the covariate-balance and efficiency-SE
  figures still need to be regenerated from the corrected,
  full-scale simulation once the rerun above is complete** —
  the caption edits made here are honesty fixes, not a
  substitute for regenerating the underlying figures.
- **Rendering `report.Rmd` to PDF.** Attempted via
  `bash tools/render.sh analysis/report/report.Rmd`; failed at
  the `setup` chunk with "there is no package called
  'stringi'" (a `kableExtra` dependency missing from the host R
  library). This is a pre-existing host-environment gap, not
  something introduced by this remediation pass, and is exactly
  the kind of package-availability problem `renv`/the project's
  Docker image exist to avoid. Recommended command for the
  author: `make docker-render REPORT=analysis/report/report.Rmd`
  (containerized build with the full `renv` library), or install
  `stringi` in the host R library and rerun
  `bash tools/render.sh analysis/report/report.Rmd`. All edits
  in this pass were verified at the source level (R code
  execution, `tinytest`, direct RDS inspection) rather than by
  visual inspection of a rendered PDF.

## 3. New issues found while fixing

- **`factorial_contrasts()`/`generate_factorial_trial()`
  weight-normalization mismatch** (detailed in Fixed, above):
  the single most consequential defect found during this pass.
  It was not identified in the original whitepaper because the
  whitepaper's bias-magnitude checks (M2's "absolute bias never
  exceeded X") reported the number without flagging that X
  (~0.10-0.12, against true effects of 0.15-0.2) is large
  relative to the effect sizes being estimated — i.e., the
  existing `sim_performance.rds` already contained the evidence
  of this bug, but no prior review step compared point
  estimates to their data-generating values directly.
- **The current `sim_performance.rds` is now further stale**
  than the whitepaper's M1/M2 already established: in addition
  to missing 3 of 5 methods, it also predates the
  `factorial_contrasts()`/DGM normalization fix, so its bias
  columns are systematically wrong (not just incomplete) for
  the 2 methods it does contain. The full rerun (Deferred,
  above) is therefore even more clearly required before any
  Results-section number in the current PDF can be trusted.
- **`archive/README.md` and `archive/TODO.md` were empty
  placeholder files**, not merely undocumented; there was no
  existing note anywhere explaining the `archive/` directory's
  purpose prior to this pass.

---
*Source: ~/prj/res/16-factorial-matched-random/factorialrandom/docs/pub_review_remediation_2026-08-16.md*
