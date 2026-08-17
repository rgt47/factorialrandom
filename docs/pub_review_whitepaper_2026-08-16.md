# Referee Review: Matched-Tuple Randomization for Factorial Clinical Trial Designs

*Review date: 2026-08-16 16:06 PDT*

Repository: `~/prj/res/16-factorial-matched-random` (inner package
repo `factorialrandom`, git commit `3839f48` at time of review).
Manuscript reviewed: `analysis/report/report.Rmd` (rendered to
`analysis/report/report.tex` and staged as
`share/report-2026-08-15-1643-3839f48.pdf`). Supporting documents:
`docs/analysis-plan.md`, `docs/literature-review.md`,
`docs/morris-audit-2026-04-17.md`. Package source under `R/`,
simulation driver `analysis/scripts/01_run_simulation.R`, simulation
output `analysis/data/derived_data/sim_performance.rds`, test suite
`inst/tinytest/test_basic.R`.

## 1. Summary of the work under review

**`analysis/report/report.Rmd`** ("Matched-Tuple Randomization for
Factorial Clinical Trial Designs: Theory and Application to
Trial-Ready Cohorts") proposes matched-tuple randomization: subjects
are grouped into tuples of size $k$ (the number of factorial
treatment combinations) using Mahalanobis or prognostic-score
distance, and each tuple member is randomly assigned a distinct
treatment combination. The manuscript states difference-in-means and
regression-adjusted estimators, an exact/conservative variance
formula, a relative-efficiency result versus simple randomization,
and a degrees-of-freedom caveat. It reports a $2\times2$ factorial
simulation study comparing five randomization strategies (simple,
stratified, minimization is described in the analysis plan but
dropped from the manuscript's comparator table in favor of
rerandomization, matched at pool ratio 1.0, matched at pool ratio
2.0) across 72 data-generating scenarios, and sketches an application
to an AHEAD-3-45-like Alzheimer's prevention trial. An R package,
`factorialrandom`, is offered as the software contribution. This is
the only manuscript-length document in the repository; no companion
paper exists yet, though the text explicitly defers the worked
application example to "a companion paper."

The two supporting planning documents (`docs/analysis-plan.md`,
`docs/literature-review.md`) are not manuscripts but internal
project-management artifacts: an aims/timeline document and a
narrative literature synthesis. They are reviewed here only insofar
as they bear on the completeness and framing of the main report,
per the instructions in this repository's review skill.

## 2. Major issues

**M1. The simulation contains no comparator data — the paper's
central empirical claim is unsupported (verified).**
`analysis/data/derived_data/sim_performance.rds` contains exactly two
values of `method`: `matched_1.0` and `matched_2.0` (1,296 rows, 72
of 72 scenarios, 2,000 replicates each — verified by loading the RDS
directly). `simple`, `stratified`, and `rerandomization` — three of
the five methods listed in the manuscript's own Table 2
(`report.Rmd:564-592`) and required for the headline comparison
promised in the Introduction and Simulation Study sections — are
entirely absent. Root cause (inspected in
`analysis/scripts/01_run_simulation.R:47-70` and
`R/match_tuples.R:194-231`): `run_one_replicate()` calls
`tuple_balance(dat, trial$tuples, ...)` unconditionally, but
`generate_factorial_trial()` (`R/dgm_factorial.R:59-102`) returns
`tuples = NULL` whenever `randomization != "matched"`. Inside
`tuple_balance()`, `matched <- !is.na(tuples$tuple_id)` on a `NULL`
tuples object silently produces `logical(0)`, `dat <- data[matched,]`
yields zero rows, `pooled_sd <- sd(numeric(0))` is `NA`, and
`if (pooled_sd < 1e-12)` throws `"missing value where TRUE/FALSE
needed"`. Verified directly: calling
`generate_factorial_trial(..., randomization = "simple")` followed by
`tuple_balance()` reproduces this exact error. Because
`run_one_replicate()` is wrapped in `tryCatch(..., error = function(e)
NULL)` inside `run_condition()`, every one of the 2,000 replicates for
`simple`, `stratified`, and `rerandomization`, across all 72
scenarios (432,000 replicate attempts), failed silently and
contributed zero rows. Every figure and every prose claim in Results
that references "randomization method" as a comparison axis
(covariate-balance figure, coverage-by-method figure, and the
Introduction/Discussion claims of efficiency gains "relative to
simple and stratified randomization") is therefore built on data that
compares matched-tuple designs only to *other* matched-tuple designs
(1x vs. 2x pool oversampling). No result in the current manuscript
demonstrates that matched-tuple randomization outperforms simple
randomization, stratified randomization, or rerandomization — the
three baselines the paper is framed against. **Remediation:** fix
`tuple_balance()` to handle `tuples = NULL` (treat all subjects as a
single "arm-only" comparison, or skip the balance diagnostic and
compute simple raw covariate SMDs across `treatment` directly), rerun
the full simulation, and rewrite every affected figure, caption, and
Results paragraph. This is a blocking defect for submission.

**M2. The manuscript's own "Provenance" paragraph misrepresents the
completeness of the simulation (verified).** `report.Rmd:611-648`
states "The full factorial sweep completed at 2,000 replicates per
cell, producing the `r n_scen_done` scenario-method-estimator rows,"
citing bias/coverage/rejection ranges as if they characterize the
full comparator set. This is only true for 2 of 5 methods (see M1);
the paragraph does not disclose that 3 of 5 methods are absent. A
code comment immediately above this text
(`report.Rmd:618-620`, "Production run hit the 30-min cap at ~25 of
72 scenarios") offers a *different*, also incorrect, explanation —
verified against the RDS, all 72 scenarios are present for the two
methods that did complete, so a 30-minute wall-clock cap stopping
mid-sweep is not consistent with the data. Two mutually inconsistent
narratives about the same run are present in the source and neither
matches the true failure mode (a code bug, not a timeout). A referee
will read this as either careless bookkeeping or an attempt to gloss
over a materially incomplete study. **Remediation:** after fixing M1
and rerunning, replace this paragraph with an accurate provenance
statement, and remove the stale/contradictory comment.

**M3. Reported "regression-adjusted" and "Neyman" point estimates are
identical by construction — the adjustment described in the text is
never applied (verified).** The manuscript states
(`report.Rmd:462-476`) $\hat\tau_c^{\text{adj}} = \hat\tau_c -
\hat\beta'(\bar X_c^+ - \bar X_c^-)$. In
`R/estimate_effects.R:219-273`, `regression_adjustment()` computes
`est_adj <- sum(weights * arm_coefs)` (line 252) and then returns
`correction = est_adj - sum(weights * arm_coefs)`, which is
algebraically identical to zero for every call, regardless of data.
Consequently `est <- est - adj$correction` in
`estimate_effects()` (line 82) leaves the point estimate unchanged
from the unadjusted difference-in-means; only `se` is replaced by the
residual-based standard error. Verified empirically:
`sim_performance.rds` shows `bias` and `emp_se` for `estimator ==
"adjusted"` equal to `estimator == "neyman"` for all 144 checked
(scenario × method) combinations examined (0 discrepancies), while
`model_se` differs by up to 0.059. This means the "efficiency-emp-se"
figure (`report.Rmd:673-697`), whose caption claims "Adjusted
estimation tightens the empirical SE as $R^2$ rises," cannot show
what it claims: empirical SE (the standard deviation of the point
estimate across replicates) is identical across estimators by
construction, since the point estimate itself never changes. Any
apparent separation between estimator lines in that figure is either
an artifact of the plotted summary or reflects `model_se`
(reported SE), not the empirical SE the caption names.
**Remediation:** fix `regression_adjustment()` to actually return the
adjustment (`correction = original_diff_in_means_estimate -
est_adj`, or restructure so `est` is directly replaced by
`est_adj`), rerun simulations, and re-verify that the adjusted
estimator's point estimates differ from the unadjusted ones before
re-plotting.

**M4. Two of three claimed $k$-tuple formation algorithms do not
exist in the software (inspected).** Methods §"Tuple Formation"
(`report.Rmd:392-411`) describes three approaches for $k>2$:
sequential greedy, hierarchical clustering (Ward's method or complete
linkage, with post-processing to enforce exact tuple size), and
integer programming (binary assignment formulation, "feasible for $N
\lesssim 1000$"). `R/match_tuples.R:43-49` implements only two named
branches, `"greedy"` and `"optimal"`, and the `"optimal"` branch is a
dead alias: `result <- match_greedy(dist_mat, k)` — identical to the
greedy branch, no integer program is formulated or solved. No
hierarchical-clustering code path exists anywhere in `R/`. The
Discussion's computational-limitation claim ("Optimal $k$-tuple
formation by integer programming is feasible for $N \lesssim 1000$
but may require heuristic approaches for larger pools,"
`report.Rmd:821-826`) is therefore a claim about a method that is not
implemented and was never benchmarked in this repository.
**Remediation:** either implement the IP and hierarchical methods and
benchmark them (per Aim 1 of `docs/analysis-plan.md`), or rewrite
Methods and Discussion to describe only the greedy heuristic actually
shipped, removing unsubstantiated complexity and feasibility claims.

**M5. Empty "Application" section and empty Discussion subsection
(inspected).** `report.Rmd:759-790`, subsections "## Matching
Quality" and "## Treatment Effect Estimation," contain section
headers with no body text, tables, or figures — only a single
placeholder sentence stating that "a full worked example ... is
reserved for a companion paper." `report.Rmd:792-794`, "## Summary of
Findings" under Discussion, is likewise an empty heading. A
manuscript cannot be evaluated as submission-ready with stubbed
sections that promise content ("the design tables in the next
subsection summarise the expected matching quality, treatment-effect
estimation precision, and coverage for that scenario,"
`report.Rmd:776-786`) and then deliver none. If the intent is to
defer the worked example to a companion paper, the promissory
sentence and the two empty subsection headers should be removed from
this manuscript, not left as visible stubs; if the intent is to
retain the application as part of this paper, the tables/figures must
be produced from the (corrected) simulation output specific to the
AHEAD-3-45-like design parameters, not merely asserted to follow from
the generic simulation grid. **Remediation:** either populate these
sections with a concrete worked example (design-specific simulation
at $N_{\text{pool}}=2000$, $r=4.0$, $k=4$) or delete the stubs and
adjust the framing sentence accordingly; write the missing
Discussion summary paragraph.

**M6. Software contribution is materially short of what the abstract
and Aim 5 promise (verified/inspected).** The abstract states
"Software implementing the proposed methods is provided in the R
package `factorialrandom`." `docs/analysis-plan.md` Aim 5 promises
"Diagnostic and visualization tools" and a worked
`ri_test()` randomization-inference function
(`docs/analysis-plan.md:299`, `## 4.5 Software Design` code block).
Inspected: no `ri_test()`, no plotting/diagnostic function, and no
vignette exists anywhere in `R/`, `man/`, or `vignettes/` (the
`vignettes/` directory is empty; `DESCRIPTION` lists `rmarkdown` and
`knitr` in `Suggests` but ships nothing that uses them). The test
suite (`inst/tinytest/test_basic.R`) contains a single trivial
assertion, `expect_true(TRUE)` — verified by running
`tinytest::run_test_dir()`, which reports "1 tests OK." There is no
executable check that `match_tuples()`, `estimate_effects()`,
`tuple_balance()`, or `generate_factorial_trial()` behave correctly,
which is how bugs M1, M3, and M4 went undetected. A referee for a
software-adjacent statistics venue (or the R Journal, if that framing
is chosen — see §5) will require substantive test coverage and a
vignette walking through the workflow shown in the abstract.
**Remediation:** write unit tests covering each exported function
(including the `tuples = NULL` path that currently crashes
`tuple_balance()`), add a vignette demonstrating the
`match_tuples()` → `randomize_tuples()` → `estimate_effects()`
pipeline, and either implement `ri_test()` and diagnostic plotting or
remove them from the stated aims.

**M7. Relative-efficiency theory (Eq. "RE", `report.Rmd:478-494`) is
asserted, not derived, and cannot currently be checked against
simulation.** The manuscript states $\text{RE} \to k/(k-1)\cdot M$ as
$\rho_w \to 1$ "as a function of the within-tuple intraclass
correlation $\rho_w$," but no derivation, appendix, or supplementary
proof is provided anywhere in the repository — this is presented as
a bare assertion in the main text. Because M1 leaves no simple
randomization data, this closed-form claim is also empirically
untested even after the derivation gap is closed. **Remediation:**
add a derivation (main text or supplementary appendix) connecting
$\rho_w$ to the potential-outcomes variance decomposition already set
up in §"Matched-Pair Designs," and validate the asymptotic formula
against the corrected simulation's empirical RE.

## 3. Minor issues

**m1. Stale ADEMP compliance summary contradicts current code
(inspected).** `report.Rmd:857-874` reproduces the April 17 audit's
"Key gaps" list verbatim, including "`RNGkind('L'Ecuyer-CMRG')` not
pinned." `analysis/scripts/01_run_simulation.R:278` already calls
`RNGkind("L'Ecuyer-CMRG")` immediately before `set.seed(20260312)`
(line 279), and lines 11-15 already contain an MCSE-based
justification for `n_replicates = 2000` that the same gap list
claims is missing. Two of the three listed gaps appear to be already
resolved in code but not in the manuscript's self-reported compliance
narrative. The remaining gap (no per-replicate `.Random.seed`
snapshot) does appear to still be genuinely open — verified, no such
capture exists in the script. **Remediation:** re-run the Morris
audit against current code before the next submission and update
this section; do not carry forward a stale audit as if current.

**m2. MCSE columns computed but never surfaced (verified, previously
flagged, still open).** `compute_performance()`
(`01_run_simulation.R:227-233`) computes `mcse_bias`,
`mcse_coverage`, `mcse_rejection` and stores them in
`sim_performance.rds` (confirmed present as columns), but no table or
figure in `report.Rmd` displays them; the Results section reports
point ranges of bias/coverage/rejection without their Monte Carlo
uncertainty. This was already identified in
`docs/morris-audit-2026-04-17.md` and remains unaddressed four months
later per the file timestamps.

**m3. Comparator table (Table 2, `report.Rmd:564-592`) describes
minimization dropped from the actual comparator set without
comment.** `docs/analysis-plan.md` §4.4.2 lists "Minimization —
Pocock-Simon dynamic balancing" as comparator method 3 of 5; the
manuscript's Table 2 replaces it with rerandomization only, silently
dropping minimization with no note explaining the change from the
plan. `docs/literature-review.md` §2.1 also discusses minimization at
length as foundational context, making its absence from the
comparator set read as an unexplained gap between planning documents
and the manuscript.

**m4. `n_factors` and design generality overstated relative to what
was tested.** The Methods section is written generally for $2^F$
designs and even for fractional factorials (Future Work), but the
entire simulation and application are restricted to $F=2$ ($k=4$).
The Discussion's practical guidance ("$k \le 8$," `report.Rmd:806`)
is not supported by any simulation at $k=8$; it is an unverified
extrapolation.

**m5. Degrees-of-freedom claim (`report.Rmd:495-504`) cites
Martin et al. (1993) for a specific numeric threshold ("$M > 30$",
`report.Rmd:819`) without reproducing or deriving the threshold in
this manuscript.** Inspected only; this reviewer did not verify the
cited source's numeric claim independently. Should be either derived
or given a direct quotation/page reference.

**m6. Figure captions state interpretive conclusions ("consistently
outperforms," `report.Rmd:652`) that belong in prose discussion, not
captions, and pre-empt the reader's own reading of the graphic** — a
presentation style point common in referee reports for statistics
journals (e.g., JASA/Biometrics figure-caption conventions favor
neutral, descriptive captions).

**m7. No data/code availability statement, and no explicit statement
of computational environment (R version, package versions, hardware,
wall-clock time) for the simulation**, which the Morris et al. (2019)
ADEMP framework the paper explicitly invokes requires for
reproducibility reporting.

**m8. `archive/` directory at the repository root contains unrelated
LaTeX scratch files** (`x22.tex`, `tab1.tex`, `sdout.tmp`,
`sdout2.tmp`) with no evident connection to the matched-tuple
manuscript (inspected `archive/README.md`, `archive/TODO.md`; both
empty/uninformative). This is repository hygiene, not a manuscript
defect, but should be resolved (removed or documented) before the
repository is shared with collaborators or a journal's reproducibility
reviewers.

## 4. What remains to be done

**(a) Required for correctness**

1. Fix `tuple_balance()` to handle `tuples = NULL` without crashing
   (M1), and remove the `tryCatch(..., error = function(e) NULL)`
   masking in `run_one_replicate()`/`run_condition()`, or at minimum
   log and surface caught errors so silent total failure of a method
   cannot recur undetected.
2. Fix `regression_adjustment()` so the point-estimate correction is
   not identically zero (M3).
3. Rerun the full simulation (all 5 methods × 72 scenarios × 2,000
   reps) after fixes 1-2, and regenerate every figure/table in
   Results and the Application section from the corrected output.
4. Correct or remove the "Provenance" paragraph and the stale
   30-minute-cap code comment (M2) once the corrected run's true
   provenance is known.
5. Either implement hierarchical-clustering and integer-programming
   tuple formation and benchmark them, or rewrite Methods/Discussion
   to describe only the shipped greedy heuristic (M4).

**(b) Required for acceptance**

6. Populate the Application section (Matching Quality, Treatment
   Effect Estimation subsections) with a concrete worked example
   using the corrected simulation machinery at the stated
   AHEAD-3-45-like design parameters, or delete the stubs (M5).
7. Write the missing Discussion "Summary of Findings" paragraph
   (M5).
8. Add a derivation (main text or appendix) for the relative
   efficiency asymptote $\text{RE} \to k/(k-1)\cdot M$ (M7), and
   validate it against corrected simulation output.
9. Wire the already-computed MCSE columns into the Results tables
   (m2), consistent with the manuscript's own stated Morris et al.
   compliance goal.
10. Add non-trivial unit tests for every exported function, including
    edge cases (`tuples = NULL`, unmatched remainders, degenerate
    covariance) (M6).
11. Add a package vignette demonstrating the full
    match/randomize/estimate workflow shown in the abstract, or
    scope the abstract down to match what is shipped (M6).
12. Re-run the Morris et al. (2019) ADEMP audit against current code
    and replace the stale April 17 summary embedded in the
    manuscript (m1).
13. Reconcile the comparator method set with the analysis plan —
    either add minimization or explicitly justify its removal (m3).
14. Add a data/code availability statement and a computational
    environment disclosure (R version, package versions, hardware,
    wall-clock time for the full sweep) (m7).

**(c) Desirable polish**

15. Rewrite figure captions to be descriptive rather than
    interpretive (m6).
16. Either run and report a $k=8$ ($2^3$) scenario to support the
    Discussion's stated operating range, or soften the claim to match
    the tested range ($k=4$ only) (m4).
17. Resolve or remove the unrelated `archive/` scratch files (m8).
18. Add a page/section reference for the Martin et al. (1993)
    degrees-of-freedom threshold citation (m5).

## 5. Recommended framing

**Paper A: the matched-tuple randomization manuscript itself.**

*(a) Plausible framings.*

1. **New methodology paper** (algorithms + inference theory +
   simulation), targeting a general statistical methodology journal
   (Biometrics, Statistics in Medicine, JCGS).
2. **Applied/translational paper for the AD trials community**,
   emphasizing the trial-ready-cohort application and targeting a
   clinical-trials-methods or AD-specific venue (e.g., Alzheimer's &
   Dementia, Clinical Trials).
3. **Software/tools paper** for the R Journal or JSS, centered on the
   `factorialrandom` package, with theory condensed to a supporting
   role.

*(b) Recommendation: framing 1, a methodology paper, but scoped
honestly around what is actually established.* The literature review
(`docs/literature-review.md` §5) correctly identifies the gap:
matched-pair theory is well developed for $k=2$ (Bai 2022, Imai 2008,
Fogarty 2018), rerandomization is well developed but not
tuple-structured (Morgan and Rubin 2012), and no published method
combines matched-tuple formation with factorial treatment structure
and valid inference. That gap claim is credible and, as far as this
review's literature check went (citations all resolve to bibliography
entries; no independent verification of each cited paper's content
was performed — inferred from the literature review's own synthesis,
not independently re-derived). However, the manuscript as it
currently stands has not yet demonstrated the gap is filled: the
empirical evidence that would show matched-tuple randomization beats
the standard $k=2$-derived alternatives (simple, stratified,
rerandomization) does not exist in the current run (M1). A referee
cannot recommend acceptance of a methodology paper whose central
comparison is missing. Framing 2 (translational/AD paper) is
premature for the same reason — worse, because the entire Application
section is currently unwritten (M5) — and additionally because the
AD-specific material (TRC-PAD, AHEAD 3-45, DIAN-TU background) is
already well covered by existing review literature the manuscript
itself cites; a clinical-trials-methods reviewer will want the
methodological contribution demonstrated first, with the AD
application as illustration, not the reverse. Framing 3 (software
paper) is premature given M6: the package currently lacks tests,
vignettes, and one of its three advertised algorithms, and a
software-paper reviewer (R Journal/JSS) will check exactly these
things first.

*(c) Implications of the recommended framing.* Title and abstract
should keep the current "Theory and Application" framing only after
the application section is either completed or removed — if removed,
retitle to something like "Matched-Tuple Randomization for Factorial
Clinical Trial Designs" without the "Application to Trial-Ready
Cohorts" subtitle, and move the TRC/AD material to a motivating
example in the Introduction rather than a dedicated Application
section with empty subsections. The Introduction should keep the
$k=2$-literature-gap argument (it is the paper's strongest asset) but
must not claim efficiency gains "relative to simple and stratified
randomization" (as `report.Rmd:96-98` currently implies) until the
comparison exists. Comparators for the corrected simulation should
remain simple, stratified, rerandomization, and matched-tuple at two
pool ratios; minimization should either be added back (per the
analysis plan) or explicitly excluded with a stated reason (e.g., "not
applicable to designs where the full pool is available at once"),
since its omission is otherwise conspicuous given how much space
`docs/literature-review.md` gives it. Target journal: Biometrics or
Statistics in Medicine given the clinical-trials framing and the
existing $k=2$ optimality literature (Bai 2022) published in
economics/statistics venues that a methods-focused biostatistics
journal would want directly engaged.

*(d) What to emphasize/de-emphasize/move.* Emphasize: the variance
theory, the degrees-of-freedom tradeoff (a genuinely underexplored
point relative to existing $k=2$ literature), and a corrected,
complete simulation comparison. De-emphasize or move to supplementary
material: the extended TRC/AD background (§"Trial-Ready Cohorts,"
currently four full subsections, `report.Rmd:270-337`) — most of this
narrates well-known program history (TRC-PAD, AHEAD 3-45, DIAN-TU,
EPAD) that is not itself a new contribution and is already
well-covered in the cited program-description papers; a condensed
paragraph plus a supplementary table of program characteristics would
serve a methods-journal reviewer better than four subsections. Move
the fractional-factorial and unequal-allocation extensions currently
in "Future Work" to a short paragraph unless they are implemented;
as written they read as scope padding.

## 6. Assessment

**Verdict: major revision, and in its current state the manuscript
is not ready for the review process at all** — several of the Major
issues above (M1, M3, M4) are the kind of defect that would normally
be caught and corrected before initial submission, not surfaced by a
referee. A referee at a statistics journal encountering M1 (the
paper's central comparator data does not exist) would most likely
recommend rejection with encouragement to resubmit once the study is
actually complete, rather than major revision, because the paper's
claimed contribution cannot currently be evaluated at all — there is
no comparison to simple, stratified, or rerandomization to assess.
The theoretical contribution (variance formulas, degrees-of-freedom
tradeoff, connection to matched-pair optimality results) is a
reasonable methodological seed and the literature positioning is
sound (inspected literature review and citation list; all citations
resolve). But between the incomplete simulation, the empty
Application section, the two point-estimator bugs (M3, and the
matching-algorithm gap M4), and the near-absence of package test
coverage, this is fairly characterized as a partially assembled draft
rather than a submission-ready manuscript. Priority for the authors:
fix M1-M4 and rerun the simulation before any further framing or
polish decisions are made, since several Results-section conclusions
may change once real comparator data exists.

## 7. Revision history

- 2026-08-16: Initial review. Established that the simulation's
  comparator data (simple, stratified, rerandomization) is entirely
  absent due to a `tuple_balance()` crash on `NULL` tuples (M1);
  identified a zero-valued point-estimate correction bug in
  `regression_adjustment()` making "adjusted" and "Neyman" point
  estimates identical by construction (M3); confirmed two of three
  claimed tuple-formation algorithms (hierarchical clustering,
  integer programming) are unimplemented, with `"optimal"` aliasing
  to `"greedy"` (M4); found the manuscript's own provenance narrative
  self-contradictory and inconsistent with the underlying data (M2);
  found the Application section and a Discussion subsection stubbed
  with no content (M5); found package test coverage limited to a
  single trivial assertion and no vignette despite abstract claims of
  a complete software contribution (M6); found the manuscript's
  embedded Morris et al. (2019) compliance summary stale relative to
  already-fixed code (m1). No prior review exists for this repository.
