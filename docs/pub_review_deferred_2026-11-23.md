# Deferred-Item Pass: factorialrandom
*2026-11-23 12:50 PST*

Scope: the seven Deferred bullets in
`docs/pub_review_remediation_2026-08-16.md`, which answers
`docs/pub_review_whitepaper_2026-08-16.md`. Current state was
established from git history (commits through `dee2dc1`,
2026-09-07) and from the uncommitted working tree, which holds a
second remediation pass recorded in
`analysis/report/referee-remediation-log-2026-11-19.md`
(answering `analysis/report/referee-report-2026-11-18.md`). That
log is used here as evidence of current state, not as a second
inventory. Nothing was committed or staged. A pre-edit copy of
`report.Rmd` is at
`/private/tmp/claude-501/-Users-zenn-prj-res-01-fisher-exact-rx2-fisherexacttestrx2-analysis-report/be4b228f-9df5-41a4-9039-7499d4f23926/scratchpad/ws16/report.Rmd.pre-deferred`
(session scratch; may be cleared).

## 1. Inventory

| # | Source log | Item | Class |
|---|---|---|---|
| 1 | 2026-08-16 | Full rerun, 5 methods x 72 scenarios x 2,000 reps | RESOLVED-SINCE |
| 2 | 2026-08-16 | Implement hierarchical-clustering and IP tuple formation (M4) | AUTHOR-DECISION |
| 3 | 2026-08-16 | Derive and validate RE asymptote `k/(k-1) M` (M7) | NO-LONGER-APPLIES |
| 4 | 2026-08-16 | Implement `ri_test()` and diagnostics (analysis plan Aim 5) | AUTHOR-DECISION |
| 5 | 2026-08-16 | Page locator in Martin et al. (1993) for the `M > 30` heuristic (m5) | NO-LONGER-APPLIES |
| 6 | 2026-08-16 | Regenerate balance and efficiency figures from the full run (m6) | RESOLVED-SINCE |
| 7 | 2026-08-16 | Render `report.Rmd` (host lacked `stringi`) | RESOLVED-SINCE |

Evidence per item:

1. `analysis/data/derived_data/sim_performance.rds` (mtime
   2026-11-19 20:56:01) has 2,592 rows, 5 methods, 72 scenarios,
   and `n_reps == 2000` in every row. No log records this run,
   so its provenance was checked directly: two conditions
   (condition 1, simple, scenario 1, first 100 replicates;
   condition 4, matched_1.0, scenario 1, first 50 replicates)
   were regenerated from the current code with the script's
   per-condition seed, in a scratch script that wrote nothing
   to the repository. Point estimates and SEs matched
   `sim_raw_results.rds` exactly (max absolute difference 0).
   [verified, for those two conditions; the remaining 358 are
   inferred from the shared seeding scheme.] The results are
   uncommitted; they exist only in the working tree.
2. `R/match_tuples.R` still has only the greedy heuristic;
   `method = 'optimal'` still aliases to it. The manuscript
   (Methods, Tuple formation; Discussion, Limitations) already
   states this. `docs/analysis-plan.md` lines 176-178 still
   describe both methods.
3. The `k/(k-1) M` limit was withdrawn on 2026-11-19 (abstract,
   Contributions, Methods, Discussion) and replaced by a labeled
   heuristic bound `RE <= 1/(1 - R^2)` compared against
   empirical RE in a new table. `grep` finds no remaining
   instance.
4. No `ri_test()` in `R/` or `NAMESPACE`.
   `docs/analysis-plan.md` lines 300 and 310 still list it. The
   abstract and Contributions no longer claim it.
5. The Degrees of freedom subsection now cites
   `@martinMatchingPower1993` only for the qualitative
   degrees-of-freedom tradeoff and says 'We do not propose a
   threshold for $M$.' No `M > 30` text remains.
6. `report_files/figure-latex/balance-asmd-1.pdf`,
   `coverage-by-method-n-1.pdf` and
   `oversampling-relative-efficiency-1.pdf` were rewritten at
   20:56:02-03 on 2026-11-19, one to two seconds after the
   2,000-replicate RDS. The stale `efficiency-emp-se` chunk no
   longer exists. Captions are descriptive and no longer
   mention pre-fix data. The prose interpreting the figures was
   not updated at that time; see section 5, N1.
7. `analysis/report/share/MANIFEST.md` records renders through
   `tools/render.sh` at 2026-11-19 13:57 and 20:56. The render
   in this pass also succeeded on the host.

## 2. Done now

No inventory item required new work beyond what had already
been done. The work in this pass is the fix for new issue N1
(section 5): the 2,000-replicate rerun replaced the data, but
the prose written for the earlier 60-replicate run still
described a reduced run and made ordering claims that the new
data partly contradict. All edits are in
`analysis/report/report.Rmd`.

- **Replicate-count text.** Removed 'using a reduced number of
  replicates' (abstract), 'not the planned 2,000', 'Single-scenario
  values are therefore noisy' and the **TODO** (Simulation
  scale), and the 'not the planned 2,000' sentence in
  Limitations. The rendered PDF previously read '2,000
  replicates per scenario and method, not the planned 2,000'.
  Method and scenario counts in Simulation scale are now inline
  R, and the 'planned number' wording is conditional on
  `n_reps_sim == 2000`. [verified: rendered text reads
  '2,000 replicates per scenario and method, the planned
  number. All 5 methods and all 72 scenarios'.]
- **Relative-efficiency paragraph.** The sentence 'Adjusted
  simple randomization reached 1.73. That equals or exceeds
  both unadjusted matched-tuple designs' was false: unadjusted
  matched tuples with 2x oversampling have RE 1.76 at
  $R^2 = 0.5$. Likewise 'its advantage over Lin-adjusted simple
  randomization appeared only when it was also adjusted' was
  marginally false for the oversampled design (1.76 vs 1.73).
  Rewritten: adjusted simple exceeds unadjusted matched tuples
  without oversampling and is close to the oversampled design;
  matched tuples and rerandomization were similarly efficient;
  without oversampling, unadjusted matched tuples were the less
  efficient of the two (1.62 vs 1.70); matched tuples clearly
  outperformed adjusted simple only when also adjusted (1.95 vs
  1.73). Added the adjusted oversampled RE as an inline value.
  The Monte Carlo caution is now scoped to differences of a
  few hundredths, citing the 1.73 vs 1.76 comparison, rather
  than to the whole table.
  Added a `stopifnot()` guard in the `headline-16` chunk so a
  future rerun that reverses any of these orderings stops the
  render instead of leaving false prose. [verified: values from
  the RDS; render passes the guard.]
- **Abstract and Discussion summary.** 'They did not outperform
  rerandomization with a calibrated threshold' changed to
  'Their efficiency was similar to that of rerandomization with
  a calibrated threshold', because the oversampled design's
  unadjusted RE (1.76) exceeds rerandomization's (1.70) while
  the non-oversampled one's (1.62) falls below. [verified
  against the RDS.]
- **Covariate balance paragraph.** 'With $p = 8$,
  rerandomization balanced as well as or better than the
  matched designs' was false at $N = 500$ (matched_2.0 0.154
  vs rerandomization 0.168 mean max ASMD). Replaced by a new
  `balance-best` chunk that names the lowest-ASMD design per
  $(p, N)$ inline, with a `stopifnot()` guard on the
  '$p = 3$: oversampled matched tuples best at every $N$'
  claim, and an inline maximum gap between rerandomization and
  the oversampled matched design at $p = 8$ (0.041). Matched
  tuples without oversampling are not described as close to
  rerandomization at $p = 8$: at $N = 100$ they are 0.482 vs
  0.375. [verified: rendered text names rerandomization at
  $N = 100, 200$ and oversampled matched tuples at $N = 500$.]
- **Oversampling paragraph.** 'small and inconsistent at this
  replicate count' replaced by inline ranges from a new
  `oversampling-ratio` chunk: Neyman SE ratio 0.99 to 1.06,
  increasing with $R^2$ at each $N$; Lin-adjusted 0.99 to 1.01.
  [verified: monotonicity in $R^2$ checked at each $N$ from the
  RDS before writing.]
- **ADEMP bullet.** Now says 'replicates per condition'.
  [verified, rendered.]

Claims re-checked and left unchanged, because the 2,000-replicate
data support them [verified against the RDS]:

- 'in line with chance' for the maximum bias: 4 of 2,592 rows
  have |bias/MCSE| > 3 (about 7 expected under no bias) and 97
  have |z| > 2 (about 118 expected); per-method mean z ranges
  -0.09 to 0.02. No clustered excess.
- At $R^2 = 0.1$, three table rows exceed the 1.111 bound:
  oversampled matched tuples, Neyman and paired (they share a
  point estimate, 1.114), and rerandomization, Neyman (1.112).
  The excess is at most 0.003. No Monte Carlo SE was computed
  for RE, so 'within Monte Carlo error' is inferred, not
  verified; the existing 'consistent with noise' sentence was
  left as is.
- Paired coverage 0.950 overall and 0.950 at $N = 100$;
  Lin-adjusted 0.947 overall, 0.945 at $N = 100$; Neyman 0.964,
  conservative for restricted designs (model SE 12% to 16%
  above empirical). 'Close to 0.95, including at $M = 25$'
  stands.
- Application: the 4x oversampled design has the smallest F1
  empirical SE (0.071), and 'about 4%' relative uncertainty
  equals $1/\sqrt{2 \cdot 299}$.

Verification run:

- `tinytest::run_test_dir('inst/tinytest')` after
  `pkgload::load_all()`: All ok, 34 results.
- `bash tools/render.sh analysis/report/report.Rmd` (from the
  package root), three times: all succeeded. The final render
  produced `analysis/report/report.pdf` and staged
  `analysis/report/share/report-2026-11-23-1253-dee2dc1-wip.pdf`;
  the final wording was checked in `report.tex`. The first two
  renders fell in the same minute, so the second overwrote
  `report-2026-11-23-1250-dee2dc1-wip.pdf` and `MANIFEST.md`
  has a duplicate row for that file. The duplicate was left as
  written by the wrapper.

## 3. Long runs

None. The only long run in the inventory (item 1) has been
done. Its wall-clock time was not logged; the 2026-11-19 log
estimated about 2.5 hours at 8 cores. If it ever needs
repeating, from the package root:
`Rscript analysis/scripts/01_run_simulation.R 2000 8`, then
`bash tools/render.sh analysis/report/report.Rmd`. The
`stopifnot()` guards added in this pass will stop the render if
the RE or balance orderings change.

## 4. Author decisions

- **Exact or hierarchical tuple formation (item 2).** Options:
  (a) keep greedy only, as the manuscript now says; (b) add an
  exact solver for small $N$ to measure the greedy optimality
  gap; (c) remove `method = 'optimal'` from `match_tuples()`,
  since it silently aliases to greedy. Evidence: the paper's
  conclusions do not depend on the partition algorithm, but a
  user calling `'optimal'` gets greedy without warning.
  Recommendation: (c) now, either by erroring or by warning on
  `'optimal'`, and (b) only if a referee asks for the gap.
- **Randomization inference and diagnostics (item 4).**
  Options: (a) implement `ri_test()` (permutations within
  tuples); (b) drop it from `docs/analysis-plan.md`. Evidence:
  the manuscript no longer claims it, and the paired $t_{M-1}$
  intervals achieve 0.950 coverage at $M = 25$, so the
  inferential case for it in this paper is weak.
  Recommendation: (b), and note it as future work if wanted.
- **Commit the 2,000-replicate results.** The RDS files
  (`sim_raw_results.rds` is 177 MB) and the two remediation
  passes are uncommitted. `sim_raw_results.rds` is already
  tracked (it shows as modified), so committing the new
  version adds about 177 MB to history;
  `sim_performance.rds` (119 KB) is what the manuscript reads.
  Recommendation: commit `sim_performance.rds` and the code
  and prose changes; decide separately whether to untrack the
  raw file or move it to LFS before committing it.

## 5. New issues found

- **N1 (fixed, see section 2). Stale prose after the
  2,000-replicate rerun.** The rendered PDF of 2026-11-19 20:56
  said '2,000 replicates per scenario and method, not the
  planned 2,000', kept a TODO asking for that same rerun, and
  contained three ordering claims that the 2,000-replicate data
  contradict (RE of adjusted simple vs unadjusted oversampled
  matched tuples; $p = 8$ balance at $N = 500$; abstract and
  Discussion 'did not outperform rerandomization').
- **N2. Wall-clock time and hardware for the full sweep are
  not reported** (whitepaper item 14 asked for them). The
  timing was not logged, so no number was added. The author
  should record it on the next full run.
- **N3. Title.** The 2026-11-19 log already flagged that 'Theory
  and Application' oversells a paper whose theory is now short
  finite-population arguments with asymptotics credited to Bai,
  Liu and Tabord-Meehan (2024). Unchanged here; an author
  decision.

---
*Source: ~/prj/res/16-factorial-matched-random/factorialrandom/docs/pub_review_deferred_2026-11-23.md*
