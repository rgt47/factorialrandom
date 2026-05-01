# Morris et al. (2019) ADEMP Audit: 16-factorial-matched-random
*2026-04-17 09:02 PDT*

## Scope

Files audited:

- `analysis/scripts/01_run_simulation.R`
- `analysis/report/report.Rmd`
- DGM / utility functions referenced therein

## ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Partial | goals stated; no explicit ADEMP header |
| DGMs documented | Met | factorial grid explicit in script |
| Factors varied factorially | Met | 90-cell factorial design |
| Estimand defined with true value | Met | treatment effect input to DGM |
| Methods justified | Met | matched-random vs factorial methods |
| Performance measures justified | Partial | listed but not mapped explicitly to aims |
| n_sim stated | Met | `n_replicates = 2000` at `01_run_simulation.R` near line 270 |
| n_sim justified via MCSE | Not met | no MCSE-based derivation |
| MCSE reported per metric | Met | `mcse_bias`, `mcse_coverage`, `mcse_rejection` computed at `01_run_simulation.R:222-228` |
| Seed set once | Met | `set.seed(20260312)` called once at line 270 |
| RNG states stored | Not met | no per-replicate `.Random.seed` capture |
| Paired comparisons | Met | same data fed to methods per rep |
| Reproducibility | Partial | seed once; RNGkind not pinned; Rmd results tables have TODO placeholders |

## Overall verdict

**Mostly compliant.**

## Gaps

- `n_replicates = 2000` not justified via a Monte Carlo SE target.
- Report has TODO placeholders around `report.Rmd:605-649` — the
  compliance narrative for MCSE is already computed in the script but
  not surfaced in the rendered tables.
- `RNGkind("L'Ecuyer-CMRG")` not pinned; exact reproducibility across R
  versions is fragile.
- No per-replicate `.Random.seed` snapshot stored.

## Remediation plan

1. Add an n_sim justification block at the top of
   `01_run_simulation.R` deriving `n_replicates` from a target MCSE.
   For rejection-rate MCSE ≤ 0.01 at p = 0.5, need n ≥ 2500. Document.
2. Wire the MCSE columns (`mcse_bias`, `mcse_coverage`, `mcse_rejection`)
   into the result tables in `report.Rmd:605-649`, replacing TODOs.
3. Pin `RNGkind("L'Ecuyer-CMRG")` immediately before the `set.seed()`
   call in `01_run_simulation.R`.
4. Store `.Random.seed` per replicate in an RDS under
   `analysis/data/derived_data/`.
5. Add explicit ADEMP headings to `report.Rmd` Methods section and
   cite Morris Table 6 for the MCSE formulae.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate
statistical methods. Stat Med 2019;38:2074-2102. doi:10.1002/sim.8086

---
*Source: ~/prj/res/16-factorial-matched-random/factorialrandom/docs/morris-audit-2026-04-17.md*
