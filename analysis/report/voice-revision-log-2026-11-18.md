# Voice and Humanizer Revision Log: report.Rmd
*2026-11-18 16:58 PST*

Target: `analysis/report/report.Rmd` (no supplement). Baseline:
`git show HEAD:analysis/report/report.Rmd` (HEAD `dee2dc1`; the
working tree was clean before this pass). Pass 2 applied the
voice in `~/.claude/commands/use_rgtvoice.md`. Pass 3 applied
the humanizer skill in file mode. Both passes were done in a
single set of edits, so this log covers them together. Nothing
was rendered, staged, or committed.

## Scope of edits

The count comes from the `git diff -U0` hunks. There are 52
prose hunks, each changing one to a few sentences, and 42
heading-only hunks. One of the prose hunks also contains the
ADEMP heading. The diff totals 227 insertions and 220
deletions; most of these lines are re-wraps.

## Kinds of changes

1. Long sentences split. A scanner flagged 18 prose sentences
   over 25 words in the original. Each was split into two or
   three sentences, with a connective at the head of the new
   sentence ('However,', 'That is', 'This', 'It'). After the
   pass the scanner reports one flag (line 123). It is a
   scanner artifact: the scanner does not split before a
   sentence that opens with an `@citation` key. The actual
   sentences are under 25 words.
2. Subordinate clauses removed. A grep of the prose for
   ', which', 'because', 'so that', 'although', 'though',
   'only if', 'where', and sentence-initial 'When' was run.
   Clauses split into their own sentences include 'which is
   solvable in polynomial time', 'which uses factorial
   experiments' (MOST), ', where $k$ is the number of
   treatment combinations,', 'so that each tuple member
   receives', 'because it is a sequential-enrollment
   procedure', 'because eligibility requires ...', and the
   participial riders 'creating a pool ...', 'permitting
   global ...', and 'providing excess subjects'. Some clauses
   were kept deliberately. These are the short conditional
   'When ...' sentences (lines 147, 244, 430, 506, and 518),
   'accepted only if' (line 267, which leads into its
   equation), and 'where ...' definitions after display
   equations, which the voice file prescribes. Relative
   clauses introduced by 'that' were not scanned
   exhaustively.
3. Staged run-ups removed (humanizer section 4). These were
   'We begin with three motivating examples.', 'The root cause
   is as follows.', and 'We begin with a summary of what this
   manuscript develops.' In each case the content of the next
   sentence was kept and the announcement was dropped or
   folded in. The count 'three' was kept ('Three examples
   are ...').
4. 'We note that' thinned. The original used it 11 times,
   which is heavy even for this voice. Eight occurrences were
   removed, with 'However,' kept where the original had a
   contrast. Three were kept as signposts on scope caveats.
   These are the k > 2 demand (Trial-ready cohorts), the
   greedy-only scope (Tuple formation), and the
   matched-only ranges and the 'central claimed contribution'
   (Provenance and Summary of findings). That makes four
   places in all; three use 'We note that' and one uses
   'We note, however,'.
5. Stock or inflated words. 'represents' became 'is', 'key
   efficiency result' became 'main efficiency result',
   'demonstrates' became 'shows', 'approximately' became
   'about', and 'rich baseline phenotyping' became 'detailed'.
   The intensifier 'essential' was dropped before
   'prerequisites'.
6. Empty connectives. 'It is clear that', 'On the other hand'
   (no contrast present) and 'in fact' were removed.
7. A grammar fix: 'neither factor, factor A alone, factor B
   alone, or both' became '... and both'.
8. Headings were changed to sentence case only, keeping proper
   nouns (Alzheimer's, TRC-PAD, AHEAD 3-45, DIAN-TU, EPAD,
   ADEMP). Prose cross-mentions were recased to match:
   '(Tuple formation)' and '(Efficiency analysis)'. No
   section was renamed, added, removed, or reordered. The
   nesting under '# References' was left as it was.

## Passages deliberately left alone

- The abstract (YAML), all code chunks, chunk options, figure
  and table captions (these sit inside chunk headers or
  `caption =` arguments), display and inline math, inline R,
  citation keys, and `\@ref()` cross-references.
- Run-in bold labels were kept. These are **Mahalanobis
  distance.**, **Prognostic score distance.**, **Sequential
  greedy.**, **Provenance.**, **Monte Carlo standard
  error.**, **TODO (before submission):**, **Verdict:**, and
  the Design parameters and Future work list labels. They
  label content rather than decorate it. The instructions
  also forbid converting the lists to prose.
- The Provenance note, the TODO, and the remediation narration
  in Results, Application, and ADEMP. Humanizer would cut much
  of this as writing about the document (section 25). However,
  it carries disclosures and hedges, and the no-delete rule
  takes precedence. These passages were reworded and split
  only.
- Claims that the referee report flags were reworded at most,
  and none was strengthened. These include 'four separate
  experiments', 'bipartite graph', 'unbiased', the RE limit
  $k/(k-1) \cdot M$, 'a fifth comparator', and
  'truncated normal'. The citation of
  `01_run_simulation.R:278-279` was left as is, although the
  lines are actually 297-298. The Provenance note still names
  only the `tuple_balance()` and `regression_adjustment()`
  fixes; that is referee comment M1, left for the remediation
  step.
- The hedges 'It appears that', 'plausible', 'is expected
  to', 'may be small', and 'what we regard as' were kept.
- No sentences were merged on grounds of repetition, and no
  sentence carrying information was deleted.

## Verification

Script:
`/private/tmp/claude-501/-Users-zenn-prj-res-01-fisher-exact-rx2-fisherexacttestrx2-analysis-report/be4b228f-9df5-41a4-9039-7499d4f23926/scratchpad/16-factorial/verify.py`.
It compares the working file against `git show HEAD:` for the
items below. The sentence-length scanner is `longsent.py` in
the same directory.

After a first draft, a meaning check over
`git diff --word-diff` found four places where a split had
shifted claim strength. All four were restored:

- 'so as to minimize' had become a flat assertion of
  optimality. It now reads 'The aim is to minimize'.
- 'offered as a heuristic consistent with' now reads 'It is
  offered as consistent with'.
- 'rather than as an established empirical finding' now reads
  'It should not be read as ...'.
- 'emerging' is again attached to the treatment structures.

| Item | Count | Result |
|---|---|---|
| YAML header, including the abstract | 1 | identical |
| Code-chunk lines | 304 | identical |
| Numbers (whole file) | 367 | identical |
| `@` citation keys | 36 | identical |
| Inline R expressions | 14 | identical |
| Display equations | 13 | identical |
| Inline math spans | 135 | identical |
| `\@ref()` tokens | 3 | identical |
| Other backticked spans | 55 | identical |
| Spelled-out number words | 32 | identical |
| Headings, case-insensitive | 58 | identical |

All of these are multisets. No em dashes, en dashes, or curly
quotes remain in the prose. Result: PASS.

Not checked: the document was not rendered, so LaTeX and knitr
behavior after the re-wrap is unverified. The re-wrapped lines
around inline R (Provenance) keep each inline expression on its
own line.
