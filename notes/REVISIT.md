# Revisit list — deferred work waiting on upstream (tplyr2 / clinify) or a merge

Things intentionally left with a workaround, to be simplified once the dependency lands.
Check this file when a tplyr2/clinify release or PR merges.

## PR outcomes (assessed 2026-07-24 — tplyr2 PR #21 merged to main; clinify PR #96 on gh_issue_95)

### tplyr2 #18/#28/#30/#31/#32 — native group_shift lab-shift tables  [RESOLVED & APPLIED]
Fully closed. `shift_denom="column"` (#21) + by-group scoping (#30/#28) + `zero_count_display` honored
by group_shift (#32/#31) => **14-6.04/.05/.06 now use native `group_shift(c(row=, column=),
shift_denom="column")`**, dropping the crossed `cols=c("TRTP","BNRIND")` count + `denoms_by` workaround.
Installed tplyr2 from PR #32 branch `issue-31-shift-zero-display` (open, targets main — repin to main
once merged). Verified: 14-6.04 difflines=9, 14-6.05=0, 14-6.06=2 — all baseline, full 30-sweep clean.
(#32 also makes pct_lt/pct_gt work on shift, for free.)

### tplyr2 #29/#30 — negative-zero in mean formatting  [RESOLVED & applied]
PR #30 normalizes a rounded `-0` to `0`. Removed the local `fix_neg_zero()` gsub from 14-6.01;
rebuilt, content-identical (difflines=8, unchanged).

### tplyr2 #19 — `denoms_by` replaces (not augments)  [RESOLVED, docs-only — no code change]
PR #21 documented the replace-not-augment behavior (no behavior change, per the issue's preferred
resolution). Our `denoms_by = c("PARAM","VISIT","TRTP","BNRIND")` calls remain correct as-is.

### tplyr2 #20 — `stats_as_columns=TRUE` + `by`  [RESOLVED & validated — no program change warranted]
PR #21 fixed it (verified: with `by=c("PARAM","AVISITN")` it now keeps by-groups as rows and emits
one `"<arm> | <stat>"` column per treatment×stat). BUT the tables that "worked around" it don't
benefit from switching: 14-6.01 uses single-stat builds (no reshape), and 14-7.01/.02/.03 put
TREATMENT in a `by` var (dummy `cols="ALL"`), so `pivot_wider(names_from=stat)` stays cleaner than
parsing `"x | stat"` labels. Left as-is; fix is validated for future use.

## Filed FRs — outcomes

### clinify #97 — regulatory-fidelity vertical rendering  [clin_row_height ADOPTED; clin_header_pad NOT adoptable]
Both verbs shipped (`clin_row_height`, `clin_header_pad`) on `development`.
- **`clin_row_height` (item 1) ADOPTED (2026-07-25):** removed the per-styler `hrule("atleast") +
  height_all()` from all three stylers (`cdisc_table_default`/`titles`/`footnotes`); one central
  `clin_row_height(body=15.35, title=11.4, footnote=11.4, rule="atleast", unit="pt")` in
  `add_titles_footnotes()`. All 30 rebuilt, ZERO change vs baseline. It also reaches the title/footnote
  blocks natively (which the option stylers already did for us, but this is the intended path).
- **`clin_header_pad` (items 2/3) NOT adoptable — reverted.** Its `above=` pads only the TOP header
  edge, but our multi-row spanned headers (14-1.03, 14-3.10, 14-4.01, 14-6.05, shift tables) rely on
  `padding.top` on EVERY header row. Adopting it regressed all of those — **14-1.03 went 2.408→2.523%
  (over the 2.5% gate)**, 14-3.10 0.336→2.001%, 14-4.01 1.122→1.679%, 14-6.05 content +3. So we KEEP
  the `padding.top=18`-on-all-header-rows buffer (works for all 30). **Finding for the maintainer:**
  `clin_header_pad` only spaces the 3 outer edges of the header block; it can't reproduce per-header-row
  spacing that multi-row spanned headers need. Report on #97.

### clinify #98 — table alignment + per-line title alignment  [title-align APPLIED]
PR #99 fixed the `set_table_properties` clobber AND added `clin_table_align()` + per-line title
`align`. **Applied:** R/titles.R now uses `clin_add_titles(align=)` / `clin_add_footnotes(align=)`
(a per-line `"left"/NA` vector) instead of the duplicate-text trick — dropped `.line_vec`'s
`c(t1,t1)` branch; all 30 rebuilt, titles render identically. Table alignment: our explicit
`set_table_properties(align="center")` now sticks post-fix (was coincidentally preserved by the
clobber-to-center before), left as-is — `clin_table_align()` is the cleaner API if we ever want a
non-center table, no output change.

### clinify #95 — header-cell auto-merge per-row control  [RESOLVED & applied]
clinify PR #96 (v0.4.0, branch gh_issue_95; merged to `development`) added
`clin_column_headers(merge=)`. Applied `merge = "spanners"` in **14-6.05 and 14-6.06** — removed the
`merge_none(part="header")` + manual `merge_at()` blocks. Both re-verified (14-6.05 content-identical;
14-6.06 difflines=2 = the intended TRANSHY p only).

### clinify — titles/footnotes from a spreadsheet  [NOT filed — deferred]
`notes/feature-requests/clinify-regulatory-fidelity.md` item 2 (our `R/titles.R`). Not filed (not in
the batch I filed); offer to file if the maintainer wants it.

## Internal (our own code / harness — no upstream dependency)

### R/efficacy.R — generalize build_efficacy_table  [DONE 2026-07-25]
Added optional params `extra_filter` (14-3.07 completers), `anl01`, `derive_chg`, `blocks`
(list of list(var=,label=,visit=)), and `use_base`. 14-3.07 and 14-3.12 are now thin
build_efficacy_table() callers (inlined bodies removed). Verified: all 8 existing callers
(14-3.01..09) unchanged pixel, 14-3.07 content-identical, 14-3.12 pixel 1.074 — zero regression.

### verify/fidelity.py — registration fragility on uniform grids  [DONE 2026-07-25]
`register()`'s `best_shift` now, among shifts within `tol` (0.5%) of the best cross-correlation
score, prefers the SMALLEST displacement — so a periodic grid's full-row-off secondary peak no
longer wins over the true near-zero shift. No-op for tables with a single clear peak. Re-swept all
30: 15/16 pixel identical, 14-1.03 +0.01% (2.398->2.408, still PASS), 14/14 content unchanged.

## Housekeeping / version pins
- **CURRENT INSTALLS (2026-07-25 pm):** tplyr2 `atorus-research/tplyr2@issue-31-shift-zero-display`
  (PR #32 branch = main + #31/#32 fix; **repin to `@main` once #32 merges**); clinify
  `atorus-research/clinify@gh_issue_98` (v0.4.0, has #95/#96 + #99/#100 = merge control, align,
  clin_table_align, clin_row_height, per-line title align). Reinstall from released tags once cut.
- **tplyr2 (earlier note):** was installed from **GitHub `atorus-research/tplyr2@main`**
  (renv record is a clean github ref, not a local path). main now contains everything we need:
  PR #21 (#18/#19/#20), #23 (group_desc by-group order), #25 (group_count by-group order +
  total/missing with by, #24), and **#27 (byfactor + nested target ordering, #16 — supersedes the
  stale #17)**. The old local combined branch is DELETED (no longer needed). Full 30-table
  rebuild+sweep on this main = zero change vs baseline.
- **clinify (installed 2026-07-24):** v0.4.0 from local clone on `gh_issue_95` (PR #96, base
  `development`, still OPEN). Reinstall from a released tag/main once #96 merges.
- renv.lock still out of sync (not snapshotted); pin/snapshot once clinify #96 lands on a release.

## Resolved (kept for provenance — no action)
- tplyr2 #13/#14 (count column order + display quirks) — merged in PR #15, incorporated;
  see `notes/feature-requests/tplyr2-count-column-order.md`, `...count-display-quirks.md`.
- tplyr2 #16 (byfactor + nested target ordering) — fixed in **PR #27** (v2 on main; the old #17
  went stale and was superseded). Verified: all 8 byfactor tables order correctly on main.
- tplyr2 #20 (stats_as_columns) — fixed in #21; the by-group ROW-ordering follow-up landed in
  **#23** (group_desc) and **#25/#24** (group_count + total/missing with a by var).

## Optional cleanup enabled by #23/#25/#27 (not done — low value, regression risk)
The by-group ROW ordering is now native (factor levels). Our wide-desc / count tables still
re-sort rows manually in assembly (e.g. `arrange(match(AVISITN, VISN))` in t_14_6_01.R, and the
sort steps in t_14_7_0{1,2,3}.R). These are now REDUNDANT but harmless (re-sorting already-ordered
data). Could be removed to simplify, but they gate no correctness and removing them risks a subtle
regression on passing tables. All 30 verified unchanged on latest main (2026-07-25) with the manual
sorts left in place.
