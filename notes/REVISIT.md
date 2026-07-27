# Revisit list — upstream (tplyr2 / clinify) follow-ups

**STATUS (2026-07-25): CLOSED OUT.** Every tplyr2 finding is resolved and applied (native
group_shift, negative-zero, byfactor/ordering, stats_as_columns, denoms_by). Every clinify finding
is resolved and either applied (merge control, per-line title align, clin_row_height) or consciously
NOT adopted with a documented reason (clin_header_pad — see #97). All 30 tables build clean and match
the references (bar documented deliberate divergences); the rebuild is in **draft PR #8**
(`nextgen → master`). The only work remaining is **version-pin housekeeping** (see the bottom):
repin tplyr2/clinify to release tags once they cut to `main`, then `renv::snapshot`.

Originally: a log of things left with a workaround, to be simplified once each dependency landed.

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

### clinify #97 — regulatory-fidelity vertical rendering  [clin_row_height ADOPTED; clin_header_pad NOT adopted — FINAL]
Both verbs shipped on `development` (`clin_row_height` in #100; `clin_header_pad` in #100→#102).
- **`clin_row_height` ADOPTED (2026-07-25):** removed the per-styler `hrule("atleast") + height_all()`
  from all three stylers; one central `clin_row_height(body=15.35, title=11.4, footnote=11.4,
  rule="atleast", unit="pt")` in `add_titles_footnotes()`. All 30 rebuilt, ZERO change vs baseline.
- **`clin_header_pad` NOT adopted (tested #100 AND the #102 redesign; both reverted).** #102 made
  `above`/`below` space EVERY header row (fixing the uniform multi-row case — 28/30 then matched). But
  two tables need NON-uniform per-row header padding that clin_header_pad's uniform args can't express:
  **14-1.03** (34pt above row 2 only, spanner→label gap) and **14-3.10** (21pt above row 1 only, to land
  the header rule at the reference y). Those fall to `flextable::padding(i=)`, but `clin_header_pad` is
  per-table config applied at `finish_table_()` *after* the styler, so a central call CLOBBERS the
  per-row `flextable::padding(i=)` (14-1.03 → 2.761%, 14-3.10 → 2.443%). So the pilot KEEPS the uniform
  buffer as `flextable::padding(padding.top=18, padding.bottom=4, part="header")` in the house styler
  (composes with the per-row exceptions). Reported both findings on #97 (comments); no further action
  our side — clin_row_height is the win, header padding stays styler-level.

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

## Pagination / page-break layout review (2026-07-25)
Reviewed all 11 multi-page tables (fan-out workflow) for defects the pagination-agnostic content
check misses. 14-6.03 came back clean (validated the earlier fix + the review). Three classes found:
- **① orphaned spacer gap** — 14-6.03 FIXED (5-blank excess → 1). Others (14-7.01/.02/.04, 14-5.01)
  are single-blank spacer ORPHANS (already 1 row) that land at a page top — no "reduce-but-keep" fix
  like 14-6.03; only remove-the-spacer-entirely (tested on 14-7.02: content-identical, gap gone, but
  DENSER — diverges from the reference's block spacing) or clinify page-control. NOT applied (more
  aggressive than the 14-6.03 precedent; the vitals blocks are still readable via their Position/
  Measure labels, but AE/CM spacers are load-bearing group separators). Bundle with class ③.
- **② spurious empty trailing page** — 14-6.01 page 18 (header+footnote, no data) FIXED: strip
  trailing all-blank rows (18→17 pages, content unchanged). Same class as 14-2.01's dropped page-4.
- **③ split-block / header-widow** — DEFERRED (user's call): a param/SOC/shift block splits across a
  page with continuation rows or the section/param label stranded unlabeled. Pervasive in 14-6.04
  (~15 transitions), plus 14-6.01 (PROTEIN/ERY-MCH headers), 14-5.01 (INVESTIGATIONS SOC + a torn
  wrapped PT label), 14-6.02 (HEMATOLOGY), 14-6.05, 14-7.x, 14-3.13. This is fundamental table
  pagination (rows break where they fall; group labels don't repeat) — needs a clinify page-control
  feature (keep-with-next / repeat-group-label / suppress-top-of-page-blank; cf. clin_group_by /
  clin_auto_page, issues #90/#91) and LIKELY the legacy RTFs paginate the same way (so it may be
  reference-consistent, not a divergence). Full per-table findings: workflow wf_28a53ca3-1e5 journal.

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

## Version pins / renv snapshot — DONE (2026-07-26)
- **tplyr2** repinned to `atorus-research/tplyr2@main` (**0.2.0**, SHA 8e94e12) — carries #33/#34/#35/#36 (PR #38)
  and #37 (PR #39) plus the earlier shift/ordering/byfactor fixes.
- **clinify** repinned to `atorus-research/clinify@main` (**0.4.0**, SHA 5721abf), the release line. NOTE: main
  and development had DIVERGED but only in NEWS.md / cran-comments.md / .gitignore — the R/ source is identical
  (main got #109/#110 via the release merge #103), so main is authoritative. development is 2 doc-commits ahead
  (a NEWS tidy not forward-merged) — a clinify-side housekeeping nit, not a code difference.
- **renv.lock** snapshotted against the verified library. (renv still reports 5 legacy packages "used but not
  installed" — huxtable/pharmaRTF/vcdExtra/plyr/assertthat — but none are referenced by the nextgen code; it's
  scan-noise from the legacy `programs/` files, pre-existing and harmless to the rebuild.)

## Upstream fixes adopted into the pilot (2026-07-26)
All 8 filed issues were addressed upstream (tplyr2 0.2.0 = PR #38 for #33–#36 + PR #39 for #37; clinify 0.4.0 =
PR #109 for #104 + PR #110 for #105; #106 closed won't-do — Word owns pagination here so page-relative blanking
is unimplementable). Features were adopted only where the native call reproduces the verified output
BYTE-IDENTICALLY (fanned-out verify-and-revert; independent full-30 sweep afterward = 0/30 changed):

**Adopted**
- **as_display() (#36)** ×12 — replaced the `grep("^res")` recast in efficacy.R, ae.R, 14-1.01, 14-2.01, 14-3.13,
  14-4.01, 14-6.01/.02/.03/.04/.05/.06.
- **coerce_character (#104)** ×16 — replaced the `final[] <- lapply(...)` NA/character boilerplate.
- **clin_spanner_rule (#105)** ×7 — replaced the hand-rolled header `hline` (dashed in 14-6.01) across 14-4.01
  and 14-6.01/.02/.03/.04/.05/.06.
- **total_group()** ×3 — 14-1.01, 14-1.02, 14-2.01 (dropped the `bind_rows(…TRT01P="Total")` duplication).
- **n_records (#34)** ×2 — 14-6.01, 14-7.02 (14-7.01/.03 pre-filter to non-missing, so n/a there).
- **missing_count (#33)** ×1 — 14-1.02 (dropped the hardcoded all-zero Missing row).
- **assoc_test (#37)** ×2 — the Fisher lab tables 14-6.02, 14-6.03.
- **14-7.04 CM rewrite** — now native `group_count(c("CMCLAS","CMDECOD"), distinct_by=, pop_data=, total_row=)`.
  **14-7.04 now uses tplyr2 → only 3 tables remain non-tplyr2** (14-1.03 crosstab, 14-3.11 MMRM, 14-3.10 wide).

**Later adoptions (2026-07-27) — after further upstream fixes**
- **AE tables 14-5.01/.02 → native pairwise assoc_test; `fisher_ae` RETIRED.** The blockers I flagged were
  fixed upstream in sequence: pairwise mode + denominator access (#40/PR#43), character-return display for the
  `*`/`>.99`/blank formatting (#47/PR#48), nested-count-layer support (#49/PR#50), and an empty/sparse-reference
  patch (PR#50 follow-up). `R/ae.R` now runs `group_count(..., assoc_test = assoc_test(fn = ae_p, reference,
  comparisons, total_row = TRUE))`; `ae_p` supplies the display verbatim. Removed `fisher_ae`, its manual 2x2,
  and the per-p `N` plumbing. Verified byte-identical (14-5.01 16pp, 14-5.02 sparse serious AEs) + full 30-sweep.
- **num_fmt → apply_formats (#41/PR#42):** once `na=`/`width=`/`pad=` shipped, `num_fmt()` became a thin wrapper
  over `apply_formats()` (byte-identical, 88/88 across call-site shapes); all 31 call sites unchanged.
- tplyr2 repinned to `@main` (post-#50) + `renv::snapshot` each round.

**CMH / Fisher adoption round (2026-07-27)** — the omnibus `assoc_test` fn receives the by-group's FULL source
subset (stratum columns included) and can return a preformatted string, so these became byte-identical:
- **14-3.13** CIBIC+ CMH, **14-6.05** lab shift CMH, **14-6.06** Hy's-law shift CMH — `cmh_pval()` retired;
  `coin::cmh_test(... | SITEGR1/BNRIND/BASE)` now runs inside `assoc_test`, guards and formats preserved.
- **14-1.02** completion-status Fisher p — now from `assoc_test` on the COMP_STAT count layer.
- `cmh_pval` and `fisher_ae` are gone; `assoc_test` supplies p-values in **7** programs. All 30 byte-identical.

**Still kept manual**
- **14-2.01** (demographics) — attempted and reverted. Two blockers: (a) its 9 CONTINUOUS `aov_p_str` p-values
  can't move (`assoc_test` has no `group_desc` support → filed **tplyr2 #51**); (b) placement — see the bug below.
  So the bespoke single-p column + `stamp()` machinery survive regardless, and `chi_p_str`/`aov_p_str` stay.
- **denom_row (#35)** for 14-6.04/.05/.06 — emitted denominator row didn't match the pilot's exact format.
- `fish_p_str` (2 calls in 14-1.02: AE / lack-of-efficacy reason) — separate single tests, not a count layer.

### Two tplyr2 bugs found during the 14-2.01 attempt — NOT yet filed (both silent, wrong numbers)
1. **`total_group()` rows leak into `assoc_test`'s fn `.data`.** The subset carries the synthetic "Total" arm
   duplicates (508 rows vs 254), so a test computed over `.data` double-counts and returns a WRONG p with no
   error (AGEGR1: 0.3347 vs the correct 0.1439). Caller must filter `.data$TRT01P != "Total"` — undocumented.
2. **Omnibus p is placed before the sort.** `merge_assoc_column()` does `wide[1L, pval1 := ...]` on the PRE-sort
   (alphabetical/dcast) frame, and `order_count_method="byfactor"` reorders afterwards — so with any
   non-alphabetical ordering the p lands on an arbitrary category row (AGEGR1 → "65-80 yrs" instead of
   "<65 yrs"; BMIBLGR1 → "25-<30" instead of "<25"). Pixel diff barely moves (0.07%); only the content check
   caught it. Should place on the first row of the FINAL display order.

## Post-cleanup redundancy review + upstream issues filed (2026-07-25)
The final code review surfaced recurring hand-written boilerplate. Each candidate was VERIFIED against the
INSTALLED packages before filing (drafting agents ran reprexes), which debunked three as already-supported.
Net: **8 issues filed** (each with a runnable reprex), **2 comments** on existing issues, **3 adopt-in-pilot**.
NONE block the pilot — all current code is verified-correct; the adopt-in-pilot items are optional simplifications.

### Filed — genuine gaps
tplyr2:
- **#33** `group_count`: `missing_count` drops the Missing row when the count is 0 (and blanks cells for columns
  with no missings) — blocks a clean 14-1.02 Missing row.
- **#34** `group_desc`: option for a record-count `n` (assessed) vs the current non-missing `n` — 14-7.0x, 14-6.01.
- **#35** `group_shift`: emit the per-baseline-group denominator ("n") row for `shift_denom="column"` — 14-6.04/.05/.06.
- **#36** add a display-ready extraction helper (`as_display()` / drop the internal `ord_*` cols) — ~11 tables.
- **#37** (scope question) association-test p-value column (Fisher/CMH) for count & shift layers — 14-5.01, 14-6.0x.
clinify:
- **#104** `clintable()`: render columns verbatim (accept non-character cols). NA-blank already works via flextable.
- **#105** `clin_column_headers`: option to draw a rule (underline) beneath a spanner header row — 14-6.0x, 14-4.01.
- **#106** in-body row-group stub suppression + spacer rows that repeat the label on continuation pages
  (distinct from `clin_group_by`'s banner+page-break model) — vitals + AE; **also resolves pagination class-③**.

### Commented on existing (not duplicated)
- clinify **#15** (auto-insert header N's) — added the pilot's hand-rolled big-N reprex.
- clinify **#7** (titles/footnotes from a data source) — added the `read_titles` / titles.xlsx reprex (also covers #6).

### Already supported — NOT filed (ADOPTED 2026-07-26 where byte-identical; see "Upstream fixes adopted" above)
- **`f_str` unpadded percent** — the "(8%)" blocker was a FALSE premise: a single-char token
  `f_str("xxx (x%)", "distinct_n", "distinct_pct")` already emits unpadded/variable-width percents. So **14-7.04
  (CM) CAN use `group_count(c("CMCLAS","CMDECOD"), distinct_by="USUBJID", ...)`** — the "4 tables don't use
  tplyr2" count is really **3 structural** (14-1.03 crosstab, 14-3.11 MMRM, 14-3.10 wide) **+ 14-7.04 adoptable**.
- **`total_group()`** — `tplyr_spec(cols="TRT01P", total_groups = list(total_group("TRT01P")))` replaces the manual
  `bind_rows(adsl0, mutate(adsl0, TRT01P="Total"))` in 14-1.01/.02 and 14-2.01 (applies to pop_data too).
- **`apply_formats()`** — the exported standalone f_str formatter (handles multi-var + lt/gt thresholds) replaces
  R/helpers.R `num_fmt`/`fish_p_str`/`chi_p_str` for external model/p-value formatting.

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
