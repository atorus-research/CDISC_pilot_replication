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

**14-2.01 ADOPTED (2026-07-27)** — the last table with hand-rolled statistics. Unblocked by **#51/PR #52**
(`assoc_test` on `group_desc`) together with the two correctness fixes that PR carried, **#53** (`total_group`
duplicates leaked into the fn's `.data`, doubling every subject) and **#54** (omnibus p placed before the
byfactor reorder). `desc_block()`/`count_block()` now each carry an `assoc_test`; `stamp()` lost its `p_at`
argument (the p travels with the block, so `bind_rows` places both p-values of a composite characteristic);
the separate `adsl3` 3-arm frame is gone. All 12 p-values verified unchanged (5 chi-square, 7 ANOVA) — including
AGEGR1 0.1439 and SEX 0.1409, the two #53/#54 had been silently corrupting. **`aov_p_str` and `chi_p_str`
retired from R/helpers.R.**

**Header (N=) from `tplyr_header_n()` (2026-07-27)** — `R/ae.R` and 14-7.04 now read the header N off the build
instead of counting ADSL a second time, so the displayed denominator is provably the one used for the percentages
(the two-sources-of-truth risk behind the old hardcoded `c(86,84,84)`). Only these two: the other 16 N-deriving
programs don't pass `pop_data` (they scope with `denoms_by`/layer totals) and adding it just to read the N back
would change their denominators; 14-1.02's builds live inside a per-block helper, so hoisting would be worse code.

**Still kept manual**
- `fish_p_str` (2 calls in 14-1.02: AE / lack-of-efficacy reason) — single 2×2 tests, not count layers.
- `arm_label()` character-width header wrapping (13 files) — filed as **clinify #112**.
- House-style header padding via `flextable::padding(part="header")` rather than `clin_header_pad` — the verb
  clobbers the per-row exceptions 14-1.03/14-3.10 need. Filed as **clinify #113** (verified reprex).
- Hand-rolled shift `n` rows in 14-6.04/.05/.06 — `denom_row` renders literal `"NA"` for an absent baseline group
  and inherits the `n_counts` format width. Diagnosed and filed as **tplyr2 #55**.

## Unmerged upstream PRs — all evaluated 2026-07-27 (pilot verified byte-identical on each)
- **tplyr2 #56** (→#55 denom_row): verified at SHA d7eb9fc. Absent baseline group now `0` (was literal `"NA"`),
  and `denom_row_format = f_str("xx","n")` gives the reference's 2-char field. **All three shift tables adopt it
  byte-identically** (14-6.04 28pp, 14-6.05 5pp, 14-6.06 1pp incl. the corrected TRANSHY p) — held as
  `scratchpad/denom_row_adoption.patch`. **➜ apply on merge**; retires the last hand-rolled n-row.
- **clinify #114 (→#107), #115 (→#101), #116 (→#113)**: tested as a merged set (throwaway worktree off
  `development`). **Full 30-sweep byte-identical** — notably #115 restores 9pt header padding after
  `clin_column_headers()`, but the house styler's `flextable::padding(part="header")` runs later and wins, so
  the pilot is unaffected. #107 needs no pilot change (no label-stripping remains). **#116 unblocks adopting
  `clin_header_pad` centrally** — verified `above=c(18,34)` and `above=18, rows=1` both preserve the per-row
  exception. **➜ optional adoption on merge**: replace the styler's raw `flextable::padding` and the two
  per-row hacks (14-1.03 row 2 = 34pt, 14-3.10 rows 1/2 = 21/2pt) with the native verb.

## clinify #112 (header wrapping) — CLOSED won't-do, and the reasoning corrects our premise
Not a guess at the rendered width: `str_wrap(width = 10)` is a **deliberate override** of it. The arm spanner
merges across three 0.45in columns (~1.29in ≈ 15 chars at Courier New 10pt), so Word legitimately breaks
`Xanomeline Low | Dose`; the reference breaks `Xanomeline | Low Dose`. **A width-derived mechanism cannot produce
a wrap narrower than the width**, so neither proposed shape could work. The constant is load-bearing and stays
in the program next to its comment. (`arm_label()` keeps the explicit width; only `sprintf` boilerplate could
ever be factored out, and that's a project helper, not a clinify change.)

## Harness gap closed (2026-07-27) — from the #112 investigation
A visibly wrong header wrap **passed every gate**: pixel moved 2.408% → 2.417% (still PASS) and
`page_text_set()` strips all whitespace, so a re-wrap is invisible to it. Added `header_compare()` +
`compare_table.py --header`, validated by reproducing that exact experiment (fails unwrapped, passes restored),
and wired into the build-vs-baseline sweep (`hdr=ok` for all 30). Also fixed `header_rule_pt()`: the rule
threshold was page-relative (>0.9 × page width) but these tables span only ~82% of a landscape page, so the
rule was never detected and the cut fell back to a fixed page fraction, landing mid-header — now measured
against the table's own ink extent. Against the reference RTFs 21/30 match exactly; the other 9 carry
pre-existing multi-line-header stacking divergences, so reference-mode is informational, not a gate.

### clinify #117 — header pitch in clin_row_height() (filed 2026-07-27)
The unshipped third of #97: row pitch → `clin_row_height`, header buffer/gap → `clin_header_pad` (+#113), but
**header pitch has no clinify lever**. Verified `clin_row_height(body=15.35,…)` moves body rowheights
0.25→0.213 while the header stays at `line_spacing=1`; only a raw `flextable::line_spacing(part="header")`
touches it (14-4.01 uses `space=0.75`). Motivated by the new header check: **9/30 tables stack multi-line
header cells differently from the reference** (e.g. 14-1.02's `p-value` on header line 5 vs the reference's 6) —
a class that was invisible to the pixel and text-set gates. Left the height-vs-leading design choice to the
maintainer (leading is what the multi-line case actually needs).

### clinify #117 RESOLVED in PR #118 (merged to development) — adopted, held for release
Both levers shipped: `clin_row_height(header =, header_leading =)`. The maintainer measured the design question
and confirmed a header *height* with `rule="atleast"` is only a floor (a 3-line cell grows past it) — only
*leading* reaches the multi-line case, so both ship with docs pointing at the right one. Verified
`header_leading=0.75` survives the pilot's styler `line_spacing(space=1, part="all")` (applied at
`finish_table_()`, after the option styler). **14-4.01 adopted it byte-identically**: its raw
`flextable::line_spacing(part="header")` is gone and `spanned_default()` is now a plain passthrough. Held as
`scratchpad/header_leading_adoption.patch` — clinify `main` is 13 commits behind `development`, so the pilot's
`@main` pin has none of #114/#115/#116/#118 yet. **➜ apply once development ships to main.**
NOTE: `header_leading` moves leading, NOT where header lines break — it does NOT address the 9/30 stacking
divergence (filed separately as #120); citing that on #117 was the wrong evidence for this lever.

### Filed 2026-07-27 from the #118 adoption
- **clinify #119** — `clin_row_height()` REPLACES on a second call instead of merging unspecified args.
  Reproduced on pinned main (predates #118): chaining `clin_row_height(title=)` after
  `clin_row_height(body=15.35)` silently reverts body to the flextable default. Cost a 3.189% pixel diff with
  no error; the pilot now threads `header_leading` through the single central call in `add_titles_footnotes()`.
- **clinify #120** — multi-line header cells stack one line higher than the reference in **9/30** tables
  (14-1.02, 14-5.01/.02, 14-6.02/.04/.05, 14-7.01/.02/.04). Filed at the maintainer's invitation on #118.
  Cosmetic and pre-existing; may end up documented as a deliberate divergence.

### clinify #120 RESOLVED — and 7 of the 9 were OUR harness lying (2026-07-27)
The maintainer's investigation found **no clinify defect**. Only outcome upstream was a docs change
(PR #122: the `merge` argument now says what merging a row of repeated labels actually looks like).
Everything else was on our side. Header gate: **9 failing → 1**, and the one remaining is documented.

**`header_compare()` had three defects — all fixed in `verify/fidelity.py`:**
1. **Shared absolute cut.** `cut = min(header_rule_pt(ref), header_rule_pt(cand))` cut both documents at
   one absolute y. The references start content at y=36/54/**72** depending on the legacy program; our
   DOCX is always 36. So one cut sliced the two at different *depths*. For 14-7.01 it landed **above the
   reference's entire column-header block** → ref 3 lines vs cand 6, a pure artifact. Now each document is
   cut at its own rule. This alone fixed 7 of the 9. See divergences.md on the y-origin.
2. **`zip` + tail slices, not a diff.** `only_ref = ref[len(cand):]` / `only_cand = cand[len(ref):]` meant
   one extra line near the top shifted every later pair, so identical lines reported as changed and lines
   present in BOTH were named "only in CAND". **This is what hid the real 14-6.04 defect** behind a wrap
   complaint. Now uses `difflib.SequenceMatcher`.
3. **Silent fallback in `header_rule_pt`.** Thresholded on summed ink per row against *all* page ink — so
   the edge-to-edge `Protocol: … Page n of m` band set the scale and the table's own rule (~82% of page
   width) never qualified, falling through to `h * 0.18`, a fixed fraction of page height landing
   mid-header. 14-7.04's candidate returned exactly 110.16 = that fallback, with nothing to say so. Now
   detects a rule as a long **contiguous** run (text never is: the longest unbroken run in a line of
   Courier is one glyph) and **raises `NoHeaderRule`** rather than guessing. Contiguity also subsumes the
   spanner-underline worry that motivated defect 1 in the first place.

**Two real differences, both our table code:**
- **14-1.02** `t_14_1_02.R` had `p = "p-value\n[1]"`; legacy `programs/t-14-1-02.R:124` has
  `"p-value [1]"` with a **space**, and the reference RTF cell carries no `\line` while its neighbours do.
  One-character fix → header gate PASSES, pixel 2.022%. The issue's premise was inverted: the reference
  cell holds ONE line and ours held two, not the same lines distributed differently.
- **14-6.04** `t_14_6_04.R` passed no `merge`, so the default merged the bottom header row's six identical
  `"Baseline"` cells into one `gridSpan=6` — **five of the six labels were absent from the output**, a
  worse defect than the wrap that was reported, and invisible because of harness defect 2. Added
  `merge = "spanners"` (what `t_14_6_05.R` already did, which is why .05 passed). Six labels restored;
  content still at its documented `difflines=9` baseline (4 of those are the pre-existing
  `ERY. MEAN CORPUSCULAR` label wrap).
- **14-6.04's `"Shift to"` wrap remains, deliberately.** The SHIFT column is 0.588in in *both* — legacy
  `\cellx` = 3508/4234/5081 — so no transcription error. The reference fits it because **its header
  renders in a proportional font**: measured pitch varies per word (`Shift` 4.104, `Baseline` 4.438,
  `Normal` 5.452) while its body is a fixed 6.000 and ours is 6.000 throughout. 8 chars is ~33-38pt in
  Times vs 48pt in Courier, against 38.34pt of usable width. Widening the column would trade faithful
  geometry for a cosmetic line break, so it is recorded in divergences.md instead. **This affects every
  header-wrap comparison** — reference header text is 10-30% narrower than ours.

**The issue text also mischaracterised 14-7.01**: it gives the header as `"Planned\nRelative\nTime"`, but
`t_14_7_01.R:80` is `PRT = "Planned Relative Time"`, a plain string, and the RTF cell is
`{Planned Relative Time}`. Both wrap naturally to three lines and match exactly. `valign = "bottom"` is
not implicated in any of the nine — it correctly pins each cell's last line to the row's bottom baseline
in both renderers.

## Open upstream (as of 2026-07-27)
- **clinify #112** header wrapping to rendered column width (the wrapping half of #15, maintainer-invited).
- **clinify #113** `clin_header_pad()` overwrites per-row `flextable::padding(i=)` instead of composing.
- **clinify #101** `as_clintable()` header padding reset by `clin_column_headers()`; **#107** haven value-label
  partial match — both have pilot workarounds.
- **tplyr2 #55** `denom_row` NA/width (the last thing keeping a hand-rolled n-row).
- Declined upstream, no action: **clinify #106** in-body row grouping (Word owns pagination here), **#15** N
  injection (mapping + population choice aren't derivable — `tplyr_header_n()` covers the number instead).
- **denom_row (#35)** for 14-6.04/.05/.06 — emitted denominator row didn't match the pilot's exact format.
- `fish_p_str` (2 calls in 14-1.02: AE / lack-of-efficacy reason) — separate single tests, not a count layer.

### Two tplyr2 bugs found during the 14-2.01 attempt — FILED as #53/#54, both FIXED in PR #52 (verified)
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
