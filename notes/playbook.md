# CDISC Pilot Next-Gen — Build Playbook & Gotchas

Practical, verified guidance for building the tables. Full API refs are in
`notes/research/`. This is the "what actually mattered" digest.

## Pipeline
1. Summarize with **tplyr2** (`tplyr_spec` + `tplyr_layers` + `tplyr_build`).
2. For cross-arm model output (ANCOVA/MMRM/CMH LS-means, diffs, p-values) fit the
   model **once** on the full data and **row-bind** formatted rows beneath the
   tplyr2 build. `group_analyze()` sees only one column's slice at a time, so it
   **cannot** do cross-arm contrasts — do not try.
3. Render with **clinify**: `clintable(df, use_labels = FALSE)` → `clin_column_headers`
   → `flextable::*` (align/valign/width) → `add_titles_footnotes()` → `write_clindoc()`.
4. Verify: `python verify/compare_table.py <id>` (renders orig RTF + new DOCX to
   PDF via the local Gotenberg/LibreOffice service, scores pixel/gridline/geom/text).

## House style (R/clinify_defaults.R) — already tuned; reuse as-is
- Letter landscape, 1" margins, Courier New 10pt, header bottom-rule only,
  italic footnotes (no top rule), bold+italic titles.
- **Exact row heights** are the key to page fidelity: body 15.35pt, title/footnote
  11.4pt (flextable adds ~17% leading vs pharmaRTF → drift/overflow without this).
- Header buffer via header `padding.top = 18`; 4pt horizontal cell padding so
  monospace strings don't wrap; `top_spacer()` blank row for the header-rule-to-body gap.

## Verified gotchas
- **Result-column order:** count layers order `res*` alphabetically by `cols` value,
  ignoring factor levels (desc layers respect factor order). Always reorder `res*`
  by the `(N=)` label — and capture labels BEFORE any `b[order(...), ]` row-subset,
  which strips column attributes. See `feature-requests/tplyr2-count-column-order.md`.
- **Header Ns:** tplyr2's `res` labels count records, not distinct subjects — set the
  `(N=)` headers yourself from a distinct-subject count.
- **Wrapping:** LibreOffice ignores `\clNoWrap`; keep cell padding tight (4pt) so
  fixed-width strings fit. Column widths taken from the reference RTF `\cellx` twips.
- **collapse_row_labels()** gives the stub-header-then-indented-stat layout; indent
  stats by naming the `format_strings` entries with leading spaces ("  n").
- **Rounding:** keep tplyr2 `IBMRounding = FALSE` (banker's) to match the reference
  R RTFs (NOT SAS round-half-up).

## Fidelity thresholds (verify/fidelity.py)
Pass = per page: pixel-diff ≤ 2.5% area, gridline shift ≤ 2px, geometry ≤ 1pt,
page-count match. Residual ~1.5–2% is sub-pixel glyph rendering between the two
engine paths — genuinely "visually identical". text-set differences are reported
but not gating (dynamic tokens / display quirks).

## Sample status
- **14-3.01** (ANCOVA) — PASS (1.7%).
- **14-2.01** (demographics, desc+count+p, 4pp) — PASS.
- **14-5.01** (nested AE) — next; see status notes.
