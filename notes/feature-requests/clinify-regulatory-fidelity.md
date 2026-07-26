# clinify FRs — regulatory-fidelity helpers

Discovered while reproducing the CDISC pilot RTFs as pixel-identical .docx. Each
of these we currently solve in a project-level `R/clinify_defaults.R`, but they
are generic enough to belong in clinify.

## 1. Exact / configurable row height ("regulatory row pitch")   **highest value**
flextable→docx→LibreOffice renders 10pt Courier body rows at ~18pt, vs the
pharmaRTF reference's ~15.35pt — a ~17% taller row. Over a page this drifts
content down and overflows pages. We fix it by forcing
`flextable::hrule(rule="exact") + height_all()` at a measured pitch (body 15.35pt,
titles/footnotes 11.4pt). A clinify option like `clin_row_height(body=, title=,
footnote=)` (or a `clinify.row_pitch` default) would make this a one-liner and
guarantee page-for-page fidelity, which matters for validated outputs.

## 2. Spreadsheet → titles/footnotes reader
clinify has no way to source titles/footnotes from a spreadsheet, yet TLF suites
almost always carry a titles dataset. We wrote `R/titles.R` (reads a titles.xlsx
with `table_number/index/type/text1/text2/align`, resolves `PAGE_FORMAT`,
`FILE_PATH`, `DATE_FORMAT` tokens, applies clinify alignment rules). Proposal: a
`clin_titles_from_data(df, table_id, ...)` with a small pluggable token/alignment
mapping so shops can point it at their own titles metadata.

## 3. Header "buffer" row + first-body spacer as first-class options
The reference leaves ~1 blank row above the column header
(pharmaRTF `set_column_header_buffer`) and a gap between the header rule and the
first body row. We reproduce these with header `padding.top` and a manual leading
blank row. A `clin_column_headers(..., buffer = TRUE)` / `clin_header_gap()`
option would be cleaner and portable.

## 4. Left-aligned single-line title without the duplicate-text trick
To get a single left-aligned title line, clinify requires passing the text twice
(`c("x","x")`) so it merges. A `justify=` per line (or auto-detect a 1-element
line that should be left) would be friendlier.

## 5. (Consider) native decimal/zero-count display options
See `tplyr2-count-display-quirks.md` — a "<N%" percent format and a zero-count
display option would remove per-project post-processing. Could live in tplyr2
(formatting) rather than clinify.
