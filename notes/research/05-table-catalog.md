I now have complete coverage. Here is the full deliverable.

---

# CDISC Pilot Replication — Table Output Catalog & 3-Table Build Recommendation

Repo root: `/Users/mstackhouse/Documents/repos/CDISC_pilot_replication`

Scope: 30 RTF outputs (`14-1.01` … `14-7.04`), the driver programs in `programs/`, the shared helpers, and the title/footnote metadata. The existing code is **hand-rolled tidyverse + huxtable + pharmaRTF** — there is no Tplyr in it today. This catalog maps each output to the `tplyr2` layer vocabulary (count / desc / shift / analyze) you asked for, so the rebuild can be scoped.

---

## STEP 1 — `data/titles.xlsx` structure

Single sheet `Sheet1`, 224 data rows. One row = one title **or** footnote line for one output. Columns:

| Column | Meaning | Observed values |
|---|---|---|
| `table_number` | Output id key, e.g. `14-1.01`. Also has a `14-1` entry (a KM figure, see note). | matches RTF basename |
| `index` | Line order **within** its type (1-based), resets per type | 1..6 |
| `type` | `title` or `footnote` | 2 values |
| `text1` | Primary text of the line (left segment when split) | free text / tokens |
| `text2` | Right segment for `align=split` lines only | tokens / blank |
| `align` | `split` \| `left` \| `center` | 3 values |
| `bold` | logical | titles `TRUE`, footnotes `FALSE` |
| `italic` | logical | **all** rows `TRUE` |
| `font` | font override | empty everywhere (falls back to Courier New) |

The reader `example_custom_reader()` (in `funcs.R`) forces these 9 column types: `c('text','numeric','text','text','text','text','logical','logical','text')`, then returns `df[df$table_number == table_number, !names(df)=='table_number']`.

**How an output id maps to lines** — every output follows the same title block:

- title idx 1 — `align=split`: left `text1="Protocol: CDISCPILOT01"`, right `text2="PAGE_FORMAT: Page %s of %s"`
- title idx 2 — `align=left`: `text1="Population: <pop>"` (All Subjects / Intent-to-Treat / Efficacy / Safety / Completers)
- title idx 3 — `align=center`: `text1="Table 14-x.xx"`
- title idx 4 — `align=center`: `text1=<the descriptive title>`

Footnotes vary in count (1–5 numbered/notes lines) and always end with:

- footnote (last idx) — `align=split`: left `text1="FILE_PATH: Source: %s"`, right `text2="DATE_FORMAT: %H:%M %A, %B %d, %Y"`

The tokens `PAGE_FORMAT:`, `FILE_PATH:`, `DATE_FORMAT:` are pharmaRTF directives: `%s` is filled at render time with page number / source path / timestamp. `split` alignment = left segment flush-left, right segment flush-right via a right tab.

The exact descriptive titles + footnotes per output are folded into the catalog (STEP 3) and the "Population" line is listed there.

**Note:** `titles.xlsx` also contains `table_number = "14-1"` → *"Figure 14-1 / Time to Dermatologic Event by Treatment Group"* (a Kaplan-Meier figure sourced from `adtte.xpt`). There is **no program and no RTF** for it in this repo, so it is not one of the 30 tables. `adtte.xpt` is present but unused by the 30.

---

## STEP 2 — Program architecture

### `programs/config.R`
```r
sdtm_lib <- "data/sdtm"
adam_lib <- "data/adam"
options(huxtable.add_colnames = FALSE)   # revert huxtable v5 default
```
> **Runnability caveat:** the actual `.xpt` files live **directly in `data/`** (e.g. `data/adsl.xpt`), not in `data/adam/` or `data/sdtm/`. As written, the `read_xpt(glue("{adam_lib}/adsl.xpt"))` calls point at `data/adam/…`. Either the data must be moved into `data/adam` + `data/sdtm`, or `config.R` must be repointed. Also, **`14-7.04` needs `data/sdtm/cm.xpt`, which is NOT present** in `data/` at all — that table cannot run until SDTM CM is supplied.

### `programs/funcs.R` — shared helpers (verified signatures)

Formatting / building blocks:
- `num_fmt(var, digits=0, size=10, int_len=3)` — rounds, formats to `int_len+digits` width, right-pads to `size`. `Vectorize`d. This is the workhorse behind decimal alignment.
- `n_pct(n, pct, n_width=3, pct_width=3, digits=0, mark_lt=TRUE)` — `"{n} ({pct}%)"`; renders `<1` when ratio `< .01` and `mark_lt`. `Vectorize`d. (Several programs **locally override** `n_pct` with a variant that prints `" 0      "` for zero — 14-6.02/03/04/05/06, 14-7.03, 14-7.04.)
- `pad_row(.data, n=1)` — append `n` blank rows (a second, index-based `pad_row(df, r)` is redefined inside the lab/VS programs).
- `get_header_n(.data, trtp=TRT01P, trtpn=TRT01PN)` — big-N per arm; builds `labels = str_wrap(glue('{trtp} (N={N})'), width=10)` with `\line` breaks (this is the `(N=xx)` column header).
- `desc_stats(.data, var, group=TRT01PN, na.rm=TRUE, int_len=3, size=10, include=c('n','Mean','SD','Median','Min','Max'))` — the **desc** engine.
- `sum_subgrp(.data, subgroup_var, order_var=NULL, include.n=TRUE, pad.row=TRUE, header_n=header_n)` — the **count** engine (n(%) categorical rows, transposed by arm).
- `get_meta`, `invert.list` — utilities.

Statistics (external models — the "analyze" content):
- `aov_p(.data, forumula)` — ANOVA p (`aov`). *(arg is misspelled `forumula` in source.)*
- `chi_p(data, results, categories)` — Pearson chi-square p.
- `fish_p(data, results, categories, width=10)` — Fisher exact p.
- `fisher_test_ae(.data)` — row-wise 2×2 Fisher for AE tables; appends `*` if `p<.15`, `>.99` cap.
- `ae_counts(.data, ..., N_counts=NULL, sort=FALSE)` — subject + event counts by arm with Fisher p_low/p_high; drives 14-5.01/02.
- `cmh_p(.data, formula, alternate=c('rmeans','cmeans','general','cor'))` — CMH via **`vcdExtra::CMHtest`** (requires the author's patched fork of vcdExtra).
- `attach_p(.data, p_value, digits=4)` — puts a p-value on row 1 of a block.
- `summary_data(data, var, week, stub_header)` — efficacy desc block (n / Mean(SD) / Median(Range)) at a visit.
- `efficacy_models(data, var=NULL, wk=NULL, model_type='ancova')` — the efficacy model engine:
  - `ancova`: `contr.sum` contrasts (Type III); `lm(CHG ~ TRTPN + SITEGR1 + BASE)` for dose-response via **`car::Anova(type=3)`**; `lm(… TRTPCD_F …)` + **`emmeans::lsmeans(weights='proportional')`** → pairwise diffs, 95% CI, p (unadjusted).
  - `repeated` (MMRM): **`lme4::lmer(CHG ~ TRTPCD_F + SITEGR1 + AWEEKC + TRTPCD_F:AWEEKC + BASE + BASE:AWEEKC + (AVISITN | USUBJID))`** + `emmeans` with `lmer.df='kenward-roger'`.

Every program ends with the same pharmaRTF chain:
```r
doc <- rtf_doc(ht, header_rows=…) %>%
  titles_and_footnotes_from_df(from.file='./data/titles.xlsx',
                               reader=example_custom_reader,
                               table_number='14-x.xx') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)   # (+ set_header_height/set_footer_height on wide tables)
write_rtf(doc, file='./outputs/14-x.xx.rtf')
```

**Filename separator quirk (confirmed):** `14-1.01`…`14-5.02` use hyphen programs (`t-14-1-01.R`), while `14-6.02`…`14-7.04` use dot programs (`t-14-6.02.R`). `t-14-6-01.R` is the odd one (hyphen). All RTF outputs use dots.

---

## STEP 3 — Full catalog (30 outputs)

Layer legend: **count** = categorical n/n(%) (Tplyr `group_count`), **desc** = summary stats (Tplyr `group_desc`), **shift** = from→to crosstab (Tplyr `group_shift`), **analyze** = injected external model/test result (no native Tplyr layer — supplied as data). "Ext. model" = needs a statistical model/test outside the tabulation.

| ID | Population / Title | Source ADaM | Type | tplyr2 layers | Ext. model | Cx |
|---|---|---|---|---|---|---|
| **14-1.01** | All Subjects · Summary of Populations | ADSL | Disposition / population counts | count | No | **Low** |
| **14-1.02** | ITT · Summary of End of Study Data | ADSL | Disposition (completion status + reason for early termination) | count | Yes — Fisher (`fish_p`) | Med |
| **14-1.03** | All Subjects · Summary of Number of Subjects By Site | ADSL | Disposition crosstab: Site × Arm × (ITT/Eff/Com) | count (multi-dim, nested spanners) | No | **High** (bespoke pivot + 4 arm×3 flag spanners) |
| **14-2.01** | ITT · Summary of Demographic and Baseline Characteristics | ADSL | Demographics / baseline characteristics | **desc + count** | Yes — ANOVA (`aov_p`) + Pearson chi-sq (`chi_p`) | Med |
| **14-3.01** | Efficacy · Primary Endpoint: ADAS-Cog(11) Chg→Wk24 LOCF | ADADAS | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA (car TypeIII + emmeans) | High |
| **14-3.02** | Efficacy · Primary Endpoint: CIBIC+ at Wk24 LOCF | ADCIBC | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA | High |
| **14-3.03** | Efficacy · ADAS-Cog Chg→Wk8 LOCF | ADADAS | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA | High |
| **14-3.04** | Efficacy · CIBIC+ at Wk8 LOCF | ADCIBC | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA | High |
| **14-3.05** | Efficacy · ADAS-Cog Chg→Wk16 LOCF | ADADAS | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA | High |
| **14-3.06** | Efficacy · CIBIC+ at Wk16 LOCF | ADCIBC | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA | High |
| **14-3.07** | Completers · ADAS-Cog Chg→Wk24, OC, Windowed | ADADAS | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA | High · ⚠ **known discrepancy** (baseline n off by 1) |
| **14-3.08** | Efficacy · ADAS-Cog Chg→Wk24, **Male** LOCF | ADADAS | Efficacy ANCOVA (sex subgroup) | desc + analyze | Yes — ANCOVA | High |
| **14-3.09** | Efficacy · ADAS-Cog Chg→Wk24, **Female** LOCF | ADADAS | Efficacy ANCOVA (sex subgroup) | desc + analyze | Yes — ANCOVA | High |
| **14-3.10** | Efficacy · ADAS-Cog Mean & Mean-Change over Time | ADADAS | Efficacy descriptive over time (LOCF + Windowed sets) | desc (by visit) | No | High (very wide, "Change from baseline" spanner) |
| **14-3.11** | Efficacy · ADAS-Cog Repeated-Measures (MMRM) Chg→Wk24 | ADADAS | Efficacy **MMRM** | analyze | Yes — lme4 `lmer` + emmeans (KR) | High · ⚠ **known discrepancy** (p-values) |
| **14-3.12** | Efficacy · Mean NPI-X Total Score Wk4–Wk24 Windowed | ADNPIX | Efficacy ANCOVA | desc + analyze | Yes — ANCOVA | High · ⚠ **known discrepancy** (counts) |
| **14-3.13** | Efficacy · CIBIC+ Categorical Analysis LOCF | ADCIBC | Categorical counts by visit + CMH | count + analyze | Yes — CMH general assoc (vcdExtra fork) | High |
| **14-4.01** | Safety · Summary of Planned Exposure to Study Drug | ADSL | Exposure descriptive (avg daily / cumulative dose) | desc | No | Med (dual Completers/Safety spanner) |
| **14-5.01** | Safety · Incidence of TEAEs by Treatment Group | ADAE (+ADSL) | **AE nested by SOC → Preferred Term** | count (nested) + analyze | Yes — Fisher per row | High |
| **14-5.02** | Safety · Incidence of TE **Serious** AEs | ADAE (+ADSL) | AE nested by SOC → PT (serious subset) | count (nested) + analyze | Yes — Fisher per row | High |
| **14-6.01** | Safety · Summary Stats for Continuous Lab Values | ADLBC + ADLBH | Lab descriptive by visit (chem + heme) | desc (param × visit) | No | High (very wide, per-arm spanners, ~30 params × visits) |
| **14-6.02** | Safety · Freq Normal/Abnormal (Beyond Normal Range) Labs | ADLBC + ADLBH | Lab L/N/H frequency | count / shift-like | Yes — Fisher 3×3 | High |
| **14-6.03** | Safety · Freq Normal/Abnormal (Clinically Sig Chg from Prev Visit) | **ADLBCPV + ADLBHPV** | Lab L/N/H frequency (previous-visit ref) | count / shift-like | Yes — Fisher 3×3 | High |
| **14-6.04** | Safety · Shifts of Lab Values by Visit (threshold ranges) | ADLBC + ADLBH | **Lab shift** baseline→post, by visit | shift | No | High (~540 rows, largest RTF) |
| **14-6.05** | Safety · Shifts of Lab Values (threshold ranges) | ADLBC + ADLBH | Lab shift + CMH | shift + analyze | Yes — CMH `rmeans` (vcdExtra fork) | High · ⚠ **known discrepancy** |
| **14-6.06** | Safety · Shifts of Hy's Law Values | ADLBHY | Hy's-law shift + CMH | shift + analyze | Yes — `mantelhaen.test` | High · ⚠ **known discrepancy** |
| **14-7.01** | Safety · Summary of Vital Signs at Baseline & EOT | ADVS (+ADSL) | Vital signs descriptive | desc | No | Med-High · ⚠ **known discrepancy** (EOT flag) |
| **14-7.02** | Safety · Summary of Vital Signs Change from Baseline at EOT | ADVS (+ADSL) | Vital signs change descriptive | desc | No | Med-High · ⚠ **known discrepancy** |
| **14-7.03** | Safety · Summary of Weight Change from Baseline at EOT | ADVS (+ADSL) | Weight change descriptive | desc | No | Med-High · ⚠ **known discrepancy** |
| **14-7.04** | Safety · Summary of Concomitant Medications (# Subjects) | **SDTM CM** (+ADSL) | ConMeds nested by therapeutic class → medication | count (nested) | No | Med-High · ⚠ **needs `data/sdtm/cm.xpt` (missing)** |

Datasets used across the set: `adsl, adadas, adcibc, adnpix, adae, adlbc, adlbh, adlbcpv, adlbhpv, adlbhy, advs` (all present) + `sdtm/cm.xpt` (absent). `adtte.xpt` present but unused.

---

## STEP 4 — Visual specification ("visually identical")

I opened three RTFs spanning the layout families — **`14-1.01.rtf`** (demographics-family / counts), **`14-3.01.rtf`** (efficacy), **`14-5.01.rtf`** (AE nested). The page frame, fonts, and title/footnote machinery are **identical across all three** (all emitted by the same pharmaRTF chain). Below is the exact, verified spec.

### Page frame (identical in all files)
```
{\rtf1\ansi\deff1                       default font = f1
{\fonttbl {\f0 Times;} {\f1 Courier New;}}{\colortbl;;}
\paperw15840\paperh12240                US Letter: 15840 twips = 11.0", 12240 twips = 8.5"
\lndscpsxn                              LANDSCAPE section  → page is 11" wide × 8.5" tall
\margl1440\margr1440\margt1440\margb1440   ALL margins = 1440 twips = 1.0 inch
\headery720\footery720                  header/footer band 720 twips = 0.5" from page edge
\fs20                                   default size 20 half-pt = 10 pt
\widowctrl\ftnbj\fet0\sectd\linex0
```
- **Paper:** US Letter, **landscape**, 1-inch margins all sides.
- **Font family:** everything is **Courier New** (monospace). `\deff1` makes f1 the default, so body cells with no explicit font inherit Courier New; titles/footnotes name `\f1` explicitly. Times (`f0`) is declared but never used. The monospace font is load-bearing — decimal/column alignment is achieved purely by fixed-width, space-padded strings from `num_fmt`/`n_pct`, **not** by decimal tabs.
- **Font size:** **10 pt** (`\fs20`) everywhere — titles, column headers, body, footnotes.

### Title block (RTF `{\header …}`) — 4 lines, Courier New 10 pt, **bold + italic** (`\b\i`)
1. Split line: left `Protocol: CDISCPILOT01`; right `Page {PAGE} of {NUMPAGES}` (live RTF page fields). Right-aligned via tab stops `\tx7245 \tqr\tx12960` and pharmaRTF `\pmartabqr`.
2. Left-aligned (`\ql`): `Population: <pop>`.
3. Centered (`\qc`): `Table 14-x.xx`.
4. Centered (`\qc`): descriptive title.

### Footnote block (RTF `{\footer …}`) — Courier New 10 pt, **italic, NOT bold** (`\i` only)
- Each numbered/`Note:` footnote on its own left-aligned (`\ql`) line. Count varies 1→5 (e.g. 14-1.01 has 1 note; 14-5.01 has 5; 14-3.01 has 3 numbered `[1][2][3]`).
- Final line is split: left `Source: programs/t-14-x-xx.R`; right timestamp `HH:MM Weekday, Month DD, YYYY` (from `DATE_FORMAT`).

### Column-header structure
- A **spacer/buffer row** is emitted first (`set_column_header_buffer(top=1)`) — an all-empty bordered row above the labels.
- The label row: **bold, centered (`\qc`), bottom-valigned (`\clvertalb`)**, with a **single 1-pt bottom border** (`\clbrdrb\brdrs\brdrw20`; 20 twips = 1 pt).
- Big-N labels wrap at width 10 via `\line`, e.g. `Placebo\line (N=86)`, `Xanomeline\line Low Dose\line (N=84)`.
- **Spanners** use merged cells (`\clmgf` first / `\clmrg` continuation). 14-5.01 has a 2-row header: row-1 arm spanners each spanning the `n(%)`+`[AEs]` pair (`merge_cells(1, 2:3)`, `4:5`, `6:7`) plus a `Fisher's Exact\line p-values` spanner over `8:9`; row-2 sub-labels `n(%) / [AEs] / n(%) / [AEs] / n(%) / [AEs] / Placebo\line vs.\line Low Dose / Placebo\line vs.\line High Dose`, with the stub header `System Organ Class/\line Preferred Term`. Lab/exposure tables similarly merge per-arm spanners with dashed or solid bottom borders (e.g. 14-6.01 uses `set_bottom_border_style('dashed')` under arm spanners).

### Body
- Cells left-aligned (`\ql`), Courier New 10 pt, not bold. Numeric alignment = pre-formatted monospace strings (e.g. `" 86 (100%)"`, `" 60 ( 70%)"`).
- **Row padding:** `set_ignore_cell_padding(TRUE)` + `top_padding=0 / bottom_padding=0` → tight single-line rows. Cell padding attributes in RTF are `\clpadl80…\clpadr80` (80 twips) but suppressed by ignore-padding.
- **Cell widths** are absolute twip boundaries per column (`\cellxN`), derived from `set_width(<ratio of text area>)` × `set_col_width(<vector of proportions>)`. Example 14-1.01: `set_width(1.1)`, `col_width=c(.4,.15,.15,.15,.15)` → boundaries `3802/5228/6653/8079/9504`. 14-3.01: `set_width(1.2)`, `c(.5,1/6,1/6,1/6)` → `5184/6912/8640/10368`.
- **Sub-row indentation** is literal leading whitespace inside the label string, not an RTF indent:
  - AE preferred terms: **two leading spaces** — `{  DIARRHOEA}` under `{GASTROINTESTINAL DISORDERS}` under `{ANY BODY SYSTEM}` (verified). SOC flush-left, PT indented 2 spaces.
  - Efficacy stat sub-labels: two leading spaces — `  n`, `  Mean (SD)`, `  Median (Range)` under a flush-left stub header.
  - Disposition / conmed rows: a leading tab `\t` in the source string.
- **Blank spacer rows** between blocks are real empty table rows (`pad_row`).

That combination — Letter/landscape, 1" margins, Courier New 10 pt, bold-italic 4-line title, italic footer with live page/date fields, 1-pt underlined bold centered headers with `(N=)` line-wraps, tight zero-padding monospace body with space-based decimal alignment and 2-space PT indents — **is the definition of "visually identical"** for this suite.

---

## STEP 5 — Recommended first 3 tables

Chosen to exercise the full breadth of `tplyr2` (desc / count / nested-count / injected analysis) and `clinify` (grouping, indentation, pagination, spanners, RTF frame), while staying on **clean, no-known-discrepancy** outputs.

### Pick 1 (desc + count, simple layout) → **`14-2.01` — Summary of Demographic and Baseline Characteristics**
- **Why:** The canonical demographics table and the only common table that exercises **both** a `desc` layer (Age, MMSE, disease duration, education, weight, height, BMI) **and** a `count` layer (Sex, Race, Age groups, Duration groups, BMI groups) in one display. Layout is a straightforward single-header, 4-arm (Placebo/Low/High/Total) + p-value column — ideal for standing up the `tplyr2`→`clinify` pipeline and the shared title/footnote frame.
- **Exercises:** desc + count layers; Total-column derivation; row-label grouping/indentation; the standard 4-line title / 3-line footnote frame.
- **Stats:** lightweight and well-behaved — ANOVA p (`aov_p`) for continuous, Pearson chi-square p (`chi_p`) for categorical. Good first look at injecting a p-value column without a heavy model. No known discrepancy.
- Simpler warm-up if you want zero stats first: **`14-1.01`** (Summary of Populations, pure `count`, cleanest layout in the suite) — but it does not touch `desc`, so `14-2.01` is the better single pick.

### Pick 2 (nested count) → **`14-5.01` — Incidence of TEAEs by SOC / Preferred Term**
- **Why:** The flagship **nested count** — subjects counted at SOC and Preferred-Term level, PTs indented 2 spaces under their SOC, `ANY BODY SYSTEM` overall row, descending-frequency sort within SOC. This is exactly the `tplyr2` nesting + `clinify` grouping / indentation / **multi-page pagination** (this RTF is 425 KB, many pages) and per-arm **spanner** exercise.
- **Exercises:** nested `count` layer; the `[AEs]` event-count companion column (a nuance beyond a plain nested count); paired-cell arm spanners; page breaking with repeated headers.
- **Stats:** per-row Fisher exact (Placebo vs Low, Placebo vs High) with `*` for p<.15 — a clean `analyze`-style injected column. No known discrepancy.
- **Iterate-first tip:** build/validate the identical structure on **`14-5.02`** (Serious AEs) first — same layout, far fewer rows, so you can eyeball correctness before scaling to 14-5.01.

### Pick 3 (statistics-driven efficacy) → **`14-3.01` — Primary Endpoint: ADAS-Cog(11) Change from Baseline to Week 24, LOCF (ANCOVA)**
- **Why:** The primary efficacy endpoint and the cleanest **model-driven** table. Combines a `desc` block (`summary_data`: n / Mean(SD) / Median(Range) at Baseline, Week 24, Change) with an **injected ANCOVA** result block (`efficacy_models`, `model_type='ancova'`): dose-response p (car Type III), LS-Means pairwise Diff(SE), 95% CI, unadjusted p (emmeans). This is the reference implementation for the "`analyze` layer + externally computed model results fed into the table" pattern, and it **matches the original numbers** (no discrepancy flagged).
- **Exercises:** desc + analyze; footnote references `[1][2][3]`; mixed stat-string formatting/alignment.
- Preferred over the MMRM table because MMRM (`14-3.11`) has a **known p-value discrepancy** — see avoid-list.

### Tables to AVOID for the initial sample (README "Notes on Data" known discrepancies)
Do not use these to judge fidelity — the reference numbers themselves diverge from the original CDISC displays for documented data/tooling reasons:

- **`14-3.07`** — ADAS-Cog Completers/Windowed: baseline n off by 1.
- **`14-3.11`** — ADAS-Cog MMRM: p-values slightly off (lme4 covariance-structure ambiguity).
- **`14-3.12`** — Mean NPI-X: counts differ (source subset no longer reproducible).
- **`14-6.05`** — Lab threshold shifts: Bilirubin/Monocytes values differ (ref-range/shift vars taken as-is).
- **`14-6.06`** — Hy's Law shifts: transaminase/bilirubin values differ (one p-value hard-blanked in code).
- **`14-7.01`, `14-7.02`, `14-7.03`** — Vital signs / weight at End-of-Treatment: EOT flag differences propagate through.

Also avoid for a *first* build for practical (not accuracy) reasons: **`14-7.04`** (requires `data/sdtm/cm.xpt`, currently **missing**); **`14-6.04`** (~540 rows, largest output — good later stress test for pagination, not a starter); and note the **`data/adam` / `data/sdtm` path mismatch** in `config.R` must be resolved before any program runs.

The recommended trio — **`14-2.01`, `14-5.01`, `14-3.01`** — is discrepancy-free, spans desc + count + nested-count + injected-analysis, and collectively stresses every part of the shared RTF visual frame.