# Intentional / documented divergences from the legacy RTFs

Per the project policy (use correct modern values; document where they differ from
the 2020 R RTFs). Verification for these tables: the LAYOUT/structure is matched
(and checked), but specific values differ for the reasons below.

## 14-3.11 - ADAS-Cog(11) Repeated Measures (MMRM)  [deliberate; VALIDATED vs SAS]
Fit with **mmrm** (unstructured within-subject covariance + Kenward-Roger) = SAS
`PROC MIXED TYPE=UN`. LS means at **Week 24** via `emmeans(fit, ~TRTPCD_F, at=Week 24)`
with **equal** weights over SITEGR1 (SAS `lsmeans`, no OM). **Matches the SAS output
in atorus GitHub issue #6 to the digit:**

| stat | Placebo | Low | High |
|---|---|---|---|
| ours (mmrm, Wk24) | 2.3 (0.69) | 1.7 (0.77) | 1.5 (0.84) |
| SAS (issue #6)    | 2.3291     | 1.7352     | 1.5009     |

Contrasts: High-PBO 0.8282 p=0.4403, Low-PBO 0.5939 p=0.5600 - exact vs SAS.

**Divergence from the legacy RTF is intended and explained by issue #6:** the original
report showed LS means **averaged over weeks 8/16/24** (main-effect `lsmeans trtpcd_f`)
despite the "Week 24" title - a documented title/number mismatch - and used lme4
(which also didn't match SAS). We produce the correct Week-24 MMRM values.

Two gotchas that mattered:
1. Use the `trtpcd_f*aweekc` interaction at Week 24, not the averaged main effect.
2. Equal weights, NOT `weights="proportional"`. Proportional/OM is correct for the
   PROC GLM ANCOVA tables (14-3.01 etc.), but this MMRM `lsmeans` is unweighted. The
   treatment *contrasts* are weight-invariant (so they matched either way), but the
   marginal LS means were ~0.13 off until this was fixed.

## Row pitch standardized (option A) - affects pixel match, not content
The legacy tables use INCONSISTENT body row pitches (14-3.01/14-2.01 ~15.35pt,
AE 14-5.01/.02 ~11.6pt, 14-3.13 ~17.35pt). We render every table at ONE house
pitch (~15.35pt via hrule="atleast" floor). Tables whose legacy pitch happened to
be 15.35 pixel-match; the others (AE, 14-3.13, ...) drift vs the legacy by design
and are verified by the pagination/pitch-agnostic CONTENT check instead. This is
the intended "clean & consistent" outcome, not an error.

## 14-3.13 - CIBIC+ Categorical (CMH)  [VALIDATED via coin]
Row-mean-scores CMH via **coin::cmh_test** (ordinal AVAL scored, treatment groups,
stratified by SITEGR1). Reproduces the legacy (forked-vcdExtra) / SAS values EXACTLY:
Wk8 0.2727, Wk16 0.4003, Wk24 0.6180. **coin is a clean CRAN replacement - the
vcdExtra fork is no longer needed.** Content identical; pixel diverges only via the
row-pitch standardization above (legacy 14-3.13 pitch was ~17.35pt).

## 14-2.01 - Duration of disease, Placebo Mean  [R-version rounding]
Raw mean is exactly 42.65. R changed `round()` in R 4.0.0, so modern R gives 42.7
(legacy R 3.6 gave 42.6). One cell; unavoidable and correct-for-modern-R.

## 14-6.01 - Summary Statistics for Continuous Lab Values  [wide transposed desc]
Content IDENTICAL to the legacy RTF except **4 cells** out of 370, all half-way (`x.x5`)
values hit by the same R 4.0.0 `round()` change as 14-2.01 (e.g. AVAL mean 94.55 -> 94.5
legacy / 94.6 modern; a Wk24 change -0.05 -> -0.1 legacy / 0.0 modern). 115 AVAL + 60 CHG
means land on `.x5` boundaries; only these 4 fall on the float representation where old/new
`round()` disagree. Unavoidable and correct-for-modern-R.
- **N is the record count at the visit** (subjects assessed), matching the reference's
  `n()`. tplyr2's built-in `n` statistic counts non-missing AVAL instead, which is 1-2 lower
  for the handful of Bilirubin/Glucose visits with an assessed-but-missing value. We compute
  N directly (records) so the displayed N matches the reference; mean/SD are tplyr2
  (na.rm, so over non-missing) - same convention as the legacy.
- tplyr2 renders a rounds-to-negative-zero mean as `-0.0`; base R `format()` (legacy) drops
  the sign to `0.0`. Normalized locally (`fix_neg_zero`). Minor tplyr2 papercut (see REVISIT).
- 18 pages vs reference 17 - the house row-pitch (option A) pagination consequence.

## Column headers are NOT bolded (option A, suite-wide)
The pilot programs are inconsistent (22/30 call `set_bold`), but the RENDERED references are
mostly regular weight. Measured: a bold house default WORSENED every pixel-gated single-page
table (14-1.01 0.907->0.985, 14-3.01 1.174->1.239, 14-1.03 2.744->2.889 nobold->bold) and
its only clearly-bold reference, 14-6.01, is content-verified (multi-page) so header weight
is invisible to its gate. So non-bold is the consistent house choice. **Exception:** 14-6.01's
reference headers ARE visibly bold; our non-bold rendering is a documented cosmetic divergence
there, immaterial to its content verification.

## 14-1.03 - Number of Subjects By Site  [content-identical, PASS]
Bespoke site x arm x flag crosstab; content and values IDENTICAL, pixel 2.40% PASS. One
build detail: clinify emits a compact 2-row header (spanner + labels) whereas the legacy left
a blank spacer row between them, so the body sat ~16pt too high and mis-registered. Added
`padding.top` to the sub-label header row to reproduce the gap and align the body.

## No pilot table needs horizontal pagination (clin_alt_pages)
Checked 14-6.01/14-3.10/14-6.02/14-6.03 (the widest candidates): every one is a single
landscape section (`\lndscpsxn`, zero `\sect`, no "continued") that fits the 9" usable
width. The "wide / alternating-page" type flagged as unproven does not actually occur in
this suite; clinify `clin_alt_pages` is not required. Wide tables (14-6.01 etc.) are just
tall single-width tables verified by CONTENT.

## Batch build 2026-07-24 (parallel agents): remaining lab / vitals / efficacy tables

### 14-6.02 / 14-6.03 - Lab value frequency (Normal/Abnormal; clin. signif. change)  [value-identical]
Every cell + all 30 Fisher's-exact p-values match the reference exactly. The `--content` line-set
check FAILs ONLY because the legacy set the table to 1.4x page width, so LibreOffice wraps every
"16( 19%)" cell onto two physical lines -> the reference is 4 pages of `) ) )` fragments; our clean
widths give 2 pages of the same values (the CAND-only full rows carry every correct value). Two
legacy artifacts: the p-value column header carries a dangling "[1]" marker with NO corresponding
footnote (legacy never wrote the Fisher footnote) - reproduced to match the reference, not invented;
if a "[1] Fisher's exact test" footnote is ever wanted it must be added to data/titles.xlsx. Legacy
`complete(nesting(...))` errors under modern tidyr (and was a no-op) -> replaced with factorial
zero-fill.

### 14-6.05 - Lab shift (threshold ranges) + CMH  [content-identical]
Shift counts/percentages identical; CMH via coin::cmh_test reproduces all six shown legacy/SAS
p-values to the digit. Degenerate strata handled to MATCH the reference blanks: a guard blanks any
test with <2 baseline strata (undefined "controlling for baseline status"), reproducing the legacy
blanks. The "ref-range/shift vars as-is" discrepancy produced NO value differences here.

### 14-6.06 - Hy's Law shift + CMH  [one corrected p-value; keep-modern]
Content identical except ONE cell: the TRANSHY (Transaminase 1.5xULN) CMH p. Legacy showed **0.015**
from a genuine BUG - `mantelhaen.test(array(unlist(N), dim=c(2,3,2)))` on a TRTP-major count vector,
which scrambles treatment/baseline across the array. Correct value **0.405**, verified two ways
(coin::cmh_test AND a correctly-shaped base mantelhaen.test). We keep the correct modern value.
HYLAW p blank in both (empty High-at-baseline stratum; legacy had it commented out "FIXME").

### 14-7.01 / 14-7.02 / 14-7.03 - Vital signs descriptive  [content/pixel-identical]
All pass (14-7.01/.02 content-identical; 14-7.03 pixel 2.461%). The documented "EOT flag"
discrepancy is baked into the shipped ADVS End-of-Treatment records (a difference vs the ORIGINAL
CDISC analysis data, not vs our reference RTF, which is built from the same ADVS), so it affects
reference and candidate identically -> NO differing cells. 14-7.02: N = record count (n()) as in the
reference; tplyr2's `n` = non-missing is 1 lower for two Xan-High EOT groups (ref 73 vs 72) that
carry an assessed-but-missing value - same record-count-N convention as 14-6.01. 14-7.03 needed a
per-table docx override reproducing the legacy `set_header_height(1)` 1-inch header band (vs the
0.5" house default) so the single-page body lands inside the harness registration window.

### 14-3.07 - ADAS ANCOVA (completers)  [identical; "off by 1" is real source data]
Content identical (30/30); model rows match to the digit (ancova_block reused as-is). The legacy
"baseline n off by 1" comment is NOT a build divergence - it is a real property of the source data
(baseline 60/28/30 vs Week-24 completer 59/27/30: one Placebo and one Xan-Low completer have a
baseline record but no observed Week-24 assessment). Both legacy and our build display it identically.

### 14-3.12 - NPI-X ANCOVA  [identical]
Pixel 1.074%, content identical (26/26). The documented NPI-X "subset counts" discrepancy is against
the UNAVAILABLE original SAS analysis dataset; both the legacy RTF and our build derive NPTOTMN from
the shipped adnpix and agree exactly. ancova_block reused unchanged.

### 14-3.10 - ADAS wide desc (Windowed + LOCF)  [identical]
Pixel 0.336%, content identical (27/27), every value byte-identical (num_fmt reused). No spanner
underline (unlike 14-6.01); header not bolded (accepted cosmetic divergence - reference is bold).
Exposed a fidelity-harness registration fragility on uniform numeric grids (see notes/REVISIT.md).

### 14-7.04 - Summary of Concomitant Medications  [content-identical]
Built once data/cm.xpt became available (2026-07-24). Nested distinct-subject count: therapeutic
class (CMCLAS) -> coded medication (CMDECOD), + a "patients receiving at least one" top row;
denominator = arm N (86/84/84). Content IDENTICAL (51/51 lines), 3 pages matching. Computed
directly (tplyr2 group_count(distinct_by=) is the natural engine, cf. AE, but its fixed-width
f_str percent can't produce the reference's UNPADDED integer percent "7 (8%)"). Medications ordered
within class by descending Placebo count then alphabetical (legacy radix order). Body valign="top"
so arm counts sit on the first line of the wrapped "at least one" label (same idiom as 14-7.x).

## The reference RTFs render their column headers in a PROPORTIONAL font (affects every header-wrap comparison)

Measured character pitch, `14-6.04` page 1, from `pdftotext -bbox` glyph widths:

| word | where | reference | ours |
|---|---|---|---|
| `Shift` | header | **4.104** | 6.000 |
| `Baseline` | header | **4.438** | 6.000 |
| `Normal` | header | **5.452** | 6.000 |
| `Week` | header | **5.960** | 6.000 |
| `CHEMISTRY` | body | 6.000 | 6.000 |
| `ALANINE` | body | 6.000 | 6.000 |

The reference's header pitch **varies per word**, which only a proportional font does;
its body is a fixed 6.000. The RTF fonttbl is `{\f0 Times}{\f1 Courier New}` and the
header cell paragraphs are `\pard\intbl\qc{\fs20 \b {Shift to}` with no `\f1`, so they
fall back to `\f0` Times while the body renders Courier. Almost certainly a
pharmaRTF/huxtable quirk rather than a decision - the legacy output is internally
inconsistent with itself. We render Courier New throughout, consistent with the legacy
*body* and with every other table in the suite.

Consequence, and the reason this is worth recording on its own: **header text in the
reference is 10-30% narrower than the same string in our output.** Any "their header
fits on one line and ours wraps" observation has to account for this before it can be
attributed to anything else. The reference's header cells also carry `\clNoWrap`, so
the legacy would have overflowed rather than wrapped even where text did not fit.

### 14-6.04 - `"Shift to"` wraps to two lines  [unavoidable consequence of the above]
The SHIFT column is **0.588in in both** - the legacy `\cellx` values are 3508/4234/5081/
then 6×0.840in, which the pilot transcribes exactly - so this is not a width transcription
error. `"Shift to"` is 8 characters: ~33-38pt in the reference's proportional header font,
which fits the 38.34pt of usable cell width, but 48pt in Courier New 10, which does not.

Given faithful column widths and a monospace header, the wrap cannot be avoided.
**Deliberately NOT "fixed" by widening the column**, which would trade a faithful
geometry for a cosmetic line break and shift all six arm columns right. Everything else
about this header now matches: the six `Baseline` sub-labels are present and in their own
columns (see below).

Note also that the legacy header was never levelled the way ours is - `outputs/14-6.04.rtf`
is a single header row of whole strings (`{Week}`, `{Shift to}`, `{Normal at Baseline}`)
stacked purely by natural wrapping. We express it as explicit levels, which is
deterministic rather than dependent on the renderer's font metrics, and is what makes the
`merge = "spanners"` argument necessary.

## Reference RTFs do not share a y origin with our DOCX (harness-relevant, not a divergence)
The references start their page content at y=**36, 54 or 72** points depending on which
legacy program wrote them; our DOCX is always 36. The pixel gate absorbs this through
`register()`, but any check that compares absolute y between the two documents will be
wrong by up to 36pt. `header_compare()` did exactly that and produced 7 false failures
before it was changed to cut each document at its own header rule.

## COMPLETE: all 30/30 tables built + verified (2026-07-24)
