# tplyr2 FR: group_shift needs a column-wise (per-shift-column) denominator option

**Package:** tplyr2 · found building CDISC pilot 14-6.04 (lab shift tables)
**Filed:** atorus-research/tplyr2 **#18** (the column-wise denominator FR) and **#19**
(the `denoms_by` replaces-not-augments docs/DX clarification), 2026-07-24.
**REVISIT when #18 lands:** rewrite `tables/t_14_6_04.R` to use native `group_shift`
with the new column-wise denominator instead of the crossed `cols=c("TRTP","BNRIND")`
count + manual `n`-row assembly. Same applies to the upcoming shift+CMH tables
**14-6.05 / 14-6.06** — build them on the crossed-cols pattern now, migrate later.
See notes/REVISIT.md.

## Context
14-6.04 is a *transposed* shift: baseline status (Normal/High) is a column dimension
crossed with treatment; shift-to (Normal/High) are rows. Each cell's percent is
**column-wise** — out of that baseline group's N (e.g. of the 81 subjects Normal at
baseline on Placebo, how many shifted to Normal/High).

## What group_shift does
`group_shift(c(row="ANRIND", column="BNRIND"))` produces the right layout and the
right **counts**, but its denominator is the **cols (treatment) total only** — it does
not offer a per-shift-column (per-baseline-group) denominator. `denoms_by` on a shift
layer did not yield the column-wise denominator either (gave odd values). So the
percentages come out wrong for this common shift display.

## Workaround used
Express the shift as a **crossed-cols count** instead:
`tplyr_spec(cols = c("TRTP","BNRIND"), group_count("ANRIND", by = c("PARAM","VISIT")))`.
The crossed column IS the baseline group, so the default per-column total is exactly the
column-wise denominator. Two subtleties worth documenting regardless of any FR:
- **`denoms_by` REPLACES the default cols grouping** (it doesn't add to it). To scope a
  crossed-cols count's denominator per by-group you must list *all* of them:
  `denoms_by = c("PARAM","VISIT","TRTP","BNRIND")`. Just `c("PARAM","VISIT")` collapses
  the denominator to the whole param×visit (wrong), and omitting them entirely uses the
  column total across *all* by-groups (also wrong).

## Ask
A `group_shift` denominator option for column-wise (per-shift-column) percentages —
i.e. "% within baseline group" — so the shift layer covers this standard display without
falling back to a crossed-cols count. Optionally, clarify in docs that `denoms_by`
replaces (not augments) the cols grouping.
