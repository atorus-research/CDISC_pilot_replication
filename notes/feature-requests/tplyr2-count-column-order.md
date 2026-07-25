> **RESOLVED** — filed as tplyr2 issues #13/#14 and fixed in merged PR #15
> (atorus-research/tplyr2). This project now uses the native behaviour/settings.

# tplyr2 FR/bug: count layers ignore `cols` factor levels for result-column order

**Severity:** medium (correctness footgun) · **Package:** tplyr2 0.1.0.9000

## What happens
When `cols` is a **factor**, `group_desc` orders the `res*` columns by the factor
levels, but `group_count` orders them **alphabetically by the `cols` value**,
ignoring the factor levels. The two layer types disagree within the same spec.

## Repro
```r
ARMS <- c("Placebo","Xanomeline Low Dose","Xanomeline High Dose","Total")
d <- dplyr::bind_rows(adsl, dplyr::mutate(adsl, TRT01P = "Total")) |>
  dplyr::mutate(TRT01P = factor(TRT01P, levels = ARMS),
                BMIBLGR1 = factor(BMIBLGR1, levels = c("<25","25-<30",">=30")))

b <- tplyr_build(
  tplyr_spec(cols = "TRT01P",
             layers = tplyr_layers(group_count("BMIBLGR1"))), d)

sapply(grep("^res", names(b), value = TRUE), \(c) attr(b[[c]], "label"))
#> res1 "Placebo (N=86)"  res2 "Total (N=254)"
#> res3 "Xanomeline High Dose (N=84)"  res4 "Xanomeline Low Dose (N=84)"
```
Expected `res1..res4` = Placebo / Xan Low / Xan High / Total (factor order), as
`group_desc` produces.

## Impact
Any table that mixes desc and count layers, or that assumes `res1` == first arm,
silently scrambles columns. We now defensively reorder `res*` by parsing each
column's `(N=)` label — but that only works because count layers *do* set a
`label` attribute (desc layers don't, so we can't cross-check them the same way).

## Suggested fix
Make `group_count` honor `cols` factor levels for result-column ordering,
matching `group_desc`. Secondarily, set the `label` attribute consistently on
`res*` columns for **all** layer types so renderers have one reliable way to map
a result column to its `cols` level.

## Related
Row-subsetting a build (`b[order(b$ord_layer_1), ]`) strips the per-column
`label` attributes (base-R `[.data.frame` behavior) — capture labels before
reordering rows. Not a tplyr2 bug, but a sharp edge worth a note in the docs.
