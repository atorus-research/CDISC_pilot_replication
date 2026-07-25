I now have a complete, verified picture of the package. Here is the exhaustive developer reference.

---

# `tplyr2` — Exhaustive Developer Reference

> Next-generation successor to `Tplyr`. A declarative, spec-based grammar of clinical summary tables built on **data.table**. Version `0.1.0.9000`, license MIT, `Depends: R (>= 4.1.0)`, `Imports: data.table, jsonlite, purrr, rlang, stringr`, `Suggests: testthat, withr, yaml, knitr, rmarkdown`.
>
> Every fact below was verified against `R/*.R`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md`, and `vignettes/*.Rmd`. Where source and prose disagree, source wins and the discrepancy is flagged.

---

## 0. Package facts

**Exported API (from `NAMESPACE`):**

Construction: `tplyr_spec`, `tplyr_layers`, `group_count`, `group_desc`, `group_shift`, `group_analyze`, `layer_settings`, `f_str`, `label`, `pop_data`, `total_group`, `custom_group`.
Build: `tplyr_build`.
Accessors: `tplyr_header_n`, `tplyr_numeric_data`, `tplyr_stats_data`, `get_data_labels`, `tplyr_layers` (also predicate helpers).
Predicates: `is_tplyr_spec`, `is_tplyr_layer`, `is_pop_data`.
Options: `tplyr2_options`.
Metadata: `tplyr_meta`, `generate_row_ids`, `tplyr_meta_result`, `tplyr_meta_subset`.
ARD: `tplyr_to_ard`, `tplyr_from_ard`.
Serialization: `tplyr_write_spec`, `tplyr_read_spec`.
Post-processing / string helpers: `apply_row_masks`, `collapse_row_labels`, `apply_conditional_format`, `apply_formats`, `str_extract_num`, `str_indent_wrap`, `replace_leading_whitespace`.
S3 `print` methods for: `tplyr_spec`, `tplyr_f_str`, `tplyr_meta`, `tplyr_pop_data`, `tplyr_total_group`, `tplyr_custom_group`.

**Bundled datasets** (CDISC ADaM samples from PHUSE Test Data Factory): `tplyr_adsl` (ADSL), `tplyr_adae` (ADAE), `tplyr_adlb` (ADLB).

---

## 1. Object model, build, and output shape

### 1.1 Workflow

Two steps: **declare** a spec (pure config, no data), then **build** against data.

```r
library(tplyr2)

spec <- tplyr_spec(
  cols = "TRT01P",
  layers = tplyr_layers(
    group_count("SEX", by = "Sex n (%)"),
    group_desc("AGE", by = "Age (Years)")
  )
)

result <- tplyr_build(spec, tplyr_adsl)
```

### 1.2 `tplyr_spec()` — exact signature

```r
tplyr_spec(
  cols,                       # character vector of column (usually treatment) variables
  where = NULL,               # BARE expression, captured via rlang::enexpr (e.g. SAFFL == "Y")
  pop_data = NULL,            # a pop_data() config object
  total_groups = NULL,        # list of total_group() objects
  custom_groups = NULL,       # list of custom_group() objects
  layers = tplyr_layers(),    # list of layer objects
  settings = NULL             # SEE GOTCHA §11 — currently unused by tplyr_build()
)
```

- `cols` accepts **multiple** variables → the output columns are the **cross** of all combinations; column labels join the parts with `" | "` and each combination gets its own `N`.
- `where` is a **bare** (unquoted) expression — the *only* place bare expressions are used. It filters all data before any layer runs. It does **not** apply to `pop_data`.
- Returns an S3 object of class `tplyr_spec`. `print.tplyr_spec` summarizes cols/where/layers.

### 1.3 `tplyr_build()` — exact signature

```r
tplyr_build(spec, data, pop_data = NULL, metadata = FALSE, ...)
```

- `spec` may be a `tplyr_spec` **or a file path** to a `.json`/`.yaml`/`.yml` spec (auto-read via `tplyr_read_spec()`).
- `pop_data` here is the actual population **data.frame** (overrides / satisfies the spec's `pop_data()` config).
- `metadata = TRUE` attaches cell-level traceability (adds a `row_id` column and a `tplyr_meta` attribute).
- `...` are **build-time spec overrides** (see §1.6).
- During the build, `scipen` is temporarily set to `getOption("tplyr2.scipen", 9999L)` to suppress scientific notation, and restored on exit.

### 1.4 What `build()` returns — the output data.frame shape

`tplyr_build()` returns a **plain data.frame** (converted from an internal data.table). Columns fall into four families, always in this left-to-right order after harmonization:

| Family | Naming | Meaning |
|---|---|---|
| Row labels | `rowlabel1`, `rowlabel2`, … | One column per `by` label, per `by` data-variable, then the target variable / stat label. |
| Results | `res1`, `res2`, … | One formatted result column per `cols` level (or per `cols`×stat when `stat_columns` is used, or per `cols`×shift-column level for shift). Each carries a `label` **attribute**. |
| Risk diff | `rdiff1`, `rdiff2`, … | One per risk-difference comparison (count layers only). Each carries a `label` attribute like `"Trt vs Ref"`. |
| Ordering | `ord_layer_index`, `ord_layer_1`, `ord_layer_2`, … | Sort keys (see §5). |

Row-label composition (order): `[by text labels...] → [by data variables...] → [target var / (multi-desc) variable name → stat label]`.

Column (`res*`) `label` attribute conventions:
- Single `cols`: `"<level> (N=<n>)"` when population N is available, else `"<level>"`.
- Multiple `cols`: `"<lvl1> | <lvl2> (N=<n>)"`.
- Shift: `"<cols level> | <shift-column level> (N=<n>)"`.
- `stat_columns`: `"<column group> (N=<n>) | <stat name>"` — the `(N=)` attaches to the column-group segment; renderers split on `" | "` to build a two-level header.

**Attributes attached to the returned data.frame:**
- `attr(result, "header_n")` — data.frame of `cols` levels + `.n` (only when `pop_data` used); retrieve via `tplyr_header_n(result)`.
- `attr(result, "numeric_data")` — named list (keyed by layer index string) of pre-format raw numeric tables; retrieve via `tplyr_numeric_data(result, layer=)`.
- `attr(result, "tplyr_meta")` — named list of `tplyr_meta` objects keyed `"row_id||column"` (only when `metadata = TRUE`).

### 1.5 How a renderer consumes the output

```r
# 1. Sort by layer, then within-layer order
result <- result[order(result$ord_layer_index, result$ord_layer_1), ]

# 2. Read column headers from res* label attributes
sapply(grep("^res", names(result), value = TRUE),
       function(c) attr(result[[c]], "label"))

# 3. (optional) collapse/mask row labels, then drop ord_* columns
display <- result[, !grepl("^ord", names(result))]
```

The build already sorts rows by `ord_layer_index` then the remaining `ord*` columns. Renderers typically: (a) collapse `rowlabel*` into one indented stub column (`collapse_row_labels()`), (b) split the two-level `res*` label attributes, (c) drop the `ord_*` columns.

### 1.6 Build-time overrides (`...`)

`apply_overrides()` merges `...` into a copy of the spec:
- `where = "STRING"` → parsed via `rlang::parse_expr()` (so overrides pass strings, not bare expr).
- `pop_data = list(...)` → merged into the existing `pop_data` config sub-fields rather than replaced.
- Any other name matching a spec field replaces it (e.g. `cols = "TRT01A"`).

```r
tplyr_build(spec, adsl, cols = "TRT01A", where = "SAFFL == 'Y'")
```

### 1.7 Total & custom groups execution order

`tplyr_build()` applies **custom groups first, then total groups** to both analysis and population data (`apply_custom_groups()` then `apply_total_groups()`). Both work by **duplicating** matching rows with the column variable overwritten:
- `total_group("TRT01P", label = "Total")` duplicates *all* rows with `TRT01P := "Total"`.
- Because totals run after customs, the Total column **includes** custom-group rows.

---

## 2. Layer types — constructors and full settings

All four constructors share the same first four positional arguments. `where` is a **bare** expression (captured by `rlang::enexpr`). `target_var`/`by`/`distinct_by`/etc. are **quoted strings**.

```r
group_count  (target_var, by = NULL, where = NULL, settings = layer_settings())
group_desc   (target_var, by = NULL, where = NULL, settings = layer_settings())
group_shift  (target_var, by = NULL, where = NULL, settings = layer_settings())
group_analyze(target_var, by = NULL, where = NULL, analyze_fn, settings = layer_settings())
```

- `target_var`: string, or vector of strings (nested counts / multi-var desc). For **shift** it must be a length-2 **named** vector `c(row = "...", column = "...")` (validated in the constructor).
- `by`: strings matching data columns become **grouping variables**; non-matching strings become **text labels**. Wrap in `label("x")` to force a string to be treated as a text label even if it collides with a column name. `by` may mix both.
- `tplyr_layers(...)`: validates all args inherit `tplyr_layer` and returns them as a list.

### 2.1 `layer_settings()` — full signature (verified defaults)

```r
layer_settings(
  format_strings = NULL,              # named list of f_str()
  stat_columns = NULL,                # count only: named list of f_str() -> one res col per stat per group
  denoms_by = NULL,                   # character vector of denominator grouping vars
  denom_where = NULL,                 # expression filtering the denominator data
  denom_ignore = NULL,                # target-var values excluded from the denominator
  distinct_by = NULL,                 # subject id var for distinct counting
  total_row = FALSE,
  total_row_label = "Total",
  total_row_count_missings = TRUE,
  missing_count = NULL,               # list config (see §2.2.4)
  missing_subjects = FALSE,
  missing_subjects_label = "Missing",
  keep_levels = NULL,                 # character vector of target-var levels to keep
  limit_data_by = NULL,               # restrict completion grid to observed combinations
  custom_summaries = NULL,            # desc: named list of quote() expressions using .var
  stats_as_columns = FALSE,           # desc: transpose stats to columns
  precision_by = NULL,                # desc: precision grouping vars
  precision_on = NULL,                # desc: variable scanned for precision (default = target)
  precision_data = NULL,              # desc: external data.frame with max_int / max_dec
  precision_cap = NULL,               # desc: c(int=, dec=)
  order_count_method = NULL,          # count: "byfactor" | "byvarn" | "bycount"
  ordering_cols = NULL,               # count: which cols level drives bycount
  result_order_var = NULL,            # count: which stat drives bycount (default "n")
  outer_sort_position = NULL,         # count: "asc"/"desc" — SEE GOTCHA §11 (not implemented)
  risk_diff = NULL,                   # count: list(comparisons=, ci=, format=)
  name = NULL                         # layer name (identification only)
)
```

**Applicability by layer type** (from the `layer_settings` roxygen; settings for an inapplicable type are silently ignored):

| Setting | Count | Desc | Shift | Analyze |
|---|:-:|:-:|:-:|:-:|
| `format_strings` | ✅ | ✅ | ✅ | ✅ |
| `stat_columns` | ✅ | | | |
| `denoms_by` | ✅ | ✅ | ✅ | |
| `denom_where` | ✅ | ✅ | ✅ | |
| `denom_ignore` | ✅ | | ✅ | |
| `distinct_by` | ✅ | | ✅ | |
| `total_row` / `total_row_label` / `total_row_count_missings` | ✅ | | | |
| `missing_count` | ✅ | | | |
| `missing_subjects` / `missing_subjects_label` | ✅ | | | |
| `keep_levels` / `limit_data_by` | ✅ | | | |
| `custom_summaries` / `stats_as_columns` | | ✅ | | |
| `precision_by` / `precision_on` / `precision_data` / `precision_cap` | | ✅ | | |
| `order_count_method` / `ordering_cols` / `result_order_var` / `outer_sort_position` | ✅ | | | |
| `risk_diff` | ✅ | | | |
| `name` | ✅ | ✅ | ✅ | ✅ |

---

### 2.2 Count layers (`group_count`)

#### 2.2.1 Basic count and defaults

Default format is `f_str("xx (xx.x%)", "n", "pct")`; the denominator is the column total. Zero-count cells are always filled (data completion via cross-join of all factor levels / observed values).

```r
tplyr_spec(cols = "TRT01P", layers = tplyr_layers(group_count("RACE")))
```

**Count statistics** available in format strings: `n`, `pct`, `total` (denominator), `distinct_n`, `distinct_pct`, `distinct_total`. Count-layer format strings use the list key **`n_counts`**:

```r
group_count("DCDECOD",
  settings = layer_settings(
    format_strings = list(n_counts = f_str("xxx (xxx.x%)", "n", "pct"))
  ))
```

#### 2.2.2 Distinct vs event counts

```r
group_count("AEDECOD",
  settings = layer_settings(
    distinct_by = "USUBJID",
    format_strings = list(
      n_counts = f_str("xxx (xx.x%) [xxx]", "distinct_n", "distinct_pct", "n")
    )))
```
`n` = event (row) count; `distinct_n` = unique-subject count (`uniqueN(USUBJID)`). The common AE cell format is `xxx (xx.x%) [xxx]` = distinct subjects (pct) [events].

#### 2.2.3 Total rows, missing counts, missing subjects, keep_levels

```r
layer_settings(total_row = TRUE, total_row_label = "Overall Total",
               total_row_count_missings = TRUE)
```
- **Total row** sums `n` across levels (and recomputes `distinct_n` from raw data to avoid double-counting). `total_row_count_missings = FALSE` excludes NA/`missing_values` from the total.
- **`missing_count`** — a list controlling a "Missing" row:
  ```r
  missing_count = list(
    missing_values = character(0),  # extra values treated as missing besides NA
    denom_exclude  = FALSE,
    sort_value     = Inf,           # sort position of the missing row
    label          = "Missing",
    format         = NULL           # optional per-row f_str (serialized specially)
  )
  ```
  Minimal usage: `missing_count = list(label = "Missing")`. Missing is detected as `is.na(target)` OR `target %in% missing_values`.
- **`missing_subjects = TRUE`** (requires `pop_data`) — adds a row counting subjects present in the population but **absent** from the analysis data (subject-level `fsetdiff` anti-join when `distinct_by` is set). Label via `missing_subjects_label` (default `"Missing"`; common alt `"Not reported"` / `"No events reported"`).
- **`keep_levels`** — after completion, filter to only these target-var levels.

#### 2.2.4 `stat_columns` (new in 0.1.0.9000, issue #10)

Named list of `f_str()`; each entry produces its **own** result column per `cols` group. Result columns interleave **column-group-major** (each arm's stat columns adjacent). Takes precedence over `format_strings` (a warning is raised if both are set).

```r
group_count("AEDECOD",
  settings = layer_settings(
    distinct_by = "USUBJID",
    stat_columns = list(
      "n (%)" = f_str("xxx (xx.x%)", "distinct_n", "distinct_pct"),
      "E"     = f_str("xxx", "n")
    )))
# res1/res2 = arm1's "n (%)" and "E"; res3/res4 = arm2's; ...
# attr(res1,"label") == "Placebo (N=86) | n (%)"
```
Constraints (enforced by `validate_stat_columns` / `validate_stat_columns_alignment`): every entry must be named; names unique; names may **not** contain `" | "` or `"(N="`; when **any** layer in a spec uses `stat_columns`, **all** layers must be count layers with the **same** stat names (else error). Works with nesting, by-vars, total/missing rows, risk_diff, metadata (populates `tplyr_meta$statistic`), serialization, and ARD.

#### 2.2.5 Nested counts (2-level)

Pass a vector of two variable names. Outer level → `rowlabel1`, inner → `rowlabel2`; outer rows have blank inner label. `ord_layer_2` encodes nesting depth (1 = outer, 2 = inner, 0 = total/missing).

```r
group_count(c("AEBODSYS", "AEDECOD"),
  settings = layer_settings(
    distinct_by = "USUBJID",
    format_strings = list(n_counts = f_str("xxx (xx.x%)", "distinct_n", "distinct_pct")),
    limit_data_by = c("AEBODSYS", "AEDECOD"),   # PTs only appear under their real SOC
    total_row = TRUE, total_row_label = "Any adverse event"
  ))
```
`limit_data_by` restricts the completion grid to observed combinations (essential for AE tables). Nested denominators: `denoms_by` may be a plain character vector (same for all levels) **or a list** (per-level). *Effective nesting depth is 2* — see §11.

#### 2.2.6 Risk difference

```r
risk_diff = list(
  comparisons = list(c("Xanomeline High Dose", "Placebo"),   # c(treatment, reference)
                     c("Xanomeline Low Dose",  "Placebo")),
  ci = 0.95,                                                  # default 0.95
  format = f_str("xx.x (xx.x, xx.x)", "rdiff", "lower", "upper")  # default shown
)
```
Uses `stats::prop.test(x, n, conf.level = ci, correct = FALSE)` (Wald, no continuity correction). Variables: `rdiff` (percentage-point diff = (p1−p2)·100), `lower`, `upper` (CI·100), `p_value`. Each comparison → one `rdiff<i>` column with `label` = `"trt vs ref"`. Computed on `distinct_*` proportions when `distinct_by` is set. **Computed before** total/missing rows are appended, so those rows show empty risk diff. Requires ≥1 `cols` variable.

---

### 2.3 Descriptive statistics layers (`group_desc`)

**Built-in statistics** (always computed): `n` (non-missing count), `mean`, `sd`, `median`, `var`, `min`, `max`, `iqr`, `q1`, `q3`, `missing` (count of NA). Also merged in: `total` (denominator) and `pct`.
- All use `na.rm = TRUE`; `min`/`max` operate on **finite** values (all-NA → `NA_real_`, never `±Inf`).
- `sd`/`var` require n>1 else `NA`.
- `iqr`/`q1`/`q3` use `getOption("tplyr2.quantile_type", 7)` (7 = R default; 3 ≈ SAS `PROC UNIVARIATE`).

Format strings are a **named list**; each **name is the row label**, each value an `f_str`. Default (when `format_strings` omitted):

```r
list(
  "n"         = f_str("xxx", "n"),
  "Mean (SD)" = f_str("xx.x (xx.xx)", "mean", "sd"),
  "Median"    = f_str("xx.x", "median"),
  "Q1, Q3"    = f_str("xx.x, xx.x", "q1", "q3"),
  "Min, Max"  = f_str("xx, xx", "min", "max"),
  "Missing"   = f_str("xxx", "missing")
)
```

**Selecting which stats appear / per-stat format & multi-line cells:** you fully control the set and order by the `format_strings` list; ordering follows list order (`ord_layer_1` = stat position). Each row is one cell format; combine multiple stats per line as needed (`"Mean (SD)"` puts mean and sd in one cell).

**Multiple target variables** — pass a vector; each variable gets its own block, the variable name becomes an extra `rowlabel` column. `ord_layer_1 = var_index*100 + stat_order`.

```r
group_desc(c("AGE","HEIGHTBL","WEIGHTBL"), settings = layer_settings(format_strings = list(...)))
```

**Custom summaries** — named `quote()` expressions using `.var` as the group's value vector; errors caught → `NA_real_`. Precedence: **layer-level > session-level (`tplyr2.custom_summaries`) > built-in** (a custom `mean` overrides the built-in mean).

```r
custom_summaries = list(geo_mean = quote(exp(mean(log(.var[.var > 0]), na.rm = TRUE))))
```

**`stats_as_columns = TRUE`** — transposes so treatment groups become rows and stat names become columns (opposite of `stat_columns`).

**Auto-precision** — see §3.4. Controlled by `precision_on`, `precision_by`, `precision_cap`, `precision_data`.

---

### 2.4 Shift layers (`group_shift`)

Cross-tabulation: `target_var = c(row = "BNRIND", column = "ANRIND")`. The shift `column` variable becomes an additional column dimension folded into `res*` alongside `cols`. `rowlabel(last)` holds the row (baseline) variable value; result columns are treatment×post-baseline combinations, labels `"<trt> | <ANRIND level> (N=n)"`.

```r
group_shift(
  c(row = "BNRIND", column = "ANRIND"),
  by = c("PARAM", "VISIT"),
  settings = layer_settings(format_strings = list(n_counts = f_str("xx (xxx.x%)", "n", "pct")))
)
```
- Default format = `f_str("xx (xx.x%)", "n", "pct")` (via `get_count_format`).
- **Use factors** on row/column vars to control completeness (zero cells appear) and order.
- Denominator default = spec `cols` only (not the shift column); supports `denoms_by`, `denom_where`, `denom_ignore`, `distinct_by`.
- **Limitations:** no nesting, no total/missing rows, no risk difference, no result-based sorting, **ignores `stat_columns`**. For text-style shifts ("Low to High"), derive a variable and use `group_count`.

---

### 2.5 Analyze layers (`group_analyze`) — custom / injected statistics

The escape hatch for arbitrary computations, incl. **feeding external model results** (ANCOVA/MMRM/CMH p-values, LSMEANS) into a display. `analyze_fn` is called **once per group** (each combination of `cols` + `by` data-variables):

```r
analyze_fn = function(.data, .target_var) { ... }
```
`.data` = data.frame subset for the group; `.target_var` = the target variable name(s). Two modes:

**Mode A — format-strings mode.** Return a **single-row** data.frame of named numeric columns; supply `format_strings`. Each format-string entry becomes an output row (name = row label), referencing the returned column names.

```r
geo_fn <- function(.data, .target_var) {
  v <- .data[[.target_var]]; v <- v[!is.na(v) & v > 0]
  data.frame(geo_mean = exp(mean(log(v))), geo_sd = exp(sd(log(v))))
}
group_analyze("AVAL", by = "Urate (umol/L)", where = AVISIT == "Baseline",
  analyze_fn = geo_fn,
  settings = layer_settings(format_strings = list(
    "Geometric Mean" = f_str("xxx.xx", "geo_mean"),
    "Geometric SD"   = f_str("xxx.xx", "geo_sd"))))
```

**Mode B — pre-formatted mode.** Omit `format_strings`; return a data.frame with `row_label` and `formatted` (character) columns. tplyr2 does **no** padding/rounding — you own alignment. Ideal for injecting externally-computed strings:

```r
# Inject pre-computed ANCOVA LSMEANS diff + p-value keyed by the group columns
ancova_fn <- function(.data, .target_var) {
  key <- unique(.data$TRTP)                # group identity from the subset
  m <- lsmeans_results[lsmeans_results$TRTP == key, ]
  data.frame(
    row_label = c("LS Mean (SE)", "Diff vs PBO (95% CI)", "p-value"),
    formatted = c(sprintf("%.1f (%.2f)", m$lsmean, m$se),
                  sprintf("%.1f (%.1f, %.1f)", m$diff, m$lcl, m$ucl),
                  format.pval(m$pval, digits = 3))
  )
}
group_analyze("CHG", by = "AVISIT", analyze_fn = ancova_fn)
```
Because `analyze_fn` receives each group's rows, join in your external model table by the group's key columns (read them off `.data`). Multi-layer specs stack analyze blocks by `ord_layer_index`. Errors surface immediately (build stops) — build error handling into the function.

---

## 3. `f_str()` — the formatting language

```r
f_str(format_string, ..., empty = NULL)
```
- `format_string`: template with **format groups** (runs of `x`/`X`/`a`/`A`, optional `.decimal`, optional `+N`) separated by **literal** text.
- `...`: quoted statistic names, one per format group (count of groups must equal count of names, else error at construction).
- `empty`: replacement when *all* group values are NA — only the `.overall` key is honored: `empty = c(.overall = "---")`.

Standalone use via `apply_formats(fmt, ..., precision = NULL)`.

### 3.1 Format codes

| Code | Meaning |
|---|---|
| `x` | Fixed-width digit slot. Each `x` reserves one char; numbers narrower are **left-padded with spaces**; wider numbers print in full (never truncated). |
| `.`  | Separates integer part from decimal part. Chars after `.` set decimal precision. |
| `X` | Like `x` but **parenthesis hugging**: leading pad-spaces shift to *after* the trailing literal so the delimiter hugs the number. Total width preserved. |
| `a` | **Auto-precision** width from data (like `x`, but width computed at build). |
| `A` | Auto-precision **with** hugging (uppercase `a`). |
| `+N` | Suffix on int or dec part: add N to the auto-determined width, e.g. `a.a+1`, `a+2.a`. |

Everything between groups is literal (`(`, `)`, `[`, `]`, `%`, `,`, spaces, etc.). Parsing regex: `[xXaA]+(\+\d+)?(\.[xXaA]+(\+\d+)?)?`.

```r
f_str("xx.xx (xx.xx)", "mean", "sd")     # two decimals each, parens literal
f_str("xxx", "n")                        # integer count
f_str("xx (xx.x%)", "n", "pct")          # count with percent (default count fmt)
f_str("xxx (xx.x%) [xxx]", "distinct_n", "distinct_pct", "n")
```

### 3.2 Variables usable per layer

- **Count / Shift:** `n`, `pct`, `total`, `distinct_n`, `distinct_pct`, `distinct_total`.
- **Desc:** `n`, `mean`, `sd`, `median`, `var`, `min`, `max`, `iqr`, `q1`, `q3`, `missing` (+ any custom-summary names, + `total`, `pct`).
- **Risk diff:** `rdiff`, `lower`, `upper`, `p_value`.
- **Analyze:** any column name returned by `analyze_fn`.

Unknown stat names → a **warning** (not error) at build (custom summaries may add names).

### 3.3 Rounding

`format_number_vec()` uses `formatC(..., format="f"/"d", width=...)`. Rounding respects `getOption("tplyr2.IBMRounding", FALSE)`: `FALSE` → R banker's rounding; `TRUE` → round-half-away-from-zero (`sign(x)*floor(abs(x)*10^d + 0.5)/10^d`, SAS-like). NA groups render as blank space of the full field width.

### 3.4 Auto-precision (desc)

Set at least one `a`/`A` in a format. tplyr2 scans `precision_on` (default = target var) to get `max_int` (integer digits) and `max_dec` (meaningful decimal places, trailing zeros stripped). `a` → `max_int`/`max_dec`; `+N` adds N. Standard practice: mean `a.a+1`, sd `a.a+2`.
- `precision_by` computes precision independently per group (e.g. per `PARAMCD`).
- `precision_cap = c(int=, dec=)` caps resolved widths (layer-level overrides global `tplyr2.precision_cap`).
- `precision_data` supplies precision directly (data.frame with `max_int`, `max_dec` + any `precision_by` cols) — bypasses scanning.

```r
group_desc("AVAL", by = "PARAMCD",
  settings = layer_settings(
    precision_by = "PARAMCD", precision_on = "AVAL",
    precision_cap = c(int = 4, dec = 3),
    format_strings = list(
      "Mean (SD)" = f_str("a.a+1 (A.A+2)", "mean", "sd", empty = c(.overall = "")),
      "Median"    = f_str("a.a+1", "median", empty = c(.overall = "NE")),
      "Min, Max"  = f_str("a.a, a.a", "min", "max"))))
```

---

## 4. Descriptive statistics — controlling appearance

Covered by §2.3 and §3.4. Key points:
- **Which stats appear + order** = the `format_strings` list contents/order.
- **Per-stat precision** = each entry's `f_str` template (auto or fixed).
- **Multi-line cells** = multiple stats per format group line (e.g. `"Mean (SD)"`, `"Q1, Q3"`, `"Min, Max"`).
- **Quantile algorithm** = `tplyr2_options(quantile_type = 3)` for SAS parity.
- **Empty (all-NA) cells** = `empty = c(.overall = "NE")` per `f_str`.

---

## 5. Sorting / ordering

### 5.1 Ordering columns generated

Internally `ordindx`/`ord1`/`ord2` are renamed at the end of build (`rename_ord_columns`):
- `ord_layer_index` — integer layer position (layer 1, 2, …).
- `ord_layer_1` — within-layer sort key. Count: target-var sort position; Desc: stat position (list order, `var_index*100+stat_order` for multi-var); Shift: row order.
- `ord_layer_2` — nested counts only: nesting depth (1 outer, 2 inner, 0 total/missing).

Build sorts rows by `c("ord_layer_index", <remaining ord cols>)`. When ≥10 result columns exist, `sort_by_numeric_suffix()` orders `res*`/`rdiff*` by numeric suffix (fixes `res10` before `res2` — a documented bug fix in NEWS).

### 5.2 Count-layer sort control

`order_count_method` (in `layer_settings`):

| Method | Behavior |
|---|---|
| `NULL` (default) | Auto-detect: factor levels → VARN companion (`<VAR>N`) → alphabetical. |
| `"byfactor"` | Factor level order of the target variable. |
| `"byvarn"` | Numeric companion column `<VAR>N` (e.g. `RACEN` for `RACE`). |
| `"bycount"` | Descending by count. `ordering_cols` picks which `cols` level drives it; `result_order_var` picks which stat (default `"n"`). |

`compute_var_order()` implements the priority. Nested outer direction is via `outer_sort_position` (**declared but not wired in — see §11**).

```r
group_count("RACE", settings = layer_settings(order_count_method = "byvarn"))
group_count("DCDECOD", settings = layer_settings(
  order_count_method = "bycount", ordering_cols = "Placebo", result_order_var = "n"))
```

The standard workflow is **build → sort by `ord_*` → drop `ord_*`**. Factor levels on the target/`by`/shift variables also control completion and display order globally.

---

## 6. Denominators

All in `layer_settings` (count/shift; desc supports `denoms_by`/`denom_where`).

| Control | Effect | Default |
|---|---|---|
| (default) | Column-total denominator (one per `cols` level). | `cols` |
| `denoms_by` | Explicit denominator grouping vars (any data vars). | `cols` |
| `denom_where` | Bare/`quote()`d expression filtering the *denominator* data only (numerator unchanged; can push pct >100%). | none |
| `denom_ignore` | Target-var **values** removed from the denominator (row still shown). | none |
| `distinct_by` | Enables subject-level `distinct_n`/`distinct_pct`/`distinct_total`, with their own denominators. | NULL |
| `pop_data()` + build `pop_data=` | Denominators & header N come from a separate population dataset. | analysis data |
| `missing_subjects` | Row for pop subjects absent from analysis data (needs `pop_data`). | FALSE |

**Header N & population data:**
```r
tplyr_spec(
  cols = "TRTA",
  pop_data = pop_data(cols = c("TRTA" = "TRT01A"),   # c(analysis_name = pop_name)
                      where = SAFFL == "Y"),          # pop_data's OWN where (spec$where does NOT apply)
  layers = tplyr_layers(group_count("AEBODSYS",
    settings = layer_settings(distinct_by = "USUBJID"))))
result <- tplyr_build(spec, tplyr_adae, pop_data = tplyr_adsl)
tplyr_header_n(result)   # data.frame of cols levels + .n
```
- `pop_data(cols, where = NULL)`: `cols` named → `c(spec_col = pop_col)`; unnamed → positional. `where` is a bare expression.
- Build renames the pop columns to the spec names, then computes `header_n` from the (filtered, group-augmented) population.
- **The spec-level `where` never filters pop_data** — pop filtering is done exclusively via `pop_data(where=)`.

---

## 7. Options and project-level config

`tplyr2_options(...)` — get all (no args) or set (named args; auto-prefixed with `tplyr2.`). Stored as ordinary R options, so `options(tplyr2.X = )` works too.

| Option | Type | Default | Purpose |
|---|---|---|---|
| `tplyr2.IBMRounding` | logical | `FALSE` | Round-half-away-from-zero (SAS) vs banker's. |
| `tplyr2.quantile_type` | integer | `7L` | `quantile()` type for q1/q3/iqr (3 ≈ SAS). |
| `tplyr2.precision_cap` | `c(int=,dec=)` | `NULL` | Global auto-precision cap (layer-level overrides). |
| `tplyr2.custom_summaries` | named list of `quote()` | `NULL` | Session-wide custom desc stats. |
| `tplyr2.scipen` | integer | `9999L` | scipen used during build to suppress sci notation. |

**Project config file mapping.** These are session options, so a project config is idiomatically the top of an analysis script or `.Rprofile`:
```r
# project setup
tplyr2_options(
  IBMRounding      = TRUE,
  quantile_type    = 3L,
  precision_cap    = c(int = 5, dec = 3),
  custom_summaries = list(
    geo_mean = quote(exp(mean(log(.var[.var > 0]), na.rm = TRUE))),
    cv       = quote(sd(.var, na.rm = TRUE) / mean(.var, na.rm = TRUE) * 100))
)
```
Table structure itself is externalized via **spec serialization** (§8.3), the natural "config file" for table definitions (`specs/table_14_1.json`). `get_data_labels(data)` extracts the `label` attribute of each column (Haven/CDISC) as a named vector for header building.

---

## 8. Metadata / ARD / traceability

### 8.1 Cell-level metadata

`tplyr_build(spec, data, metadata = TRUE)` adds:
- a `row_id` column (`generate_row_ids()`: `<layer>_<rowlabel1>_<rowlabel2>...`, e.g. `1_F`, `2_Mean (SD)`); deterministic.
- `attr(result, "tplyr_meta")`: named list keyed `"row_id||column"`.

```r
result <- tplyr_build(spec, tplyr_adsl, metadata = TRUE)
m <- tplyr_meta_result(result, "1_F", "res1")   # a tplyr_meta object
src <- tplyr_meta_subset(result, "1_F", "res1", tplyr_adsl)  # actual source rows
```

`tplyr_meta(names, filters, layer_index, anti_join = NULL, statistic = NULL)`:
- `filters` = list of R call expressions (combined with `&`) reproducing the source subset.
- `names` = variables referenced.
- `statistic` = which stat the cell shows (populated for `stat_columns`).
- `anti_join` = for missing-subjects rows: `tplyr_meta_anti_join(join_meta, on = distinct_by)` — needs `pop_data` passed to `tplyr_meta_subset(..., pop_data=)`.

Filters incorporate: column-level value (with total/custom-group translation — totals → no filter, custom → `%in%`), row target/by values, spec `where`, layer `where`. Row types classified: `normal`/`total`/`missing`/`missing_subjects`, each with distinct filter logic (e.g. total omits the target filter). Common use: QC cell verification and Shiny drill-down.

### 8.2 Numeric data & ARD

- `tplyr_numeric_data(result, layer = NULL)` → raw pre-format stats (list keyed by layer, or one layer's data.frame). `tplyr_stats_data(result, layer, statistic)` → the layer's data if it contains that stat.
- `tplyr_to_ard(result)` → long CDISC-style ARD: `analysis_id` (layer), grouping columns, `stat_name`, `stat_value`. Internal columns (`.sort_*`, `formatted*`, `rowlabel*`, etc.) are dropped; captures **every** computed stat (not just displayed).
- `tplyr_from_ard(ard, spec)` → reconstruct a formatted table by re-applying the spec's format strings (round-trips; supports `stat_columns` for count, ignores it for shift). The spec acts as a reusable formatting template — the same ARD can be re-rendered with a different spec.

```r
ard <- tplyr_to_ard(result)
rebuilt <- tplyr_from_ard(ard, spec)
```

### 8.3 Serialization

`tplyr_write_spec(spec, path)` / `tplyr_read_spec(path)` — format by extension (`.json` via jsonlite; `.yaml`/`.yml` via the suggested `yaml` package). What is preserved:
- Expressions (`where`, `denom_where`, custom summaries) → deparsed to `{"_expr": "..."}`, re-parsed with `rlang::parse_expr()`.
- `f_str` → `{format_string, vars, _class:"tplyr_f_str"[, empty]}`.
- `label()` in `by` → `{value, _type:"label"}`.
- `analyze_fn` → `{"_fn": "<deparsed function>"}` (closures capturing environment vars may not survive — write self-contained functions).
- `precision_data` → list of columns → data.frame.
- Shift `target_var` names (`row`/`column`) restored on read if lost.
- `tplyr_build()` accepts a spec **path** directly.

---

## 9. Migration from `Tplyr` v1 (the legacy repo idioms)

Three shifts: **declare not pipe**, **quote variable names**, **collect settings in one object**.

| Tplyr v1 | tplyr2 |
|---|---|
| `tplyr_table(data, treat_var)` | `tplyr_spec(cols = "treat_var")` — data at build time |
| `%>% add_layer(...)` | `layers = tplyr_layers(...)` inside the spec |
| `group_count(SEX)` (bare symbol) | `group_count("SEX")` (quoted string) |
| `group_desc(AGE)` | `group_desc("AGE")` |
| `group_count(vars(AEBODSYS, AEDECOD))` | `group_count(c("AEBODSYS", "AEDECOD"))` |
| `group_shift(vars(...))` | `group_shift(c(row = "BNRIND", column = "ANRIND"))` |
| `set_format_strings(f_str("xx (xx.x%)", n, pct))` | `format_strings = list(n_counts = f_str("xx (xx.x%)", "n", "pct"))` in `layer_settings()` |
| `set_distinct_by(USUBJID)` | `distinct_by = "USUBJID"` |
| `set_denoms_by(TRT01P)` | `denoms_by = "TRT01P"` |
| `set_where(...)` | `where =` in layer or spec (bare expr) |
| `set_order_count_method("bycount")` | `order_count_method = "bycount"` |
| `set_ordering_cols("X")` | `ordering_cols = "X"` |
| `add_total_group(...)` | `total_group()` in spec `total_groups` |
| `set_pop_data(...)` | `pop_data()` config in spec + `pop_data=` at build |
| `set_nest_count(TRUE)` display | `collapse_row_labels(result, nest = TRUE)` |
| `build()` | `tplyr_build(spec, data)` |
| `f_str("xx (xx.x%)", n, pct)` | `f_str("xx (xx.x%)", "n", "pct")` (quoted names) |

**Behavioral / structural changes to flag when porting legacy code:**
- Variable names are **strings everywhere** except `where` (bare). This is required for JSON/YAML serialization.
- Treatment variable is absorbed into `cols` (no separate `treat_var`); `cols` can be multi-variable (crossed columns).
- Count format-string list key is **`n_counts`**; desc keys are the row labels.
- Nesting depth is nominally recursive but **effectively 2-level** (v1 was fixed 2-level).
- The **output structure is compatible**: `rowlabel*`, `res*` (with `label` attrs), `rdiff*`, `ord_layer_*` — downstream renderers should need minimal changes (v1 used `row_label*`/`ord_layer_*`; tplyr2 uses `rowlabel*`).
- New in tplyr2: `stat_columns` (multi-column-per-arm), `group_analyze()`, spec serialization, metadata, ARD conversion, `tplyr_numeric_data()`. data.table backend gives ~4× median speedup (shift up to ~180×+).

---

## 10. Analyze layer for injected external statistics

See §2.5 for the contract. For **feeding external model output** (ANCOVA/MMRM/CMH/LSMEANS) into a display:
- Prefer **pre-formatted mode** when the model output is already display-ready or needs bespoke p-value formatting: return `row_label` + `formatted`.
- Prefer **format-strings mode** when you want tplyr2 to align/round numeric model outputs: return a single-row numeric data.frame with named columns and map them via `format_strings`.
- `analyze_fn` runs once per (`cols` × `by`-data-var) group; read the group's key columns off `.data` to join your external results table.
- No `distinct_by`/denominators/risk_diff apply to analyze layers; the layer only supports `format_strings`, `by`, `where`, `name`.
- Analyze layers integrate in multi-layer specs; each gets its own `ord_layer_index`. Errors abort the build (wrap edge cases yourself). `NA_real_` renders as blank of correct width in format-strings mode.

---

## 11. Limitations, gotchas, and TODOs

Verified from source/NEWS/scratch:

1. **`spec$settings` is inert.** The `settings` parameter of `tplyr_spec()` is stored but **never read** by `tplyr_build()`. No `spec$settings` reference exists in `R/build.R`.
2. **`outer_sort_position` is not implemented.** It is a `layer_settings()` parameter and is serialized, but there is **no build/ordering code that uses it** (the sort vignette describes it aspirationally). Grep confirms no usage in `build_count.R`/`ordering.R`.
3. **`order_count_method = "bycount"` may not actually reorder (single-variable count path).** `compute_count_sort_keys()` computes `.ord_tv`/`.ord_by_*` keys, but those columns are **not** carried into `cast_to_wide()` (only `rowlabel*` + value columns are), so they are dropped. After casting, `ord1` is `seq_len(.N)` (dcast order). The post-cast override sets `ord1` via `compute_var_order()` for `"byfactor"`/`"byvarn"` only; the `"bycount"` branch is an explicit **no-op** whose comment references a "pre-sort below" that does not exist in that path. **Net effect (from reading, not execution): `byfactor`/`byvarn` reorder; `bycount` likely leaves alphabetical order despite the sort vignette's claim of negated counts.** *Verify empirically before relying on bycount.*
4. **Nesting is effectively 2-level.** The counting loop iterates all `target_vars`, but the interleave/sort (`sort_nested`) and total/missing handling reference `target_vars[1]` (outer) and `target_vars[2]` (inner) only. DESIGN.md aspires to recursive N-level; requirements/vignettes treat it as 2-level. Deeper nesting is not covered by tests or docs.
5. **Shift layer restrictions:** no nesting, no total/missing rows, no risk difference, no result-based sorting, and **`stat_columns` is ignored** (shift uses only `get_count_format()`, i.e. a single format).
6. **`stat_columns` cross-layer rule:** if any layer uses `stat_columns`, all layers must be count layers with identical stat names, or `tplyr_build()` errors (positional res-column alignment constraint). Stat names cannot contain `" | "` or `"(N="`.
7. **`empty` supports only `.overall`.** `apply_formats()` only checks the `.overall` key; per-group empty keys are not honored.
8. **`denom_where` can produce percentages > 100%** by design (shifts denominator frame without changing numerator). Same with `denom_ignore`.
9. **`missing_subjects` requires `pop_data`**; returns NULL otherwise. Its metadata anti-join is only built when `distinct_by` is set.
10. **`stats_as_columns` (desc) is heuristic** for extra non-stat label columns (only handles a single constant non-stat label, shifting it to `rowlabel1`); complex `by` layouts may not transpose cleanly.
11. **Lifecycle: experimental.** Version `0.1.0.9000`. README badge marks the package experimental.
12. **Known NEWS fix already in place:** result/rdiff columns now sorted by numeric suffix (was lexicographic `res10` before `res2`).
13. **Deferred follow-ups noted in `scratch/issue-10-implementation-plan.md`** (not implemented): `distinct = TRUE` risk_diff option, an exported header-frame parser helper, relaxing the mixed-layer `stat_columns` guard, and a `split_cell_columns()` convenience.

**Key file map** for further digging: build orchestration `R/build.R`; count `R/build_count.R` (incl. `cast_to_wide`, `build_col_labels`, nested helpers); desc `R/build_desc.R`; shift `R/build_shift.R`; analyze `R/build_analyze.R`; formatting `R/f_str.R`; precision `R/precision.R`; ordering `R/ordering.R`; risk diff `R/risk_diff.R`; metadata `R/metadata.R`; ARD `R/ard.R`; serialization `R/serialize.R`; validation `R/validate.R`; settings `R/layer_settings.R`; constructors `R/layers.R`; groups/pop `R/pop_data.R`; options `R/options.R`; post-processing `R/post_process.R`.