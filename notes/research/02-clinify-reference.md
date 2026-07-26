I have everything I need. I've read all source files, all vignettes, metadata, and verified the numeric defaults by running R against the installed package. Here is the complete reference.

---

# `clinify` — Exhaustive Developer Reference

`clinify` (Atorus Research; version `0.3.1` in source, Apache-2.0) extends **flextable** and **officer** to produce regulatory clinical tables/listings/figures with Microsoft Word `.docx` as the first-class output. Lifecycle is *experimental*; it is on CRAN.

**Verification note:** All signatures, defaults, and option names below were read directly from `R/*.R`, the `man/*.Rd`, the vignettes, and `NAMESPACE`. Numeric page defaults were confirmed by running R against the installed package (`clin_default_table_width()` returns `9.725`). Two discrepancies between the docs and the shipped code are flagged explicitly in the *Limitations / Gotchas* section — do not trust the older vignette prose over the source.

---

## 0. Mental model (read this first)

The single most important architectural fact:

> **A `clintable` IS a `flextable`** with one extra element attached: `x$clinify_config` (a plain list). Almost every `clin_*()` verb does nothing except *write a field into `clinify_config` and return the object.* No pagination, slicing, styling, borders, or page numbers are computed until a **terminal action** runs: `print()`, `knit_print()`, `clindoc()` / `as_clindoc()`, or `write_clindoc()`.

Consequences:
- Native flextable/officer functions (`flextable::align()`, `bold()`, `color()`, `bg()`, `merge_v()`, `officer::*`) work directly on a `clintable` at any point in the pipe. This is used heavily (e.g. `flextable::align(align="center", part="header")`).
- Order of `clin_*()` verbs mostly doesn't matter (they just set config), **except** styling you apply by hand can be overridden by the default styling functions that run at render time (see §2, §8).

### Object model (all inherit from flextable/officer)

| Object | Inherits | Built by | File |
|---|---|---|---|
| `clintable` | `flextable` | `clintable(df, …)`, `as_clintable(ft)` | `R/clintable.R` |
| `clinpage` | `flextable` | `slice_clintable()` (internal; one rendered page = row/col subset) | `R/clinpage.R`, `R/slice_clintable.R` |
| `clindoc` | `officer::rdocx` | `clindoc(…)`, `as_clindoc(x)` | `R/clindoc.R` |

### `clinify_config` fields (the deferred state)

Written by the verbs, read by the renderers:

`pagination_method` (`"default"` / `"custom"`), `page_by`, `max_rows`, `group_by`, `caption_by`, `group_when`, `key_cols`, `col_groups`, `auto_page_var`, `auto_page_when`, `auto_page_drop`, `titles`, `footnotes`, `footnote_page`, and the computed `pagination_idx`.

Any of `clin_page_by()`, `clin_group_by()`, `clin_alt_pages()` flips `pagination_method` to `"custom"`. Absent those, it stays `"default"`.

### Exported API (from `NAMESPACE`)

Verbs: `clintable`, `as_clintable`, `clindoc`, `as_clindoc`, `clin_page_by`, `clin_group_by`, `clin_alt_pages`, `clin_auto_page`, `clin_col_widths`, `clin_column_headers`, `clin_group_pad`, `clin_add_titles`, `clin_add_footnotes`, `clin_add_footnote_page`, `new_title_footnote`, `clin_replace_pagenums`, `make_grouped_pagenums`, `write_clindoc`.
Default-style functions (overridable options): `clinify_docx_default`, `clinify_titles_default`, `clinify_footnotes_default`, `clinify_table_default`, `clinify_caption_default`, `clinify_grouplabel_default`, plus helper `clin_default_table_width`.
S3: `print.clintable`, `knit_print.clintable`.

---

## 1. Core objects: `clintable` / `as_clintable`

### `clintable(x, page_by = NULL, group_by = NULL, use_labels = TRUE, ...)`

```r
clintable <- function(x, page_by = NULL, group_by = NULL, use_labels = TRUE, ...)
```

- `x` — **a data.frame** (any data.frame). It passes `x` and `...` straight to `flextable::flextable(x, ...)`, then converts via `as_clintable()`.
- **tplyr2 / Tplyr input:** there is no special adapter. You feed the *final presentation data.frame* — e.g. the result of `Tplyr::build()` after sorting, row-masking (`apply_row_masks()`) and `select()`. clinify treats every column you hand it as a display column. (The `end_to_end` vignette does exactly this with Tplyr.)
- `use_labels = TRUE` (default) → calls `headers_from_labels_()`, turning each column's `label` attribute into a header, splitting nested header levels on the literal delimiter `"||"` (see §4).
- `page_by`, `group_by` here are convenience pass-throughs to `as_clintable()` (equivalent to calling `clin_page_by()` / `clin_group_by()` later).

### `as_clintable(x, page_by = NULL, group_by = NULL)`

```r
as_clintable <- function(x, page_by = NULL, group_by = NULL)
```
- `x` **must be a flextable** (`stopifnot(inherits(x, "flextable"))`). Use this to wrap a flextable you built yourself (or from rtables/gtsummary → flextable).
- Sets `pagination_method` to `"custom"` if either `page_by`/`group_by` is supplied, else `"default"`.
- **Sets body/header padding at construction time** (this is *not* in the render-time table default): body `padding.top = padding.bottom = 0.1`; header first row `padding.top = 9`; header last row `padding.bottom = 9`.
- Sets `class(x) <- c("clintable", "flextable")`.

### Column → display-column mapping
Columns map 1:1 to `x$col_keys` (the flextable column keys = the data.frame column names). Widths (§3), headers (§4), and alignment are all addressed by these column-name keys. Some verbs *remove* columns from display after using them for logic (the `page_by`/`group_by`/`caption_by` vars are sliced out at render time in `prep_pagination_()`; `clin_auto_page(drop=TRUE)` and `clin_group_pad(drop=TRUE)` drop their driving variable).

---

## 2. `clinify_defaults` — project-level configuration (CRITICAL)

This is the mechanism for a **reusable, project- or organization-level clinify config**. There is no single `clinify_defaults()` function; instead clinify uses **six R `options()`**, populated in `.onLoad` (`R/zzz.R`), each read at render time via `getOption(...)`. You override a project's look by assigning your own values to these options (typically in a project startup file — `.Rprofile`, a `config.R` you `source()`, or a package `.onLoad`).

### The six options

| Option | Type | Controls | Read by |
|---|---|---|---|
| `clinify_docx_default` | `officer::prop_section` object | **Page size, orientation, margins** (and thus usable table width). Landscape by default. | `write_clindoc()`, `clin_default_table_width()` |
| `clinify_titles_default` | `function(x, ...)` → flextable | Styling of the **title block** (header of the docx) | title rendering |
| `clinify_footnotes_default` | `function(x, ...)` → flextable | Styling of the **footnote block** (footer of the docx) and footnote pages | footnote rendering |
| `clinify_table_default` | `function(x, ...)` → flextable | Styling of the **table body + header** (borders, font, size, layout) | `add_clintable_()`, print |
| `clinify_caption_default` | `function(x, ...)` → flextable | Styling of in-body **captions** from `caption_by` | `get_table_()`, print |
| `clinify_grouplabel_default` | `function(x, ...)` → flextable | Styling of **group labels** placed above the header (from `group_by`) | `get_table_()`, print |

`.onLoad` sets all six to the package's built-in functions, but preserves any user values already present (it saves `op <- options()` and re-applies them after setting defaults) — so **if you set these options in `.Rprofile` before the package loads, or after, your versions win.** The safe, explicit pattern is to set them after `library(clinify)`.

### Enumerating every controllable setting

There is **no flat list of scalar options** (no `clinify.font_size = 9` style knobs). Every visual property is expressed as flextable/officer calls *inside* these functions. The full inventory of what the shipped defaults control (verified against `R/styles.R`):

**Page geometry — `clinify_docx_default()`** (returns a `prop_section`):
- `page_size = officer::page_size(orient = "landscape")` → **landscape**. Default page is officer's A4: **11.693 in wide × 8.268 in tall** in landscape (portrait A4 is 8.264 × 11.694).
- `type = "continuous"`.
- `page_margins` carried from the officer default docx template via `officer::docx_dim(officer::read_docx())$margins` with `gutter = 0` appended: **top/bottom/left/right ≈ 0.984 in, header/footer ≈ 0.492 in, gutter = 0.** (Margins must be present so `clin_col_widths()` can compute usable width.)
- To go **portrait, letter, custom margins**, replace this function (see snippet below). Orientation/size/margins are the only geometry knobs and they all live here.

**Titles — `clinify_titles_default(x, ...)`:**
- All borders removed; font color black; **font `"Courier New"`, size `9`**; not bold, not italic.
- Width set to `clin_default_table_width() / 2` (title flextable has two columns: left/right).
- `line_spacing = 1` (all parts); cell padding top/bottom `= 0`; `set_table_properties(layout = "fixed")`.
- Calls `clin_replace_pagenums(x)` — swaps `{PAGE}`/`{NUMPAGES}` for Word fields.

**Footnotes — `clinify_footnotes_default(x, ...)`:**
- All borders removed, then a **single top border (black, width 1)** over the top body row (the rule above footnotes).
- Same font/size/spacing/padding/layout as titles; width `= clin_default_table_width() / 2`; calls `clin_replace_pagenums(x)`.
- Also used to style dedicated **footnote pages**.

**Table body + header — `clinify_table_default(x, ...)`:**
- Borders: removed everywhere, then horizontal lines drawn **only in the header** — an `hline` between header levels (default `fp_border()`), `hline_top` at the very top, and `hline_bottom` under the header. Blank interior header cells have their bottom borders explicitly set to `style="none", width=0` so spanner gaps don't show a stray rule (the bottom-most header row is exempted). **The body has no borders.**
- Font `"Courier New"`, size `9`, both parts; `set_table_properties(layout = "fixed")`.
- (Padding is *not* re-set here — it was set once in `as_clintable()`.)

**Captions — `clinify_caption_default(x, ...)`:** font `"Courier New"`, size `9` on the `"footer"` part (must run *after* the footer is added, hence a separate option).

**Group labels — `clinify_grouplabel_default(x, ...)`:** removes the top line above the group label, shifts the top rule down to row 1, and re-asserts `"Courier New"`/size 9 (font can be lost when a header is deleted, e.g. figures). This is what moves the group label *above* the header line (NEWS 0.3.0).

### How to set them once for a project and have them propagate

Because they are session options, set them once at session start and every `print()`/`write_clindoc()` for the rest of the session picks them up. Put this in a project `defaults.R` (there's a starter at `inst/defaults_template.R`, retrievable with `system.file("defaults_template.R", package = "clinify")`):

```r
# ---- project_clinify_defaults.R : source() this after library(clinify) ----
library(clinify)
library(flextable)
library(officer)

# 1) Page geometry: PORTRAIT US-LETTER with 1" margins (override landscape A4)
my_docx_default <- function() {
  officer::prop_section(
    page_size    = officer::page_size(width = 8.5, height = 11, orient = "portrait"),
    type         = "continuous",
    page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1,
                                     header = 0.5, footer = 0.5, gutter = 0)
  )
}

# 2) House table style: Times New Roman 10pt, thin header rules
my_table_default <- function(x, ...) {
  x <- flextable::border_remove(x)
  x <- flextable::hline_top(x, part = "header",
                            border = officer::fp_border(color = "black", width = 1))
  x <- flextable::hline_bottom(x, part = "header",
                               border = officer::fp_border(color = "black", width = 1))
  x <- flextable::font(x, part = "all", fontname = "Times New Roman")
  x <- flextable::fontsize(x, part = "all", size = 10)
  x <- flextable::set_table_properties(x, layout = "fixed")
  x
}

# 3) Titles/footnotes MUST keep the page-number substitution + fixed layout,
#    and width should track the (now portrait) usable width.
my_titles_default <- function(x, ...) {
  x <- flextable::border_remove(x)
  x <- flextable::font(x, fontname = "Times New Roman")
  x <- flextable::fontsize(x, size = 10)
  x <- flextable::width(x, width = clin_default_table_width() / 2)
  x <- flextable::line_spacing(x, space = 1, part = "all")
  x <- flextable::set_table_properties(x, layout = "fixed")
  clin_replace_pagenums(x)          # <-- do not drop this if you want {PAGE}
}

options(
  clinify_docx_default       = my_docx_default(),   # NOTE: a prop_section OBJECT, not a function
  clinify_titles_default     = my_titles_default,   # NOTE: the FUNCTION itself
  clinify_footnotes_default  = clinify_footnotes_default,  # keep shipped
  clinify_table_default      = my_table_default,
  clinify_caption_default    = clinify_caption_default,
  clinify_grouplabel_default = clinify_grouplabel_default
)
```

**Critical asymmetry to remember:** `clinify_docx_default` is stored as the **evaluated `prop_section` object** (`= my_docx_default()`), while the five styling entries are stored as **the functions themselves** (`= my_titles_default`, no parentheses). This mirrors `.onLoad` and the vignette.

If you override the page geometry, keep margins populated (the landscape default deliberately re-attaches them) — `clin_default_table_width()` errors if `clinify_docx_default` is unset and returns wrong widths if margins are missing.

To scope changes temporarily, use `withr::with_options(list(clinify_table_default = ...), { ... })`.

---

## 3. `clin_col_widths()` — percent-based / auto column widths

```r
clin_col_widths <- function(x, ...)
```
- `...` — **named numeric args**, each a **proportion of the usable page width** (0–1), keyed by column name: `clin_col_widths(mpg = .2, cyl = .2, disp = .15)`. Non-numeric → error. Any value `> 1` → error (with a hint to call `clin_alt_pages()` first if you meant per-page proportions). Unknown column names → error.
- Usable width comes from `clin_default_table_width()` = `page_width − (left + right margins)` (default **9.725 in** in the shipped landscape config).
- **Leftover columns** (those you didn't name) split the remaining width **evenly**: `lo_wd = ((1 - sum(named)) / n_leftovers) * usable_width`. If you name columns summing to > 1 for a standard table → error.
- **There is no true "auto-fit" call** — flextable's own `flextable::autofit()`/`dim_pretty()` still work on the object, but the clinify idiom is: name the important columns, let the rest divide the remainder. To let *every* column auto-size, just don't call `clin_col_widths()` (but note `layout="fixed"` is set by the table default).
- **Alternating pages (`clin_alt_pages`) interaction:** when `col_groups` is set, key-column proportions are applied to (and repeated on) *every* page; each `col_group`'s named widths are applied, and any unnamed columns in a group absorb the remaining page width so each page fills the full width. Guard: `sum(key_cols + one col_group)` must be ≤ 1 per page or it errors. Widths must be given as page proportions and `clin_alt_pages()` must be called first.

```r
clintable(mtcars) |>
  clin_col_widths(mpg = .2, cyl = .2, disp = .15, vs = .15)   # rest split evenly
```

---

## 4. `clin_column_headers()` — multi-level headers / spanners / `(N=x)`

```r
clin_column_headers <- function(x, ...)
```
- `...` — **named character vectors**, keyed by column name. Each **element of the vector is one header row**, applied **top → bottom** (first element = top level, last = bottom). All args must be character; all names must be existing columns.
- **Uneven depth:** the deepest arg sets the number of header rows; shorter vectors **bind to the bottom** (a 2-element header for a column in a 3-row layout occupies rows 2–3). Implementation fills a `depth × ncol` matrix bottom-up.
- **Spanners:** cells sharing identical text are **merged horizontally** (via `merge_h(part="header")`) — and by repeating the same text vertically across adjacent columns you get grouped spanners. (The shipped `clinify_table_default` additionally suppresses borders under blank spanner-gap cells.)
- **`(N=x)` treatment-arm headers and line breaks:** put `\n` inside a string to force a new visual line within a header cell; interpolate the N with `sprintf`. Nested arm example (from `end_to_end.Rmd`, Ns pulled from `Tplyr`'s `t$header_n`):

```r
header_n <- t$header_n$n
names(header_n) <- t$header_n$TRT01P

clintable(dat) |>
  clin_column_headers(
    row_label1 = "",
    row_label2 = "",
    var1_Placebo = sprintf("Placebo\n(N=%s)", header_n["Placebo"]),
    `var1_Xanomeline Low Dose`  = c("Xanomeline", sprintf("Low Dose\n(N=%s)",  header_n["Xanomeline Low Dose"])),
    `var1_Xanomeline High Dose` = c("Xanomeline", sprintf("High Dose\n(N=%s)", header_n["Xanomeline High Dose"]))
  )
# "Xanomeline" repeated across the two dose columns => a merged spanner over both,
# with "Low Dose (N=..)" / "High Dose (N=..)" beneath.
```

- **Labels alternative (`use_labels = TRUE`):** instead of calling `clin_column_headers()`, set each column's `label` attribute and split levels with the literal `"||"` delimiter; `clintable()` calls `headers_from_labels_()` → `clin_column_headers()` for you. Repeated tokens across columns merge into spanners:

```r
iris2 <- iris
attr(iris2$Sepal.Length, "label") <- "Flower||Sepal||Length"
attr(iris2$Sepal.Width,  "label") <- "Flower||Sepal||Width"
attr(iris2$Petal.Length, "label") <- "Flower||Petal||Length"
attr(iris2$Petal.Width,  "label") <- "Flower||Petal||Width"
clintable(iris2)   # "Flower" spans all four; "Sepal"/"Petal" span their pairs
```

Header alignment is not set here — apply `flextable::align(part = "header")` yourself.

---

## 5. Grouping, indentation & pagination verbs

All defaults verified from source. Note the shared `when` argument (`match.arg`) — **its default is `"change"`** for `clin_group_by`, `clin_auto_page`, and `clin_group_pad` (the first element of `c("change","notempty")`).

### `clin_page_by(x, page_by, max_rows = 10)`
Data-driven page breaks. Break inserted wherever `page_by` **changes** value (uses `find_split_inds(..., "change")`; not a sort). If `page_by` is omitted, splits every `max_rows` rows (default 10) instead. Sets `pagination_method = "custom"`. The page variable is dropped from the displayed table at render time.

```r
clintable(dat) |> clin_page_by("page")
clintable(mtcars) |> clin_page_by(max_rows = 10)   # or clin_page_by(max_rows = 5)
```

### `clin_group_by(x, group_by, caption_by = NULL, when = c("change","notempty"))`
- `group_by` — character vector of one or more variables rendered as **label lines ABOVE the column header** (and pulled out of the body). Multiple vars → multiple stacked label lines.
- `caption_by` — single variable rendered as an **in-body caption** (a footer line below the table body, above the page footer).
- `when` — `"change"` (new group when value changes) or `"notempty"` (new group when value is populated). Values carried forward with `zoo::na.locf`.
- When combined with `clin_page_by()`, a page break occurs when **any** of page/group/caption vars change. Sets `pagination_method = "custom"`.

```r
clintable(dat2) |>
  clin_page_by("page") |>
  clin_group_by(c("groups1","groups2"), caption_by = "captions", when = "change")
```

### `clin_auto_page(x, group_var, when = c("change","notempty"), drop = FALSE)`
A **separate strategy** from pre-slicing: it does not compute page breaks. It flags rows with flextable's `keep_with_next()` so **Word itself** avoids splitting a group across a page boundary (the row before each new group is marked *not* keep-with-next to allow a break there). Effect is invisible in HTML preview — only in the `.docx`. `drop = TRUE` removes `group_var` from the table. Least-effort pagination; recommended when it suffices.

```r
clintable(dat) |> clin_auto_page("ord_layer_index", when = "change", drop = TRUE)
```

### `clin_group_pad(x, pad_by, size = 9, when = c("change","notempty"), drop = FALSE)`
Vertical **indentation/whitespace between groups**: adds `padding.top = size + 5` (default `9 + 5 = 14`) to the first row of each group (the very first table row is never padded: `splits[1] <- FALSE`). `drop = TRUE` removes `pad_by`. This is how group separation is achieved; horizontal indentation of row labels is expected to come pre-baked in the incoming data (leading spaces), not from a clinify function.

```r
clintable(dat) |> clin_group_pad("row_label1", when = "notempty")           # size 14
clintable(mtcars) |> clin_group_pad("gear", size = 15)
```

### `clin_alt_pages(x, key_cols, col_groups)`
Splits a **too-wide** table across pages horizontally. `key_cols` (character vector) repeat on every page; `col_groups` (list of character vectors) — one group appended per page. With `col_groups` of length 3 you get 3 column-page variants, each = `key_cols + one group`. **Should be paired with `clin_page_by()`** so row chunking is defined; if no row pagination is set, clinify emits a `NOTE` and defaults to **20 rows per page**. Sets `pagination_method = "custom"`. Output order: for each row-chunk, emit all col-group variants before advancing to the next row-chunk.

```r
clintable(dat2) |>
  clin_page_by("page") |>
  clin_alt_pages(
    key_cols  = c("mpg","cyl","hp"),
    col_groups = list(c("disp","drat","wt"), c("qsec","vs","am"), c("gear","carb"))
  )
```

### `make_grouped_pagenums(var, rows)` (helper)
Builds a page-index vector for alt-page tables: assigns sequential page numbers, starting a new page every `rows` rows **or** whenever `var` changes (input must be pre-sorted so identical values are contiguous). Useful to synthesize a `page_by` variable when using `clin_group_by()` with alternating pages.

```r
library(dplyr)
head(ToothGrowth, 20) |> mutate(page = make_grouped_pagenums(dose, 7))
```

---

## 6. Titles / footnotes model + page numbers

### The model
- `clin_add_titles(x, ls = NULL, ft = NULL)`, `clin_add_footnotes(...)`, `clin_add_footnote_page(...)` — all thin wrappers over `add_titles_footnotes_()`. **Exactly one** of `ls` or `ft` must be supplied (both or neither → error). They accept a `clintable` **or** a `clindoc`.
- `ls` — a **list of character vectors**, each vector = one line; **max 2 elements per vector** (`new_title_footnote()` errors on length > 2 — see gotchas).
- Alignment rules (per line, by element count):
  - **1 element:** center in a **title**, left in a **footnote**.
  - **2 elements:** split — left element left-aligned, right element right-aligned.
  - **Single left-aligned title:** pass a 2-element vector with the *same text twice* (`c("Left","Left")`) — merges and keeps left alignment.
- `ft` — supply your **own flextable** as the title/footnote block for full control (`clin_add_titles(ft = flextable::flextable(...))`).
- **How many lines:** unlimited lines (list length); each line ≤ 2 columns.

### `new_title_footnote(x, sect = c("titles","footnotes","footnote_page"))`
Builds the underlying two-column flextable (columns `Left`/`Right`, header/footer parts deleted) from the list, applying the alignment/merge rules above. Exposed so you can build, then further style, a block before attaching it.

### Sourcing titles from a spreadsheet (titles.xlsx)?
**No built-in support.** A repo-wide search found no `readxl`/`openxlsx`/`xlsx`/`read_excel`/`.csv` usage anywhere in `R/`, vignettes, or docs. Titles/footnotes come only from (a) an R list of character vectors, or (b) a pre-built flextable via `ft`. To drive from a `titles.xlsx`, you must read it yourself and assemble the `ls` list, e.g.:

```r
# DIY spreadsheet-driven titles (clinify has NO native reader):
tt <- readxl::read_excel("titles.xlsx")                 # your own columns/schema
ls <- split(tt, tt$line) |> lapply(\(d) as.character(d$text))
clintable(dat) |> clin_add_titles(ls)
```

### Footnote page — `clin_add_footnote_page()`
Same interface. Produces a dedicated **preceding page** of footnotes (styled with `clinify_footnotes_default`) inserted before the table, then a page break. Use when there are too many footnotes for the footer.

### Page numbers — `clin_replace_pagenums(x)`
Finds cells containing `{PAGE}` or `{NUMPAGES}` and replaces them with real Word field objects (`officer` `as_word_field`) for current / total pages. It is **called inside `clinify_titles_default` and `clinify_footnotes_default`**, so simply putting `"Page {PAGE} of {NUMPAGES}"` in a title/footnote line "just works" in Word. In HTML preview those placeholders render blank. If you write custom title/footnote default functions, **keep the `clin_replace_pagenums()` call** or you lose page numbering.

```r
clin_add_titles(list(
  c("Protocol: CDISCPILOT01", "Page {PAGE} of {NUMPAGES}"),
  c("Table 14-2.01"),
  c("Summary of Demographic and Baseline Characteristics")
))
```

`make_grouped_pagenums()` (see §5) is the complementary helper for *data-driven* page indexing (distinct from the Word field substitution).

---

## 7. `clindoc()` / `as_clindoc()` and `write_clindoc()`

### Building the document
```r
clindoc(...)        # one clintable, several clintables, OR a single list of clintables
as_clindoc(x)       # a single clintable -> clindoc
```
- `clindoc(ct)` (single) → calls `as_clindoc()`.
- `clindoc(ct1, ct2, ...)` **or** `clindoc(list(ct1, ct2))` → multiple tables in one document, **separated by page breaks** (a page break after every table except the last). All must be clintables.
- **Multiple-table titles/footnotes:** apply them to the **`clindoc`** (not the individual tables) — the `clin_add_*()` functions work on a `clindoc` identically. Title/footnote config on individual clintables in a multi-table doc is **ignored**. Rationale (per `write` vignette): Word applies header/footer once to the section and repeats them on every page, so specifying per-table would be redundant/clunky.
- `as_clindoc()` extracts a single table's `titles`/`footnotes` into the doc config, prepends any `footnote_page`, then adds the table.
- A `clindoc` inherits `officer::rdocx`, so any `officer` body/section function works on it.

```r
ct <- clintable(head(mtcars,10))
doc <- clindoc(ct)                       # or clindoc(ct1, ct2) / clindoc(list(...))

tables <- lapply(list(mtcars, iris), \(x) clintable(head(x, 10)))
doc <- clindoc(tables) |>
  clin_add_titles(list(c("Left","Page {PAGE} of {NUMPAGES}"), c("Just the middle"))) |>
  clin_add_footnotes(list(c("Source: prog.R", format(Sys.time(), "%H:%M %A, %B %d, %Y"))))
```

### Writing — `write_clindoc(x, file)`
```r
write_clindoc <- function(x, file)
```
- `x` may be a `clintable` (auto-converted via `as_clindoc()`) **or** a `clindoc`.
- Pulls `settings <- getOption("clinify_docx_default")` (the `prop_section`).
- If titles present → styled with `clinify_titles_default` → set as `settings$header_default` (a `block_list`) ⇒ **repeats on every page**.
- If footnotes present → styled with `clinify_footnotes_default` → `settings$footer_default` ⇒ **repeats on every page**.
- A footnote page (if any) is inserted at the document start, followed by a page break.
- `settings` applied via `officer::body_set_default_section()`; then `print(doc, target = file)` writes the `.docx`.
- Page setup (size/orientation/margins), page headers/footers, and page numbering are therefore **all governed by `clinify_docx_default` + the title/footnote default functions** — set them at project level (§2).

```r
write_clindoc(ct,  file.path(tempdir(), "demo.docx"))   # clintable
write_clindoc(doc, file.path(tempdir(), "demo.docx"))   # clindoc
```

**Custom pagination path:** for `pagination_method == "custom"`, `add_clintable_()` runs `prep_pagination_()` (builds `pagination_idx`), then `doc_alternating_()` emits each page's slice (`slice_clintable()`) with `run_pagebreak()` between pages, using officer's `body_append_*_context` API. `default` method just does `flextable::body_add_flextable()`. `clin_auto_page()` runs `auto_page_()` (adds `keep_with_next` flags) before adding.

### HTML preview — `print.clintable(x, n = 3, nrows = 15, apply_defaults = TRUE, ...)`
- `n` — pages to render (default 3, custom pagination only); `nrows` — rows per page when rows aren't set by pagination (default 15); `apply_defaults = TRUE` applies the default styling functions for the preview. Renders styled HTML with a small JS page-switcher; titles/footnotes shown as separate flextables above/below the body. `{PAGE}`/`{NUMPAGES}` appear blank in HTML. `knit_print.clintable()` is the knitr equivalent.

---

## 8. Default clinical-table look (regulatory)

Verified from `R/styles.R` + `as_clintable()`:

| Property | Default |
|---|---|
| **Font** | `Courier New` (monospace) — header, body, titles, footnotes, captions, group labels |
| **Font size** | `9 pt` everywhere |
| **Orientation / page** | **Landscape**, A4 (11.693 × 8.268 in) |
| **Margins** | ≈ 0.984 in top/bottom/left/right; ≈ 0.492 in header/footer; gutter 0 |
| **Usable table width** | 9.725 in (`clin_default_table_width()`) |
| **Table layout** | `fixed` (widths honored, no auto-stretch) |
| **Table borders** | Header only: top rule, rule(s) between header levels, bottom rule; blank spanner-gap cells have their bottom border suppressed. **Body: no borders.** |
| **Footnote border** | Single black top rule (width 1) above the footnote block |
| **Title borders** | None |
| **Line spacing** | Titles/footnotes: `1`. Body: flextable default. |
| **Body padding** | top/bottom `0.1` (set in `as_clintable`); header first row `padding.top = 9`, last header row `padding.bottom = 9` |
| **Group-label** | Placed above the header line, `Courier New`/9, top rule shifted down |

**Decimal alignment:** clinify does **not** auto-align on the decimal point. Numeric alignment is expected to be handled either by (a) the incoming data being pre-formatted as fixed-width strings (typical Tplyr output — `Courier New` monospacing then visually aligns them), or (b) manual `flextable::align()` calls. There is no `align_decimal`-style helper.

**Monospacing / indentation:** monospacing is achieved purely via the `Courier New` font. Horizontal indentation of nested row labels is expected to be baked into the data as leading spaces — `print.clintable` even post-processes the HTML (`gsub` to `&nbsp;`) so leading/among spaces survive in the preview; in Word the fixed-width font preserves them. Vertical separation between groups is done with `clin_group_pad()`.

**Note:** `inst/defaults_template.R` is a *starter* template and is **not identical** to the shipped functions — its `clinify_table_default` calls `merge_v`/`merge_h` and sets padding and uses `width=0.2` header borders, whereas the shipped `R/styles.R` version does merging elsewhere (in `clin_column_headers`) and sets padding in `as_clintable()`. Treat `inst/defaults_template.R` as a copy-paste starting point, not a spec of current defaults.

---

## 9. Full runnable end-to-end example (data.frame → styled clintable → docx)

Self-contained (uses `mtcars`; no external packages beyond clinify). Mirrors the vignette patterns.

```r
library(clinify)
library(flextable)   # clin_* objects accept native flextable verbs
library(officer)

## ---- (optional) project-level defaults: set ONCE, propagate everywhere ----
## Landscape A4 / Courier New 9 are already the shipped defaults; override only if needed.
## options(clinify_docx_default = my_docx_default(),
##         clinify_table_default = my_table_default, ...)   # see section 2

## 1) A presentation-ready data.frame (in practice: Tplyr::build() |> arrange() |>
##    apply_row_masks() |> select(...); here we mock it)
dat <- mtcars
dat$page   <- c(rep(1,10), rep(2,10), rep(3,10), 4,4)   # drives page breaks
dat$grp    <- ifelse(dat$page <= 2, "Group A", "Group B")# drives group label + auto-page
dat$model  <- rownames(mtcars)
rownames(dat) <- NULL
dat <- dat[, c("grp","model","mpg","cyl","hp","disp","drat","wt","page")]

## 2) Build + style the clintable (every clin_* just records config; flextable verbs work inline)
ct <- clintable(dat) |>
  # data-driven page breaks
  clin_page_by("page") |>
  # label line ABOVE the header, from data
  clin_group_by("grp") |>
  # whitespace between model groups (vertical separation)
  clin_group_pad("cyl", when = "change", size = 12) |>
  # multi-level / spanner headers, with \n line breaks and (N=)
  clin_column_headers(
    model = c("", "Car"),
    mpg   = c("Efficiency", "MPG"),
    cyl   = c("Engine", "Cyl"),
    hp    = c("Engine", "HP"),
    disp  = c("Engine", "Disp"),
    drat  = c("Drivetrain", "Rear\naxle"),
    wt    = c("Drivetrain", "Weight\n(1000 lb)")
  ) |>
  # native flextable alignment on the clintable
  flextable::align(align = "center", part = "header") |>
  flextable::align(j = c("mpg","cyl","hp","disp","drat","wt"), align = "right") |>
  # column widths as proportions of usable width (rest split evenly)
  clin_col_widths(model = .28, mpg = .10, cyl = .10, hp = .10) |>
  # titles: 2-element line splits L/R; {PAGE}/{NUMPAGES} become Word fields
  clin_add_titles(list(
    c("Protocol: CDISCPILOT01", "Page {PAGE} of {NUMPAGES}"),
    c("Table 14-1.01"),
    c("Demonstration Table")
  )) |>
  # footnotes: 2-element line splits L/R; single top rule drawn above
  clin_add_footnotes(list(
    c("Source: /path/to/prog.R", format(Sys.time(), "%H:%M %A, %B %d, %Y"))
  ))

## 3) Preview in the IDE / knit (styled HTML, page-switcher, 3 pages by default)
print(ct)                    # or print(ct, n = 5, nrows = 20)

## 4a) Write straight to .docx (clintable -> clindoc happens automatically)
write_clindoc(ct, file.path(tempdir(), "demo_table.docx"))

## 4b) Or build the officer::rdocx yourself first (e.g. to add more officer content,
##     or to combine several tables in one document)
doc <- clindoc(ct)
write_clindoc(doc, file.path(tempdir(), "demo_table.docx"))

## Multi-table document (titles/footnotes go on the DOC, not the tables):
tables <- lapply(list(head(mtcars,10), head(iris,10)), clintable)
doc2 <- clindoc(tables) |>
  clin_add_titles(list(c("Study X","Page {PAGE} of {NUMPAGES}"), c("Combined Output")))
write_clindoc(doc2, file.path(tempdir(), "combined.docx"))
```

A wide-table (alternating pages) variant:

```r
clintable(dat) |>
  clin_page_by("page") |>                                   # required for row chunking
  clin_alt_pages(
    key_cols  = c("model","mpg"),                            # repeated on every page
    col_groups = list(c("cyl","hp"), c("disp","drat"), c("wt")) # one group per page
  ) |>
  clin_col_widths(model = .3, mpg = .2, cyl = .25) |>        # key proportions carry per page
  write_clindoc(file.path(tempdir(), "wide.docx"))
```

---

## 10. Limitations / gotchas

**Breaking changes (NEWS.md):**
- **0.3.0:** titles/footnotes now split into **at most 2 parts** (was 3) to avoid line wrapping — *breaking*. `new_title_footnote()` errors on any line with > 2 elements.
- **0.3.0:** `write_clintable()` **renamed to `write_clindoc()`** — *breaking*.
- **0.3.0:** `clin_auto_page()` gained `drop`; **does not drop the page var by default** now — *breaking*.
- **0.3.1:** minimum `officer (>= 0.7.2)`; group-label styling retains configured font size; fixes to single-column slicing, spanner-gap handling across page splits, and column-vector simplification.

**Documentation vs. code discrepancies (trust the code):**
- The **Getting Started vignette (`clinify.Rmd`)** still says *"If three elements are provided, they align left, right, and center."* This is **stale** — the code (`new_title_footnote`) hard-errors on > 2 elements ("All sublists must have length <= 2"). Use at most 2 per line.
- `inst/defaults_template.R` differs from the shipped `R/styles.R` defaults (does `merge_v`/`merge_h` in the table default, sets padding, uses 0.2pt header borders). It's a starter, not the current spec.
- The `defaults` vignette narrates the styling functions but omits `clin_replace_pagenums()` in its `clinify_table_default` listing — the real title/footnote defaults are where page-number substitution lives; don't drop it when overriding.

**Version-fragility / internals:**
- `slice_clintable()` (used for all custom/alt pagination and column-dropping) reaches directly into flextable/officer internal structures (`$body$dataset`, `complex_tabpart`, `fpstruct`, `chunkset_struct`, `$spans$rows`, `$styles$cells/pars/text`). Git history shows repeated breakage on flextable/officer upgrades — hence the `officer (>= 0.7.2)` pin. **Treat any flextable/officer bump as high-risk** for slicing and the (numerous) snapshot tests.
- Font-metric-based slicing estimates are platform-dependent (tests were hardened for r-devel-fedora-gcc).

**Behavioral gotchas:**
- `clin_default_table_width()` **errors** if `clinify_docx_default` is unset; if you override page geometry, keep margins populated or column-width math is wrong.
- `clin_alt_pages()` without `clin_page_by()` emits a `NOTE` and silently defaults to **20 rows/page**.
- `page_by`/`group_by`/`caption_by` variables are **removed from the displayed table** at render time; `clin_auto_page(drop=)` and `clin_group_pad(drop=)` default to **keeping** their driving variable (`drop = FALSE`).
- `clin_auto_page()` and `{PAGE}` fields have **no visible effect in the HTML preview** — only in Word.
- Multi-table `clindoc()`: title/footnote config on the individual clintables is **ignored**; put them on the `clindoc`.
- The shared `when` argument defaults to **`"change"`** (first `match.arg` value) for `clin_group_by`, `clin_auto_page`, `clin_group_pad`.
- **No spreadsheet/xlsx/csv ingestion** anywhere — titles/footnotes are R lists or flextables only.
- No decimal-alignment or auto-fit helper — rely on monospaced pre-formatted strings and manual `flextable::align()`; `layout="fixed"` means widths are respected literally.
- Installed build in this environment reports `0.3.0`; source tree is `0.3.1` (numbers above verified against the installed package where run).

**Key file map:** objects `R/clintable.R`, `R/clinpage.R`, `R/clindoc.R`; verbs `R/pagination.R`, `R/col_width.R`, `R/column_headers.R`, `R/add_titles_footnotes.R`, `R/group_pad.R`; defaults/options `R/styles.R`, `R/zzz.R`, `inst/defaults_template.R`; render `R/print.R`, `R/write.R`, `R/slice_clintable.R`; page numbers `R/pagenums.R`; helpers `R/helpers.R`.