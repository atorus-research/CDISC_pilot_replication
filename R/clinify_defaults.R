# R/clinify_defaults.R
# -----------------------------------------------------------------------------
# Project-level clinify house style for the CDISC Pilot next-gen outputs.
#
# clinify exposes its look through six session options (see ?clinify::clinify_defaults).
# Sourcing this file overrides them to match the pilot's visual spec, verified
# against the original RTFs:
#   * US-Letter LANDSCAPE, 1" margins all sides, 0.5" header/footer band
#   * Courier New, 10pt everywhere (clinify ships 9pt)
#   * Titles: bold + italic, no borders
#   * Footnotes: italic (not bold), NO rule above them
#   * Table: header carries a single 1pt bottom rule under the last header row;
#           no top rule, no inter-level rule; body has no borders
#   * layout = "fixed" so column widths are honoured literally (monospace align)
#
# Set once per session (see R/setup.R). Because clinify reads these options at
# render time, every write_clindoc()/print() thereafter inherits the style.
# -----------------------------------------------------------------------------

library(flextable)
library(officer)

# Reference body row pitch (~15.35pt), measured from the original RTFs rendered
# through LibreOffice. Used to force an exact row height (see cdisc_table_default).
# House-standard body row pitch for the whole suite: one deliberate, compact
# single-line pitch for 10pt Courier (~15.4pt / 0.213in). Applied consistently
# to every table (the legacy outputs used inconsistent per-table heights; we do
# not reproduce that). A `getOption("cdisc.body_row_in")` escape hatch exists but
# is intentionally unused.
CDISC_ROW_IN   <- 15.35 / 72
# Title/footnote single-line pitch (~11.4pt).
CDISC_TITLE_IN <- 11.4 / 72

# 1) Page geometry: prop_section OBJECT (note: stored evaluated, not as a fn) ---
cdisc_docx_default <- function() {
  officer::prop_section(
    page_size    = officer::page_size(width = 8.5, height = 11, orient = "landscape"),
    type         = "continuous",
    page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1,
                                     header = 0.5, footer = 0.5, gutter = 0)
  )
}

.cdisc_rule <- function(width = 1) officer::fp_border(color = "black", width = width, style = "solid")

# 2) Table body + header: Courier New 10, header bottom rule only -------------
cdisc_table_default <- function(x, ...) {
  x <- flextable::border_remove(x)
  # single 1pt rule under the bottom-most header row (matches \clbrdrb\brdrw20)
  x <- flextable::hline_bottom(x, part = "header", border = .cdisc_rule(1))
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  x <- flextable::fontsize(x, part = "all", size = 10)
  # NOTE: headers are NOT bolded. Although 22/30 legacy programs call set_bold, the
  # RENDERED references are mostly regular weight (a bold house default measurably WORSENED
  # the pixel-gated single-page tables 14-1.01/14-3.01/14-1.03, and its only clear
  # beneficiary, 14-6.01, is content-verified so header weight is invisible to its gate).
  # Non-bold is the consistent option-A choice; 14-6.01's visibly-bold reference is the
  # documented cosmetic exception (see notes/divergences.md).
  # Match the reference cells' 4pt (80 twip) horizontal padding. The reference
  # relied on \clNoWrap (which LibreOffice ignores), so tight padding is what
  # actually lets fixed-width strings like "21.0 (  5;61)" fit without wrapping.
  x <- flextable::padding(x, padding.left = 2, padding.right = 2, part = "all")
  # Tight rows to match the reference (top/bottom padding 0, single line spacing).
  # clinify's as_clintable() sets header padding.top/bottom = 9; trim it so the
  # header block matches the reference's compact buffer + label rows.
  x <- flextable::line_spacing(x, space = 1, part = "all")
  x <- flextable::padding(x, padding.top = 0, padding.bottom = 0, part = "body")
  # Header buffer: ~15pt of space above the column labels (reference's set_column_header_buffer),
  # applied to every header row via flextable padding. NOT via clinify's clin_header_pad(): #102
  # made above/below uniform per-row (good), but two tables need NON-uniform per-row header padding
  # for fidelity -- 14-1.03 (34pt on row 2, spanner->label gap) and 14-3.10 (21pt on row 1, to align
  # the header rule at the reference y) -- and clin_header_pad's config precedence clobbers their
  # per-table flextable::padding(i=). So we keep the styler padding, which composes with those. (#97)
  x <- flextable::padding(x, padding.top = 18, padding.bottom = 4, part = "header")
  # House row pitch (body/title/footnote) is set once, centrally, via clinify's
  # clin_row_height() in add_titles_footnotes() (clinify #97) - it replaces the former
  # per-styler hrule("atleast") + height_all() and reaches the title/footnote blocks too.
  x <- flextable::set_table_properties(x, layout = "fixed")
  x
}

# 3) Titles: Courier New 10, BOLD + ITALIC, no borders, keep page-number fields
cdisc_titles_default <- function(x, ...) {
  x <- flextable::border_remove(x)
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  x <- flextable::fontsize(x, part = "all", size = 10)
  x <- flextable::bold(x, part = "all")
  x <- flextable::italic(x, part = "all")
  x <- flextable::width(x, width = clinify::clin_default_table_width() / 2)
  x <- flextable::line_spacing(x, space = 1, part = "all")
  x <- flextable::padding(x, padding.top = 0, padding.bottom = 0,
                          padding.left = 0, padding.right = 0, part = "all")
  # Title pitch (~11.4pt) now comes from clin_row_height(title=) (see add_titles_footnotes).
  x <- flextable::set_table_properties(x, layout = "fixed")
  clinify::clin_replace_pagenums(x)   # {PAGE}/{NUMPAGES} -> Word fields
}

# 4) Footnotes: Courier New 10, ITALIC only, NO top rule ----------------------
cdisc_footnotes_default <- function(x, ...) {
  x <- flextable::border_remove(x)            # (drops clinify's default top rule)
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  x <- flextable::fontsize(x, part = "all", size = 10)
  x <- flextable::italic(x, part = "all")
  x <- flextable::width(x, width = clinify::clin_default_table_width() / 2)
  x <- flextable::line_spacing(x, space = 1, part = "all")
  x <- flextable::padding(x, padding.top = 0, padding.bottom = 0,
                          padding.left = 0, padding.right = 0, part = "all")
  # Footnote pitch (~11.4pt) now comes from clin_row_height(footnote=) (see add_titles_footnotes).
  x <- flextable::set_table_properties(x, layout = "fixed")
  clinify::clin_replace_pagenums(x)
}

# 5) Captions & 6) group labels: keep Courier New 10 --------------------------
cdisc_caption_default <- function(x, ...) {
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  flextable::fontsize(x, part = "all", size = 10)
}
cdisc_grouplabel_default <- function(x, ...) {
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  flextable::fontsize(x, part = "all", size = 10)
}

# Apply everything ------------------------------------------------------------
apply_cdisc_clinify_defaults <- function() {
  options(
    clinify_docx_default       = cdisc_docx_default(),   # evaluated object
    clinify_titles_default     = cdisc_titles_default,   # functions (no parens)
    clinify_footnotes_default  = cdisc_footnotes_default,
    clinify_table_default      = cdisc_table_default,
    clinify_caption_default    = cdisc_caption_default,
    clinify_grouplabel_default = cdisc_grouplabel_default
  )
  invisible(TRUE)
}

apply_cdisc_clinify_defaults()
