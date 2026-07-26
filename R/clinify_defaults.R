# R/clinify_defaults.R — project clinify house style (fonts, rules, page geometry) for the pilot outputs

library(flextable)
library(officer)

#' Build the page-geometry section for the pilot outputs
#'
#' US-Letter landscape, 1in margins, 0.5in header/footer band.
#'
#' @return An evaluated officer `prop_section` (stored as an object, not a function).
cdisc_docx_default <- function() {
  officer::prop_section(
    page_size    = officer::page_size(width = 8.5, height = 11, orient = "landscape"),
    type         = "continuous",
    page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1,
                                     header = 0.5, footer = 0.5, gutter = 0)
  )
}

#' Create a solid black border of the given width
#'
#' @param width Border width in points.
#' @return An officer `fp_border`.
.cdisc_rule <- function(width = 1) officer::fp_border(color = "black", width = width, style = "solid")

#' Style the table body and header (Courier New 10, header bottom rule)
#'
#' @param x A flextable.
#' @param ... Unused; present for the clinify styler interface.
#' @return The styled flextable.
cdisc_table_default <- function(x, ...) {
  x <- flextable::border_remove(x)
  # single 1pt rule under the bottom-most header row (matches the reference RTF border)
  x <- flextable::hline_bottom(x, part = "header", border = .cdisc_rule(1))
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  x <- flextable::fontsize(x, part = "all", size = 10)
  # header labels are deliberately not bolded (matches the rendered references)
  # tight horizontal padding keeps fixed-width strings from wrapping
  x <- flextable::padding(x, padding.left = 2, padding.right = 2, part = "all")
  # compact rows: single line spacing, no top/bottom padding in the body
  x <- flextable::line_spacing(x, space = 1, part = "all")
  x <- flextable::padding(x, padding.top = 0, padding.bottom = 0, part = "body")
  # ~18pt of buffer above the column labels, 4pt below
  x <- flextable::padding(x, padding.top = 18, padding.bottom = 4, part = "header")
  # row pitch (body/title/footnote) is set centrally in add_titles_footnotes()
  x <- flextable::set_table_properties(x, layout = "fixed")
  x
}

#' Style the title block (Courier New 10, bold + italic, page-number fields)
#'
#' @param x A flextable.
#' @param ... Unused; present for the clinify styler interface.
#' @return The styled flextable with `{PAGE}`/`{NUMPAGES}` replaced by Word fields.
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
  x <- flextable::set_table_properties(x, layout = "fixed")
  clinify::clin_replace_pagenums(x)   # {PAGE}/{NUMPAGES} -> Word fields
}

#' Style the footnote block (Courier New 10, italic, no top rule)
#'
#' @param x A flextable.
#' @param ... Unused; present for the clinify styler interface.
#' @return The styled flextable with `{PAGE}`/`{NUMPAGES}` replaced by Word fields.
cdisc_footnotes_default <- function(x, ...) {
  x <- flextable::border_remove(x)            # drops clinify's default top rule
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  x <- flextable::fontsize(x, part = "all", size = 10)
  x <- flextable::italic(x, part = "all")
  x <- flextable::width(x, width = clinify::clin_default_table_width() / 2)
  x <- flextable::line_spacing(x, space = 1, part = "all")
  x <- flextable::padding(x, padding.top = 0, padding.bottom = 0,
                          padding.left = 0, padding.right = 0, part = "all")
  x <- flextable::set_table_properties(x, layout = "fixed")
  clinify::clin_replace_pagenums(x)
}

#' Style the caption block (Courier New 10)
#'
#' @param x A flextable.
#' @param ... Unused; present for the clinify styler interface.
#' @return The styled flextable.
cdisc_caption_default <- function(x, ...) {
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  flextable::fontsize(x, part = "all", size = 10)
}

#' Style the group-label block (Courier New 10)
#'
#' @param x A flextable.
#' @param ... Unused; present for the clinify styler interface.
#' @return The styled flextable.
cdisc_grouplabel_default <- function(x, ...) {
  x <- flextable::font(x, part = "all", fontname = "Courier New")
  flextable::fontsize(x, part = "all", size = 10)
}

#' Register the CDISC house-style stylers as clinify session options
#'
#' @return `TRUE`, invisibly.
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
