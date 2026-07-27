# R/helpers.R — shared formatters for row-bound model output and (N=) column labels

library(stringr)

#' Format a number as a fixed-width, right-padded string
#'
#' Thin wrapper over `tplyr2::apply_formats()`: builds an `f_str` of `int_len`
#' integer digits (plus `digits` decimals) and right-pads the result to `size`.
#' NA renders as "".
#'
#' @param var Numeric value(s) to format
#' @param digits Number of decimal places
#' @param size Total field width, right-padded with spaces
#' @param int_len Minimum width of the integer part
#' @return Character, "" when `var` is NA
num_fmt <- function(var, digits = 0, size = 10, int_len = 3) {
  fmt <- if (digits > 0) paste0(strrep("x", int_len), ".", strrep("x", digits)) else strrep("x", int_len)
  tplyr2::apply_formats(tplyr2::f_str(fmt, "x"), x = var, na = "", width = size, pad = "right")
}

#' Build a wrapped "<Arm> (N=n)" column label
#' @param name Arm name
#' @param n Population count
#' @param width Wrap width in characters
#' @return Character scalar with "\n" line breaks
arm_label <- function(name, n, width = 10) {
  stringr::str_wrap(sprintf("%s (N=%s)", name, n), width = width)
}

#' Append blank rows to a data frame
#' @param .data Data frame to pad
#' @param n Number of blank rows to append
#' @return `.data` with `n` blank rows appended
pad_row <- function(.data, n = 1) {
  .data[(nrow(.data) + 1):(nrow(.data) + n), ] <- ""
  .data
}

#' Prepend a single blank row matching a data frame's columns
#' @param df Data frame to prepend a spacer row to
#' @return `df` with one leading blank row (house-style header/body gap)
top_spacer <- function(df) {
  blank <- df[1, , drop = FALSE]
  blank[] <- ""
  dplyr::bind_rows(blank, df)
}

#' Compute a Fisher's exact p-value, formatted with a "<.0001" floor
#' @param res Response factor values
#' @param cats Grouping factor values
#' @param width Field width for non-floored values
#' @return Character scalar
fish_p_str <- function(res, cats, width = 10) {
  p <- suppressWarnings(fisher.test(factor(res), factor(cats))$p.value)
  if (round(p, 4) == 0) return("<.0001")
  format(round(p, 4), width = width, nsmall = 4)
}
