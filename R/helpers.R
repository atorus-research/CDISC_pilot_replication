# R/helpers.R
# -----------------------------------------------------------------------------
# Minimal shared helpers. Kept intentionally small: summarization/formatting is
# tplyr2's job. These cover only the pieces tplyr2 does not own — namely
# formatting externally-fit MODEL output (LS-means diffs, CIs, p-values) that is
# row-bound beneath a tplyr2 build, and building wrapped (N=) column labels.
# -----------------------------------------------------------------------------

library(stringr)

# Fixed-width monospace number formatter, byte-for-byte compatible with the
# reference programs' num_fmt() so row-bound model strings align identically.
num_fmt <- function(var, digits = 0, size = 10, int_len = 3) {
  if (is.na(var)) return("")
  nsmall <- digits
  if (digits > 0) digits <- digits + 1
  stringr::str_pad(
    format(round(var, nsmall), width = (int_len + digits), nsmall = nsmall),
    side = "right", width = size
  )
}
num_fmt <- Vectorize(num_fmt)

# Wrapped "<Arm> (N=n)" column label, matching get_header_n()'s str_wrap(width=10).
# Returns a single string with "\n" line breaks (clinify renders \n as a break).
arm_label <- function(name, n, width = 10) {
  stringr::str_wrap(sprintf("%s (N=%s)", name, n), width = width)
}

# Append n blank rows to a data.frame (reference pad_row()).
pad_row <- function(.data, n = 1) {
  .data[(nrow(.data) + 1):(nrow(.data) + n), ] <- ""
  .data
}

# Prepend a single blank row matching `df`'s columns. Reproduces the reference's
# gap between the header rule and the first body row (a house-style spacer).
top_spacer <- function(df) {
  blank <- df[1, , drop = FALSE]
  blank[] <- ""
  dplyr::bind_rows(blank, df)
}

# One-way ANOVA p-value vs treatment, formatted like the reference (width 10, 4dp).
aov_p_str <- function(data, var, trt = "TRT01P") {
  f <- stats::as.formula(paste(var, "~", trt))
  p <- summary(stats::aov(f, data, na.action = stats::na.omit))[[1]][["Pr(>F)"]][1]
  format(round(p, 4), width = 10, nsmall = 4)
}

# Fisher's exact p-value, reference format ("<.0001" floor; width controls the
# right-justified field for non-floored values).
fish_p_str <- function(res, cats, width = 10) {
  p <- suppressWarnings(fisher.test(factor(res), factor(cats))$p.value)
  if (round(p, 4) == 0) return("<.0001")
  format(round(p, 4), width = width, nsmall = 4)
}

# Pearson chi-square p-value vs treatment, reference format ("<.0001" floor).
chi_p_str <- function(data, var, trt = "TRT01P") {
  p <- suppressWarnings(
    stats::chisq.test(factor(data[[var]]), factor(data[[trt]]))$p.value)
  if (round(p, 4) == 0) return("<.0001")
  format(round(p, 4), width = 10, nsmall = 4)
}

# NOTE: the former fix_count_quirks() (regex post-processing for "<1%" and bare-0
# cells) was removed once tplyr2 PR #15 added native `pct_lt` / `zero_count_display`
# layer settings (tplyr2 issues #13/#14). Those are used directly in the tables now.
