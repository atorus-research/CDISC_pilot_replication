# R/tplyr2_defaults.R
# -----------------------------------------------------------------------------
# Project-level tplyr2 options.
#
# The reference outputs were produced with base-R rounding (round(), i.e.
# round-half-to-even / "banker's"), so we KEEP tplyr2's default
# IBMRounding = FALSE to match them (NOT SAS round-half-up). quantile_type is
# left at 7 (base-R default) since the sample tables use median/min/max only.
# Set explicitly here to document intent and make it easy to flip for SAS parity.
# -----------------------------------------------------------------------------

library(tplyr2)

apply_cdisc_tplyr2_defaults <- function() {
  tplyr2::tplyr2_options(
    IBMRounding   = FALSE,   # match the reference R RTFs (banker's rounding)
    quantile_type = 7L       # base-R default (SAS PROC UNIVARIATE ~= 3)
  )
  invisible(TRUE)
}

apply_cdisc_tplyr2_defaults()
