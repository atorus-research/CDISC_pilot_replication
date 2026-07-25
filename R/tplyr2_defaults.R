# R/tplyr2_defaults.R — project tplyr2 options (rounding, quantile type) for reference parity

library(tplyr2)

#' Apply the project-level tplyr2 options
#'
#' Sets rounding and quantile behaviour to match the reference R RTFs.
#'
#' @return `TRUE`, invisibly.
apply_cdisc_tplyr2_defaults <- function() {
  tplyr2::tplyr2_options(
    IBMRounding   = FALSE,   # match the reference R RTFs (banker's rounding)
    quantile_type = 7L       # base-R default (SAS PROC UNIVARIATE ~= 3)
  )
  invisible(TRUE)
}

apply_cdisc_tplyr2_defaults()
