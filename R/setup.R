# R/setup.R — sourced by every table program; loads the stack and applies project defaults

suppressPackageStartupMessages({
  library(tidyverse)
  library(haven)
  library(tplyr2)
  library(clinify)
  library(flextable)
  library(officer)
})

# Resolve the project root regardless of the caller's working directory.
.PROJ <- tryCatch(
  {
    # walk up until we find renv.lock
    d <- getwd()
    while (!file.exists(file.path(d, "renv.lock")) && dirname(d) != d) d <- dirname(d)
    d
  },
  error = function(e) getwd()
)

DATA_DIR    <- file.path(.PROJ, "data")
OUTPUT_DIR  <- file.path(.PROJ, "outputs")
TITLES_XLSX <- file.path(DATA_DIR, "titles.xlsx")
# Pristine master worktree holding the reference RTFs (for fidelity comparison).
REF_DIR     <- file.path(dirname(.PROJ), "CDISC_pilot_replication-original", "outputs")

source(file.path(.PROJ, "R", "tplyr2_defaults.R"))
source(file.path(.PROJ, "R", "clinify_defaults.R"))
source(file.path(.PROJ, "R", "titles.R"))

#' Read an ADaM dataset from `data/` by name
#'
#' @param name Dataset name without extension (e.g. `"adsl"`).
#' @return A tibble read from `data/<name>.xpt`.
read_adam <- function(name) {
  haven::read_xpt(file.path(DATA_DIR, paste0(name, ".xpt")))
}

# Footer timestamp: NULL -> each table reproduces its own reference-RTF timestamp (fidelity mode).
FIDELITY_DATE <- NULL
