# R/setup.R
# -----------------------------------------------------------------------------
# Central entry point for the next-gen CDISC Pilot build.
# Each table program starts with: source(here::here("R", "setup.R"))
# (or source("R/setup.R") when run from the project root).
#
# Loads the core stack (tplyr2, clinify, tidyverse, haven), applies the
# project-level tplyr2 options and clinify house style, and exposes the
# titles.xlsx reader.
# -----------------------------------------------------------------------------

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
    # rprojroot-free: walk up until we find the .Rproj / renv.lock
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

# Small helper: read an ADaM xpt by name from data/
read_adam <- function(name) {
  haven::read_xpt(file.path(DATA_DIR, paste0(name, ".xpt")))
}

# The reference RTFs baked a per-table timestamp into the footer. For the
# fidelity comparison we reproduce each table's own timestamp verbatim (extracted
# from its reference RTF); in production this would be Sys.time(). NULL means
# "look it up per table" (see add_titles_footnotes / ref_timestamp in R/titles.R).
FIDELITY_DATE <- NULL
