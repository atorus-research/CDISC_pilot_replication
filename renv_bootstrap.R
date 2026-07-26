# renv_bootstrap.R — one-time environment setup for the next-gen build.
# Seeds the project library from existing installs (fast, no recompile), then
# installs the atorus dev packages from their local source clones.
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

setwd("/Users/mstackhouse/Documents/repos/CDISC_pilot_replication")

# bare = no code scanning (legacy programs reference packages we are dropping)
renv::init(bare = TRUE, restart = FALSE)

# Copy already-installed CRAN packages (and their deps) into the project library
pkgs <- c(
  "tidyverse", "dplyr", "tidyr", "purrr", "stringr", "tibble", "readr",
  "ggplot2", "forcats", "glue",
  "haven", "readxl",
  "flextable", "officer",
  "emmeans", "car", "mmrm", "broom",
  "data.table", "jsonlite", "rlang", "zoo", "knitr", "htmltools",
  "tidyselect", "magrittr"
)
renv::hydrate(packages = pkgs)

# Install the atorus packages from the local maintained source clones (latest dev)
renv::install("/Users/mstackhouse/Documents/repos/tplyr2")
renv::install("/Users/mstackhouse/Documents/repos/clinify")

# Capture everything in the lockfile
renv::snapshot(type = "all", prompt = FALSE)

cat("\nRENV_BOOTSTRAP_DONE\n")
cat("tplyr2:", as.character(utils::packageVersion("tplyr2")), "\n")
cat("clinify:", as.character(utils::packageVersion("clinify")), "\n")
cat("mmrm:", as.character(utils::packageVersion("mmrm")), "\n")
