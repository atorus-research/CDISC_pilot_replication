# t_14_3_02.R
# Table 14-3.02: Primary Endpoint Analysis: CIBIC+ - Summary at Week 24 - LOCF  (Population: Efficacy)
# Produces: outputs/14-3.02.docx
# Source: adcibc (PARAMCD CIBICVAL); tplyr2 descriptive summary + ANCOVA without baseline covariate, via build_efficacy_table().

source("R/setup.R")
source("R/helpers.R")
source("R/efficacy.R")

build_efficacy_table("14-3.02", source_path = "programs/t-14-3-02.R",
                     dataset = "adcibc", paramcd = "CIBICVAL", week = 24, endpoint = "CIBIC")
