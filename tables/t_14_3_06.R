# t_14_3_06.R
# Table 14-3.06: CIBIC+ - Summary at Week 16 - LOCF  (Population: Efficacy)
# Produces: outputs/14-3.06.docx
# Source: adcibc (PARAMCD CIBICVAL); tplyr2 descriptive summary + ANCOVA without baseline covariate, via build_efficacy_table().

source("R/setup.R")
source("R/helpers.R")
source("R/efficacy.R")

build_efficacy_table("14-3.06", source_path = "programs/t-14-3-06.R",
                     dataset = "adcibc", paramcd = "CIBICVAL", week = 16, endpoint = "CIBIC")
