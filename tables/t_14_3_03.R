# t_14_3_03.R
# Table 14-3.03: ADAS Cog (11) - Change from Baseline to Week 8 - LOCF  (Population: Efficacy)
# Produces: outputs/14-3.03.docx
# Source: adadas (PARAMCD ACTOT); tplyr2 descriptive summary + ANCOVA with baseline covariate, via build_efficacy_table().

source("R/setup.R")
source("R/helpers.R")
source("R/efficacy.R")

build_efficacy_table("14-3.03", source_path = "programs/t-14-3-03.R",
                     dataset = "adadas", paramcd = "ACTOT", week = 8, endpoint = "ADAS")
