# tables/t_14_3_01.R
# Table 14-3.01 — Primary Endpoint: ADAS-Cog(11) Change from Baseline to Week 24, LOCF
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.01", source_path = "programs/t-14-3-01.R",
                     dataset = "adadas", paramcd = "ACTOT", week = 24, endpoint = "ADAS")
