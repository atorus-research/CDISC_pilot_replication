# Table 14-3.05 — ADAS-Cog(11) Change from Baseline to Week 16, LOCF
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.05", source_path = "programs/t-14-3-05.R",
                     dataset = "adadas", paramcd = "ACTOT", week = 16, endpoint = "ADAS")
