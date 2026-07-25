# Table 14-3.08 — ADAS-Cog(11) Change from Baseline to Week 24, Male, LOCF
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.08", source_path = "programs/t-14-3-08.R",
                     dataset = "adadas", paramcd = "ACTOT", week = 24, endpoint = "ADAS", sex = "M")
