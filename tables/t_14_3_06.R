# Table 14-3.06 — CIBIC+ Summary at Week 16, LOCF
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.06", source_path = "programs/t-14-3-06.R",
                     dataset = "adcibc", paramcd = "CIBICVAL", week = 16, endpoint = "CIBIC")
