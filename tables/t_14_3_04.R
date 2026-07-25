# Table 14-3.04 — CIBIC+ Summary at Week 8, LOCF
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.04", source_path = "programs/t-14-3-04.R",
                     dataset = "adcibc", paramcd = "CIBICVAL", week = 8, endpoint = "CIBIC")
