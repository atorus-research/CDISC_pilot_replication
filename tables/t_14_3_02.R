# Table 14-3.02 — Primary Endpoint: CIBIC+ Summary at Week 24, LOCF
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.02", source_path = "programs/t-14-3-02.R",
                     dataset = "adcibc", paramcd = "CIBICVAL", week = 24, endpoint = "CIBIC")
