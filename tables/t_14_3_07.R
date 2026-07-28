# t_14_3_07.R
# Table 14-3.07: ADAS Cog (11) - Change from Baseline to Week 24 - Completers at Wk 24-Observed Cases-Windowed   (Population: Completers)
# Produces: outputs/14-3.07.docx
# ADAS ANCOVA shape (adadas / PARAMCD ACTOT) restricted to Week 24 completers, via build_efficacy_table().
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")

build_efficacy_table("14-3.07", "programs/t-14-3-07.R", "adadas", "ACTOT", 24,
                     endpoint = "ADAS",
                     # completers: a subject-level population criterion, so it also scopes
                     # the header N via the population spec; DTYPE selects records only.
                     extra_filter = COMP24FL == "Y" & DTYPE != "LOCF",
                     pop_filter = COMP24FL == "Y")
