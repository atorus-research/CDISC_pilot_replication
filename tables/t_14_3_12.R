# t_14_3_12.R
# Table 14-3.12: Mean NPI-X Total Score from Week 4 through Week 24 - Windowed   (Population: Efficacy)
# Produces: outputs/14-3.12.docx
# NPI-X ANCOVA (adnpix / PARAMCD NPTOTMN) on derived change at the Week 4-24 window, via build_efficacy_table().
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")

build_efficacy_table("14-3.12", "programs/t-14-3-12.R", "adnpix", "NPTOTMN", 98,
                     anl01 = FALSE, derive_chg = TRUE, use_base = TRUE,
                     blocks = list(
                       list(var = "AVAL", label = "Baseline",           visit = 0),
                       list(var = "AVAL", label = "Mean of Weeks 4-24", visit = 98)))
