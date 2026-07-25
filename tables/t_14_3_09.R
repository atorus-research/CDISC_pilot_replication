# t_14_3_09.R
# Table 14-3.09: ADAS Cog (11) - Change from Baseline to Week 24 in Female Subjects - LOCF   (Population: Efficacy)
# Produces: outputs/14-3.09.docx
# ADAS ANCOVA (adadas / PARAMCD ACTOT), Week 24 LOCF, female subjects, via build_efficacy_table().
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.09", source_path = "programs/t-14-3-09.R",
                     dataset = "adadas", paramcd = "ACTOT", week = 24, endpoint = "ADAS", sex = "F")
