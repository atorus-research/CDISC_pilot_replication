# t_14_3_08.R
# Table 14-3.08: ADAS Cog (11) - Change from Baseline to Week 24 in Male Subjects - LOCF   (Population: Efficacy)
# Produces: outputs/14-3.08.docx
# ADAS ANCOVA (adadas / PARAMCD ACTOT), Week 24 LOCF, male subjects, via build_efficacy_table().
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")
build_efficacy_table("14-3.08", source_path = "programs/t-14-3-08.R",
                     dataset = "adadas", paramcd = "ACTOT", week = 24, endpoint = "ADAS", sex = "M")
