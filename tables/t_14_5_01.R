# t_14_5_01.R
# Table 14-5.01: Incidence of Treatment Emergent Adverse Events by Treatment Group   (Population: Safety)
# Produces: outputs/14-5.01.docx
# Source: ADAE + ADSL; TEAE counts by SOC/PT with Fisher's exact p-values (build_ae_table).
source("R/setup.R"); source("R/helpers.R"); source("R/ae.R")
build_ae_table("14-5.01", source_path = "programs/t-14-5-01.R", serious = FALSE)
