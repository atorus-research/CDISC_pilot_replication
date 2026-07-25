# t_14_5_02.R
# Table 14-5.02: Incidence of Treatment Emergent Serious Adverse Events by Treatment Group   (Population: Safety)
# Produces: outputs/14-5.02.docx
# Source: ADAE (serious only) + ADSL; TEAE counts by SOC/PT with Fisher's exact p-values (build_ae_table).
source("R/setup.R"); source("R/helpers.R"); source("R/ae.R")
build_ae_table("14-5.02", source_path = "programs/t-14-5-02.R", serious = TRUE)
