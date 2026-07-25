# tables/t_14_3_12.R
# Table 14-3.12 — Mean NPI-X Total Score from Week 4 through Week 24 (Windowed), ANCOVA
#
# An ADAS-style ANCOVA (change-from-baseline response, baseline covariate) with a
# bespoke shape, expressed through build_efficacy_table()'s optional parameters:
#   * dataset adnpix / PARAMCD == "NPTOTMN"; the base filter drops ANL01FL
#     (anl01 = FALSE) and CHG is derived (derive_chg = TRUE).
#   * TWO descriptive blocks, both AVAL: "Baseline" (AVISITN 0) and
#     "Mean of Weeks 4-24" (AVISITN 98) — supplied via `blocks`.
#   * the ANCOVA model is fit on CHG at AVISITN 98 with the baseline covariate
#     (week = 98, use_base = TRUE).
#
# DOCUMENTED DISCREPANCY: programs/t-14-3-12.R notes the "Mean of Weeks 4-24"
# counts cannot be reproduced exactly because the derived NPTOTMN windowing
# depends on an original analysis dataset NOT shipped in the CDISC pilot package.
# Both the legacy RTF and this build derive NPTOTMN from the SAME shipped adnpix,
# so they agree with each other; the discrepancy is against the unavailable
# original SAS analysis data.
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")

build_efficacy_table("14-3.12", "programs/t-14-3-12.R", "adnpix", "NPTOTMN", 98,
                     anl01 = FALSE, derive_chg = TRUE, use_base = TRUE,
                     blocks = list(
                       list(var = "AVAL", label = "Baseline",           visit = 0),
                       list(var = "AVAL", label = "Mean of Weeks 4-24", visit = 98)))
