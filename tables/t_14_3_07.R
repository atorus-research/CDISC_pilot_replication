# tables/t_14_3_07.R
# Table 14-3.07 — ADAS-Cog(11) Change from Baseline to Week 24, COMPLETERS
#   (Observed Cases, Windowed)
#
# Structurally identical to 14-3.01 (the ADAS ANCOVA shape); the only difference
# is the COMPLETERS subset, passed via build_efficacy_table()'s `extra_filter`.
# Completers filter (verbatim from programs/t-14-3-07.R):
#   COMP24FL=="Y" & EFFFL=="Y" & PARAMCD=="ACTOT" & ANL01FL=="Y" & DTYPE!="LOCF"
#
# DATA NOTE (documented discrepancy, reproduced faithfully):
#   The legacy program flags the baseline n on Placebo/Xan-Low as "off by 1"
#   (60 and 28) relative to the completer week-24 n (59 and 27). This is inherent
#   to the source data: one Placebo and one Xan-Low completer have a baseline
#   record but no windowed week-24 observed-case assessment, so they count in
#   Baseline but not Week 24 / Change. Our build reproduces the reference RTF
#   exactly (Baseline 60/28/30; Week 24 & Change 59/27/30) — no divergence.
source("R/setup.R"); source("R/helpers.R"); source("R/efficacy.R")

build_efficacy_table("14-3.07", "programs/t-14-3-07.R", "adadas", "ACTOT", 24,
                     endpoint = "ADAS",
                     extra_filter = COMP24FL == "Y" & DTYPE != "LOCF")
