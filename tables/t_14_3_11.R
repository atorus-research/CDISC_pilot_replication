# tables/t_14_3_11.R
# Table 14-3.11 — ADAS-Cog(11) Repeated Measures Analysis of Change from Baseline
#   to Week 24 (Efficacy). Model-only display (LS Means + pairwise contrasts).
#
# Fit with mmrm() using an UNSTRUCTURED within-subject covariance + Kenward-Roger
# (matches SAS PROC MIXED / TYPE=UN). This deliberately differs from the legacy
# RTF, which used lme4 with a random-slope term and was documented as not
# matching SAS. Per the project's "correct modern values" policy, mmrm is used.
source("R/setup.R"); source("R/helpers.R")
suppressPackageStartupMessages({library(mmrm); library(emmeans)})

TABLE <- "14-3.11"; SOURCE <- "programs/t-14-3-11.R"

adas <- read_adam("adadas") |>
  filter(EFFFL == "Y", PARAMCD == "ACTOT", ANL01FL == "Y", DTYPE != "LOCF", AVISITN > 0) |>
  mutate(TRTPCD_F = factor(case_when(TRTPN == 0 ~ "Pbo", TRTPN == 54 ~ "Xan_Lo",
                                     TRTPN == 81 ~ "Xan_Hi"),
                           levels = c("Xan_Hi", "Xan_Lo", "Pbo")),
         AWEEKC  = factor(AVISIT), SITEGR1 = factor(SITEGR1), USUBJID = factor(USUBJID))
hn <- adas |> distinct(USUBJID, TRTP) |>
  mutate(TRTP = factor(TRTP, levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))) |>
  count(TRTP) |> deframe()

fit <- mmrm(CHG ~ SITEGR1 + TRTPCD_F * AWEEKC + BASE * AWEEKC + us(AWEEKC | USUBJID),
            data = adas, reml = TRUE, method = "Kenward-Roger", vcov = "Kenward-Roger-Linear")

# Week-24 LS means (the title's endpoint; the original report averaged over weeks
# via the main-effect lsmeans, a documented title/number mismatch - see issue #6).
# DEFAULT (equal) weights over SITEGR1 to match SAS `lsmeans` (no OM option);
# validated to the digit against the SAS output in issue #6.
emm <- emmeans(fit, ~ TRTPCD_F, at = list(AWEEKC = "Week 24"))
lsm <- as.data.frame(emm)
ls_of <- function(lvl) {
  r <- lsm[lsm$TRTPCD_F == lvl, ]
  sprintf("%s (%s)", num_fmt(r$emmean, int_len = 1, digits = 1, size = 3),
                     num_fmt(r$SE,     int_len = 1, digits = 2, size = 4))
}
pw <- as.data.frame(summary(contrast(emm, method = list(
  "Xan_Lo - Pbo"    = c(0, 1, -1), "Xan_Hi - Pbo" = c(1, 0, -1),
  "Xan_Hi - Xan_Lo" = c(1, -1, 0)), adjust = NULL), infer = c(TRUE, TRUE)))
p_of  <- function(cn) num_fmt(pw$p.value[pw$contrast == cn],  int_len = 4, digits = 3, size = 12)
d_of  <- function(cn) { r <- pw[pw$contrast == cn, ]
  sprintf("%s (%s)", num_fmt(r$estimate, int_len = 2, digits = 1, size = 4),
                     num_fmt(r$SE,       int_len = 1, digits = 2, size = 4)) }
ci_of <- function(cn) { r <- pw[pw$contrast == cn, ]
  sprintf("(%s;%s)", num_fmt(r$lower.CL, int_len = 2, digits = 1, size = 4),
                     num_fmt(r$upper.CL, int_len = 1, digits = 1, size = 3)) }

blank <- function() tibble(row_label = "", res1 = "", res2 = "", res3 = "")
final <- top_spacer(bind_rows(
  tibble(row_label = "LS Means (SE)", res1 = ls_of("Pbo"), res2 = ls_of("Xan_Lo"), res3 = ls_of("Xan_Hi")),
  blank(),
  tibble(row_label = c("p-value(Xan - Placebo)", "  Diff of LS Means (SE)", "  95% CI"),
         res1 = "",
         res2 = c(p_of("Xan_Lo - Pbo"), d_of("Xan_Lo - Pbo"), ci_of("Xan_Lo - Pbo")),
         res3 = c(p_of("Xan_Hi - Pbo"), d_of("Xan_Hi - Pbo"), ci_of("Xan_Hi - Pbo"))),
  blank(),
  tibble(row_label = c("p-value(Xan High - Xan Low)", "  Diff of LS Means (SE)", "  95% CI"),
         res1 = "", res2 = "",
         res3 = c(p_of("Xan_Hi - Xan_Lo"), d_of("Xan_Hi - Xan_Lo"), ci_of("Xan_Hi - Xan_Lo")))))
final[] <- lapply(final, function(x) { attr(x, "label") <- NULL; as.character(x) })
cat("\n===== 14-3.11 (mmrm) =====\n"); print(as.data.frame(final))

ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    row_label = "",
    res1 = arm_label("Placebo", hn[["Placebo"]]),
    res2 = arm_label("Xanomeline Low Dose",  hn[["Xanomeline Low Dose"]]),
    res3 = arm_label("Xanomeline High Dose", hn[["Xanomeline High Dose"]])) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "row_label", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "row_label", width = 3.6) |>
  flextable::width(j = c("res1", "res2", "res3"), width = 1.2) |>
  flextable::set_table_properties(align = "center")
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
