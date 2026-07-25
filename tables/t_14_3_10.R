# tables/t_14_3_10.R
# Table 14-3.10 — ADAS-Cog(11): Mean and Mean Change from Baseline over Time
#
# WIDE descriptive over Windowed + LOCF visits. Arms are ROW BLOCKS (Placebo,
# Xan.Low, Xan.High), each an 8-slot visit ladder (Baseline, 3 Windowed weeks,
# 3 LOCF weeks, spacer). Columns: AVAL n/Mean/Std/Med./Min./Max., baseline
# Mean/Std, then a spanned "Change from baseline" block (CHG Mean/Std/Med./Min./
# Max.). Baseline rows carry AVAL stats only (no change stats).
#
# The summary is a plain grouped dplyr build formatted with the shared byte-exact
# num_fmt() — identical stat engine and sizes to the legacy programs/t-14-3-10.R,
# so displayed values match the reference character-for-character. Only the
# rendering is re-expressed with clinify (multi-row header + spanner) in place of
# the legacy huxtable/pharmaRTF.
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-3.10"; SOURCE <- "programs/t-14-3-10.R"

# --- data (reproduces programs/t-14-3-10.R exactly) --------------------------
adas <- read_adam("adadas") |>
  filter(EFFFL == "Y", PARAMCD == "ACTOT", ITTFL == "Y",
         AVISITN %in% c(0, 8, 16, 24), ANL01FL == "Y") |>
  mutate(SET = "LOCF") |>
  select(TRTPN, TRTP, AVISIT, AVISITN, AVAL, BASE, CHG, DTYPE, SET)

# Display order of the 8 per-arm visit slots (slot 8 is a blank spacer row).
visits <- tibble(
  ORD    = rep(1:8, 3),
  AVISIT = rep(c("Baseline", "Week 8 (Windowed)", "Week 16 (Windowed)",
                 "Week 24 (Windowed)", "Week 8 LOCF", "Week 16 LOCF",
                 "Week 24 LOCF", ""), 3),
  TRTPN  = c(rep(0, 8), rep(54, 8), rep(81, 8)))

# LOCF set + observed-post-baseline WINDOWED set; format AVISIT and TRTP labels.
step1 <- adas |>
  bind_rows(adas |> filter(AVISITN != 0, DTYPE != "LOCF") |> mutate(SET = "WINDOWED")) |>
  mutate(
    AVISIT = case_when(
      SET == "WINDOWED" & AVISITN != 0 ~ paste(AVISIT, "(Windowed)"),
      SET == "LOCF"     & AVISITN != 0 ~ paste(AVISIT, "LOCF"),
      TRUE ~ AVISIT),
    TRTP = case_when(TRTPN == 0 ~ "Placebo", TRTPN == 54 ~ "Xan.Low",
                     TRTPN == 81 ~ "Xan.High"))

# AVAL descriptive stats (all visit slots) — joined to the 8-slot ladder.
aval <- step1 |>
  group_by(TRTPN, TRTP, AVISITN, AVISIT, SET) |>
  summarize(
    n    = num_fmt(n(),           int_len = 2, size = 2),
    mean = num_fmt(mean(AVAL),    digits = 1, int_len = 2, size = 4),
    sd   = num_fmt(sd(AVAL),      digits = 2, int_len = 2, size = 5),
    md   = num_fmt(median(AVAL),  digits = 1, int_len = 2, size = 4),
    mn   = num_fmt(min(AVAL),     int_len = 2, size = 4),
    mx   = num_fmt(max(AVAL),     int_len = 2, size = 4),
    .groups = "drop") |>
  full_join(visits, by = c("TRTPN", "AVISIT"))

# Baseline mean/std + change-from-baseline stats (non-baseline visits only).
chg <- step1 |>
  group_by(TRTPN, TRTP, AVISITN, AVISIT, SET) |>
  filter(AVISITN != 0) |>
  summarize(
    meanc = num_fmt(mean(CHG),   digits = 1, int_len = 1, size = 4),
    sdc   = num_fmt(sd(CHG),     digits = 2, int_len = 1, size = 4),
    mdc   = num_fmt(median(CHG), digits = 1, int_len = 1, size = 4),
    mnc   = num_fmt(min(CHG),    int_len = 3, size = 4),
    mxc   = num_fmt(max(CHG),    int_len = 2, size = 4),
    bmn   = num_fmt(mean(BASE),  digits = 1, int_len = 2, size = 4),
    bsd   = num_fmt(sd(BASE),    digits = 2, int_len = 2, size = 5),
    .groups = "drop")

final <- left_join(aval, chg, by = c("TRTPN", "TRTP", "AVISITN", "AVISIT", "SET")) |>
  arrange(TRTPN, ORD) |>
  ungroup() |>
  mutate(TRTP = ifelse(ORD == 1, TRTP, "")) |>   # arm label only on the Baseline row
  select(TRTP, AVISIT, n, mean, sd, md, mn, mx, bmn, bsd, meanc, sdc, mdc, mnc, mxc)
final[] <- lapply(final, function(x) { x[is.na(x)] <- ""; as.character(x) })

# --- render: two-row header, spanned "Change from baseline" block ------------
# Header row 1 = spanner (only over the 5 CHG columns; blank elsewhere); row 2 =
# labels, with "Bsln\nMean"/"Bsln\nStd" as stacked two-line labels. valign=bottom
# so the single-line labels sit on the baseline and "Bsln" rides one line above.
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    TRTP   = c("", ""),
    AVISIT = c("", ""),
    n      = c("", "nc"),
    mean   = c("", "Mean"),
    sd     = c("", "Std"),
    md     = c("", "Med."),
    mn     = c("", "Min."),
    mx     = c("", "Max."),
    bmn    = c("", "Bsln\nMean"),
    bsd    = c("", "Bsln\nStd"),
    meanc  = c("---Change from baseline---", "Mean"),
    sdc    = c("---Change from baseline---", "Std"),
    mdc    = c("---Change from baseline---", "Med."),
    mnc    = c("---Change from baseline---", "Min."),
    mxc    = c("---Change from baseline---", "Max.")) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(part = "body", align = "center") |>            # numeric cells centered (legacy)
  flextable::align(j = c("TRTP", "AVISIT"), part = "body", align = "left")

# Column widths reproduce the reference's per-column ratios exactly (the legacy
# set_col_width() proportions) scaled to the 9.0" landscape table width, so every
# centered numeric cell lands at the reference x-position. Courier New 10pt
# (~0.083"/char) fits with no wrapping; the 5 CHG columns (2.43") carry the
# 26-char spanner on one line.
wv <- 9.0 * c(TRTP = .09, AVISIT = .19, n = .05,
              mean = .06, sd = .06, md = .06, mn = .05, mx = .05,
              bmn = .06, bsd = .06,
              meanc = .06, sdc = .05, mdc = .05, mnc = .05, mxc = .06)
# The reference table renders ~2pt narrower than ratio*9" (huxtable does not lay
# the ratios out to exactly ratio*9): its numeric grid sits 2pt left of ours with
# the stub flush. Trim AVISIT by 2pt so the whole numeric block slides onto the
# reference x-grid (stub left edge is unchanged), matching the reference's ~8.97"
# total and aligning every centred cell.
wv["AVISIT"] <- wv["AVISIT"] - 2/72
ct <- flextable::width(ct, j = names(wv), width = unname(wv))

ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)

# The house header buffer (padding.top = 18) is applied to every header row,
# which floats the single-line spanner row far above the label rows. The
# reference instead clusters the three header lines (spanner / "Bsln" / labels)
# tightly and carries the buffer above the spanner. Redistribute the header
# padding per-table so the three lines sit at the reference y-positions while the
# body (already row-pitch aligned) stays put: buffer above the spanner row, tight
# gap to the label row, and the buffer's slack moved below the labels.
sd <- function(x, ...) {
  x <- cdisc_table_default(x)
  # The reference used set_ignore_cell_padding(TRUE) (zero cell padding). The
  # house default's 2pt L/R padding shifts the whole body grid ~2pt right of the
  # reference (and off from the centred title block). Match the reference's zero
  # padding: the body grid sits flush at the left margin and every centred
  # numeric cell lands on the reference x-position.
  x <- flextable::padding(x, padding.left = 0, padding.right = 0, part = "all")
  x <- flextable::padding(x, i = 1, padding.top = 21, padding.bottom = 0, part = "header")
  # Row2 bottom padding is kept tiny so the header rule sits tight under the label
  # row (as in the reference, ~144pt) rather than being pushed down against the
  # body; the rule->body gap is instead carried by the first body row's top
  # padding. (Aligning the strong header rule is what stops the fidelity harness's
  # row-projection registration from mis-locking by one row on this 24-row grid.)
  x <- flextable::padding(x, i = 2, padding.top = 2,  padding.bottom = 1, part = "header")
  x <- flextable::padding(x, i = 1, padding.top = 18, padding.bottom = 2, part = "body")
  # The reference table is flush-left at the margin; flextable defaults to a
  # page-centred table, which (against the centred title block) leaves the body
  # grid a couple of points off. Pin it left so the grid starts at the margin.
  x <- flextable::set_table_properties(x, layout = "fixed", align = "left")
  x
}
old <- options(clinify_table_default = sd)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), "  width:", sum(wv), "\n")
