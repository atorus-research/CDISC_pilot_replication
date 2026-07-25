# tables/t_14_2_01.R
# Table 14-2.01 — Summary of Demographic and Baseline Characteristics (ITT)
#
# Continuous characteristics: tplyr2 group_desc (n / Mean / SD / Median / Min / Max).
# Categorical characteristics: tplyr2 group_count (n (%)) with reference display
#   quirks (<1%, bare 0) reproduced. Total column via a bound "Total" arm.
# p-value column: one-way ANOVA (continuous) / Pearson chi-square (categorical),
#   computed on the 3 real arms and placed on the appropriate row.
# -----------------------------------------------------------------------------
source("R/setup.R")
source("R/helpers.R")

TABLE  <- "14-2.01"
SOURCE <- "programs/t-14-2-01.R"

ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose", "Total")

# ---- Data + display derivations (match reference) ---------------------------
adsl0 <- read_adam("adsl") |>
  filter(ITTFL == "Y") |>
  mutate(
    RACE_DISPLAY = case_when(
      ETHNIC == "HISPANIC OR LATINO"            ~ "Hispanic",
      RACE == "WHITE"                           ~ "Caucasian",
      RACE == "BLACK OR AFRICAN AMERICAN"       ~ "African Descent",
      RACE == "AMERICAN INDIAN OR ALASKA NATIVE" ~ "Other"),
    SEX      = if_else(SEX == "M", "Male", "Female"),
    DURDSGR1 = paste(DURDSGR1, "months"),
    AGEGR1   = paste(AGEGR1, "yrs")
  )

# 3-arm frame for p-values; +Total frame for display
adsl3 <- adsl0
adslT <- bind_rows(adsl0, mutate(adsl0, TRT01P = "Total")) |>
  mutate(
    TRT01P       = factor(TRT01P, levels = ARMS),
    AGEGR1       = factor(AGEGR1,   levels = c("<65 yrs", "65-80 yrs", ">80 yrs")),
    SEX          = factor(SEX,      levels = c("Male", "Female")),
    RACE_DISPLAY = factor(RACE_DISPLAY, levels = c("Caucasian", "African Descent", "Hispanic", "Other")),
    DURDSGR1     = factor(DURDSGR1, levels = c("<12 months", ">=12 months")),
    BMIBLGR1     = factor(BMIBLGR1, levels = c("<25", "25-<30", ">=30"))
  )

# ---- Block builders via tplyr2 ----------------------------------------------
desc_fs <- list(
  "n"      = f_str("xxx",    "n"),
  "Mean"   = f_str("xxx.x",  "mean"),
  "SD"     = f_str("xxx.xx", "sd"),
  "Median" = f_str("xxx.x",  "median"),
  "Min"    = f_str("xxx.x",  "min"),
  "Max"    = f_str("xxx.x",  "max")
)

# Reorder tplyr2 result columns into ARMS order by their (N=) label. tplyr2 count
# layers order res columns alphabetically by the cols value (ignoring factor
# levels), which differs from desc layers - so never assume res1==first arm.
reslist <- function(b) {
  # tplyr2 (>= PR #15) orders count/desc/shift result columns by the cols factor
  # levels, so res1..resN are already in ARMS order (Placebo/Low/High/Total).
  b <- b[order(b$ord_layer_1), , drop = FALSE]
  rescols <- grep("^res", names(b), value = TRUE)
  out <- tibble(rowlbl2 = as.character(b$rowlabel1))
  for (i in seq_along(rescols)) out[[paste0("res", i)]] <- as.character(b[[rescols[i]]])
  out
}

desc_block <- function(var) {
  s <- tplyr_spec(cols = "TRT01P",
                  layers = tplyr_layers(group_desc(var,
                    settings = layer_settings(format_strings = desc_fs))))
  pad_row(reslist(tplyr_build(s, adslT)))
}

count_block <- function(var, incl_n = FALSE) {
  # Native tplyr2 display conventions (PR #15): pct_lt=1 renders sub-1% as "<1",
  # zero_count_display="count_only" renders zero-count cells as a bare count.
  # Category ROWS ordered by the target's factor levels via order_count_method
  # ="byfactor" (tplyr2 PR #17 fixed byfactor to honor factor levels).
  s <- tplyr_spec(cols = "TRT01P",
                  layers = tplyr_layers(group_count(var,
                    settings = layer_settings(
                      format_strings = list(n_counts = f_str("xxx (xxx%)", "n", "pct")),
                      order_count_method = "byfactor",
                      pct_lt = 1, zero_count_display = "count_only"))))
  blk <- reslist(tplyr_build(s, adslT))
  if (incl_n) {
    nrow_df <- adslT |> filter(!is.na(.data[[var]])) |> count(TRT01P, .drop = FALSE)
    getn <- function(a) sprintf("%3d", nrow_df$n[nrow_df$TRT01P == a])
    blk <- bind_rows(
      tibble(rowlbl2 = "n", res1 = getn(ARMS[1]), res2 = getn(ARMS[2]),
             res3 = getn(ARMS[3]), res4 = getn(ARMS[4])),
      blk)
  }
  pad_row(blk)
}

# stamp rowlbl1 (first row) and p-values (by row index) onto a block
stamp <- function(block, name, p_at = list()) {
  block$rowlbl1 <- ""; block$rowlbl1[1] <- name
  block$pcol <- ""
  for (i in names(p_at)) block$pcol[as.integer(i)] <- p_at[[i]]
  select(block, rowlbl1, rowlbl2, res1, res2, res3, res4, pcol)
}

# ---- Characteristics --------------------------------------------------------
age <- stamp(bind_rows(desc_block("AGE"), count_block("AGEGR1")), "Age (y)",
             list("1" = aov_p_str(adsl3, "AGE"), "8" = chi_p_str(adsl3, "AGEGR1")))
sex <- stamp(count_block("SEX", incl_n = TRUE), "Sex",
             list("1" = chi_p_str(adsl3, "SEX")))
race <- stamp(count_block("RACE_DISPLAY", incl_n = TRUE), "Race (Origin)",
              list("1" = chi_p_str(adsl3, "RACE_DISPLAY")))
mmse <- stamp(desc_block("MMSETOT"), "MMSE", list("1" = aov_p_str(adsl3, "MMSETOT")))
dur  <- stamp(bind_rows(desc_block("DURDIS"), count_block("DURDSGR1")), "Duration of disease ",
              list("1" = aov_p_str(adsl3, "DURDIS"), "8" = chi_p_str(adsl3, "DURDSGR1")))
educ <- stamp(desc_block("EDUCLVL"), "Years of education", list("1" = aov_p_str(adsl3, "EDUCLVL")))
wt   <- stamp(desc_block("WEIGHTBL"), "Baseline weight(kg)", list("1" = aov_p_str(adsl3, "WEIGHTBL")))
ht   <- stamp(desc_block("HEIGHTBL"), "Baseline height(cm)", list("1" = aov_p_str(adsl3, "HEIGHTBL")))
bmi  <- stamp(bind_rows(desc_block("BMIBL"), count_block("BMIBLGR1")), "Baseline BMI",
              list("1" = aov_p_str(adsl3, "BMIBL"), "8" = chi_p_str(adsl3, "BMIBLGR1")))

final <- top_spacer(bind_rows(age, sex, race, mmse, dur, educ, wt, ht, bmi))

# ---- clinify render ---------------------------------------------------------
Ns <- adsl0 |> count(TRT01P) |> deframe()
Ntot <- sum(Ns)
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    rowlbl1 = "", rowlbl2 = "",
    res1 = arm_label("Placebo", Ns[["Placebo"]]),
    res2 = arm_label("Xanomeline Low Dose",  Ns[["Xanomeline Low Dose"]]),
    res3 = arm_label("Xanomeline High Dose", Ns[["Xanomeline High Dose"]]),
    res4 = arm_label("Total", Ntot),
    pcol = "p-value\n[1]"
  ) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = c("rowlbl1", "rowlbl2"), part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "rowlbl1", width = 1.8) |>
  flextable::width(j = "rowlbl2", width = 1.8) |>
  flextable::width(j = "res1", width = 1.081) |>
  flextable::width(j = c("res2", "res3", "res4", "pcol"), width = 1.08)

ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
cat("wrote", file.path(OUTPUT_DIR, paste0(TABLE, ".docx")), "\n")
