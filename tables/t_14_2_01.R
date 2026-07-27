# t_14_2_01.R
# Table 14-2.01: Summary of Demographic and Baseline Characteristics   (Population: Intent-to-Treat)
# Produces: outputs/14-2.01.docx
# Source: ADSL (ITTFL == "Y"); tplyr2 group_desc/group_count with ANOVA and chi-square p-values.
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R")
source("R/helpers.R")

TABLE  <- "14-2.01"
SOURCE <- "programs/t-14-2-01.R"

ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose", "Total")

# ---- Data + display derivations ---------------------------------------------
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

# Display frame carries factor levels for row/column ordering. The "Total" column is
# emitted natively via total_group() in the block builders below, and the p-values come
# from each layer's assoc_test (which sees the real arms only, not the Total duplicates).
adslG <- adsl0 |>
  mutate(
    TRT01P       = factor(TRT01P, levels = ARMS[1:3]),
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

# p-value column format; the test functions below return the finished string, which
# assoc_test passes through verbatim, so this format is only a placeholder.
P_FMT <- f_str("xxxxxxxxxx", "p")

#' Build the ANOVA test function for a continuous characteristic
#' @param var Name of the continuous variable
#' @return A function of the layer's source subset returning a formatted p-value
aov_fn <- function(var) {
  function(.data) {
    f <- stats::as.formula(paste(var, "~ TRT01P"))
    p <- summary(stats::aov(f, .data, na.action = stats::na.omit))[[1]][["Pr(>F)"]][1]
    format(round(p, 4), width = 10, nsmall = 4)
  }
}

#' Build the chi-square test function for a categorical characteristic
#' @param var Name of the categorical variable
#' @return A function of the layer's source subset returning a formatted p-value
#'   ("<.0001" below the display floor)
chi_fn <- function(var) {
  function(.data) {
    p <- suppressWarnings(
      stats::chisq.test(factor(.data[[var]]), factor(.data$TRT01P))$p.value)
    if (round(p, 4) == 0) return("<.0001")
    format(round(p, 4), width = 10, nsmall = 4)
  }
}

#' Reshape a tplyr2 build into rowlbl2 plus res columns in ARMS order
#' @param b A tplyr2 build result
#' @return Tibble with rowlbl2, res1..resN (Placebo/Low/High/Total) and the
#'   assoc_test p-value column `pcol`
reslist <- function(b) {
  # as_display() returns the display-ready, already-ordered frame (rowlabel*/res*/pval*),
  # dropping the internal ord*/row_id columns — no manual re-sort needed.
  b <- as_display(b)
  rescols <- grep("^res", names(b), value = TRUE)
  out <- tibble(rowlbl2 = as.character(b$rowlabel1))
  for (i in seq_along(rescols)) out[[paste0("res", i)]] <- as.character(b[[rescols[i]]])
  # assoc_test places the p on the block's first display row; blank elsewhere.
  out$pcol <- if ("pval1" %in% names(b)) coalesce(as.character(b$pval1), "") else ""
  out
}

#' Build a descriptive-statistics block (n/Mean/SD/Median/Min/Max) for a variable
#' @param var Name of the continuous variable to summarize
#' @return Padded tibble with rowlbl2 and res1..res4 columns
desc_block <- function(var) {
  s <- tplyr_spec(cols = "TRT01P",
                  total_groups = list(total_group("TRT01P")),
                  layers = tplyr_layers(group_desc(var,
                    settings = layer_settings(format_strings = desc_fs,
                      # ANOVA p across arms; lands on this block's first row ("n").
                      assoc_test = assoc_test(fn = aov_fn(var), format = P_FMT)))))
  pad_row(reslist(tplyr_build(s, adslG)))
}

#' Build an n(%) count block for a categorical variable
#' @param var Name of the categorical variable to count
#' @param incl_n Whether to prepend an "n" row of non-missing counts per arm
#' @return Padded tibble with rowlbl2 and res1..res4 columns
count_block <- function(var, incl_n = FALSE) {
  # pct_lt = 1 shows sub-1% as "<1"; zero_count_display = "count_only" shows zero cells as a bare count;
  # order_count_method = "byfactor" orders category rows by the target's factor levels.
  s <- tplyr_spec(cols = "TRT01P",
                  total_groups = list(total_group("TRT01P")),
                  layers = tplyr_layers(group_count(var,
                    settings = layer_settings(
                      format_strings = list(n_counts = f_str("xxx (xxx%)", "n", "pct")),
                      order_count_method = "byfactor",
                      pct_lt = 1, zero_count_display = "count_only",
                      # Chi-square p across arms; lands on this block's first category row.
                      assoc_test = assoc_test(fn = chi_fn(var), format = P_FMT)))))
  blk <- reslist(tplyr_build(s, adslG))
  if (incl_n) {
    # adslG holds only the 3 real arms; the Total column count is their sum.
    nrow_df <- adslG |> filter(!is.na(.data[[var]])) |> count(TRT01P, .drop = FALSE)
    n_by <- setNames(nrow_df$n, as.character(nrow_df$TRT01P))
    getn <- function(a) sprintf("%3d", if (identical(a, "Total")) sum(nrow_df$n) else n_by[[a]])
    # The prepended "n" row carries the p, so move it up off the first category row.
    n_row <- tibble(rowlbl2 = "n", res1 = getn(ARMS[1]), res2 = getn(ARMS[2]),
                    res3 = getn(ARMS[3]), res4 = getn(ARMS[4]), pcol = blk$pcol[1])
    blk$pcol[1] <- ""
    blk <- bind_rows(n_row, blk)
  }
  pad_row(blk)
}

#' Stamp the characteristic label onto a block's first row, then select display columns
#' @param block Block tibble with rowlbl2, res and pcol columns
#' @param name Row label for the first row (rowlbl1)
#' @return Tibble with columns rowlbl1, rowlbl2, res1..res4, pcol
stamp <- function(block, name) {
  block$rowlbl1 <- ""; block$rowlbl1[1] <- name
  select(block, rowlbl1, rowlbl2, res1, res2, res3, res4, pcol)
}

# ---- Characteristics --------------------------------------------------------
age  <- stamp(bind_rows(desc_block("AGE"), count_block("AGEGR1")), "Age (y)")
sex  <- stamp(count_block("SEX", incl_n = TRUE), "Sex")
race <- stamp(count_block("RACE_DISPLAY", incl_n = TRUE), "Race (Origin)")
mmse <- stamp(desc_block("MMSETOT"), "MMSE")
dur  <- stamp(bind_rows(desc_block("DURDIS"), count_block("DURDSGR1")), "Duration of disease ")
educ <- stamp(desc_block("EDUCLVL"), "Years of education")
wt   <- stamp(desc_block("WEIGHTBL"), "Baseline weight(kg)")
ht   <- stamp(desc_block("HEIGHTBL"), "Baseline height(cm)")
bmi  <- stamp(bind_rows(desc_block("BMIBL"), count_block("BMIBLGR1")), "Baseline BMI")

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
