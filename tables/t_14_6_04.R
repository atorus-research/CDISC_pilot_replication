# t_14_6_04.R
# Table 14-6.04: Shifts of Laboratory Values During Treatment, Categorized Based on
#                Threshold Ranges, by Visit   (Population: Safety)
# Produces: outputs/14-6.04.docx
# Source: adlbc + adlbh (adsl for arm Ns); tplyr2 group_shift column-% shift-to ANRIND
#   by baseline BNRIND (Normal/High), per PARAM x visit.
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R")
source("R/helpers.R")

TABLE <- "14-6.04"
SOURCE <- "programs/t-14-6.04.R"
CHEM <- c("ALANINE AMINOTRANSFERASE","ALBUMIN","ALKALINE PHOSPHATASE","ASPARTATE AMINOTRANSFERASE",
  "BILIRUBIN","CALCIUM","CHLORIDE","CHOLESTEROL","CREATINE KINASE","CREATININE",
  "GAMMA GLUTAMYL TRANSFERASE","GLUCOSE","PHOSPHATE","POTASSIUM","PROTEIN","SODIUM","URATE","UREA NITROGEN")
HEM <- c("BASOPHILS","EOSINOPHILS","ERY. MEAN CORPUSCULAR HB CONCENTRATION","ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "ERY. MEAN CORPUSCULAR VOLUME","ERYTHROCYTES","HEMATOCRIT","HEMOGLOBIN","LEUKOCYTES","LYMPHOCYTES",
  "MONOCYTES","PLATELET")
PARAM_ORDER <- c(CHEM, HEM)
RECODE <- c(
  "Alanine Aminotransferase (U/L)"="ALANINE AMINOTRANSFERASE","Albumin (g/L)"="ALBUMIN",
  "Alkaline Phosphatase (U/L)"="ALKALINE PHOSPHATASE","Aspartate Aminotransferase (U/L)"="ASPARTATE AMINOTRANSFERASE",
  "Bilirubin (umol/L)"="BILIRUBIN","Calcium (mmol/L)"="CALCIUM","Chloride (mmol/L)"="CHLORIDE",
  "Cholesterol (mmol/L)"="CHOLESTEROL","Creatine Kinase (U/L)"="CREATINE KINASE","Creatinine (umol/L)"="CREATININE",
  "Gamma Glutamyl Transferase (U/L)"="GAMMA GLUTAMYL TRANSFERASE","Glucose (mmol/L)"="GLUCOSE",
  "Phosphate (mmol/L)"="PHOSPHATE","Potassium (mmol/L)"="POTASSIUM","Protein (g/L)"="PROTEIN",
  "Sodium (mmol/L)"="SODIUM","Urate (umol/L)"="URATE","Blood Urea Nitrogen (mmol/L)"="UREA NITROGEN",
  "Basophils (GI/L)"="BASOPHILS","Eosinophils (GI/L)"="EOSINOPHILS",
  "Ery. Mean Corpuscular HGB Concentration (mmol/L)"="ERY. MEAN CORPUSCULAR HB CONCENTRATION",
  "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))"="ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "Ery. Mean Corpuscular Volume (fL)"="ERY. MEAN CORPUSCULAR VOLUME","Erythrocytes (TI/L)"="ERYTHROCYTES",
  "Hematocrit"="HEMATOCRIT","Hemoglobin (mmol/L)"="HEMOGLOBIN","Leukocytes (GI/L)"="LEUKOCYTES",
  "Lymphocytes (GI/L)"="LYMPHOCYTES","Monocytes (GI/L)"="MONOCYTES","Platelet (GI/L)"="PLATELET")
VIS <- c("WEEK 2","WEEK 4","WEEK 6","WEEK 8","WEEK 12","WEEK 16","WEEK 20","WEEK 24","WEEK 26")
COLS <- c("Placebo N","Placebo H","Xanomeline Low Dose N","Xanomeline Low Dose H",
          "Xanomeline High Dose N","Xanomeline High Dose H")

comb <- bind_rows(read_adam("adlbc"), read_adam("adlbh")) |>
  filter(SAFFL == "Y", AVISITN != 99) |>
  mutate(PARAM = recode(PARAM, !!!RECODE)) |>
  filter(!is.na(VISIT), !is.na(TRTP), !is.na(PARAM), PARAM %in% PARAM_ORDER, VISIT %in% VIS,
         # threshold shift uses Normal/High only; drop Low so the n-row count and
         # the tplyr2 shift count agree (reference factors to N/H then drops the rest)
         ANRIND %in% c("N", "H"), BNRIND %in% c("N", "H")) |>
  mutate(PARAM = factor(PARAM, levels = PARAM_ORDER),
         VISIT = factor(VISIT, levels = VIS),
         TRTP  = factor(TRTP, levels = c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")),
         ANRIND = factor(ANRIND, levels = c("N", "H")),
         BNRIND = factor(BNRIND, levels = c("N", "H")))

# shift cells: n(%) shifting to each ANRIND, column-wise per (TRTP x BNRIND) baseline group.
# group_shift with shift_denom = "column"; res columns come out TRTP x BNRIND
# (Placebo|N, Placebo|H, Low|N, ...) in the same order as COLS.
# as_display() returns the ordered, display-ready frame (rowlabel*/res* as
# character), dropping the internal ord*/row_id cols.
b <- as_display(tplyr_build(tplyr_spec(cols = "TRTP",
       layers = tplyr_layers(group_shift(c(row = "ANRIND", column = "BNRIND"), by = c("PARAM", "VISIT"),
         settings = layer_settings(
           shift_denom = "column",
           format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
           order_count_method = "byfactor", zero_count_display = "count_only",
           # the layer's own per-baseline-group denominator, as the reference's 2-char "n" row
           denom_row = TRUE, denom_row_label = "n",
           denom_row_format = f_str("xx", "n"))))), comb))
rl <- grep("^rowlabel", names(b), value = TRUE)
rc <- grep("^res", names(b), value = TRUE)   # PARAM,VISIT,ANRIND | 6
shift <- tibble(PARAM = b[[rl[1]]], VISIT = b[[rl[2]]], SHIFT = b[[rl[3]]])
for (i in seq_along(rc)) shift[[COLS[i]]] <- b[[rc[i]]]
shift$SHIFT <- recode(shift$SHIFT, "N" = "Normal", "H" = "High")

# n rows: emitted by the shift layer itself (denom_row), so the displayed denominator is
# the one the percentages were computed against.
nrows <- shift |> filter(SHIFT == "n")

# assemble per PARAM x VISIT: n, Normal, High (drop High if all zero)

#' Test whether every shift cell is zero (blank or a formatted 0)
#' @param v Character vector of formatted cells
#' @return TRUE when all cells are " 0" / "0" / "  0"
zero_cell <- function(v) all(v %in% c(" 0", "0", "  0"))
acc <- list()
for (pv in PARAM_ORDER) {
  for (vs in VIS) {
    nr <- nrows |> filter(PARAM == pv, VISIT == vs)
    if (nrow(nr) == 0 || all(as.integer(sub("\\s", "", unlist(nr[COLS]))) == 0)) next
    sh <- shift |> filter(PARAM == pv, VISIT == vs)
    normal <- sh |> filter(SHIFT == "Normal")
    high <- sh |> filter(SHIFT == "High")
    blk <- bind_rows(nr[, c("PARAM", "VISIT", "SHIFT", COLS)],
                     normal[, c("PARAM", "VISIT", "SHIFT", COLS)])
    if (nrow(high) && !zero_cell(unlist(high[COLS]))) blk <- bind_rows(blk, high[, c("PARAM", "VISIT", "SHIFT", COLS)])
    blk$WEEK <- ""
    blk$WEEK[1] <- sub("WEEK ", "", vs)
    blk$LBL <- ""
    if (vs == "WEEK 2") blk$LBL[1] <- pv
    acc[[length(acc) + 1]] <- blk
  }
}
body <- bind_rows(acc)

#' Build a two-line section-header block spanning the full stub width
#' @param t1 First header line (e.g. "CHEMISTRY")
#' @param t2 Second header line (e.g. "----------")
#' @return A two-row tibble with blank value cells
sec <- function(t1, t2) tibble(LBL = c(t1, t2), WEEK = "", SHIFT = "",
  !!!setNames(rep(list(""), 6), COLS))
first_hem <- min(which(body$LBL %in% HEM))
final <- bind_rows(
  tibble(LBL = "", WEEK = "", SHIFT = "", !!!setNames(rep(list(""), 6), COLS)),
  sec("CHEMISTRY", "----------"),
  body[1:(first_hem - 1), ],
  sec("HEMATOLOGY", "----------"),
  body[first_hem:nrow(body), ]
) |> select(LBL, WEEK, SHIFT, all_of(COLS))

# render: 3-row header (arm spanner x baseline), spanner underlines
Ns <- read_adam("adsl") |> filter(ARM != "Screen Failure") |> count(TRT01P) |> deframe()

#' Format an arm spanner label with its N
#' @param short Short arm label
#' @param arm Full arm name for the N lookup
#' @return "short (N=n)"
arm_hdr <- function(short, arm) sprintf("%s (N=%s)", short, Ns[[arm]])
ct <- clintable(final, use_labels = FALSE, coerce_character = TRUE) |>
  clin_column_headers(
    LBL = "", WEEK = c("", "", "Week"), SHIFT = c("", "", "Shift to"),
    `Placebo N`              = c(arm_hdr("Placebo", "Placebo"), "Normal at", "Baseline"),
    `Placebo H`              = c(arm_hdr("Placebo", "Placebo"), "High at", "Baseline"),
    `Xanomeline Low Dose N`  = c(arm_hdr("Xan. Low", "Xanomeline Low Dose"), "Normal at", "Baseline"),
    `Xanomeline Low Dose H`  = c(arm_hdr("Xan. Low", "Xanomeline Low Dose"), "High at", "Baseline"),
    `Xanomeline High Dose N` = c(arm_hdr("Xan. High", "Xanomeline High Dose"), "Normal at", "Baseline"),
    `Xanomeline High Dose H` = c(arm_hdr("Xan. High", "Xanomeline High Dose"), "High at", "Baseline"),
    # merge = "spanners" merges every header row except the bottom one. The bottom row
    # here is six columns all labelled "Baseline", and the default merges identical
    # adjacent cells on every row - which collapsed all six into one cell rendering a
    # single centred "Baseline", so five of the six labels were absent from the output.
    # Same reason t_14_6_05.R passes this argument.
    merge = "spanners") |>
  clin_spanner_rule() |>   # solid rule under each arm spanner (cols 4:5, 6:7, 8:9)
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "LBL", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "LBL", width = 2.436) |>
  flextable::width(j = "WEEK", width = 0.504) |>
  flextable::width(j = "SHIFT", width = 0.588) |>
  flextable::width(j = COLS, width = 0.84)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)

# The arm-spanner underline is drawn by clin_spanner_rule() above (survives the
# house styler's border_remove()), so the default clinify_table_default applies.
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
cat("rows:", nrow(final), "\n")
