# t_14_6_05.R
# Table 14-6.05: Shifts of Laboratory Values During Treatment, Categorized Based on
#                Threshold Ranges   (Population: Safety)
# Produces: outputs/14-6.05.docx
# Source: adlbc + adlbh (ANL01FL one record per subject per analyte, adsl for arm Ns);
#   tplyr2 group_shift column-% shift-to ANRIND by baseline (Normal/High), coin::cmh_test p per analyte.
source("R/setup.R")
source("R/helpers.R")
suppressPackageStartupMessages(library(coin))

TABLE <- "14-6.05"
SOURCE <- "programs/t-14-6.05.R"
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
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
COLS <- c("Placebo N","Placebo H","Xanomeline Low Dose N","Xanomeline Low Dose H",
          "Xanomeline High Dose N","Xanomeline High Dose H")

comb <- bind_rows(read_adam("adlbc"), read_adam("adlbh")) |>
  filter(SAFFL == "Y", ANL01FL == "Y", AVISITN != 99) |>
  mutate(PARAM = recode(PARAM, !!!RECODE)) |>
  # threshold shift analyses use Normal/High only; drop anything else so the
  # computed n-row and the tplyr2 counts agree.
  filter(!is.na(PARAM), !is.na(TRTP), !is.na(BNRIND), !is.na(ANRIND),
         PARAM %in% PARAM_ORDER, ANRIND %in% c("N", "H"), BNRIND %in% c("N", "H")) |>
  mutate(PARAM = factor(PARAM, levels = PARAM_ORDER),
         TRTP  = factor(TRTP, levels = ARMS),
         ANRIND = factor(ANRIND, levels = c("N", "H")),
         BNRIND = factor(BNRIND, levels = c("N", "H")))

# shift cells: n(%) shifting to each ANRIND, column-wise per (TRTP x BNRIND) baseline group.
# group_shift with shift_denom = "column"; res columns come out TRTP x BNRIND
# (Placebo|N, Placebo|H, Low|N, ...) matching COLS.
b <- tplyr_build(tplyr_spec(cols = "TRTP",
       layers = tplyr_layers(group_shift(c(row = "ANRIND", column = "BNRIND"), by = "PARAM",
         settings = layer_settings(
           shift_denom = "column",
           format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
           order_count_method = "byfactor", zero_count_display = "count_only")))), comb)
rl <- grep("^rowlabel", names(b), value = TRUE)
rc <- grep("^res", names(b), value = TRUE)   # PARAM,ANRIND | 6
shift <- tibble(PARAM = as.character(b[[rl[1]]]), SHIFT = as.character(b[[rl[2]]]))
for (i in seq_along(rc)) shift[[COLS[i]]] <- as.character(b[[rc[i]]])
shift$SHIFT <- recode(shift$SHIFT, "N" = "Normal", "H" = "High")

# n rows: baseline-group N per (PARAM, TRTP, BNRIND) column
nrows <- comb |> count(PARAM, TRTP, BNRIND, .drop = FALSE) |>
  mutate(col = paste(TRTP, BNRIND)) |>
  select(PARAM, col, n) |> pivot_wider(names_from = col, values_from = n, values_fill = 0)
nrows[COLS] <- lapply(nrows[COLS], function(x) sprintf("%2d", x))
nrows$SHIFT <- "n"
nrows$PARAM <- as.character(nrows$PARAM)

#' Compute the per-analyte row-mean-scores CMH p-value
#' @param pv PARAM (analyte) label
#' @return Formatted p-value (3 dp), or "" when the test is not applicable
#'
#' Ordinal shift (ANRIND scored 1/2) across treatment groups, stratified by baseline
#' status BNRIND. Blank when nobody shifted to High, or when the abnormal-at-baseline
#' stratum is empty so "controlling for baseline status" is vacuous (e.g. LYMPHOCYTES).
cmh_pval <- function(pv) {
  d <- comb |> filter(PARAM == pv)
  if (all(d$ANRIND == "N")) return("")
  if (dplyr::n_distinct(d$BNRIND) < 2) return("")
  dd <- d |> transmute(ANRIND = ordered(as.character(ANRIND), c("N", "H")),
                       TRTP   = factor(as.character(TRTP), levels = ARMS),
                       BNRIND = factor(as.character(BNRIND), levels = c("N", "H")))
  tryCatch(
    num_fmt(as.numeric(pvalue(cmh_test(ANRIND ~ TRTP | BNRIND, data = dd,
                                       scores = list(ANRIND = c(1, 2))))),
            digits = 3, int_len = 1),
    error = function(e) "")
}

# assemble per PARAM: n (+p), Normal, High (drop High if all six cells zero)

#' Test whether every shift cell is zero (blank or a formatted 0)
#' @param v Character vector of formatted cells
#' @return TRUE when all cells are "0" or ""
zero_cell <- function(v) all(sub("\\s", "", v) %in% c("0", ""))
acc <- list()
for (pv in PARAM_ORDER) {
  nr <- nrows |> filter(PARAM == pv)
  if (nrow(nr) == 0 || all(as.integer(sub("\\s", "", unlist(nr[COLS]))) == 0)) next
  nr$LBL <- pv
  nr$PVAL <- cmh_pval(pv)
  sh <- shift |> filter(PARAM == pv)
  normal <- sh |> filter(SHIFT == "Normal")
  high <- sh |> filter(SHIFT == "High")
  normal$LBL <- ""
  normal$PVAL <- ""
  blk <- bind_rows(nr[, c("LBL", "SHIFT", COLS, "PVAL")],
                   normal[, c("LBL", "SHIFT", COLS, "PVAL")])
  if (nrow(high) && !zero_cell(unlist(high[COLS]))) {
    high$LBL <- ""
    high$PVAL <- ""
    blk <- bind_rows(blk, high[, c("LBL", "SHIFT", COLS, "PVAL")])
  }
  acc[[length(acc) + 1]] <- blk
}
body <- bind_rows(acc)

#' Build a section-header / spacer row spanning the full stub width
#' @param t Header text (e.g. "CHEMISTRY", "----------", or "")
#' @return A one-row tibble with blank value and p-value cells
sec <- function(t) tibble(LBL = t, SHIFT = "", !!!setNames(rep(list(""), 6), COLS), PVAL = "")
first_hem <- min(which(body$LBL %in% HEM))
final <- bind_rows(
  sec(""),
  sec("CHEMISTRY"), sec("----------"),
  body[1:(first_hem - 1), ],
  sec("HEMATOLOGY"), sec("----------"),
  body[first_hem:nrow(body), ]
) |> select(LBL, SHIFT, all_of(COLS), PVAL)
final[] <- lapply(final, function(x) {
  x[is.na(x)] <- ""
  as.character(x)
})

# render: 3-row header (arm spanner x baseline), per-arm spanner underlines
Ns <- read_adam("adsl") |> filter(ARM != "Screen Failure") |> count(TRT01P) |> deframe()

#' Format an arm spanner label with its N
#' @param short Short arm label
#' @param arm Full arm name for the N lookup
#' @return "short (N=n)"
arm_hdr <- function(short, arm) sprintf("%s (N=%s)", short, Ns[[arm]])
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    LBL = "", SHIFT = c("", "Shift", "[1]"),
    `Placebo N`              = c(arm_hdr("Placebo", "Placebo"), "Normal at", "Baseline"),
    `Placebo H`              = c(arm_hdr("Placebo", "Placebo"), "High at", "Baseline"),
    `Xanomeline Low Dose N`  = c(arm_hdr("Xan. Low", "Xanomeline Low Dose"), "Normal at", "Baseline"),
    `Xanomeline Low Dose H`  = c(arm_hdr("Xan. Low", "Xanomeline Low Dose"), "High at", "Baseline"),
    `Xanomeline High Dose N` = c(arm_hdr("Xan. High", "Xanomeline High Dose"), "Normal at", "Baseline"),
    `Xanomeline High Dose H` = c(arm_hdr("Xan. High", "Xanomeline High Dose"), "High at", "Baseline"),
    PVAL = c("", "p-\nvalue", "[2]"),
    # merge = "spanners" merges every header row except the bottom one, so the arm
    # spanners merge while the six repeated "Baseline" sub-labels each stay in their own column.
    merge = "spanners") |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "LBL", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::align(j = "PVAL", part = "body", align = "right") |>
  flextable::width(j = "LBL", width = 2.79) |>
  flextable::width(j = "SHIFT", width = 0.81) |>
  flextable::width(j = COLS, width = 0.81) |>
  flextable::width(j = "PVAL", width = 0.54)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)

#' Apply the house default plus a 1pt underline beneath each arm spanner (header row 1)
#' @param x A flextable
#' @param ... Unused
#' @return The styled flextable
sd <- function(x, ...) {
  x <- cdisc_table_default(x)
  bd <- officer::fp_border(color = "black", width = 1)
  x <- flextable::hline(x, i = 1, j = 3:4, border = bd, part = "header")
  x <- flextable::hline(x, i = 1, j = 5:6, border = bd, part = "header")
  flextable::hline(x, i = 1, j = 7:8, border = bd, part = "header")
}
old <- options(clinify_table_default = sd)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), "\n")
