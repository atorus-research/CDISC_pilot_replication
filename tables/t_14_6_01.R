# t_14_6_01.R
# Table 14-6.01: Summary Statistics for Continuous Laboratory Values   (Population: Safety)
# Produces: outputs/14-6.01.docx
# Source: adlbc/adlbh (adsl for arm Ns); tplyr2 group_desc mean(SD) of AVAL & CHG by
#   visit, assembled here into transposed per-arm N | Mean (SD) | Change columns.
source("R/setup.R")
source("R/helpers.R")

TABLE <- "14-6.01"
SOURCE <- "programs/t-14-6-01.R"
ARM <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
VISN <- c(0, 2, 4, 6, 8, 12, 16, 20, 24, 26, 99)
VISLAB <- c("Bsln", "Wk 2", "Wk 4", "Wk 6", "Wk 8", "Wk 12", "Wk 16", "Wk 20",
            "Wk 24", "Wk 26", "End[1]")
COLS <- c("P_N","P_M","P_C","L_N","L_M","L_C","H_N","H_M","H_C")

CHEM_RC <- c(
  "Alanine Aminotransferase (U/L)"="ALANINE AMINOTRANSFERASE","Albumin (g/L)"="ALBUMIN",
  "Alkaline Phosphatase (U/L)"="ALKALINE PHOSPHATASE","Aspartate Aminotransferase (U/L)"="ASPARTATE AMINOTRANSFERASE",
  "Bilirubin (umol/L)"="BILIRUBIN","Calcium (mmol/L)"="CALCIUM","Chloride (mmol/L)"="CHLORIDE",
  "Cholesterol (mmol/L)"="CHOLESTEROL","Creatine Kinase (U/L)"="CREATINE KINASE","Creatinine (umol/L)"="CREATININE",
  "Gamma Glutamyl Transferase (U/L)"="GAMMA GLUTAMYL TRANSFERASE","Glucose (mmol/L)"="GLUCOSE",
  "Phosphate (mmol/L)"="PHOSPHATE","Potassium (mmol/L)"="POTASSIUM","Protein (g/L)"="PROTEIN",
  "Sodium (mmol/L)"="SODIUM","Urate (umol/L)"="URATE","Blood Urea Nitrogen (mmol/L)"="UREA NITROGEN")
HEM_RC <- c(
  "Basophils (GI/L)"="BASOPHILS","Eosinophils (GI/L)"="EOSINOPHILS",
  "Ery. Mean Corpuscular HGB Concentration (mmol/L)"="ERY. MEAN CORPUSCULAR HB CONCENTRATION",
  "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))"="ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "Ery. Mean Corpuscular Volume (fL)"="ERY. MEAN CORPUSCULAR VOLUME","Erythrocytes (TI/L)"="ERYTHROCYTES",
  "Hematocrit"="HEMATOCRIT","Hemoglobin (mmol/L)"="HEMOGLOBIN","Leukocytes (GI/L)"="LEUKOCYTES",
  "Lymphocytes (GI/L)"="LYMPHOCYTES","Monocytes (GI/L)"="MONOCYTES","Platelet (GI/L)"="PLATELET")
HEM_DROP <- c("Anisocytes","Poikilocytes","Microcytes","Macrocytes","Polychromasia")

#' Read and filter a lab ADaM to the safety population and analysis visits
#' @param name ADaM dataset name (without extension)
#' @param recode_map Named vector remapping PARAM to short display labels
#' @param drop Character vector of PARAM values to exclude
#' @return A filtered, recoded lab data frame
read_lab <- function(name, recode_map, drop = character()) {
  read_adam(name) |>
    filter(SAFFL == "Y", (AVISITN != 99 | (AVISITN == 99 & AENTMTFL == "Y")),
           AVISIT != "UNSCHEDULED", !(PARAM %in% drop), AVISITN %in% VISN) |>
    mutate(PARAM = recode(PARAM, !!!recode_map))
}

#' Extract PARAM/visit labels and the three arm stat columns from a tplyr2 build
#' @param b A tplyr2 build data frame
#' @return A tibble with PARAM, AVISITN and the c1/c2/c3 stat cells
extract <- function(b) {
  rl <- grep("^rowlabel", names(b), value = TRUE)
  rc <- grep("^res", names(b), value = TRUE)
  tibble(PARAM = as.character(b[[rl[1]]]), AVISITN = as.integer(as.character(b[[rl[2]]])),
         c1 = as.character(b[[rc[1]]]), c2 = as.character(b[[rc[2]]]), c3 = as.character(b[[rc[3]]]))
}

#' Blank out a stat cell that contains no digits (e.g. "(      )")
#' @param x Character vector of formatted stat cells
#' @return `x` with non-numeric entries replaced by ""
blank_if_empty <- function(x) ifelse(grepl("[0-9]", x), x, "")

#' Build the wide per-arm stat block for one lab section
#' @param dat Prepared lab records for the section
#' @param params Character vector of PARAM values in display order
#' @return A wide tibble of N / Mean (SD) / Change cells per arm and PARAM x visit
section_wide <- function(dat, params) {
  dat <- dat |> mutate(TRTP = factor(TRTP, levels = ARM),
                       PARAM = factor(PARAM, levels = params),
                       AVISITN = factor(AVISITN, levels = VISN))
  #' Summarise one variable to mean (SD) by PARAM and visit via tplyr2
  #' @param var Variable to summarise (AVAL or CHG)
  #' @return A tibble of formatted mean (SD) cells per arm
  msd <- function(var) extract(tplyr_build(tplyr_spec(cols = "TRTP", layers = tplyr_layers(
    group_desc(var, by = c("PARAM", "AVISITN"), settings = layer_settings(
      format_strings = list(mean_sd = f_str("xxx.x (xxx.xx)", "mean", "sd")))))), dat))
  MA <- msd("AVAL")
  MC <- msd("CHG")
  # N is the record count at each visit (reference definition), not tplyr2's non-missing-AVAL count
  Ncnt <- dat |> count(PARAM, AVISITN, TRTP, .drop = FALSE) |>
    mutate(PARAM = as.character(PARAM), AVISITN = as.integer(as.character(AVISITN)),
           TRTP = c(P = "P", L = "L", H = "H")[match(TRTP, ARM)]) |>
    pivot_wider(names_from = TRTP, values_from = n, values_fill = 0, names_glue = "{TRTP}_N")
  w <- MA |> left_join(MC, by = c("PARAM", "AVISITN"), suffix = c("_A", "_C")) |>
    left_join(Ncnt, by = c("PARAM", "AVISITN"))
  arms <- c("P", "L", "H")
  for (i in seq_along(arms)) {
    a <- arms[i]
    ci <- paste0("c", i)
    n <- w[[paste0(a, "_N")]]
    w[[paste0(a, "_N")]] <- ifelse(n == 0, "", sprintf("%2d", n))
    w[[paste0(a, "_M")]] <- ifelse(n == 0, "", blank_if_empty(w[[paste0(ci, "_A")]]))
    w[[paste0(a, "_C")]] <- ifelse(n == 0, "", blank_if_empty(w[[paste0(ci, "_C")]]))
  }
  w
}

#' Build a full-width spacer/label row with blank value cells
#' @param lbl Text for the VISIT (stub) column
#' @return A one-row tibble with `lbl` in VISIT and "" in all value columns
blank_row <- function(lbl = "") tibble(VISIT = lbl, !!!setNames(rep(list(""), 9), COLS))

#' Assemble a section into stacked header / value / spacer rows
#' @param section_label Section heading text (e.g. "CHEMISTRY")
#' @param wide Wide stat block from `section_wide()`
#' @param params Character vector of PARAM values in display order
#' @return A tibble of body rows for the section
assemble <- function(section_label, wide, params) {
  out <- list(blank_row(section_label))
  for (p in params) {
    pv <- wide |> filter(PARAM == p) |> arrange(match(AVISITN, VISN))
    if (nrow(pv) == 0) next
    vr <- tibble(VISIT = paste0("  ", VISLAB[match(pv$AVISITN, VISN)]))
    for (cc in COLS) vr[[cc]] <- pv[[cc]]
    out[[length(out) + 1]] <- blank_row(p)
    out[[length(out) + 1]] <- vr
    out[[length(out) + 1]] <- blank_row("")
  }
  bind_rows(out)
}

chem <- read_lab("adlbc", CHEM_RC)
hema <- read_lab("adlbh", HEM_RC, drop = HEM_DROP)
chem_p <- sort(unique(chem$PARAM))
hema_p <- sort(unique(hema$PARAM))

final <- bind_rows(
  assemble("CHEMISTRY",  section_wide(chem, chem_p), chem_p),
  assemble("HEMATOLOGY", section_wide(hema, hema_p), hema_p)
) |> select(VISIT, all_of(COLS))
final[] <- lapply(final, function(x) {
  x[is.na(x)] <- ""
  as.character(x)
})
# drop trailing all-blank spacer rows; keep interior spacers
final <- final[seq_len(max(which(rowSums(final != "") > 0))), , drop = FALSE]

# rows whose label spans the full width (section + param headers): unindented, no values
merge_rows <- which(final$VISIT != "" & !startsWith(final$VISIT, "  "))

# render: spanner + multi-line label header; dashed spanner underline
Ns <- read_adam("adsl") |> filter(ARM != "Screen Failure") |> count(TRT01P) |> deframe()

#' Return an arm spanner label unchanged
#' @param arm Arm label
#' @return `arm`
spn <- function(arm) arm
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    VISIT = c("", "Visit"),
    P_N = c("Placebo", "N"),        P_M = c("Placebo", "Mean (SD)"),        P_C = c("Placebo", "Change\nfrom Bsln\nMean (SD)"),
    L_N = c("Xanomeline Low", "N"), L_M = c("Xanomeline Low", "Mean (SD)"), L_C = c("Xanomeline Low", "Change\nfrom Bsln\nMean (SD)"),
    H_N = c("Xanomeline High", "N"),H_M = c("Xanomeline High", "Mean (SD)"),H_C = c("Xanomeline High", "Change\nfrom Bsln\nMean (SD)")) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "VISIT", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "VISIT", width = 0.82) |>
  flextable::width(j = c("P_N", "L_N", "H_N"), width = 0.24) |>
  flextable::width(j = c("P_M", "P_C", "L_M", "L_C", "H_M", "H_C"), width = 1.24)
for (r in merge_rows) ct <- flextable::merge_at(ct, i = r, j = 1:10, part = "body")
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)

#' Apply the house default plus a dashed spanner underline and tight cell padding
#' @param x A flextable
#' @param ... Unused
#' @return The styled flextable
#'
#' Padding is tightened so the 14-char "xxx.x (xxx.xx)" cells fit the 1.24" columns
#' without wrapping; the dashed rule underlines each arm spanner (header row 1).
sd <- function(x, ...) {
  x <- cdisc_table_default(x)
  x <- flextable::padding(x, padding.left = 1, padding.right = 1, part = "all")
  flextable::hline(x, i = 1, j = 2:10,
    border = officer::fp_border(color = "black", width = 1, style = "dashed"), part = "header")
}
old <- options(clinify_table_default = sd)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), " merge_rows:", length(merge_rows), "\n")
