# t_14_6_02.R
# Table 14-6.02: Frequency of Normal and Abnormal (Beyond Normal Range) Laboratory
#                Values During Treatment   (Population: Safety)
# Produces: outputs/14-6.02.docx
# Source: adlbc/adlbh (one on-treatment analysis record per subject per PARAM, adsl for
#   arm Ns); tplyr2 group_count subject n(%) by arm x LBNRIND category, Fisher's-exact p per PARAM.
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R")
source("R/helpers.R")

TABLE <- "14-6.02"
SOURCE <- "programs/t-14-6.02.R"
ARMS  <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

# PARAM recode (long -> short) and display order.
CHEM_RC <- c(
  "Alanine Aminotransferase (U/L)"="ALANINE AMINOTRANSFERASE","Albumin (g/L)"="ALBUMIN",
  "Alkaline Phosphatase (U/L)"="ALKALINE PHOSPHATASE","Aspartate Aminotransferase (U/L)"="ASPARTATE AMINOTRANSFERASE",
  "Bilirubin (umol/L)"="BILIRUBIN","Calcium (mmol/L)"="CALCIUM","Chloride (mmol/L)"="CHLORIDE",
  "Cholesterol (mmol/L)"="CHOLESTEROL","Creatine Kinase (U/L)"="CREATINE KINASE","Creatinine (umol/L)"="CREATININE",
  "Gamma Glutamyl Transferase (U/L)"="GAMMA GLUTAMYL TRANSFERASE","Glucose (mmol/L)"="GLUCOSE",
  "Phosphate (mmol/L)"="PHOSPHATE","Potassium (mmol/L)"="POTASSIUM","Protein (g/L)"="PROTEIN",
  "Sodium (mmol/L)"="SODIUM","Urate (umol/L)"="URATE","Blood Urea Nitrogen (mmol/L)"="UREA NITROGEN")
CHEM_ORD <- c("ALBUMIN","ALKALINE PHOSPHATASE","ALANINE AMINOTRANSFERASE","ASPARTATE AMINOTRANSFERASE",
  "BILIRUBIN","UREA NITROGEN","CALCIUM","CHOLESTEROL","CREATINE KINASE","CHLORIDE","CREATININE",
  "GAMMA GLUTAMYL TRANSFERASE","GLUCOSE","POTASSIUM","SODIUM","PHOSPHATE","PROTEIN","URATE")
HEM_RC <- c(
  "Basophils (GI/L)"="BASOPHILS","Eosinophils (GI/L)"="EOSINOPHILS",
  "Ery. Mean Corpuscular HGB Concentration (mmol/L)"="ERY. MEAN CORPUSCULAR HB CONCENTRATION",
  "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))"="ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "Ery. Mean Corpuscular Volume (fL)"="ERY. MEAN CORPUSCULAR VOLUME","Erythrocytes (TI/L)"="ERYTHROCYTES",
  "Hematocrit"="HEMATOCRIT","Hemoglobin (mmol/L)"="HEMOGLOBIN","Leukocytes (GI/L)"="LEUKOCYTES",
  "Lymphocytes (GI/L)"="LYMPHOCYTES","Monocytes (GI/L)"="MONOCYTES","Platelet (GI/L)"="PLATELET")
HEM_ORD <- c("BASOPHILS","EOSINOPHILS","HEMATOCRIT","HEMOGLOBIN","LYMPHOCYTES",
  "ERY. MEAN CORPUSCULAR HEMOGLOBIN","ERY. MEAN CORPUSCULAR HB CONCENTRATION","ERY. MEAN CORPUSCULAR VOLUME",
  "MONOCYTES","PLATELET","ERYTHROCYTES","LEUKOCYTES")

#' Read and prep a lab dataset for the count table
#' @param name ADaM dataset name (without extension)
#' @param rc Named vector remapping PARAM to short labels
#' @param ord Character vector giving the PARAM display order
#' @return On-treatment analysis records with PARAM/TRTP/RIND factored to a full 3x3 per PARAM
prep <- function(name, rc, ord) {
  d <- read_adam(name) |> filter(SAFFL == "Y", ANL01FL == "Y", AVISITN != 99)
  d$PARAM <- recode(d$PARAM, !!!rc)
  d$RIND  <- recode(d$LBNRIND, "LOW" = "L", "NORMAL" = "N", "HIGH" = "H")
  d |> filter(!is.na(PARAM), !is.na(TRTP), RIND %in% c("L", "N", "H")) |>
    mutate(PARAM = factor(PARAM, ord),
           TRTP  = factor(TRTP, ARMS),
           RIND  = factor(RIND, c("L", "N", "H")))
}

#' Build a Fisher's-exact p-value function for a per-PARAM `assoc_test()`
#' @param na_rule When TRUE, blank the p-value where every Low and High cell is 0
#' @return A function of the by-group data subset returning a scalar p-value (or NA)
#'
#' The test is over the 3x3 (arm x category) table. na_rule returns NA (rendered blank)
#' when all subjects are Normal (every Low and High cell 0), where Fisher's test is undefined.
fisher_p <- function(na_rule = FALSE) function(.data) {
  m <- matrix(as.integer(table(factor(.data$TRTP, ARMS), factor(.data$RIND, c("L", "N", "H")))),
              nrow = 3)                                       # rows = arm, cols = L/N/H
  if (na_rule && all(m[, 1] == 0) && all(m[, 3] == 0)) return(NA_real_)
  suppressWarnings(fisher.test(m)$p.value)
}

#' Build the 9-column n(%) block plus per-PARAM Fisher p-value
#' @param d Prepped lab records from `prep()`
#' @param na_rule Passed to `fisher_p()`; blanks the p where every Low and High cell is 0
#' @return A tibble with a STUB label column, res1..res9 n(%) cells, and a formatted PVAL
#'
#' Percent denominator is the per-(PARAM, arm) total, so zeros display as a bare " 0".
#' PARAM is carried as the layer's `by` so assoc_test() keys the Fisher test per analyte;
#' the constant target ("x") gives one output row per PARAM, and as_display() surfaces the
#' p-value as a trailing pval1 column (NA/error -> blank) which we relabel PVAL.
count_wide <- function(d, na_rule = FALSE) {
  d$STAT <- factor("x")
  b <- tplyr_build(tplyr_spec(cols = c("TRTP", "RIND"),
        layers = tplyr_layers(group_count("STAT", by = "PARAM",
          settings = layer_settings(
            format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
            denoms_by = c("PARAM", "TRTP"),
            order_count_method = "byfactor",
            zero_count_display = "count_only",
            assoc_test = assoc_test(fn = fisher_p(na_rule),
                                    format = f_str("x.xxx", "p"),
                                    label = "p-val"))))), d)
  # as_display() returns rowlabel1 (PARAM) + rowlabel2 (the constant "x") + the 9 res n(%)
  # cells (arm(L,N,H) x 3 arms) + pval1, already ordered by the build; drop the dummy target
  # column and relabel PARAM -> STUB, pval1 -> PVAL.
  as_display(b) |> select(-rowlabel2) |> rename(STUB = rowlabel1, PVAL = pval1)
}

chem_d <- prep("adlbc", CHEM_RC, CHEM_ORD)
hem_d  <- prep("adlbh", HEM_RC, HEM_ORD)
chem   <- count_wide(chem_d, na_rule = FALSE)
heme   <- count_wide(hem_d,  na_rule = TRUE)

COLS <- c(paste0("res", 1:9), "PVAL")

#' Build a section/blank row (STUB text, all value cells blank)
#' @param stub Text for the STUB column
#' @return A one-row tibble with `stub` in STUB and "" in all value columns
row9 <- function(stub = "") tibble(STUB = stub, !!!setNames(rep(list(""), 10), COLS))

final <- bind_rows(
  row9(""), row9("CHEMISTRY"), row9("----------"), chem,
  row9(""), row9("HEMATOLOGY"), row9("----------"), heme)

# render: 2-level header (arm spanner x Low/Normal/High) + Fisher p column.
# Arm Ns are the safety population per arm (adsl, excluding Screen Failure).
Ns  <- read_adam("adsl") |> filter(ARM != "Screen Failure") |> count(TRT01P) |> deframe()

#' Format an arm spanner label with its N
#' @param short Short arm label
#' @param arm Full arm name for the N lookup
#' @return "short (N=n)"
sp  <- function(short, arm) sprintf("%s (N=%s)", short, Ns[[arm]])
PBO <- sp("Placebo", "Placebo")
LOW <- sp("Xan. Low", "Xanomeline Low Dose")
HIGH <- sp("Xan. High", "Xanomeline High Dose")
ct <- clintable(final, use_labels = FALSE, coerce_character = TRUE) |>
  clin_column_headers(
    STUB = "",
    res1 = c(PBO, "Low"),  res2 = c(PBO, "Normal"),  res3 = c(PBO, "High"),
    res4 = c(LOW, "Low"),  res5 = c(LOW, "Normal"),  res6 = c(LOW, "High"),
    res7 = c(HIGH, "Low"), res8 = c(HIGH, "Normal"), res9 = c(HIGH, "High"),
    PVAL = c("", "p-val\n[1]")) |>
  # solid rule beneath each arm spanner, auto-derived to span exactly its 3 value
  # columns (skips the blank STUB/PVAL header cells); survives the styler's border_remove().
  clin_spanner_rule() |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "STUB", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "STUB", width = 1.8) |>
  flextable::width(j = paste0("res", 1:9), width = 0.74) |>
  flextable::width(j = "PVAL", width = 0.52)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)

# The per-spanner rules are drawn by clin_spanner_rule() above; the full-width bottom rule
# under the labels row comes from the house cdisc_table_default() styler.
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
cat("rows:", nrow(final), "\n")
