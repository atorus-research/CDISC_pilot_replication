# t_14_6_03.R
# Table 14-6.03: Frequency of Normal and Abnormal (Clinically Significant Change from
#                Previous Visit) Laboratory Values During Treatment   (Population: Safety)
# Produces: outputs/14-6.03.docx
# Source: adlbcpv/adlbhpv (one on-treatment analysis record per subject per PARAM, adsl for
#   arm Ns); tplyr2 group_count subject n(%) by arm x ANRIND category, Fisher's-exact p per PARAM.
source("R/setup.R")
source("R/helpers.R")

TABLE <- "14-6.03"
SOURCE <- "programs/t-14-6.03.R"
ARMS  <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

# PARAM recode (long -> short) and display order.
# The HGB-concentration key is truncated ("...normal rang") exactly as stored.
CHEM_RC <- c(
  "Albumin (g/L) change from previous visit, relative to normal range"="ALBUMIN",
  "Alkaline Phosphatase (U/L) change from previous visit, relative to normal range"="ALKALINE PHOSPHATASE",
  "Alanine Aminotransferase (U/L) change from previous visit, relative to normal range"="ALANINE AMINOTRANSFERASE",
  "Aspartate Aminotransferase (U/L) change from previous visit, relative to normal range"="ASPARTATE AMINOTRANSFERASE",
  "Bilirubin (umol/L) change from previous visit, relative to normal range"="BILIRUBIN",
  "Blood Urea Nitrogen (mmol/L) change from previous visit, relative to normal range"="UREA NITROGEN",
  "Calcium (mmol/L) change from previous visit, relative to normal range"="CALCIUM",
  "Cholesterol (mmol/L) change from previous visit, relative to normal range"="CHOLESTEROL",
  "Creatine Kinase (U/L) change from previous visit, relative to normal range"="CREATINE KINASE",
  "Chloride (mmol/L) change from previous visit, relative to normal range"="CHLORIDE",
  "Creatinine (umol/L) change from previous visit, relative to normal range"="CREATININE",
  "Gamma Glutamyl Transferase (U/L) change from previous visit, relative to normal range"="GAMMA GLUTAMYL TRANSFERASE",
  "Glucose (mmol/L) change from previous visit, relative to normal range"="GLUCOSE",
  "Potassium (mmol/L) change from previous visit, relative to normal range"="POTASSIUM",
  "Sodium (mmol/L) change from previous visit, relative to normal range"="SODIUM",
  "Phosphate (mmol/L) change from previous visit, relative to normal range"="PHOSPHATE",
  "Protein (g/L) change from previous visit, relative to normal range"="PROTEIN",
  "Urate (umol/L) change from previous visit, relative to normal range"="URATE")
CHEM_ORD <- c("ALBUMIN","ALKALINE PHOSPHATASE","ALANINE AMINOTRANSFERASE","ASPARTATE AMINOTRANSFERASE",
  "BILIRUBIN","UREA NITROGEN","CALCIUM","CHOLESTEROL","CREATINE KINASE","CHLORIDE","CREATININE",
  "GAMMA GLUTAMYL TRANSFERASE","GLUCOSE","POTASSIUM","SODIUM","PHOSPHATE","PROTEIN","URATE")
HEM_RC <- c(
  "Basophils (GI/L) change from previous visit, relative to normal range"="BASOPHILS",
  "Eosinophils (GI/L) change from previous visit, relative to normal range"="EOSINOPHILS",
  "Hematocrit change from previous visit, relative to normal range"="HEMATOCRIT",
  "Hemoglobin (mmol/L) change from previous visit, relative to normal range"="HEMOGLOBIN",
  "Lymphocytes (GI/L) change from previous visit, relative to normal range"="LYMPHOCYTES",
  "Ery. Mean Corpuscular Hemoglobin (fmol(Fe)) change from previous visit, relative to normal range"="ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "Ery. Mean Corpuscular HGB Concentration (mmol/L) change from previous visit, relative to normal rang"="ERY. MEAN CORPUSCULAR HB CONCENTRATION",
  "Ery. Mean Corpuscular Volume (fL) change from previous visit, relative to normal range"="ERY. MEAN CORPUSCULAR VOLUME",
  "Monocytes (GI/L) change from previous visit, relative to normal range"="MONOCYTES",
  "Platelet (GI/L) change from previous visit, relative to normal range"="PLATELET",
  "Erythrocytes (TI/L) change from previous visit, relative to normal range"="ERYTHROCYTES",
  "Leukocytes (GI/L) change from previous visit, relative to normal range"="LEUKOCYTES")
HEM_ORD <- c("BASOPHILS","EOSINOPHILS","HEMATOCRIT","HEMOGLOBIN","LYMPHOCYTES",
  "ERY. MEAN CORPUSCULAR HEMOGLOBIN","ERY. MEAN CORPUSCULAR HB CONCENTRATION","ERY. MEAN CORPUSCULAR VOLUME",
  "MONOCYTES","PLATELET","ERYTHROCYTES","LEUKOCYTES")

#' Read and prep a previous-visit lab dataset for the count table
#' @param name ADaM dataset name (without extension)
#' @param rc Named vector remapping PARAM to short labels
#' @param ord Character vector giving the PARAM display order
#' @return On-treatment analysis records with PARAM/TRTP/RIND factored to a full 3x3 per PARAM
prep <- function(name, rc, ord) {
  d <- read_adam(name) |> filter(SAFFL == "Y", ANL01FL == "Y", AVISITN != 99)
  d$PARAM <- recode(d$PARAM, !!!rc)
  d$RIND  <- as.character(d$ANRIND)
  d |> filter(!is.na(PARAM), !is.na(TRTP), RIND %in% c("L", "N", "H")) |>
    mutate(PARAM = factor(PARAM, ord),
           TRTP  = factor(TRTP, ARMS),
           RIND  = factor(RIND, c("L", "N", "H")))
}

#' Compute the Fisher's-exact p-value over the 3x3 (arm x category) table for one PARAM
#' @param .data Raw source records for a single PARAM (all TRTP x RIND levels)
#' @return Scalar Fisher's-exact p-value
#'
#' Used as the count layer's assoc_test function; a value is produced for every PARAM
#' (no all-Normal guard). Rows of the matrix are arms, columns are the L/N/H categories.
fisher_p <- function(.data) {
  cnt <- .data |> count(TRTP, RIND, .drop = FALSE) |> arrange(TRTP, RIND)
  m <- matrix(cnt$n, nrow = 3, ncol = 3, byrow = TRUE)   # rows = arm, cols = L/N/H
  suppressWarnings(fisher.test(m)$p.value)
}

#' Build the 9-column n(%) block plus the per-PARAM Fisher p-value
#' @param d Prepped lab records from `prep()`
#' @return A tibble with a STUB label column, res1..res9 n(%) cells and a PVAL column
#'
#' Percent denominator is the per-(PARAM, arm) total. PARAM is carried as the layer's
#' `by` variable over a constant target, so the assoc_test runs `fisher_p()` once per
#' PARAM and lands the formatted p-value as pval1 on that PARAM's single output row.
#' f_str("x.xxx", "p") matches the previous num_fmt(size = 5) 3-dp fixed-width string.
count_wide <- function(d) {
  b <- tplyr_build(tplyr_spec(cols = c("TRTP", "RIND"),
        layers = tplyr_layers(group_count("DUMMY", by = "PARAM",
          settings = layer_settings(
            format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
            denoms_by = c("PARAM", "TRTP"),
            order_count_method = "byfactor",
            zero_count_display = "count_only",
            assoc_test = assoc_test(fn = fisher_p,
                                    format = f_str("x.xxx", "p"),
                                    label = "p-value [1]"))))),
        mutate(d, DUMMY = factor("x")))
  # as_display() returns the ordered display frame (rowlabel1 = PARAM, the constant
  # rowlabel2 target, res1..res9 and pval1); the build is already in ord_layer_1 order.
  out <- as_display(b)
  out$rowlabel2 <- NULL
  names(out)[names(out) == "rowlabel1"] <- "STUB"
  names(out)[names(out) == "pval1"]     <- "PVAL"
  out$PVAL <- ifelse(is.na(out$PVAL), "", out$PVAL)
  out
}

chem_d <- prep("adlbcpv", CHEM_RC, CHEM_ORD)
hem_d  <- prep("adlbhpv", HEM_RC, HEM_ORD)
chem   <- count_wide(chem_d)
heme   <- count_wide(hem_d)

COLS <- c(paste0("res", 1:9), "PVAL")

#' Build a section/blank row (STUB text, all value cells blank)
#' @param stub Text for the STUB column
#' @return A one-row tibble with `stub` in STUB and "" in all value columns
row9 <- function(stub = "") tibble(STUB = stub, !!!setNames(rep(list(""), 10), COLS))

# single blank row separates the CHEMISTRY and HEMATOLOGY sections
final <- bind_rows(
  row9(""), row9("CHEMISTRY"), row9("----------"), chem,
  row9(""),
  row9("HEMATOLOGY"), row9("----------"), heme)

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
  clin_spanner_rule() |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "STUB", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "STUB", width = 1.8) |>
  flextable::width(j = paste0("res", 1:9), width = 0.74) |>
  flextable::width(j = "PVAL", width = 0.52)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)

# clin_spanner_rule() (in the clintable pipe above) draws the solid rule beneath each
# arm spanner, spanning exactly the value columns it covers. It is applied after the
# house-style table styler, so it survives that styler's border_remove(); the full-width
# bottom rule under the labels row comes from cdisc_table_default().
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
cat("rows:", nrow(final), "\n")
