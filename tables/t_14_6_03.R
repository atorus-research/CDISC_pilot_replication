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

#' Build the 9-column n(%) block (PARAM rows x arm-category columns)
#' @param d Prepped lab records from `prep()`
#' @return A tibble with a STUB label column and res1..res9 n(%) cells
#'
#' Percent denominator is the per-(PARAM, arm) total.
count_wide <- function(d) {
  b <- tplyr_build(tplyr_spec(cols = c("TRTP", "RIND"),
        layers = tplyr_layers(group_count("PARAM",
          settings = layer_settings(
            format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
            denoms_by = c("PARAM", "TRTP"),
            order_count_method = "byfactor",
            zero_count_display = "count_only")))), d)
  b <- b[order(b$ord_layer_1), , drop = FALSE]
  rc <- grep("^res", names(b), value = TRUE)
  out <- tibble(STUB = as.character(b$rowlabel1))
  for (i in seq_along(rc)) out[[paste0("res", i)]] <- as.character(b[[rc[i]]])
  out
}

#' Compute the Fisher's-exact p-value per PARAM over the 3x3 (arm x category) table
#' @param d Prepped lab records from `prep()`
#' @param size Field width passed to `num_fmt()`
#' @return A tibble of STUB and formatted PVAL (3 dp)
#'
#' A p-value is computed for every PARAM (no all-Normal guard).
pvals <- function(d, size = 5) {
  cnt <- d |> count(PARAM, TRTP, RIND, .drop = FALSE) |> arrange(PARAM, TRTP, RIND)
  cnt |> group_by(PARAM) |> group_map(function(x, k) {
    m <- matrix(x$n, nrow = 3, ncol = 3, byrow = TRUE)   # rows = arm, cols = L/N/H
    p <- suppressWarnings(fisher.test(m)$p.value)
    tibble(STUB = as.character(k$PARAM),
           PVAL = num_fmt(p, digits = 3, int_len = 1, size = size))
  }) |> bind_rows()
}

#' Attach the p-value column to a count block (joined on the PARAM label)
#' @param block A count block from `count_wide()`
#' @param d Prepped lab records from `prep()`
#' @return `block` with a PVAL column (NA -> "")
with_p <- function(block, d) {
  block |> left_join(pvals(d), by = "STUB") |>
    mutate(PVAL = ifelse(is.na(PVAL), "", PVAL))
}

chem_d <- prep("adlbcpv", CHEM_RC, CHEM_ORD)
hem_d  <- prep("adlbhpv", HEM_RC, HEM_ORD)
chem   <- with_p(count_wide(chem_d), chem_d)
heme   <- with_p(count_wide(hem_d),  hem_d)

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
final[] <- lapply(final, function(x) {
  x[is.na(x)] <- ""
  attr(x, "label") <- NULL
  as.character(x)
})

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
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    STUB = "",
    res1 = c(PBO, "Low"),  res2 = c(PBO, "Normal"),  res3 = c(PBO, "High"),
    res4 = c(LOW, "Low"),  res5 = c(LOW, "Normal"),  res6 = c(LOW, "High"),
    res7 = c(HIGH, "Low"), res8 = c(HIGH, "Normal"), res9 = c(HIGH, "High"),
    PVAL = c("", "p-val\n[1]")) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "STUB", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "STUB", width = 1.8) |>
  flextable::width(j = paste0("res", 1:9), width = 0.74) |>
  flextable::width(j = "PVAL", width = 0.52)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)

#' Apply the house default plus a solid underline beneath each arm spanner
#' @param x A flextable
#' @param ... Unused
#' @return The styled flextable
#'
#' The rule spans the 9 value columns (cols 2:10) on header row 1; the full-width bottom
#' rule under the labels row comes from cdisc_table_default().
sd <- function(x, ...) {
  x <- cdisc_table_default(x)
  flextable::hline(x, i = 1, j = 2:10, border = officer::fp_border(color = "black", width = 1), part = "header")
}
old <- options(clinify_table_default = sd)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), "\n")
