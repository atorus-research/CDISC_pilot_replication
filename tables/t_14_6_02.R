# t_14_6_02.R
# Table 14-6.02: Frequency of Normal and Abnormal (Beyond Normal Range) Laboratory
#                Values During Treatment   (Population: Safety)
# Produces: outputs/14-6.02.docx
# Source: adlbc/adlbh (one on-treatment analysis record per subject per PARAM, adsl for
#   arm Ns); tplyr2 group_count subject n(%) by arm x LBNRIND category, Fisher's-exact p per PARAM.
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

#' Build the 9-column n(%) block (PARAM rows x arm-category columns)
#' @param d Prepped lab records from `prep()`
#' @return A tibble with a STUB label column and res1..res9 n(%) cells
#'
#' Percent denominator is the per-(PARAM, arm) total, so zeros display as a bare " 0".
count_wide <- function(d) {
  b <- tplyr_build(tplyr_spec(cols = c("TRTP", "RIND"),
        layers = tplyr_layers(group_count("PARAM",
          settings = layer_settings(
            format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
            denoms_by = c("PARAM", "TRTP"),
            order_count_method = "byfactor",
            zero_count_display = "count_only")))), d)
  b <- b[order(b$ord_layer_1), , drop = FALSE]
  rc <- grep("^res", names(b), value = TRUE)       # 9: arm(L,N,H) x 3 arms
  out <- tibble(STUB = as.character(b$rowlabel1))
  for (i in seq_along(rc)) out[[paste0("res", i)]] <- as.character(b[[rc[i]]])
  out
}

#' Compute the Fisher's-exact p-value per PARAM
#' @param d Prepped lab records from `prep()`
#' @param na_rule When TRUE, blank the p-value where every Low and High cell is 0
#' @param size Field width passed to `num_fmt()`
#' @return A tibble of STUB and formatted PVAL (NA -> "")
#'
#' The test is over the 3x3 (arm x category) table. na_rule leaves the p-value blank when
#' all subjects are Normal (every Low and High cell 0), where Fisher's test is undefined.
pvals <- function(d, na_rule = FALSE, size = 4) {
  cnt <- d |> count(PARAM, TRTP, RIND, .drop = FALSE) |> arrange(PARAM, TRTP, RIND)
  cnt |> group_by(PARAM) |> group_map(function(x, k) {
    m <- matrix(x$n, nrow = 3, ncol = 3, byrow = TRUE)   # rows = arm, cols = L/N/H
    p <- if (na_rule && all(m[, 1] == 0) && all(m[, 3] == 0)) NA_real_
         else suppressWarnings(fisher.test(m)$p.value)
    tibble(STUB = as.character(k$PARAM),
           PVAL = num_fmt(p, digits = 3, int_len = 1, size = size))  # NA -> ""
  }) |> bind_rows()
}

#' Attach the p-value column to a count block (joined on the PARAM label)
#' @param block A count block from `count_wide()`
#' @param d Prepped lab records from `prep()`
#' @param na_rule Passed through to `pvals()`
#' @return `block` with a PVAL column (NA -> "")
with_p <- function(block, d, na_rule = FALSE) {
  block |> left_join(pvals(d, na_rule = na_rule), by = "STUB") |>
    mutate(PVAL = ifelse(is.na(PVAL), "", PVAL))
}

chem_d <- prep("adlbc", CHEM_RC, CHEM_ORD)
hem_d  <- prep("adlbh", HEM_RC, HEM_ORD)
chem   <- with_p(count_wide(chem_d), chem_d, na_rule = FALSE)
heme   <- with_p(count_wide(hem_d),  hem_d,  na_rule = TRUE)

COLS <- c(paste0("res", 1:9), "PVAL")

#' Build a section/blank row (STUB text, all value cells blank)
#' @param stub Text for the STUB column
#' @return A one-row tibble with `stub` in STUB and "" in all value columns
row9 <- function(stub = "") tibble(STUB = stub, !!!setNames(rep(list(""), 10), COLS))

final <- bind_rows(
  row9(""), row9("CHEMISTRY"), row9("----------"), chem,
  row9(""), row9("HEMATOLOGY"), row9("----------"), heme)
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
