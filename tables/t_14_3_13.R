# t_14_3_13.R
# Table 14-3.13: CIBIC+ - Categorical Analysis - LOCF   (Population: Efficacy)
# Produces: outputs/14-3.13.docx
# Per-visit CIBIC+ category counts (adcibc) with a row-mean-scores CMH p-value via coin::cmh_test.
source("R/setup.R"); source("R/helpers.R")
suppressPackageStartupMessages(library(coin))

TABLE <- "14-3.13"; SOURCE <- "programs/t-14-3-13.R"
CIB <- c("Marked improvement", "Moderate improvement", "Minimal improvement", "No Change",
         "Minimal worsening", "Moderate worsening", "Marked worsening")   # AVAL 1..7
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

cbic <- read_adam("adcibc") |>
  filter(EFFFL == "Y", ITTFL == "Y", AVISITN %in% c(8, 16, 24), ANL01FL == "Y") |>
  mutate(AVALC = factor(CIB[AVAL], levels = CIB),
         TRTP  = factor(TRTP, levels = ARMS))
hn <- cbic |> distinct(USUBJID, TRTP) |> count(TRTP) |> deframe()   # header N (79/81/74)

#' Compute the row-mean-scores CMH p-value for one visit.
#' @param wk Analysis visit number (8, 16, or 24).
#' @return Character string with the formatted p-value (ordinal AVAL by treatment, stratified by site group).
cmh_pval <- function(wk) {
  d <- cbic |> filter(AVISITN == wk) |> mutate(AVAL = ordered(AVAL), SITEGR1 = factor(SITEGR1))
  p <- pvalue(cmh_test(AVAL ~ TRTP | SITEGR1, data = d, scores = list(AVAL = seq_len(nlevels(d$AVAL)))))
  num_fmt(p, digits = 4, size = 5, int_len = 1)
}

#' Build the display rows for one visit block.
#' @param wk Analysis visit number (8, 16, or 24).
#' @param lab Visit label shown on the block's first row.
#' @return A padded tibble of one n row plus the seven CIBIC category rows.
visit_block <- function(wk, lab) {
  d <- cbic |> filter(AVISITN == wk)
  b <- tplyr_build(tplyr_spec(cols = "TRTP", layers = tplyr_layers(group_count("AVALC",
         settings = layer_settings(
           format_strings = list(n_counts = f_str("xx (xxx%)", "n", "pct")),
           order_count_method = "byfactor", zero_count_display = "count_only")))), d)
  b <- b[order(b$ord_layer_1), , drop = FALSE]
  rc <- grep("^res", names(b), value = TRUE)
  cats <- tibble(AVALC = as.character(b$rowlabel1),
                 res1 = as.character(b[[rc[1]]]), res2 = as.character(b[[rc[2]]]),
                 res3 = as.character(b[[rc[3]]]), p = "")
  vn <- d |> count(TRTP, .drop = FALSE) |> deframe()
  n_row <- tibble(AVALC = "n",
                  res1 = sprintf("%2d", vn[["Placebo"]]),
                  res2 = sprintf("%2d", vn[["Xanomeline Low Dose"]]),
                  res3 = sprintf("%2d", vn[["Xanomeline High Dose"]]),
                  p = cmh_pval(wk))
  blk <- bind_rows(n_row, cats)
  blk$AVISIT <- ""; blk$AVISIT[1] <- lab
  pad_row(select(blk, AVISIT, AVALC, res1, res2, res3, p))
}

final <- top_spacer(bind_rows(visit_block(8, "Week 8"), visit_block(16, "Week 16"),
                              visit_block(24, "Week 24")))
final[] <- lapply(final, function(x) { x[is.na(x)] <- ""; as.character(x) })

ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    AVISIT = "", AVALC = "Assessment",
    res1 = arm_label("Placebo", hn[["Placebo"]]),
    res2 = arm_label("Xanomeline Low Dose",  hn[["Xanomeline Low Dose"]]),
    res3 = arm_label("Xanomeline High Dose", hn[["Xanomeline High Dose"]]),
    p    = "p-value\n[1]") |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = c("AVISIT", "AVALC"), part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "AVISIT", width = 0.9) |>
  flextable::width(j = "AVALC", width = 2.7) |>
  flextable::width(j = c("res1", "res2", "res3", "p"), width = 0.9) |>
  flextable::set_table_properties(align = "center")
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
