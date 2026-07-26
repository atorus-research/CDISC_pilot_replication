# t_14_1_01.R
# Table 14-1.01: Summary of Populations   (Population: All Subjects)
# Produces: outputs/14-1.01.docx
# Source: ADSL; per population flag, n(%) with flag = "Y" by arm + Total via tplyr2 group_count.
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-1.01"; SOURCE <- "programs/t-14-1-01.R"
ARMS  <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

adsl0 <- read_adam("adsl")
adsl <- adsl0 |>
  mutate(TRT01P = factor(TRT01P, levels = ARMS),
         COMPL  = if_else(DCDECOD == "COMPLETED", "Y", "N"))
Ns <- adsl0 |> count(TRT01P) |> deframe(); Ntot <- sum(Ns)

#' Build one display row of n(%) for a population flag
#' @param flag Name of the population flag variable, counted where its value is "Y"
#' @param label Row label shown in the stub
#' @return One-row tibble with rowlbl1 and res1..res4 (Placebo/Low/High/Total)
pop_row <- function(flag, label) {
  s <- tplyr_spec(cols = "TRT01P",
                  total_groups = list(total_group("TRT01P")),
                  layers = tplyr_layers(group_count(flag,
                    settings = layer_settings(
                      format_strings = list(n_counts = f_str("xxx (xxx%)", "n", "pct")),
                      keep_levels = "Y"))))
  b <- as_display(tplyr_build(s, adsl))
  tibble(rowlbl1 = label,
         res1 = b$res1, res2 = b$res2, res3 = b$res3, res4 = b$res4)
}

final <- top_spacer(bind_rows(
  pop_row("ITTFL",    "Intent-To-Treat (ITT)"),
  pop_row("SAFFL",    "Safety"),
  pop_row("EFFFL",    "Efficacy"),
  pop_row("COMP24FL", "Complete Week 24"),
  pop_row("COMPL",    "Complete Study")))

ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    rowlbl1 = "",
    res1 = arm_label("Placebo", Ns[["Placebo"]]),
    res2 = arm_label("Xanomeline Low Dose",  Ns[["Xanomeline Low Dose"]]),
    res3 = arm_label("Xanomeline High Dose", Ns[["Xanomeline High Dose"]]),
    res4 = arm_label("Total", Ntot)) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "rowlbl1", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "rowlbl1", width = 2.64) |>
  flextable::width(j = c("res1", "res2", "res3", "res4"), width = 0.99) |>
  flextable::set_table_properties(align = "center")
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
