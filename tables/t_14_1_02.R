# t_14_1_02.R
# Table 14-1.02: Summary of End of Study Data   (Population: Intent-to-Treat)
# Produces: outputs/14-1.02.docx
# Source: ADSL; completion status + termination reason as n(%) by arm/Total, Fisher's exact p-values.
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-1.02"; SOURCE <- "programs/t-14-1-02.R"
ARMS  <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose", "Total")

adsl <- read_adam("adsl")
# DCSREAS -> display label; map order sets the display order
reason_lab <- c(
  "Adverse Event"      = "Adverse Event",
  "Death"              = "Death",
  "Lack of Efficacy"   = "Lack of Efficacy[2]",
  "Lost to Follow-up"  = "Lost to Follow-up",
  "Withdrew Consent"   = "Subject decided to withdraw",
  "Physician Decision" = "Physician decided to withdraw subject",
  "I/E Not Met"        = "Protocol criteria not met",
  "Protocol Violation" = "Protocol violation",
  "Sponsor Decision"   = "Sponsor decision")
adsl <- adsl |>
  mutate(
    COMP_STAT = factor(if_else(COMP24FL == "Y", "Completed Week 24",
                               "Early Termination (prior to Week 24)"),
                       levels = c("Completed Week 24", "Early Termination (prior to Week 24)")),
    DCREASCD  = factor(reason_lab[DCSREAS], levels = unname(reason_lab)))

adslT <- bind_rows(adsl, mutate(adsl, TRT01P = "Total")) |>
  mutate(TRT01P = factor(TRT01P, levels = ARMS))

#' Build an n(%) count block for a categorical variable
#' @param data_with_total Analysis data including the bound "Total" arm
#' @param var Name of the categorical variable to count
#' @return Tibble with an indented rowlbl and res1..res4 (Placebo/Low/High/Total); denominators are arm N via pop_data
count_block <- function(data_with_total, var) {
  s <- tplyr_spec(cols = "TRT01P",
                  pop_data = pop_data(cols = c("TRT01P" = "TRT01P")),
                  layers = tplyr_layers(group_count(var,
                    settings = layer_settings(
                      format_strings = list(n_counts = f_str("xxx (xxx%)", "n", "pct")),
                      order_count_method = "byfactor"))))
  b <- tplyr_build(s, data_with_total, pop_data = adslT)
  b <- b[order(b$ord_layer_1), , drop = FALSE]
  rc <- grep("^res", names(b), value = TRUE)              # Placebo/Low/High/Total
  out <- tibble(rowlbl = paste0("  ", as.character(b$rowlabel1)))
  for (i in seq_along(rc)) out[[paste0("res", i)]] <- as.character(b[[rc[i]]])
  out
}
miss_row <- tibble(rowlbl = "  Missing", res1 = "  0 (  0%)", res2 = "  0 (  0%)",
                   res3 = "  0 (  0%)", res4 = "  0 (  0%)")
#' Build a section-header row with a label and empty result cells
#' @param t Section-header text
#' @return One-row tibble with rowlbl = t and empty res1..res4
sec  <- function(t) tibble(rowlbl = t, res1 = "", res2 = "", res3 = "", res4 = "")
#' Build a fully blank spacer row
#' @return One-row tibble with all cells empty
blank <- function() tibble(rowlbl = "", res1 = "", res2 = "", res3 = "", res4 = "")

# Completion status (denominator = full arm N)
comp <- count_block(adslT, "COMP_STAT")
comp_p <- fish_p_str(adsl$COMP24FL, adsl$TRT01P)

# Reason for early termination (numerator = terminated; denom = full arm N)
term_dat <- adslT |> filter(COMP24FL == "N")
term <- count_block(term_dat, "DCREASCD")
# %in% (not ==): completers have NA DCREASCD and must count as 0, not be dropped from the test.
ae_p  <- fish_p_str(as.integer(adsl$DCREASCD %in% "Adverse Event"),      adsl$TRT01P, width = 6)
loe_p <- fish_p_str(as.integer(adsl$DCREASCD %in% "Lack of Efficacy[2]"), adsl$TRT01P, width = 6)

# p-values sit on each section's first data row and the Lack-of-Efficacy row.
comp$p <- ""; comp$p[1] <- comp_p
term$p <- ""; term$p[1] <- ae_p
term$p[term$rowlbl == "  Lack of Efficacy[2]"] <- loe_p
miss_row$p <- ""

final <- bind_rows(
  blank(), sec("Completion Status:"), comp, miss_row,
  blank(), sec("Reason for Early Termination (prior to Week 24):"), term, miss_row)
final <- final |> mutate(p = ifelse(is.na(p), "", p))
final[] <- lapply(final, function(x) { attr(x, "label") <- NULL; as.character(x) })

Ns <- adsl |> count(TRT01P) |> deframe(); Ntot <- sum(Ns)
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    rowlbl = "",
    res1 = arm_label("Placebo", Ns[["Placebo"]]),
    res2 = arm_label("Xanomeline Low Dose",  Ns[["Xanomeline Low Dose"]]),
    res3 = arm_label("Xanomeline High Dose", Ns[["Xanomeline High Dose"]]),
    res4 = arm_label("Total", Ntot),
    p    = "p-value\n[1]") |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "rowlbl", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "rowlbl", width = 3.6) |>
  flextable::width(j = "res1", width = 1.081) |>
  flextable::width(j = c("res2", "res3", "res4", "p"), width = 1.08) |>
  flextable::set_table_properties(align = "center")
# Merge section-header rows across all columns so the long label doesn't wrap in the stub.
sec_rows <- which(final$rowlbl %in% c("Completion Status:",
                                      "Reason for Early Termination (prior to Week 24):"))
for (r in sec_rows) ct <- flextable::merge_at(ct, i = r, j = seq_len(ncol(final)), part = "body")
ct <- flextable::align(ct, i = sec_rows, align = "left", part = "body")
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
