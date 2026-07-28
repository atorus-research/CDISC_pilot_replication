# t_14_1_02.R
# Table 14-1.02: Summary of End of Study Data   (Population: Intent-to-Treat)
# Produces: outputs/14-1.02.docx
# Source: ADSL; completion status + termination reason as n(%) by arm/Total, Fisher's exact p-values.
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-1.02"; SOURCE <- "programs/t-14-1-02.R"
ARMS  <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

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
    DCREASCD  = factor(reason_lab[DCSREAS], levels = unname(reason_lab)),
    TRT01P    = factor(TRT01P, levels = ARMS))

#' Build an n(%) count block for a categorical variable
#' @param data Analysis data; the "Total" arm is generated natively via total_group
#' @param var Name of the categorical variable to count
#' @param assoc Optional assoc_test() spec; when supplied, its native omnibus
#'   p-value column (pval1) is carried out as a `p` column, landing on the
#'   block's first output row
#' @return Tibble with an indented rowlbl and res1..res4 (Placebo/Low/High/Total); denominators are arm N via pop_data
count_block <- function(data, var, assoc = NULL) {
  s <- tplyr_spec(cols = "TRT01P",
                  total_groups = list(total_group("TRT01P")),
                  # Intent-to-treat population, stated in the spec.
                  pop_data = pop_data(cols = "TRT01P", where = ITTFL == "Y"),
                  layers = tplyr_layers(group_count(var,
                    settings = layer_settings(
                      format_strings = list(n_counts = f_str("xxx (xxx%)", "n", "pct")),
                      order_count_method = "byfactor",
                      missing_count = list(label = "Missing"),
                      assoc_test = assoc))))
  b <- tplyr_build(s, data, pop_data = adsl)
  # order() keeps the NA-ord missing_count row last (as_display would place it first)
  b <- b[order(b$ord_layer_1), , drop = FALSE]
  rc <- grep("^res", names(b), value = TRUE)              # Placebo/Low/High/Total
  out <- tibble(rowlbl = paste0("  ", as.character(b$rowlabel1)))
  for (i in seq_along(rc)) out[[paste0("res", i)]] <- as.character(b[[rc[i]]])
  # Native omnibus p (pval1) already sits on the first row post-order()
  if (!is.null(assoc)) out$p <- as.character(b$pval1)
  out
}
#' Build a section-header row with a label and empty result cells
#' @param t Section-header text
#' @return One-row tibble with rowlbl = t and empty res1..res4
sec  <- function(t) tibble(rowlbl = t, res1 = "", res2 = "", res3 = "", res4 = "")
#' Build a fully blank spacer row
#' @return One-row tibble with all cells empty
blank <- function() tibble(rowlbl = "", res1 = "", res2 = "", res3 = "", res4 = "")

# Completion status (denominator = full arm N)
# Native omnibus assoc_test: Fisher's exact of COMP24FL vs TRT01P computed once
# over the layer's raw subset, landing on the block's first row ("Completed
# Week 24"). The subset arrives with the total_group() rows appended - a "Total"
# pseudo-arm duplicating every subject - so restrict to ARMS to compute the same
# 2x3 test the hand-rolled helper did. Character return reproduces the "<.0001"
# floor and fixed-width format verbatim.
comp <- count_block(adsl, "COMP_STAT",
  assoc = assoc_test(
    fn = function(.data) {
      d <- .data[as.character(.data$TRT01P) %in% ARMS, , drop = FALSE]
      p <- suppressWarnings(
        fisher.test(factor(d$COMP24FL), factor(as.character(d$TRT01P)))$p.value)
      if (round(p, 4) == 0) return("<.0001")
      format(round(p, 4), width = 10, nsmall = 4)
    },
    format = f_str("x.xxxx", "p")))

# Reason for early termination (numerator = terminated; denom = full arm N)
# These two p-values stay hand-rolled: omnibus assoc_test emits ONE p per
# by-group on the group's first row, but this block needs two (Adverse Event and
# Lack of Efficacy) on two non-adjacent rows of a single un-by'd layer. The
# tests also span the full ITT population, while the layer's subset is
# terminators only.
term_dat <- adsl |> filter(COMP24FL == "N")
term <- count_block(term_dat, "DCREASCD")
# %in% (not ==): completers have NA DCREASCD and must count as 0, not be dropped from the test.
ae_p  <- fish_p_str(as.integer(adsl$DCREASCD %in% "Adverse Event"),      adsl$TRT01P, width = 6)
loe_p <- fish_p_str(as.integer(adsl$DCREASCD %in% "Lack of Efficacy[2]"), adsl$TRT01P, width = 6)

# term p-values sit on the section's first data row and the Lack-of-Efficacy row.
# The trailing "Missing" row in each block is emitted natively via missing_count.
term$p <- ""; term$p[1] <- ae_p
term$p[term$rowlbl == "  Lack of Efficacy[2]"] <- loe_p

final <- bind_rows(
  blank(), sec("Completion Status:"), comp,
  blank(), sec("Reason for Early Termination (prior to Week 24):"), term)
final <- final |> mutate(p = ifelse(is.na(p), "", p))

Ns <- adsl |> count(TRT01P) |> deframe(); Ntot <- sum(Ns)
ct <- clintable(final, use_labels = FALSE, coerce_character = TRUE) |>
  clin_column_headers(
    rowlbl = "",
    res1 = arm_label("Placebo", Ns[["Placebo"]]),
    res2 = arm_label("Xanomeline Low Dose",  Ns[["Xanomeline Low Dose"]]),
    res3 = arm_label("Xanomeline High Dose", Ns[["Xanomeline High Dose"]]),
    res4 = arm_label("Total", Ntot),
    # One line, as the legacy program had it: `programs/t-14-1-02.R` builds this
    # header as "p-value [1]" with a space, and the reference RTF cell carries no
    # \line while its neighbours do. Splitting it made the cell two lines deep; with
    # the header bottom-aligned that put "p-value" a line higher than the reference.
    # 11 chars at 6pt is 66pt, inside the column's 73.76pt of usable width.
    p    = "p-value [1]") |>
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
