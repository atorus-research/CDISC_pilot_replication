# t_14_6_06.R
# Table 14-6.06: Shifts of Hy's Law Values During Treatment   (Population: Safety)
# Produces: outputs/14-6.06.docx
# Source: adlbhy (adsl for arm Ns); tplyr2 group_shift column-% shift-to AVAL by baseline
#   BASE (Normal=0/High=1) for two Hy's-law analytes (TRANSHY, HYLAW), coin::cmh_test p per analyte.
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R")
source("R/helpers.R")
suppressPackageStartupMessages(library(coin))

TABLE <- "14-6.06"
SOURCE <- "programs/t-14-6.06.R"
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
COLS <- c("Placebo N","Placebo H","Xanomeline Low Dose N","Xanomeline Low Dose H",
          "Xanomeline High Dose N","Xanomeline High Dose H")

# One record per subject per analyte: worst (max AVAL) then latest (max AVISITN),
# among post-baseline records with a non-missing baseline.
adlbhy <- read_adam("adlbhy") |>
  filter(SAFFL == "Y", PARAMCD %in% c("TRANSHY", "HYLAW"), !is.na(BASE), AVISITN > 0) |>
  group_by(USUBJID) |> filter(AVAL == max(AVAL)) |> filter(AVISITN == max(AVISITN)) |> ungroup() |>
  mutate(BASE = factor(BASE, c(0, 1)), AVAL = factor(AVAL, c(0, 1)), TRTP = factor(TRTP, ARMS))

#' Build the shift cells and n-row totals for one analyte
#' @param pc PARAMCD of the analyte ("TRANSHY" or "HYLAW")
#' @return A list with `n` (the n-row, incl. PVAL) and `shifts` (the shift-to rows)
analyte_rows <- function(pc) {
  d <- adlbhy |> filter(PARAMCD == pc)
  # group_shift with shift_denom = "column"; no `by` (one analyte per call) so the column
  # denominator is each TRTP x BASE baseline group. res columns come out TRTP x BASE.
  # as_display() returns the display-ready frame (rowlabel1/res1..res6 + pval1), dropping
  # the internal ord* columns; the build is already factor-ordered so no re-sort.
  #
  # The analyte p-value is computed natively by tplyr2's omnibus assoc_test: the supplied
  # fn runs once over this analyte's RAW source subset (all arms, both AVAL levels, incl.
  # the BASE stratum) and returns the finished display verbatim. Row-mean-scores CMH, AVAL
  # scored 1/2 across treatment groups stratified by baseline status BASE; "" guards when
  # nobody shifted to High, or when the abnormal-at-baseline stratum is empty so
  # "controlling for baseline status" is vacuous (e.g. HYLAW); tryCatch -> "" on model
  # failure. With no `by` var the result lands on the layer's first output row and surfaces
  # as `pval1`, so it is lifted off onto the n-row below where the pilot places it.
  b <- as_display(tplyr_build(tplyr_spec(cols = "TRTP",
        layers = tplyr_layers(group_shift(c(row = "AVAL", column = "BASE"),
          settings = layer_settings(
            shift_denom = "column",
            format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
            order_count_method = "byfactor", zero_count_display = "count_only",
            # the layer's own per-baseline-group denominator, as the reference's 2-char "n" row
            denom_row = TRUE, denom_row_label = "n",
            denom_row_format = f_str("xx", "n"),
            assoc_test = assoc_test(
              fn = function(.data) {
                if (all(.data$AVAL == "0")) return("")
                if (dplyr::n_distinct(.data$BASE) < 2) return("")
                dd <- .data |> transmute(AVAL = ordered(as.character(AVAL), c("0", "1")),
                                         TRTP = factor(as.character(TRTP), levels = ARMS),
                                         BASE = factor(as.character(BASE), levels = c("0", "1")))
                tryCatch(
                  num_fmt(as.numeric(pvalue(cmh_test(AVAL ~ TRTP | BASE, data = dd,
                                                     scores = list(AVAL = c(1, 2))))),
                          digits = 3, int_len = 1),
                  error = function(e) "")
              },
              format = f_str("x.xxxx", "p")))))), d))
  rc <- grep("^res", names(b), value = TRUE)
  rl <- grep("^rowlabel", names(b), value = TRUE)
  all_rows <- tibble(SHIFT = as.character(b[[rl[1]]]))
  for (i in seq_along(rc)) all_rows[[COLS[i]]] <- as.character(b[[rc[i]]])
  # The layer emits its own per-baseline-group denominator as the leading "n" row
  # (denom_row), so the displayed denominator is the one the percentages used.
  n_row <- all_rows |> filter(SHIFT == "n")
  sh    <- all_rows |> filter(SHIFT != "n")
  n_row$PVAL <- dplyr::first(as.character(b$pval1))
  sh$PVAL <- ""
  list(n = n_row, shifts = sh)
}

#' Build a full-width blank/lead row
#' @return A one-row tibble with blank stub, value and p-value cells
blank_row <- function() tibble(LBL = "", SHIFT = "", !!!setNames(rep(list(""), 6), COLS), PVAL = "")

#' Assemble one analyte block: optional wrapped-label lead row, then n / 0 / 1 rows
#' @param pc PARAMCD of the analyte
#' @param lead_label Optional first line of a wrapped stub label (or NULL)
#' @param n_label Stub label placed on the n-row
#' @return A tibble of the analyte's body rows
block <- function(pc, lead_label, n_label) {
  r <- analyte_rows(pc)
  n_row <- r$n
  n_row$LBL <- n_label
  #' Test whether every shift cell is zero (blank or a formatted 0)
  #' @param v Character vector of formatted cells
  #' @return TRUE when all cells are "0" or ""
  zero_cell <- function(v) all(sub("\\s", "", v) %in% c("0", ""))
  rows <- list(n_row[, c("LBL", "SHIFT", COLS, "PVAL")])
  for (k in seq_len(nrow(r$shifts))) {
    row <- r$shifts[k, ]
    row$LBL <- ""
    if (row$SHIFT == "1" && zero_cell(unlist(row[COLS]))) next   # drop all-zero High row
    rows[[length(rows) + 1]] <- row[, c("LBL", "SHIFT", COLS, "PVAL")]
  }
  out <- bind_rows(rows)
  if (!is.null(lead_label)) {                                    # wrapped-label lead row
    ll <- blank_row()
    ll$LBL <- lead_label
    out <- bind_rows(ll, out)
  }
  out
}

final <- bind_rows(
  blank_row(),
  block("TRANSHY", NULL, "Transaminase 1.5 x ULN"),
  blank_row(),
  block("HYLAW", "Total Bili 1.5 x ULN and", "Transaminase 1.5 x ULN")
) |> select(LBL, SHIFT, all_of(COLS), PVAL)

# render: same header as 14-6.05 (arm spanner x baseline, Shift[1], p-value[2])
Ns <- read_adam("adsl") |> filter(ARM != "Screen Failure") |> count(TRT01P) |> deframe()

#' Format an arm spanner label with its N
#' @param short Short arm label
#' @param arm Full arm name for the N lookup
#' @return "short (N=n)"
arm_hdr <- function(short, arm) sprintf("%s (N=%s)", short, Ns[[arm]])
ct <- clintable(final, use_labels = FALSE, coerce_character = TRUE) |>
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
    # spanners merge while the repeated "Baseline" leaf labels stay as separate cells.
    merge = "spanners") |>
  # 1pt rule beneath each arm spanner, across only that arm's two baseline columns
  # (auto-derived from the merged header runs); applied after the house styler, so it
  # survives border_remove(). Replaces per-spanner flextable::hline() calls.
  clin_spanner_rule() |>
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

write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
cat("rows:", nrow(final), "\n")
