# tables/t_14_6_06.R
# Table 14-6.06 — Shifts of Hy's Law Values During Treatment (Safety).
# Same transposed-shift layout as 14-6.05 for two derived Hy's-law analytes:
#   * Transaminase 1.5 x ULN                     (PARAMCD == "TRANSHY")
#   * Total Bili 1.5 x ULN and Transaminase ...  (PARAMCD == "HYLAW")
# Baseline status (Normal=0 / High=1) is a COLUMN dimension crossed with treatment;
# shift-to (n / 0 / 1) are rows, percentages column-wise per baseline group ->
# tplyr2 cols=c("TRTP","BASE") count (the 14-6.04/14-6.05 crossed-cols workaround
# for tplyr2 FR #18). Per-analyte CMH p-value (row-mean-scores, baseline-stratified)
# via coin::cmh_test. adlbhy is used (adlbc/adlbh lack the Hy's-law derivations).
source("R/setup.R"); source("R/helpers.R")
suppressPackageStartupMessages(library(coin))

TABLE <- "14-6.06"; SOURCE <- "programs/t-14-6.06.R"
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
COLS <- c("Placebo N","Placebo H","Xanomeline Low Dose N","Xanomeline Low Dose H",
          "Xanomeline High Dose N","Xanomeline High Dose H")

# One record per subject per analyte: worst (max AVAL) then latest (max AVISITN),
# among post-baseline records with a non-missing baseline (legacy selection).
adlbhy <- read_adam("adlbhy") |>
  filter(SAFFL == "Y", PARAMCD %in% c("TRANSHY", "HYLAW"), !is.na(BASE), AVISITN > 0) |>
  group_by(USUBJID) |> filter(AVAL == max(AVAL)) |> filter(AVISITN == max(AVISITN)) |> ungroup() |>
  mutate(BASE = factor(BASE, c(0, 1)), AVAL = factor(AVAL, c(0, 1)), TRTP = factor(TRTP, ARMS))

# Row-mean-scores CMH (ordinal shift AVAL scored 1/2, treatment groups, stratified
# by baseline status). coin is the clean CRAN CMH (the 14-3.13 pattern). Blank when
# nobody shifted to High, or when the abnormal-at-baseline stratum is empty so
# "controlling for baseline status" is vacuous (reproduces the legacy blank for HYLAW).
cmh_pval <- function(d) {
  if (all(d$AVAL == "0")) return("")
  if (dplyr::n_distinct(d$BASE) < 2) return("")
  dd <- d |> transmute(AVAL = ordered(as.character(AVAL), c("0", "1")),
                       TRTP = factor(as.character(TRTP), levels = ARMS),
                       BASE = factor(as.character(BASE), levels = c("0", "1")))
  tryCatch(
    num_fmt(as.numeric(pvalue(cmh_test(AVAL ~ TRTP | BASE, data = dd,
                                       scores = list(AVAL = c(1, 2))))),
            digits = 3, int_len = 1),
    error = function(e) "")
}

# shift cells (n(%) per baseline group, column-wise) + n-row totals for one analyte.
analyte_rows <- function(pc) {
  d <- adlbhy |> filter(PARAMCD == pc)
  # Native group_shift with shift_denom="column" (tplyr2 #28/#30 + #31/#32). No `by` here (one
  # analyte per call), so the column denominator is each TRTP x BASE baseline group. res -> TRTP x BASE.
  b <- tplyr_build(tplyr_spec(cols = "TRTP",
        layers = tplyr_layers(group_shift(c(row = "AVAL", column = "BASE"),
          settings = layer_settings(
            shift_denom = "column",
            format_strings = list(n_counts = f_str("xx(xxx%)", "n", "pct")),
            order_count_method = "byfactor", zero_count_display = "count_only")))), d)
  rc <- grep("^res", names(b), value = TRUE); rl <- grep("^rowlabel", names(b), value = TRUE)
  sh <- tibble(SHIFT = as.character(b[[rl[1]]]))
  for (i in seq_along(rc)) sh[[COLS[i]]] <- as.character(b[[rc[i]]])
  nr <- d |> count(TRTP, BASE, .drop = FALSE) |>
    mutate(col = paste(TRTP, recode(as.character(BASE), "0" = "N", "1" = "H"))) |>
    select(col, n) |> pivot_wider(names_from = col, values_from = n, values_fill = 0)
  nr[COLS] <- lapply(nr[COLS], function(x) sprintf("%2d", x))
  n_row <- tibble(SHIFT = "n"); for (c in COLS) n_row[[c]] <- nr[[c]]
  n_row$PVAL <- cmh_pval(d)
  sh$PVAL <- ""
  list(n = n_row, shifts = sh)
}

# assemble a block: optional label-only lead row, then n / 0 / 1 rows.
blank_row <- function() tibble(LBL = "", SHIFT = "", !!!setNames(rep(list(""), 6), COLS), PVAL = "")
block <- function(pc, lead_label, n_label) {
  r <- analyte_rows(pc)
  n_row <- r$n; n_row$LBL <- n_label
  zero_cell <- function(v) all(sub("\\s", "", v) %in% c("0", ""))
  rows <- list(n_row[, c("LBL", "SHIFT", COLS, "PVAL")])
  for (k in seq_len(nrow(r$shifts))) {
    row <- r$shifts[k, ]; row$LBL <- ""
    if (row$SHIFT == "1" && zero_cell(unlist(row[COLS]))) next   # drop all-zero High row
    rows[[length(rows) + 1]] <- row[, c("LBL", "SHIFT", COLS, "PVAL")]
  }
  out <- bind_rows(rows)
  if (!is.null(lead_label)) {                                    # wrapped-label lead row
    ll <- blank_row(); ll$LBL <- lead_label
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
final[] <- lapply(final, function(x) { x[is.na(x)] <- ""; as.character(x) })

# --- render: identical header to 14-6.05 (arm spanner x baseline, Shift[1], p-value[2])
Ns <- read_adam("adsl") |> filter(ARM != "Screen Failure") |> count(TRT01P) |> deframe()
arm_hdr <- function(short, arm) sprintf("%s (N=%s)", short, Ns[[arm]])
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    LBL = "", SHIFT = c("", "Shift", "[1]"),
    `Placebo N`              = c(arm_hdr("Placebo", "Placebo"), "Normal at", "Baseline"),
    `Placebo H`              = c(arm_hdr("Placebo", "Placebo"), "High at", "Baseline"),
    `Xanomeline Low Dose N`  = c(arm_hdr("Xan. Low", "Xanomeline Low Dose"), "Normal at", "Baseline"),
    `Xanomeline Low Dose H`  = c(arm_hdr("Xan. Low", "Xanomeline Low Dose"), "High at", "Baseline"),
    `Xanomeline High Dose N` = c(arm_hdr("Xan. High", "Xanomeline High Dose"), "Normal at", "Baseline"),
    `Xanomeline High Dose H` = c(arm_hdr("Xan. High", "Xanomeline High Dose"), "High at", "Baseline"),
    PVAL = c("", "p-\nvalue", "[2]"),
    # clinify #95: merge="spanners" merges every header row EXCEPT the bottom one, so the arm
    # spanners merge while the repeated "Baseline" leaf labels stay as separate cells. (Was a
    # merge_none() + manual merge_at() workaround before clinify gained per-row merge control.)
    merge = "spanners") |>
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

sd <- function(x, ...) { x <- cdisc_table_default(x)
  bd <- officer::fp_border(color = "black", width = 1)
  x <- flextable::hline(x, i = 1, j = 3:4, border = bd, part = "header")
  x <- flextable::hline(x, i = 1, j = 5:6, border = bd, part = "header")
  flextable::hline(x, i = 1, j = 7:8, border = bd, part = "header") }
old <- options(clinify_table_default = sd)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), "\n")
