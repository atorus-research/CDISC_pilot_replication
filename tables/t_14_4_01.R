# t_14_4_01.R
# Table 14-4.01: Summary of Planned Exposure to Study Drug, as of End of Study   (Population: Safety)
# Produces: outputs/14-4.01.docx
# Source: ADSL; descriptive stats of AVGDD and CUMDOSE by arm x population (Completers / Safety).
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-4.01"; SOURCE <- "programs/t-14-4-01.R"
GRP <- c("0_C", "54_C", "81_C", "0_S", "54_S", "81_S")

adsl <- read_adam("adsl")
compl  <- adsl |> filter(COMP24FL == "Y") |>
  transmute(AVGDD, CUMDOSE, TRTPCD = paste0(TRT01PN, "_C"))
safety <- adsl |> filter(SAFFL == "Y") |>
  transmute(AVGDD, CUMDOSE, TRTPCD = paste0(TRT01PN, "_S"))
adsl_ <- bind_rows(compl, safety) |> mutate(TRTPCD = factor(TRTPCD, levels = GRP))
Ns <- adsl_ |> count(TRTPCD, .drop = FALSE) |> deframe()

# dose values are large (cumulative up to thousands) -> int width 5
fs <- list(
  "n"      = f_str("xxxxx",    "n"),
  "Mean"   = f_str("xxxxx.x",  "mean"),
  "SD"     = f_str("xxxxx.xx", "sd"),
  "Median" = f_str("xxxxx.x",  "median"),
  "Min"    = f_str("xxxxx.x",  "min"),
  "Max"    = f_str("xxxxx.x",  "max")
)

#' Build a descriptive-statistics block for one dose variable
#' @param var Name of the numeric variable to summarize.
#' @param label Label placed on the block's first row.
#' @return A padded tibble with row labels and six result columns.
dblock <- function(var, label) {
  b <- tplyr_build(tplyr_spec(cols = "TRTPCD",
         layers = tplyr_layers(group_desc(var,
           settings = layer_settings(format_strings = fs)))), adsl_)
  d <- as_display(b)                           # display-ready + ordered; drops ord*/row_id cols
  rc <- grep("^res", names(d), value = TRUE)   # six result cols, in GRP factor order
  out <- tibble(rowlbl1 = "", rowlbl2 = as.character(d$rowlabel1))
  out$rowlbl1[1] <- label
  for (i in seq_along(rc)) out[[paste0("res", i)]] <- as.character(d[[rc[i]]])
  pad_row(out)
}

final <- top_spacer(bind_rows(
  dblock("AVGDD",   "Average daily dose (mg)"),
  dblock("CUMDOSE", "Cumulative dose at end of study [2]")))

sp_c <- "Completers at Week 24"; sp_s <- "Safety Population [1]"
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    rowlbl1 = "", rowlbl2 = "",
    res1 = c(sp_c, arm_label("Placebo", Ns[["0_C"]])),
    res2 = c(sp_c, arm_label("Xanomeline Low Dose",  Ns[["54_C"]])),
    res3 = c(sp_c, arm_label("Xanomeline High Dose", Ns[["81_C"]])),
    res4 = c(sp_s, arm_label("Placebo", Ns[["0_S"]])),
    res5 = c(sp_s, arm_label("Xanomeline Low Dose",  Ns[["54_S"]])),
    res6 = c(sp_s, arm_label("Xanomeline High Dose", Ns[["81_S"]]))) |>
  clin_spanner_rule() |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = c("rowlbl1", "rowlbl2"), part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::width(j = "rowlbl1", width = 3.24) |>
  flextable::width(j = "rowlbl2", width = 0.63) |>
  flextable::width(j = c("res1", "res4"), width = 0.9) |>
  flextable::width(j = c("res2", "res3", "res5", "res6"), width = 0.99) |>
  flextable::set_table_properties(align = "center")
# header_leading compacts this table's multi-line arm labels; clinify applies it after
# the option styler, so it survives the house line_spacing(space = 1, part = "all").
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE, date = FIDELITY_DATE,
                           header_leading = 0.75)

#' Apply the house table default and compact the header
#' @param x A flextable to style.
#' @param ... Unused; matches the `clinify_table_default` signature.
#' @return The styled flextable.
#' @details The spanner underline is drawn by `clin_spanner_rule()` in the table
#'   pipeline, which clinify applies after this styler so it survives border_remove().
spanned_default <- function(x, ...) cdisc_table_default(x)
old <- options(clinify_table_default = spanned_default)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
