# t_14_7_03.R
# Table 14-7.03: Summary of Weight Change from Baseline at End of Treatment   (Population: Safety)
# Produces: outputs/14-7.03.docx
# Source: ADVS Weight (kg) AVAL and change (CHG) at Baseline/Week 24/End of Treatment via tplyr2 group_desc; ADSL for arm N.
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-7.03"; SOURCE <- "programs/t-14-7.03.R"
ARM      <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
ARM_DISP <- c("Placebo", "Xan.Low", "Xan.High")

# Per-stat f_str widths (int_len 2, so "xx.x"/"xx.xx"); values wider than the format
# expand rather than truncate.
FS <- list("n"      = f_str("xx",   "n"),
           "Mean"   = f_str("xx.x", "mean"),
           "SD"     = f_str("xx.xx","sd"),
           "Median" = f_str("xx.x", "median"),
           "Min."   = f_str("xx.x", "min"),
           "Max."   = f_str("xx.x", "max"))

advs <- read_adam("advs") |> filter(PARAM == "Weight (kg)")
advs <- advs |>
  mutate(TRTP  = factor(TRTP, levels = ARM),
         PRTFL = case_when(ABLFL == "Y" ~ "Baseline",
                           AVISITN == 24 ~ "Week 24",
                           AVISITN == 99 ~ "End of Trt.",
                           TRUE ~ NA_character_)) |>
  filter(!is.na(PRTFL))

#' Summarize a variable into six stat columns per TRTP/PRTFL
#' @param dat ADVS rows to summarize
#' @param var Name of the analysis variable to summarize (e.g. "AVAL", "CHG")
#' @return A wide tibble with one row per group and n/Mean/SD/Median/Min./Max. columns
build_wide <- function(dat, var) {
  dat$ALL <- "x"
  b <- tplyr_build(tplyr_spec(cols = "ALL", layers = tplyr_layers(
    group_desc(var, by = c("TRTP", "PRTFL"),
               settings = layer_settings(format_strings = FS)))), dat)
  b |>
    transmute(TRTP = as.character(rowlabel1), PRTFL = as.character(rowlabel2),
              stat = as.character(rowlabel3), val = as.character(res1)) |>
    tidyr::pivot_wider(names_from = stat, values_from = val)
}

Ns <- read_adam("adsl") |> count(ARM = ARM) |> deframe()

COLS <- c("n", "Mean", "SD", "Median", "Min.", "Max.")
blank_row <- tibble(Measure = "", Treatment = "", N = "", PRT = "",
                    !!!setNames(rep(list(""), 6), COLS))

#' Build one arm's block of rows (visits down the stub) plus a trailing blank row
#' @param wide Wide stat tibble from build_wide()
#' @param arm Treatment arm value to filter on
#' @param arm_disp Display label for the arm, shown on the block's first row
#' @param visits Visit labels in display order
#' @return A tibble of the arm's rows followed by one blank spacer row
arm_rows <- function(wide, arm, arm_disp, visits) {
  pv <- wide |> filter(TRTP == arm) |> arrange(match(PRTFL, visits))
  b <- tibble(Measure = "", Treatment = "", N = "", PRT = pv$PRTFL)
  for (cc in COLS) b[[cc]] <- pv[[cc]]
  b$Treatment[1] <- arm_disp
  b$N[1] <- as.character(Ns[[arm]])
  bind_rows(b, blank_row)
}

#' Stack all arms into one measure block, labelling only its first row
#' @param wide Wide stat tibble from build_wide()
#' @param measure_label Measure label shown on the block's first row
#' @param visits Visit labels in display order
#' @return A tibble with all arms' rows for the measure
assemble <- function(wide, measure_label, visits) {
  out <- lapply(seq_along(ARM), function(i) arm_rows(wide, ARM[i], ARM_DISP[i], visits))
  res <- bind_rows(out)
  res$Measure[1] <- measure_label      # measure label only on the block's first row
  res
}

wt_wide  <- build_wide(advs, "AVAL")
chg_wide <- build_wide(advs |> filter(PRTFL != "Baseline"), "CHG")

final <- bind_rows(
  assemble(wt_wide,  "Weight (kg)",             c("Baseline", "Week 24", "End of Trt.")),
  assemble(chg_wide, "Weight Change\nfrom Baseline", c("Week 24", "End of Trt."))
) |> select(Measure, Treatment, N, PRT, all_of(COLS))
final[] <- lapply(final, function(x) { x[is.na(x)] <- ""; as.character(x) })

# render: single-row header; arms down the stub, stats centered across.
ct <- clintable(final, use_labels = FALSE) |>
  clin_column_headers(
    Measure = "Measure", Treatment = "Treatment", N = "N",
    PRT = "Planned Relative Time",
    n = "n", Mean = "Mean", SD = "SD", Median = "Median", Min. = "Min.", Max. = "Max.") |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::valign(part = "body", valign = "top") |>
  flextable::align(part = "body", align = "left") |>
  flextable::align(j = c("N", COLS), part = "body", align = "center") |>
  # PRT is 1.82in so "Planned Relative Time" fits on one header line in Courier 10pt;
  # the extra width is taken from the Measure column to keep the numeric columns' positions.
  flextable::width(j = "Measure",   width = 2.088) |>
  flextable::width(j = "Treatment", width = 1.35) |>
  flextable::width(j = "N",         width = 0.45) |>
  flextable::width(j = "PRT",       width = 1.82) |>
  flextable::width(j = "n",         width = 0.36) |>
  flextable::width(j = COLS[-1],    width = 0.585)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)

# Page geometry override: use a 1in header band for this table (restored after write).
old <- options(clinify_docx_default = officer::prop_section(
  page_size    = officer::page_size(width = 8.5, height = 11, orient = "landscape"),
  type         = "continuous",
  page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1,
                                   header = 1, footer = 0.5, gutter = 0)))
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), "\n")
