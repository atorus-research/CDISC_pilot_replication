# t_14_7_02.R
# Table 14-7.02: Summary of Vital Signs Change from Baseline at End of Treatment   (Population: Safety)
# Produces: outputs/14-7.02.docx
# Source: ADVS change-from-baseline (CHG) at Week 24/End of Treatment summarized with tplyr2 group_desc; ADSL for arm N.
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-7.02"; SOURCE <- "programs/t-14-7.02.R"
ARM      <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
ARM_DISP <- c("Placebo", "Xan.Low", "Xan.High")
PARAMS   <- c("Systolic Blood Pressure (mmHg)", "Diastolic Blood Pressure (mmHg)", "Pulse (bpm)")
VIS      <- c("Week 24", "End of Trt.")
COLS     <- c("n", "Mean", "SD", "Median", "Min.", "Max.")

FS <- list("n"      = f_str("xx",    "n_records"),
           "Mean"   = f_str("xxx.x", "mean"),
           "SD"     = f_str("xx.xx", "sd"),
           "Median" = f_str("xxx.x", "median"),
           "Min."   = f_str("xxx.x", "min"),
           "Max."   = f_str("xxx.x", "max"))

advs <- read_adam("advs") |> filter(SAFFL == "Y", !is.na(BASE))
advs <- advs |>
  mutate(EOTFL = AVISIT == "End of Treatment",
         W24FL = AVISIT == "Week 24") |>
  filter(EOTFL | W24FL,
         PARAM %in% c("Diastolic Blood Pressure (mmHg)", "Pulse Rate (beats/min)",
                      "Systolic Blood Pressure (mmHg)")) |>
  mutate(PARAM = recode(PARAM, "Pulse Rate (beats/min)" = "Pulse (bpm)"),
         PRTFL = if_else(EOTFL, "End of Trt.", "Week 24"))

#' Summarize CHG into six stat columns per PARAM/ATPT/TRTP/PRTFL
#' @param dat ADVS rows to summarize
#' @return A wide tibble with one row per group and n/Mean/SD/Median/Min./Max. columns
build_wide <- function(dat) {
  dat$ALL <- "x"
  b <- tplyr_build(tplyr_spec(cols = "ALL", layers = tplyr_layers(
    group_desc("CHG", by = c("PARAM", "ATPT", "TRTP", "PRTFL"),
               settings = layer_settings(format_strings = FS)))), dat)
  wide <- b |>
    transmute(PARAM = as.character(rowlabel1), ATPT = as.character(rowlabel2),
              TRTP = as.character(rowlabel3), PRTFL = as.character(rowlabel4),
              stat = as.character(rowlabel5), val = as.character(res1)) |>
    tidyr::pivot_wider(names_from = stat, values_from = val)
  # n is the record count at each visit (records assessed = non-missing + missing),
  # emitted directly by the n_records stat; Mean/SD/... drop NA.
  wide
}

w  <- build_wide(advs)
Ns <- read_adam("adsl") |> count(ARM = ARM) |> deframe()
ATPTS <- sort(unique(w$ATPT))

blank_row <- tibble(Measure = "", Position = "", Treatment = "", N = "", PRT = "",
                    !!!setNames(rep(list(""), 6), COLS))

out <- list()
for (p in PARAMS) for (ai in seq_along(ATPTS)) for (ti in seq_along(ARM)) {
  grp <- w |> filter(PARAM == p, ATPT == ATPTS[ai], TRTP == ARM[ti]) |>
    arrange(match(PRTFL, VIS))
  if (nrow(grp) == 0) next
  blk <- tibble(Measure = "", Position = "", Treatment = "", N = "", PRT = grp$PRTFL)
  for (cc in COLS) blk[[cc]] <- grp[[cc]]
  blk$Treatment[1] <- ARM_DISP[ti]
  blk$N[1] <- as.character(Ns[[ARM[ti]]])
  if (ti == 1)            blk$Position[1] <- ATPTS[ai]
  if (ti == 1 && ai == 1) blk$Measure[1]  <- p
  out[[length(out) + 1]] <- blk
  out[[length(out) + 1]] <- blank_row
}
final <- bind_rows(out) |> select(Measure, Position, Treatment, N, PRT, all_of(COLS))

# coerce_character coerces every column to character; NA renders blank via flextable's default na_str.
ct <- clintable(final, use_labels = FALSE, coerce_character = TRUE) |>
  clin_column_headers(
    Measure = "Measure", Position = "Position", Treatment = "Treatment", N = "N",
    PRT = "Planned Relative Time",
    n = "n", Mean = "Mean", SD = "SD", Median = "Median", Min. = "Min.", Max. = "Max.") |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = c("Measure", "Position", "PRT"), part = "header", align = "left") |>
  flextable::valign(part = "body", valign = "top") |>
  flextable::align(part = "body", align = "left") |>
  flextable::align(j = "Treatment", part = "body", align = "center") |>
  flextable::width(j = "Measure",   width = 1.45) |>
  flextable::width(j = "Position",  width = 1.27) |>
  flextable::width(j = "Treatment", width = 0.90) |>
  flextable::width(j = "N",         width = 0.35) |>
  flextable::width(j = "PRT",       width = 0.85) |>
  flextable::width(j = "n",         width = 0.35) |>
  flextable::width(j = COLS[-1],    width = 0.62)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
cat("rows:", nrow(final), "\n")
