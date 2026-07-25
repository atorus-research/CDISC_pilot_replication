# tables/t_14_7_01.R
# Table 14-7.01 — Summary of Vital Signs at Baseline and End of Treatment (Safety)
# Descriptive summary with the analysis dimensions DOWN the stub (Measure/PARAM,
# Position/ATPT, Treatment/TRTP, Planned Relative Time/visit) and the six statistics
# n/Mean/SD/Median/Min/Max ACROSS as columns. tplyr2 group_desc is the summary engine
# (by = PARAM,ATPT,TRTP,visit; each stat a single-stat f_str); group_desc emits stats
# as ROWS so they are pivoted to columns here, then repeated stub labels are blanked
# and a spacer row is inserted after each treatment block (reference pad_row). Multi-
# page -> verified by pagination-agnostic CONTENT (text-set), so stub/label wrapping is
# reproduced via column widths + top valign to match the reference's physical lines.
#
# EOT = ADVS "End of Treatment" derived visit (AVISIT == "End of Treatment", AVISITN 99)
# used as-is from the ADaM. See report re: the documented EOT-flag data discrepancy.
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-7.01"; SOURCE <- "programs/t-14-7.01.R"
ARM      <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
ARM_DISP <- c("Placebo", "Xan.Low", "Xan.High")
PARAMS   <- c("Systolic Blood Pressure (mmHg)", "Diastolic Blood Pressure (mmHg)", "Pulse (bpm)")
VIS      <- c("Baseline", "Week 24", "End of Trt.")
COLS     <- c("n", "Mean", "SD", "Median", "Min.", "Max.")

FS <- list("n"      = f_str("xx",    "n"),
           "Mean"   = f_str("xxx.x", "mean"),
           "SD"     = f_str("xx.xx", "sd"),
           "Median" = f_str("xxx.x", "median"),
           "Min."   = f_str("xxx.x", "min"),
           "Max."   = f_str("xxx.x", "max"))

advs <- read_adam("advs") |> filter(SAFFL == "Y", ANL01FL == "Y")
advs <- advs |>
  mutate(EOTFL = AVISIT == "End of Treatment",
         W24FL = AVISIT == "Week 24") |>
  filter(EOTFL | W24FL | ABLFL == "Y",
         PARAM %in% c("Diastolic Blood Pressure (mmHg)", "Pulse Rate (beats/min)",
                      "Systolic Blood Pressure (mmHg)")) |>
  mutate(PARAM = recode(PARAM, "Pulse Rate (beats/min)" = "Pulse (bpm)"),
         PRTFL = case_when(ABLFL == "Y" ~ "Baseline", W24FL ~ "Week 24", EOTFL ~ "End of Trt."))

# tplyr2 desc -> per (PARAM, ATPT, TRTP, PRTFL) six stat columns.
build_wide <- function(dat) {
  dat$ALL <- "x"
  b <- tplyr_build(tplyr_spec(cols = "ALL", layers = tplyr_layers(
    group_desc("AVAL", by = c("PARAM", "ATPT", "TRTP", "PRTFL"),
               settings = layer_settings(format_strings = FS)))), dat)
  b |>
    transmute(PARAM = as.character(rowlabel1), ATPT = as.character(rowlabel2),
              TRTP = as.character(rowlabel3), PRTFL = as.character(rowlabel4),
              stat = as.character(rowlabel5), val = as.character(res1)) |>
    tidyr::pivot_wider(names_from = stat, values_from = val)
}

w  <- build_wide(advs)
Ns <- read_adam("adsl") |> count(ARM = ARM) |> deframe()
ATPTS <- sort(unique(w$ATPT))

blank_row <- tibble(Measure = "", Position = "", Treatment = "", N = "", PRT = "",
                    !!!setNames(rep(list(""), 6), COLS))

# Assemble: stub labels shown only on the first row of their group; spacer row after
# each treatment block. valign top (at render) keeps each visit's stats on the first
# physical line while the tall Measure/Position cells wrap beneath — matching the ref.
out <- list()
for (p in PARAMS) for (ai in seq_along(ATPTS)) for (ti in seq_along(ARM)) {
  grp <- w |> filter(PARAM == p, ATPT == ATPTS[ai], TRTP == ARM[ti]) |>
    arrange(match(PRTFL, VIS))
  if (nrow(grp) == 0) next
  blk <- tibble(Measure = "", Position = "", Treatment = "", N = "", PRT = grp$PRTFL)
  for (cc in COLS) blk[[cc]] <- grp[[cc]]
  blk$Treatment[1] <- ARM_DISP[ti]
  blk$N[1] <- as.character(Ns[[ARM[ti]]])
  if (ti == 1)            blk$Position[1] <- ATPTS[ai]   # first arm of the position block
  if (ti == 1 && ai == 1) blk$Measure[1]  <- p          # first arm+position of the param block
  out[[length(out) + 1]] <- blk
  out[[length(out) + 1]] <- blank_row
}
final <- bind_rows(out) |> select(Measure, Position, Treatment, N, PRT, all_of(COLS))
final[] <- lapply(final, function(x) { x[is.na(x)] <- ""; as.character(x) })

# --- render: header labels bottom-aligned; PRT header "Planned Relative Time" wraps to
# 3 lines in its narrow column (matches ref); body stats left, treatment centered.
ct <- clintable(final, use_labels = FALSE) |>
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
