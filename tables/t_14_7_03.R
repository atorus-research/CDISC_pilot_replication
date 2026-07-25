# tables/t_14_7_03.R
# Table 14-7.03 — Summary of Weight Change from Baseline at End of Treatment (Safety)
# Two stacked descriptive blocks, arms DOWN the stub (not columns): Weight (kg) AVAL
# at Baseline/Week 24/EOT, then Weight Change (CHG) at Week 24/EOT. Stats n/Mean/SD/
# Median/Min/Max run ACROSS as columns. tplyr2 group_desc is the summary engine; because
# group_desc lays stats out as ROWS, each stat is a single-stat f_str and the six are
# pivoted to columns here. EOT = ADVS "End of Treatment" derived visit (AVISITN==99),
# used as-is from the ADaM (see report re: the documented EOT-flag data discrepancy).
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-7.03"; SOURCE <- "programs/t-14-7.03.R"
ARM      <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
ARM_DISP <- c("Placebo", "Xan.Low", "Xan.High")

# Per-stat f_str matching the reference num_fmt widths (int_len=2, so "xx.x"/"xx.xx";
# values that overflow the width, e.g. 108.0, expand rather than truncate — verified
# byte-identical to the legacy num_fmt output).
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

# tplyr2 desc -> per-(TRTP, PRTFL) six stat columns
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

# One arm's rows (visits down), then a trailing blank row (reference pad_row after EOT).
arm_rows <- function(wide, arm, arm_disp, visits) {
  pv <- wide |> filter(TRTP == arm) |> arrange(match(PRTFL, visits))
  b <- tibble(Measure = "", Treatment = "", N = "", PRT = pv$PRTFL)
  for (cc in COLS) b[[cc]] <- pv[[cc]]
  b$Treatment[1] <- arm_disp
  b$N[1] <- as.character(Ns[[arm]])
  bind_rows(b, blank_row)
}

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

# --- render: single-row header; arms down the stub, stats centered across
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
  # Reference column proportions were Measure 0.25 / PRT 0.185 of a 9in table. The
  # reference header font is narrower than Courier, so it fit "Planned Relative Time"
  # in a 1.665in PRT column on one line; in Courier 10pt that string needs ~1.82in or
  # it wraps to two header rows (which misregisters the whole body). Widen PRT to 1.82
  # and take the 0.16in from the (over-generous) Measure column so the six numeric stat
  # columns stay on the reference's x-positions (they carry the bulk of the body ink).
  flextable::width(j = "Measure",   width = 2.088) |>
  flextable::width(j = "Treatment", width = 1.35) |>
  flextable::width(j = "N",         width = 0.45) |>
  flextable::width(j = "PRT",       width = 1.82) |>
  flextable::width(j = "n",         width = 0.36) |>
  flextable::width(j = COLS[-1],    width = 0.585)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)

# Per-table page geometry: the legacy program called set_header_height(1), baking a
# 1in header band (\headery1440) into the reference RTF, so its title band + body sit
# ~0.5in lower than the suite's 0.5in-header house default. Match it here (option
# override restored after write, same pattern as t_14_6_01.R's table-default override)
# so the single-page pixel comparison registers.
old <- options(clinify_docx_default = officer::prop_section(
  page_size    = officer::page_size(width = 8.5, height = 11, orient = "landscape"),
  type         = "continuous",
  page_margins = officer::page_mar(top = 1, bottom = 1, left = 1, right = 1,
                                   header = 1, footer = 0.5, gutter = 0)))
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), "\n")
