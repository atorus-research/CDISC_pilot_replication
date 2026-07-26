# t_14_1_03.R
# Table 14-1.03: Summary of Number of Subjects By Site   (Population: All Subjects)
# Produces: outputs/14-1.03.docx
# Source: ADSL (ITTFL == "Y"); direct crosstab of site x arm subject counts for ITT/Eff/Com flags.
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-1.03"; SOURCE <- "programs/t-14-1-03.R"
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
CODE <- c(Placebo = "P", `Xanomeline Low Dose` = "L", `Xanomeline High Dose` = "H", Total = "T")
COLS <- c("P_ITT","P_Eff","P_Com","L_ITT","L_Eff","L_Com",
          "H_ITT","H_Eff","H_Com","T_ITT","T_Eff","T_Com")

a <- read_adam("adsl") |> filter(ITTFL == "Y") |>
  mutate(SITEGR1 = as.character(SITEGR1), SITEID = as.character(SITEID),
         TRT01P = as.character(TRT01P))

# counts per site x arm x flag; then a Total arm summing across the 3 arms
cnt <- a |> group_by(SITEGR1, SITEID, TRT01P) |>
  summarize(ITT = sum(ITTFL == "Y"), Eff = sum(EFFFL == "Y"),
            Com = sum(COMP24FL == "Y"), .groups = "drop")
cnt <- bind_rows(cnt, cnt |> group_by(SITEGR1, SITEID) |>
  summarize(across(c(ITT, Eff, Com), sum), .groups = "drop") |> mutate(TRT01P = "Total"))

w <- cnt |> mutate(k = CODE[TRT01P]) |>
  pivot_wider(id_cols = c(SITEGR1, SITEID), names_from = k, values_from = c(ITT, Eff, Com),
              values_fill = 0, names_glue = "{k}_{.value}") |>
  arrange(SITEGR1, SITEID)
w <- bind_rows(w, w |> summarize(across(all_of(COLS), sum)) |>
                 mutate(SITEGR1 = "TOTAL", SITEID = ""))

final <- tibble(PID = w$SITEGR1, SID = w$SITEID)
for (cc in COLS) final[[cc]] <- w[[cc]]

# header: arm spanner (str_wrap width 10) + ITT/Eff/Com sub-labels
N <- a |> count(TRT01P) |> deframe()
#' Build a wrapped arm spanner label with its N
#' @param arm Treatment-arm name
#' @param n Subject count for the arm
#' @return Spanner label string wrapped at width 10 with an "(N=n)" line
alab <- function(arm, n) paste0(gsub("\n", "\n", stringr::str_wrap(arm, width = 10)), "\n(N=", n, ")")
sp <- c(P = alab("Placebo", N[["Placebo"]]),
        L = alab("Xanomeline Low Dose", N[["Xanomeline Low Dose"]]),
        H = alab("Xanomeline High Dose", N[["Xanomeline High Dose"]]),
        T = paste0("Total\n(N=", nrow(a), ")"))
#' Build a two-line column header (arm spanner + sub-label)
#' @param code Arm code (P, L, H, or T)
#' @param sub Sub-label (ITT, Eff, or Com)
#' @return Length-2 character vector: arm spanner and sub-label
hdr <- function(code, sub) c(sp[[code]], sub)
ct <- clintable(final, use_labels = FALSE, coerce_character = TRUE) |>
  clin_column_headers(
    PID = c("", "Pooled\nId"), SID = c("", "Site\nId"),
    P_ITT = hdr("P", "ITT"), P_Eff = hdr("P", "Eff"), P_Com = hdr("P", "Com"),
    L_ITT = hdr("L", "ITT"), L_Eff = hdr("L", "Eff"), L_Com = hdr("L", "Com"),
    H_ITT = hdr("H", "ITT"), H_Eff = hdr("H", "Eff"), H_Com = hdr("H", "Com"),
    T_ITT = hdr("T", "ITT"), T_Eff = hdr("T", "Eff"), T_Com = hdr("T", "Com")) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = c("PID", "SID"), part = "body", align = "center") |>
  flextable::align(j = COLS, part = "body", align = "right") |>
  flextable::width(j = "PID", width = 0.75) |>
  flextable::width(j = "SID", width = 0.60) |>
  flextable::width(j = COLS, width = 0.45)
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)

#' Apply clinify table defaults plus the arm-spanner underline and sub-label padding
#' @param x A flextable to style
#' @param ... Additional arguments (unused)
#' @return Styled flextable with a fixed layout, spanner underline, and padded sub-label row
sd <- function(x, ...) { x <- cdisc_table_default(x)
  x <- flextable::set_table_properties(x, layout = "fixed", align = "center")
  x <- flextable::padding(x, i = 2, padding.top = 34, part = "header")
  flextable::hline(x, i = 1, j = 3:14,
    border = officer::fp_border(color = "black", width = 1), part = "header") }
old <- options(clinify_table_default = sd)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
options(old)
cat("rows:", nrow(final), "\n")
