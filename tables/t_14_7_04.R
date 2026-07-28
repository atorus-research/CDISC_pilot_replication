# t_14_7_04.R
# Table 14-7.04: Summary of Concomitant Medications (Number of Subjects)   (Population: Safety)
# Produces: outputs/14-7.04.docx
# Source: CM distinct-subject counts by therapeutic class (CMCLAS) and medication (CMDECOD); ADSL for arm N.
library(tidyverse)
library(tplyr2)
library(clinify)

source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-7.04"; SOURCE <- "programs/t-14-7.04.R"
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
TOPLBL <- "Patients receiving at least one concomitant medication"

cm   <- read_adam("cm")
adsl <- read_adam("adsl") |> filter(ARM %in% ARMS)
armmap <- adsl |> distinct(USUBJID, ARM)                # CM carries no treatment var
cm2   <- cm   |> inner_join(armmap, by = "USUBJID") |> mutate(ARM = factor(ARM, ARMS))
adslp <- adsl |> mutate(ARM = factor(ARM, ARMS))

# Native tplyr2 nested distinct-subject counts: outer therapeutic class (CMCLAS) x inner
# medication (CMDECOD), n(%) = distinct subjects / arm N. pop_data supplies the ADSL arm
# denominators (86/84/84) so an arm absent from a level renders "  0" (zero_count_display).
# f_str("xxx (x%)", ...) => width-3 n and an unpadded single-char percent, e.g. "  8 (9%)".
# total_row rolls the distinct-subject grand total into the TOPLBL "at least one" row.
b <- tplyr_build(tplyr_spec(cols = "ARM",
       pop_data = pop_data(cols = "ARM"),
       layers = tplyr_layers(group_count(c("CMCLAS", "CMDECOD"),
         settings = layer_settings(
           distinct_by = "USUBJID",
           format_strings = list(n_counts = f_str("xxx (x%)", "distinct_n", "distinct_pct")),
           zero_count_display = "count_only",
           total_row = TRUE, total_row_label = TOPLBL)))),
     cm2, pop_data = adslp)
# Header (N=) per arm from the build's own header N, so the displayed denominator is
# provably the one used for the percentages (rather than a second count of ADSL).
hn   <- tplyr_header_n(b)
Nvec <- setNames(hn$.n, as.character(hn[[1]]))

disp <- as_display(b)                                   # rowlabel1=CMCLAS, rowlabel2=CMDECOD, res1..3
names(disp)[match(c("res1", "res2", "res3"), names(disp))] <- c("P", "L", "H")

# split the display frame: the grand-total "at least one" row, the class (outer) rows, and
# the medication (inner) rows. as_display() emits meds alphabetically within each class.
top     <- disp[disp$rowlabel1 == TOPLBL, ]
classes <- disp[disp$rowlabel2 == "" & disp$rowlabel1 != TOPLBL, ]
meds    <- disp[disp$rowlabel2 != "", ]
meds$pbo <- as.integer(stringr::str_extract(meds$P, "\\d+"))   # Placebo count drives med order

#' Build one display row from a single-row slice of the display frame
#' @param label Row label (therapeutic class or indented medication)
#' @param r One-row data frame carrying P/L/H formatted count strings
#' @return A one-row tibble with LBL/P/L/H columns
row_of <- function(label, r) tibble(LBL = label, P = r$P, L = r$L, H = r$H)
blank <- tibble(LBL = "", P = "", L = "", H = "")

acc <- list(row_of(TOPLBL, top))
for (cls in sort(unique(classes$rowlabel1))) {          # classes alphabetical
  acc[[length(acc) + 1]] <- blank
  acc[[length(acc) + 1]] <- row_of(cls, classes[classes$rowlabel1 == cls, ])
  # order: descending by Placebo count, ties broken by medication name ascending
  # (stable sort over the already-alphabetical med rows)
  md <- meds[meds$rowlabel1 == cls, ]
  md <- md[order(-md$pbo), ]
  for (i in seq_len(nrow(md)))
    acc[[length(acc) + 1]] <- row_of(paste0("    ", md$rowlabel2[i]), md[i, ])
}
final <- bind_rows(acc)

# render: house defaults; centered narrow table.
#' Wrap arm_label() to produce a multi-line "Arm\n(N=n)" header label
#' @param arm Arm name
#' @param n Arm denominator
#' @return A wrapped header label string
al <- function(arm, n) arm_label(arm, n)
ct <- clintable(final, use_labels = FALSE, coerce_character = TRUE) |>
  clin_column_headers(
    LBL = "Therapeutic class, n (%)",
    P = al("Placebo", Nvec[["Placebo"]]),
    L = al("Xanomeline Low Dose", Nvec[["Xanomeline Low Dose"]]),
    H = al("Xanomeline High Dose", Nvec[["Xanomeline High Dose"]])) |>
  flextable::valign(part = "header", valign = "bottom") |>
  flextable::align(part = "header", align = "center") |>
  flextable::align(j = "LBL", part = "header", align = "left") |>
  flextable::align(part = "body", align = "left") |>
  flextable::valign(part = "body", valign = "top") |>   # counts sit on the wrapped label's first line
  flextable::width(j = "LBL", width = 4.0) |>
  flextable::width(j = c("P", "L", "H"), width = 1.1) |>
  flextable::set_table_properties(align = "center")
ct <- add_titles_footnotes(ct, TABLE, source_path = SOURCE)
write_clindoc(ct, file.path(OUTPUT_DIR, paste0(TABLE, ".docx")))
cat("rows:", nrow(final), " classes:", nrow(classes), " meds:", nrow(meds), "\n")
