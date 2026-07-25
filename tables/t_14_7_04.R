# tables/t_14_7_04.R
# Table 14-7.04 — Summary of Concomitant Medications (Number of Subjects) (Safety)
# Nested distinct-subject counts: therapeutic class (CMCLAS) -> coded medication
# (CMDECOD), plus a "patients receiving at least one" top row. Denominator = arm N.
# Bespoke display (n right-justified width 3; UNPADDED integer percent "(8%)"; bare "  0"
# for zero; medications ordered within class by descending Placebo count then alphabetical)
# is reproduced directly - tplyr2 group_count(distinct_by=) is the natural count engine
# (cf. the AE table) but its fixed-width f_str percent can't produce the unpadded "(8%)".
source("R/setup.R"); source("R/helpers.R")

TABLE <- "14-7.04"; SOURCE <- "programs/t-14-7.04.R"
ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

cm   <- read_adam("cm")
adsl <- read_adam("adsl") |> filter(ARM %in% ARMS)
Nvec <- adsl |> count(ARM) |> deframe()                 # arm denominators (86/84/84)
armmap <- adsl |> distinct(USUBJID, ARM)

# n( pct%) exactly as the legacy n_pct(): width-3 n, unpadded round() integer %, bare "  0"
fmt <- function(n, arm) ifelse(is.na(n) | n == 0, "  0",
                               sprintf("%3d (%d%%)", n, as.integer(round(100 * n / Nvec[[arm]]))))

# distinct-subject counts per arm at three levels
n_atleast <- cm |> distinct(USUBJID) |> inner_join(armmap, by = "USUBJID") |> count(ARM)
n_class   <- cm |> distinct(USUBJID, CMCLAS) |> inner_join(armmap, by = "USUBJID") |>
  count(CMCLAS, ARM)
n_med     <- cm |> distinct(USUBJID, CMCLAS, CMDECOD) |> inner_join(armmap, by = "USUBJID") |>
  count(CMCLAS, CMDECOD, ARM)

# per-arm count for a (level filtered) tibble -> named integer vector over ARMS (0 if absent)
counts <- function(df) {
  v <- setNames(integer(length(ARMS)), ARMS)
  m <- setNames(df$n, as.character(df$ARM)); v[names(m)] <- m; v
}
row_of <- function(label, cnt) tibble(LBL = label,
  P = fmt(cnt[["Placebo"]], "Placebo"),
  L = fmt(cnt[["Xanomeline Low Dose"]], "Xanomeline Low Dose"),
  H = fmt(cnt[["Xanomeline High Dose"]], "Xanomeline High Dose"))
blank <- tibble(LBL = "", P = "", L = "", H = "")

acc <- list(row_of("Patients receiving at least one concomitant medication", counts(n_atleast)))
for (cls in sort(unique(n_class$CMCLAS))) {
  acc[[length(acc) + 1]] <- blank
  acc[[length(acc) + 1]] <- row_of(cls, counts(n_class |> filter(CMCLAS == cls)))
  meds <- n_med |> filter(CMCLAS == cls)
  # order: descending by Placebo count, then medication name ascending (legacy radix order)
  ord <- meds |> group_by(CMDECOD) |>
    summarize(pbo = sum(n[ARM == "Placebo"]), .groups = "drop") |>
    arrange(desc(pbo), CMDECOD)
  for (md in ord$CMDECOD)
    acc[[length(acc) + 1]] <- row_of(paste0("    ", md), counts(meds |> filter(CMDECOD == md)))
}
final <- bind_rows(acc)

# --- render (3-page table -> content-verified). House defaults; centered narrow table.
al <- function(arm, n) arm_label(arm, n)   # str_wrap width 10 -> multi-line "Arm\n(N=n)"
ct <- clintable(final, use_labels = FALSE) |>
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
cat("rows:", nrow(final), " classes:", length(unique(n_class$CMCLAS)), " meds:", nrow(distinct(n_med, CMDECOD)), "\n")
