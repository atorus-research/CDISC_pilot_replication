# R/efficacy.R — build the ANCOVA efficacy tables (ADAS change-from-baseline; CIBIC score)

# Descriptive-statistic format strings shared by the descriptive blocks.
.eff_fs <- list(
  "  n"              = f_str("xx", "n"),
  "  Mean (SD)"      = f_str("xx.x (xx.xx)", "mean", "sd"),
  "  Median (Range)" = f_str("xx.x (xxx;xx)", "median", "min", "max")
)

#' Fit the ANCOVA model and format its result rows for the table body
#' @param dat Analysis data for one endpoint, filtered to the analysis population
#' @param week Analysis visit (AVISITN) to model
#' @param use_base Include the baseline value as a covariate
#' @return A tibble of formatted rows (row_label, res1, res2, res3)
ancova_block <- function(dat, week, use_base) {
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op))
  md <- dat |>
    filter(AVISITN == week) |>
    mutate(TRTPCD = case_when(TRTPN == 0 ~ "Pbo", TRTPN == 54 ~ "Xan_Lo",
                              TRTPN == 81 ~ "Xan_Hi"),
           TRTPCD_F = factor(TRTPCD, levels = c("Xan_Hi", "Xan_Lo", "Pbo")))
  if (use_base) {                                   # ADAS: response CHG + baseline covariate
    m1 <- lm(CHG ~ TRTPN + SITEGR1 + BASE, data = md)
    m2 <- lm(CHG ~ TRTPCD_F + SITEGR1 + BASE, data = md)
  } else {                                          # CIBIC: response AVAL, no baseline covariate
    m1 <- lm(AVAL ~ TRTPN + SITEGR1, data = md)
    m2 <- lm(AVAL ~ TRTPCD_F + SITEGR1, data = md)
  }
  dose_p <- car::Anova(m1, type = 3)[2, "Pr(>F)"]
  lsm <- emmeans::lsmeans(m2, ~ TRTPCD_F, weights = "proportional")
  pw  <- as_tibble(summary(emmeans::contrast(lsm, method = "pairwise", adjust = NULL)))
  ci  <- as_tibble(confint(emmeans::contrast(lsm, method = "pairwise", adjust = NULL)))
  pwc <- pw |> left_join(ci |> select(contrast, lower.CL, upper.CL), by = "contrast") |>
    mutate(
      p       = num_fmt(p.value,  int_len = 4, digits = 3, size = 12),
      diff_se = sprintf("%s (%s)", num_fmt(estimate, int_len = 2, digits = 1, size = 4),
                                   num_fmt(SE,       int_len = 1, digits = 2, size = 4)),
      ci      = sprintf("(%s;%s)", num_fmt(lower.CL, int_len = 2, digits = 1, size = 4),
                                   num_fmt(upper.CL, int_len = 1, digits = 1, size = 3)))
  #' Pull one formatted cell for a pairwise contrast
  #' @param cn Contrast name (e.g. "Xan_Hi - Pbo")
  #' @param col Column to extract ("p", "diff_se", or "ci")
  #' @return Character scalar
  g <- function(cn, col) pwc[[col]][pwc$contrast == cn]
  #' Build an empty ANCOVA spacer row
  #' @return A one-row tibble of blank cells
  blank <- function() tibble(row_label = "", res1 = "", res2 = "", res3 = "")
  bind_rows(
    blank(),
    tibble(row_label = "p-value(Dose Response) [1][2]", res1 = "", res2 = "",
           res3 = num_fmt(dose_p, int_len = 4, digits = 3, size = 12)),
    blank(),
    tibble(row_label = c("p-value(Xan - Placebo) [1][3]", "  Diff of LS Means (SE)", "  95% CI"),
           res1 = "",
           res2 = c(g("Xan_Lo - Pbo", "p"), g("Xan_Lo - Pbo", "diff_se"), g("Xan_Lo - Pbo", "ci")),
           res3 = c(g("Xan_Hi - Pbo", "p"), g("Xan_Hi - Pbo", "diff_se"), g("Xan_Hi - Pbo", "ci"))),
    blank(),
    tibble(row_label = c("p-value(Xan High - Xan Low) [1][3]", "  Diff of LS Means (SE)", "  95% CI"),
           res1 = "", res2 = "",
           res3 = c(g("Xan_Hi - Xan_Lo", "p"), g("Xan_Hi - Xan_Lo", "diff_se"), g("Xan_Hi - Xan_Lo", "ci")))
  )
}

#' Build and write an ANCOVA efficacy table
#' @param table Table id (e.g. "14-3.01")
#' @param source_path Source program path shown in the footer
#' @param dataset ADaM dataset name to read
#' @param paramcd PARAMCD value to subset to
#' @param week Analysis visit (AVISITN) and label week
#' @param endpoint Endpoint shape: "ADAS" (change-from-baseline) or "CIBIC" (score)
#' @param sex Optional SEX subset
#' @param extra_filter Unquoted expression AND-ed onto the base EFFFL/ITTFL/PARAMCD filter
#' @param anl01 Apply the ANL01FL == "Y" filter
#' @param derive_chg Compute CHG = AVAL - BASE when the dataset lacks it
#' @param blocks Optional descriptive-block spec: list of list(var=, label=, visit=)
#' @param use_base Include the baseline covariate in the ANCOVA; defaults to endpoint == "ADAS"
#' @return The assembled table data frame (invisibly); writes the .docx as a side effect
build_efficacy_table <- function(table, source_path, dataset, paramcd, week,
                                 endpoint = c("ADAS", "CIBIC"), sex = NULL,
                                 extra_filter = NULL, anl01 = TRUE,
                                 derive_chg = FALSE, blocks = NULL, use_base = NULL) {
  endpoint <- match.arg(endpoint)
  ef <- rlang::enquo(extra_filter)
  dat <- read_adam(dataset) |> filter(EFFFL == "Y", ITTFL == "Y", PARAMCD == paramcd)
  if (anl01)                       dat <- dat |> filter(ANL01FL == "Y")
  if (!rlang::quo_is_null(ef))     dat <- dat |> filter(!!ef)
  if (!is.null(sex))               dat <- dat |> filter(SEX == sex)
  if (derive_chg)                  dat <- dat |> mutate(CHG = AVAL - BASE)
  dat <- dat |> mutate(TRTP = factor(TRTP, levels = c("Placebo", "Xanomeline Low Dose",
                                                      "Xanomeline High Dose")))
  hn <- dat |> distinct(USUBJID, TRTP) |> count(TRTP) |> deframe()

  wk_lab <- paste("Week", week)
  # rlang::inject substitutes the visit value into the layer `where` before tplyr2
  # captures it: `where` is a bare expression, so a plain `visit` variable would not resolve.
  #' Build a descriptive-block group_desc layer from a block spec
  #' @param bk Block spec: list(var=, label=, visit=)
  #' @return A tplyr2 group_desc layer
  mk_desc <- function(bk) rlang::inject(
    group_desc(bk$var, by = label(bk$label), where = AVISITN == !!bk$visit,
               settings = layer_settings(format_strings = .eff_fs)))
  if (is.null(blocks)) {
    blocks <- if (endpoint == "ADAS")
      list(list(var = "AVAL", label = "Baseline",             visit = 0),
           list(var = "AVAL", label = wk_lab,                 visit = week),
           list(var = "CHG",  label = "Change from Baseline", visit = week))
    else
      list(list(var = "AVAL", label = wk_lab, visit = week))
  }
  if (is.null(use_base)) use_base <- (endpoint == "ADAS")
  desc <- tplyr_build(tplyr_spec(cols = "TRTP",
                                 layers = do.call(tplyr_layers, lapply(blocks, mk_desc))), dat) |>
    collapse_row_labels() |>
    select(row_label, res1, res2, res3)

  final <- top_spacer(bind_rows(desc, ancova_block(dat, week, use_base)))
  final[] <- lapply(final, function(x) { attr(x, "label") <- NULL; as.character(x) })

  ct <- clintable(final, use_labels = FALSE) |>
    clin_column_headers(
      row_label = "",
      res1 = arm_label("Placebo",              hn[["Placebo"]]),
      res2 = arm_label("Xanomeline Low Dose",  hn[["Xanomeline Low Dose"]]),
      res3 = arm_label("Xanomeline High Dose", hn[["Xanomeline High Dose"]])) |>
    flextable::valign(part = "header", valign = "bottom") |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = "row_label", part = "header", align = "left") |>
    flextable::align(part = "body", align = "left") |>
    flextable::width(j = "row_label", width = 3.6) |>
    flextable::width(j = c("res1", "res2", "res3"), width = 1.2) |>
    flextable::set_table_properties(align = "center")
  ct <- add_titles_footnotes(ct, table, source_path = source_path, date = FIDELITY_DATE)
  write_clindoc(ct, file.path(OUTPUT_DIR, paste0(table, ".docx")))
  invisible(final)
}
