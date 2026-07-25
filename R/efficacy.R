# R/efficacy.R
# -----------------------------------------------------------------------------
# Reusable builder for the ANCOVA efficacy tables (14-3.01 .. 14-3.09, ex .07).
#
# Two endpoint shapes:
#   endpoint = "ADAS"  -> change-from-baseline (CHG); descriptive blocks
#                         Baseline / Week N / Change; ANCOVA model includes the
#                         baseline value as a covariate.
#   endpoint = "CIBIC" -> single impression score (AVAL) at Week N; one
#                         descriptive block; ANCOVA model has NO baseline covariate.
#
# Descriptive block via tplyr2 group_desc; the ANCOVA (dose-response + pairwise
# LS-means) is fit once with lm + car + emmeans and row-bound beneath, formatted
# byte-for-byte via num_fmt (matches the reference display).
# -----------------------------------------------------------------------------

.eff_fs <- list(
  "  n"              = f_str("xx", "n"),
  "  Mean (SD)"      = f_str("xx.x (xx.xx)", "mean", "sd"),
  "  Median (Range)" = f_str("xx.x (xxx;xx)", "median", "min", "max")
)

# ANCOVA model block -> rows (row_label, res1/res2/res3) matching the reference.
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
  g <- function(cn, col) pwc[[col]][pwc$contrast == cn]
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

# Optional parameters (all default to the original ADAS/CIBIC behavior):
#   extra_filter : an unquoted expression AND-ed onto the base filter after
#                  EFFFL/ITTFL/PARAMCD (e.g. the 14-3.07 completers subset).
#   anl01        : apply the ANL01FL=="Y" filter (FALSE for 14-3.12 / adnpix).
#   derive_chg   : compute CHG = AVAL - BASE when the dataset lacks it (14-3.12).
#   blocks       : override the descriptive-block spec, a list of
#                  list(var=, label=, visit=); default is the endpoint's shape.
#   use_base     : ANCOVA baseline covariate; default (endpoint=="ADAS").
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
  # rlang::inject substitutes the visit value INTO the layer `where` before
  # tplyr2 captures it (layer where is a bare expression, so a plain `visit`
  # variable would not resolve).
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
