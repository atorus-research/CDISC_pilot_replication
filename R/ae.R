# R/ae.R
# -----------------------------------------------------------------------------
# Reusable builder for the AE-incidence-by-SOC/Preferred-Term tables
# (14-5.01 all TEAEs; 14-5.02 serious TEAEs).
#
# tplyr2 nested group_count + stat_columns -> paired n(%) / [AEs] columns per arm
# and the "ANY BODY SYSTEM" total row (factor-ordered columns, native bare-count
# on zero). Then: re-sort PTs within SOC by descending High-dose count, inject
# per-row Fisher's exact p-values, blank [AEs] on zero, insert per-SOC spacers,
# and render with a 2-row spanned header repeated on each page.
# -----------------------------------------------------------------------------

build_ae_table <- function(table, source_path, serious = FALSE) {
  ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
  adae <- read_adam("adae") |> filter(SAFFL == "Y", TRTEMFL == "Y")
  if (serious) adae <- adae |> filter(AESER == "Y")
  adae <- adae |> mutate(TRTA = factor(TRTA, levels = ARMS))
  adsl <- read_adam("adsl")
  N <- c(86, 84, 84)   # safety-population denominators (Placebo, Low, High)

  spec <- tplyr_spec(
    cols = "TRTA",
    pop_data = pop_data(cols = c("TRTA" = "TRT01A")),
    layers = tplyr_layers(group_count(c("AEBODSYS", "AEDECOD"),
      settings = layer_settings(
        distinct_by = "USUBJID",
        stat_columns = list("n" = f_str("xx (xx.x%)", "distinct_n", "distinct_pct"),
                            "e" = f_str("[x]", "n")),
        limit_data_by = c("AEBODSYS", "AEDECOD"),
        total_row = TRUE, total_row_label = "ANY BODY SYSTEM",
        zero_count_display = "count_only")))
  )
  b <- tplyr_build(spec, adae, pop_data = adsl)
  rc <- grep("^res", names(b), value = TRUE)   # res1..res6, factor-ordered

  lead_int <- function(s) suppressWarnings(as.integer(sub("^\\s*([0-9]+).*$", "\\1", s)))
  d <- tibble(
    soc = as.character(b$rowlabel1), pt = as.character(b$rowlabel2), depth = b$ord_layer_2,
    npct_0  = as.character(b[[rc[1]]]), e_0  = as.character(b[[rc[2]]]),
    npct_54 = as.character(b[[rc[3]]]), e_54 = as.character(b[[rc[4]]]),
    npct_81 = as.character(b[[rc[5]]]), e_81 = as.character(b[[rc[6]]])
  ) |>
    mutate(n0 = lead_int(npct_0), n54 = lead_int(npct_54), n81 = lead_int(npct_81)) |>
    mutate(soc_rank = if_else(depth == 0, 0L, 1L), soc_key = if_else(depth == 0, "", soc)) |>
    arrange(soc_rank, soc_key, depth, desc(n81), pt)

  fisher_ae <- function(n_p, n_a, N_p, N_a) {
    if ((n_p + n_a) == 0) return("")
    p <- fisher.test(matrix(c(n_p, n_a, N_p - n_p, N_a - n_a), nrow = 2))$p.value
    disp <- num_fmt(p, digits = 3, size = 5, int_len = 1)
    if (p > .99) ">.99" else if (p < .15) paste0(disp, "*") else paste0(disp, " ")
  }
  d <- d |> rowwise() |>
    mutate(p_low  = fisher_ae(n0, n54, N[1], N[2]),
           p_high = fisher_ae(n0, n81, N[1], N[3])) |>
    ungroup() |>
    mutate(
      AETERM = if_else(depth == 2, paste0("  ", pt), soc),
      c_ae_0  = if_else(n0  > 0, e_0,  ""),
      c_ae_54 = if_else(n54 > 0, e_54, ""),
      c_ae_81 = if_else(n81 > 0, e_81, ""))

  disp <- d |> transmute(
    AETERM, npct_0, cAEs_0 = c_ae_0, npct_54, cAEs_54 = c_ae_54,
    npct_81, cAEs_81 = c_ae_81, p_low, p_high,
    grp = if_else(depth == 0L, "___ANY", soc))
  # blank spacer after ANY and after each SOC block
  blankrow <- disp[1, ]; blankrow[] <- ""
  g <- disp$grp; acc <- list()
  for (i in seq_len(nrow(disp))) {
    acc[[length(acc) + 1]] <- disp[i, ]
    if (i == nrow(disp) || g[i] != g[i + 1]) acc[[length(acc) + 1]] <- blankrow
  }
  final <- top_spacer(bind_rows(acc) |> select(-grp))

  hdr <- function(arm, n) arm_label(arm, n)
  ct <- clintable(final, use_labels = FALSE) |>
    clin_column_headers(
      AETERM  = c("", "System Organ Class/\nPreferred Term"),
      npct_0  = c(hdr("Placebo", N[1]), "n(%)"),               cAEs_0  = c(hdr("Placebo", N[1]), "[AEs]"),
      npct_54 = c(hdr("Xanomeline Low Dose", N[2]), "n(%)"),   cAEs_54 = c(hdr("Xanomeline Low Dose", N[2]), "[AEs]"),
      npct_81 = c(hdr("Xanomeline High Dose", N[3]), "n(%)"),  cAEs_81 = c(hdr("Xanomeline High Dose", N[3]), "[AEs]"),
      p_low   = c("Fisher's Exact\np-values", "Placebo\nvs.\nLow Dose"),
      p_high  = c("Fisher's Exact\np-values", "Placebo\nvs.\nHigh Dose")) |>
    flextable::valign(part = "header", valign = "bottom") |>
    flextable::align(part = "header", align = "center") |>
    flextable::align(j = "AETERM", part = "header", align = "left") |>
    flextable::align(part = "body", align = "left") |>
    flextable::width(j = "AETERM", width = 2.7) |>
    flextable::width(j = c("npct_0", "npct_54", "npct_81"), width = 0.9) |>
    flextable::width(j = c("cAEs_0", "cAEs_54", "cAEs_81"), width = 0.63) |>
    flextable::width(j = "p_low", width = 0.81) |>
    flextable::width(j = "p_high", width = 0.9)
  ct <- add_titles_footnotes(ct, table, source_path = source_path, date = FIDELITY_DATE)
  write_clindoc(ct, file.path(OUTPUT_DIR, paste0(table, ".docx")))
  invisible(final)
}
