# R/ae.R — build the AE-incidence-by-SOC/preferred-term tables (all and serious TEAEs)

library(tidyverse)
library(tplyr2)
library(clinify)

#' Build and write an AE-incidence-by-SOC/preferred-term table
#' @param table Table id (e.g. "14-5.01")
#' @param source_path Source program path shown in the footer
#' @param serious Restrict to serious TEAEs (AESER == "Y")
#' @return The assembled table data frame (invisibly); writes the .docx as a side effect
build_ae_table <- function(table, source_path, serious = FALSE) {
  ARMS <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")
  adae <- read_adam("adae") |> filter(SAFFL == "Y", TRTEMFL == "Y")
  if (serious) adae <- adae |> filter(AESER == "Y")
  adae <- adae |> mutate(TRTA = factor(TRTA, levels = ARMS))
  adsl <- read_adam("adsl")

  #' Fisher's exact p-value display for one arm-comparison 2x2 (reference vs comparison)
  #' @param m 2x2 matrix c(n_ref, n_cmp, N_ref - n_ref, N_cmp - n_cmp) supplied by assoc_test
  #' @return "" when neither arm has an event, ">.99" ceiling, "*" flag for p < .15
  ae_p <- function(m) {
    if (sum(m[, 1]) == 0) return(NA_character_)
    p <- fisher.test(m)$p.value
    disp <- num_fmt(p, digits = 3, size = 5, int_len = 1)
    if (p > .99) ">.99" else if (p < .15) paste0(disp, "*") else paste0(disp, " ")
  }

  spec <- tplyr_spec(
    cols = "TRTA",
    # Safety population, stated in the spec rather than assumed from raw ADSL.
    pop_data = pop_data(cols = c("TRTA" = "TRT01A"), where = SAFFL == "Y"),
    layers = tplyr_layers(group_count(c("AEBODSYS", "AEDECOD"),
      settings = layer_settings(
        distinct_by = "USUBJID",
        stat_columns = list("n" = f_str("xx (xx.x%)", "distinct_n", "distinct_pct"),
                            "e" = f_str("[x]", "n")),
        limit_data_by = c("AEBODSYS", "AEDECOD"),
        total_row = TRUE, total_row_label = "ANY BODY SYSTEM",
        zero_count_display = "count_only",
        # Preferred terms by descending subject count within each SOC, taken from the
        # High Dose column; the SOC level keeps its own (alphabetical) order.
        order_count_method = "bycount", result_order_var = "distinct_n",
        ordering_cols = "Xanomeline High Dose",
        # Pairwise Fisher's exact p (Placebo vs each active arm) on every SOC/PT/total row;
        # ae_p supplies the "* / >.99 / blank" display verbatim (character return).
        assoc_test = assoc_test(
          fn = ae_p, format = f_str("xxxxxx", "p"),
          reference = "Placebo",
          comparisons = c("Xanomeline Low Dose", "Xanomeline High Dose"),
          total_row = TRUE))))
  )
  b <- tplyr_build(spec, adae, pop_data = adsl)
  # Header (N=) per arm, read from the build's own header N so the displayed denominator is
  # provably the one tplyr2 used for the percentages (rather than a second count of ADSL).
  hn <- tplyr_header_n(b)
  N <- unname(setNames(hn$.n, as.character(hn[[1]]))[ARMS])
  bd <- as_display(b)   # display-ready frame: rowlabel*/res1..res6, row order preserved (ord/row_id dropped)

  #' Parse the leading integer count out of an "n (xx.x%)" cell
  #'
  #' Only used to blank the `[AEs]` event column where an arm has no subjects; the row
  #' ordering is the layer's own (`order_count_method = "bycount"`), not derived from
  #' these parsed values.
  #' @param s Character vector of formatted count cells
  #' @return Integer vector of leading counts (NA when absent)
  lead_int <- function(s) suppressWarnings(as.integer(sub("^\\s*([0-9]+).*$", "\\1", s)))
  d <- tibble(
    soc = as.character(bd$rowlabel1), pt = as.character(bd$rowlabel2), depth = b$ord_layer_2,
    npct_0  = as.character(bd$res1), e_0  = as.character(bd$res2),
    npct_54 = as.character(bd$res3), e_54 = as.character(bd$res4),
    npct_81 = as.character(bd$res5), e_81 = as.character(bd$res6),
    p_low = as.character(bd$pval1), p_high = as.character(bd$pval2)   # pairwise Fisher via assoc_test
  ) |>
    mutate(n0 = lead_int(npct_0), n54 = lead_int(npct_54), n81 = lead_int(npct_81)) |>
    # The layer already orders SOCs alphabetically and preferred terms by descending
    # High-Dose count; only the grand-total row has to move to the front of the table.
    mutate(soc_rank = if_else(depth == 0, 0L, 1L)) |>
    arrange(soc_rank) |>            # stable, so the layer's own ordering survives

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

  #' Build a wrapped arm column header
  #' @param arm Arm name
  #' @param n Population count
  #' @return Character scalar column label
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
