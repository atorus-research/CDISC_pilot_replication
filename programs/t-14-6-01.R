# t-14-6-01.R
#   CDISC Pilot Table 14-4.01

library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in the ADLB datasets
adlbc <- read_xpt(glue("{adam_lib}/adlbc.xpt")) %>%
  filter(SAFFL == 'Y' & (AVISITN != 99 | (AVISITN == 99 & AENTMTFL=='Y')))

adlbc$PARAM<- recode(adlbc$PARAM,
                     "Alanine Aminotransferase (U/L)" = "ALANINE AMINOTRANSFERASE",
                     "Albumin (g/L)" = "ALBUMIN",
                     "Alkaline Phosphatase (U/L)" = "ALKALINE PHOSPHATASE",
                     "Aspartate Aminotransferase (U/L)" = "ASPARTATE AMINOTRANSFERASE",
                     "Bilirubin (umol/L)" = "BILIRUBIN",
                     "Calcium (mmol/L)" = "CALCIUM",
                     "Chloride (mmol/L)" = "CHLORIDE",
                     "Cholesterol (mmol/L)" = "CHOLESTEROL",
                     "Creatine Kinase (U/L)" = "CREATINE KINASE",
                     "Creatinine (umol/L)" = "CREATININE",
                     "Gamma Glutamyl Transferase (U/L)" = "GAMMA GLUTAMYL TRANSFERASE",
                     "Glucose (mmol/L)" = "GLUCOSE",
                     "Phosphate (mmol/L)" = "PHOSPHATE",
                     "Potassium (mmol/L)" = "POTASSIUM",
                     "Protein (g/L)" = "PROTEIN",
                     "Sodium (mmol/L)" = "SODIUM",
                     "Urate (umol/L)" = "URATE",
                     "Blood Urea Nitrogen (mmol/L)" = "UREA NITROGEN")

adlbh <- read_xpt(glue("{adam_lib}/adlbh.xpt")) %>%
  filter(SAFFL == 'Y' & !(PARAM %in% c('Anisocytes', 'Poikilocytes', 'Microcytes', 'Macrocytes', 'Polychromasia'))
         & (AVISITN != 99 | (AVISITN == 99 & AENTMTFL=='Y')))

adlbh$PARAM<- recode(adlbh$PARAM,
                     "Basophils (GI/L)" = "BASOPHILS",
                     "Eosinophils (GI/L)" = "EOSINOPHILS",
                     "Ery. Mean Corpuscular HGB Concentration (mmol/L)" = "ERY. MEAN CORPUSCULAR HB CONCENTRATION",
                     "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))" = "ERY. MEAN CORPUSCULAR HEMOGLOBIN",
                     "Ery. Mean Corpuscular Volume (fL)" = "ERY. MEAN CORPUSCULAR VOLUME",
                     "Erythrocytes (TI/L)" = "ERYTHROCYTES",
                     "Hematocrit" = "HEMATOCRIT",
                     "Hemoglobin (mmol/L)" = "HEMOGLOBIN",
                     "Leukocytes (GI/L)" = "LEUKOCYTES",
                     "Lymphocytes (GI/L)" = "LYMPHOCYTES",
                     "Monocytes (GI/L)" = "MONOCYTES",
                     "Platelet (GI/L)" = "PLATELET")

# Template for assigning display visit values
visit_names <- data.frame(
  AVISITN = c(0, 2, 4, 6, 8, 12, 16, 20, 24, 26, 99),
  AVISIT = c("  Bsln", "  Wk 2", "  Wk 4", "  Wk 6", "  Wk 8", "  Wk 12",
            "  Wk 16", "  Wk 20", "  Wk 24", "  Wk 26", "  End[1]"),
  stringsAsFactors = FALSE
)

test_summary <- function(x, df_=NULL) {
  # Build up the visit table and attach on the end visit (using flag)
  visits <- df_ %>%
    # Filter to the specified test
    filter(AVISIT != 'UNSCHEDULED' & PARAM == x)

  # Summarize results by visit and treatment
  res <- visits %>%
    group_by(AVISITN, TRTPN) %>%
    summarize(n = n(),
              mean_res = mean(AVAL, na.rm=TRUE),
              sd_res = sd(AVAL, na.rm=TRUE))

  # Summarize change from baseline by visit and treatment
  chgbl <- visits %>%
    filter(AVISITN != 1) %>%
    group_by(AVISITN, TRTPN) %>%
    summarize(mean_cbl = mean(CHG, na.rm=TRUE),
              sd_cbl = sd(CHG, na.rm=TRUE))

  # Build the display string
  df <- merge(res, chgbl, by = c('AVISITN', 'TRTPN'), all=TRUE) %>%
    rowwise() %>%
    mutate(
      N =
        ifelse(
          !is.na(n),
          num_fmt(n, size=2, int_len=2),
          ''),
      msr =
        ifelse(
          !is.na(mean_res),
          as.character(glue('{num_fmt(mean_res, size=5, digits=1, int_len=3)} ({num_fmt(sd_res, size=6, digits=2, int_len=3)})')),
          ''),
      msc =
        ifelse(
          !is.na(mean_cbl),
          as.character(glue('{num_fmt(mean_cbl, size=5, digits=1, int_len=3)} ({num_fmt(sd_cbl, size=6, digits=2, int_len=3)})')),
          '')
    ) %>%
    # Transpose the treatments out
    select(AVISITN, TRTPN, N, msr, msc) %>%
    pivot_wider(names_from = TRTPN, values_from = c(N, msr, msc)) %>%
    # Merge in the visits
    merge(visit_names, by='AVISITN') %>%
    arrange(AVISITN) %>%
    select(AVISIT, N_0, msr_0, msc_0, N_54, msr_54, msc_54, N_81, msr_81, msc_81) %>%
    pad_row()

  # Stub header
  stub_head = data.frame(AVISIT = x, stringsAsFactors = FALSE)

  final <- bind_rows(stub_head, df)
  ht <- huxtable::as_hux(final) %>%
    huxtable::merge_cells(1, 1:5)
  ht

}

add_group_head <- function(ht, group) {
  # Make a three row subset to grab names
  head_ <- ht[1:3, ]
  # Blank everything out
  head_[,] <- ''
  # First value is the group label
  head_[1, 1] <- group
  # Merge the cells
  head_ <- huxtable::merge_cells(head_, 1, 1:5)
  # Bind to the table
  rbind(head_, ht)
}

# Summarize all the chemistry data
chem <- do.call(rbind, lapply(sort(unique(adlbc$PARAM)), test_summary, df_=adlbc)) %>%
  add_group_head('CHEMISTRY')

# Summarize all the hematology data
hema <- do.call(rbind, lapply(sort(unique(adlbh$PARAM)), test_summary, df_=adlbh)) %>%
  add_group_head('HEMATOLOGY')

# Bind those two
ht <- rbind(chem, hema)

# Make the column headers
col_headers <- ht[5:6, ] # Stealing out a chunk of the table with no cell merging
col_headers[1, ] <- c('', 'Placebo', '', '', 'Xanomeline Low', '', '', 'Xanomeline High', '', '')
col_headers[2, ] <- c('Visit', 'N', 'Mean (SD)', 'Change\\line from Bsln\\line Mean (SD)',
                               'N', 'Mean (SD)', 'Change\\line from Bsln\\line Mean (SD)',
                               'N', 'Mean (SD)', 'Change\\line from Bsln\\line Mean (SD)')

# Now
col_headers <- col_headers %>%
  # Placebo spanner
  huxtable::merge_cells(1, 2:4) %>%
  huxtable::set_bottom_border(1, 2:4, 1) %>%
  huxtable::set_bottom_border_style(1, 2:4, 'dashed') %>%
  # Xanomeline Low spanner
  huxtable::merge_cells(1, 5:7) %>%
  huxtable::set_bottom_border(1, 5:7, 1) %>%
  huxtable::set_bottom_border_style(1, 5:7, 'dashed') %>%
  # Xanomeline High spanner
  huxtable::merge_cells(1, 8:10) %>%
  huxtable::set_bottom_border(1, 8:10, 1) %>%
  huxtable::set_bottom_border_style(1, 8:10, 'dashed') %>%
  # Bottom border
  huxtable::set_bottom_border(2, 1:10, value=1) %>%
  # bold it all
  huxtable::set_bold(value=TRUE) %>%
  huxtable::set_align(value='center') %>%
  huxtable::set_valign(value='bottom')

final <- rbind(col_headers, ht) %>%
  huxtable::set_width(1.5) %>%
  huxtable::set_escape_contents(FALSE) %>%
  huxtable::set_col_width(1:10, value=c(.1, .03, .14, .14, .03, .14, .14, .03, .14, .14)) %>%
  huxtable::set_bottom_padding(0) %>%
  huxtable::set_top_padding(0)

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(final, header_rows = 2) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-6.01') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-6.01.rtf')



