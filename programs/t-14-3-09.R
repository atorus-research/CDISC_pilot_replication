# t-14-3-09.R
#   CDISC Pilot Table 14-3.09

library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(huxtable)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in the ADLB datasets ----
adas <- read_xpt(glue("{adam_lib}/adadas.xpt")) %>%
  filter(EFFFL == "Y" & PARAMCD == 'ACTOT' & ANL01FL == 'Y' & SEX == "F")

# Calculate the header Ns ----
header_n <- adas %>%
  distinct(USUBJID, TRTP, TRTPN) %>%
  get_header_n(TRTP, TRTPN)

column_headers <- header_n %>%
  select(-N) %>%
  pivot_wider(names_from = TRTPN, values_from=labels) %>%
  mutate(rowlbl1 = '')

# Run each group
summary_portion <- bind_rows(summary_data(adas, AVAL, 0,  'Baseline'),
                             summary_data(adas, AVAL, 24, 'Week 24'),
                             summary_data(adas, CHG,  24, 'Change from Baseline')) %>%
  pad_row()

# Gather the model data
model_portion <- efficacy_models(adas, 'CHG', 24)

final <- bind_rows(column_headers, summary_portion, model_portion) %>%
  select(rowlbl1, `0`, `54`, `81`)

# Make the table
ht <- as_hux(final) %>%
  huxtable::set_bold(1, 1:ncol(final), TRUE) %>%
  huxtable::set_align(1, 1:ncol(final), 'center') %>%
  huxtable::set_valign(1, 1:ncol(final), 'bottom') %>%
  huxtable::set_bottom_border(1, 1:ncol(final), 1) %>%
  huxtable::set_width(1.2) %>%
  huxtable::set_escape_contents(FALSE) %>%
  huxtable::set_col_width(c(.5, 1/6, 1/6, 1/6))

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht) %>% pharmaRTF:::titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-3.09') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-3.09.rtf')



