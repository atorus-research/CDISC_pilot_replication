# t-14-3-02.R
#   CDISC Pilot Table 14-3.02

library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(huxtable)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in the ADLB datasets ----
cibc <- read_xpt(glue("{adam_lib}/adcibc.xpt")) %>%
  filter(EFFFL == "Y" & ITTFL=='Y' & PARAMCD == 'CIBICVAL' & ANL01FL == 'Y')

# Calculate the header Ns ----
header_n <- cibc %>%
  distinct(USUBJID, TRTP, TRTPN) %>%
  get_header_n(TRTP, TRTPN)

column_headers <- header_n %>%
  select(-N) %>%
  pivot_wider(names_from = TRTPN, values_from=labels) %>%
  mutate(rowlbl1 = '')

# Run each group
summary_portion <- summary_data(cibc, AVAL, 8, 'Week 8') %>%
  pad_row()

# Gather the model data
model_portion <- efficacy_models(cibc, 'AVAL', 8)

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
ht

# Write into doc object and pull titles/footnotes from excel file
## TODO: `titles_and_footnotes_from_df`` should be an exported function so remove internal reference when updated
doc <- rtf_doc(ht) %>% pharmaRTF:::titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-3.04') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-3.04.rtf')



