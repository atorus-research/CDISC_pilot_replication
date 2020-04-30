# t-14-3-12.R
#   CDISC Pilot Table 14-3.12

library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(huxtable)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in the ADLB datasets ----
npix <- read_xpt(glue("{adam_lib}/adnpix.xpt")) %>%
  filter(EFFFL == 'Y' & ITTFL == 'Y' & PARAMCD == 'NPTOTMN') %>%
  mutate(CHG = AVAL - BASE)

# Calculate the header Ns ----
header_n <- npix %>%
  distinct(USUBJID, TRTP, TRTPN) %>%
  get_header_n(TRTP, TRTPN)

column_headers <- header_n %>%
  select(-N) %>%
  pivot_wider(names_from = TRTPN, values_from=labels) %>%
  mutate(rowlbl1 = '')

# Run each group
# NOTE: The counts of Mean of Weeks 4-24 do not match the original data. This was programmed using
#       the derived NPTOTMN variable. Following the original analysis results metadata, as best as
#       I can tell, we're following the same subsetting rules. This means that the counts are a
#       discrepancy between the original analysis data, which is not available in the current
#       CDISC pilot package. The subsequent statistical summaries therefore also do not match.
summary_portion <- bind_rows(summary_data(npix, AVAL, 0 , 'Baseline'),
                             summary_data(npix, AVAL,  98, 'Mean of Weeks 4-24')) %>%
  pad_row()

# Gather the model data
model_portion <- efficacy_models(npix, 'CHG', 98)

final <- bind_rows(column_headers, summary_portion, model_portion) %>%
  select(rowlbl1, `0`, `54`, `81`)

## Create the table

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
doc <- rtf_doc(ht) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-3.12') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-3.12.rtf')



