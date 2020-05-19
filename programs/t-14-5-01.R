# t-14-5-01.R
#   CDISC Pilot Table 14-5.01

library(glue)
library(tidyverse)
library(haven)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in ADSL
adae <- read_xpt(glue("{adam_lib}/adae.xpt")) %>%
  filter(SAFFL == 'Y' & TRTEMFL == 'Y')

adsl <- read_xpt(glue("{adam_lib}/adsl.xpt"))

# Header N ----
header_n <- adsl %>%
  get_header_n()

# Overall counts
overall <- ae_counts(adae, N_counts = header_n) %>%
  mutate(AETERM = 'ANY BODY SYSTEM', AEBODSYS = 'ANY BODY SYSTEM', ord1=1, ord2=1)

# System Organ Class counts
bodsys <- ae_counts(adae, AEBODSYS, N_counts = header_n) %>%
  mutate(AETERM = AEBODSYS, ord1=2, ord2=1) %>%
  arrange(AEBODSYS)

pad <- bodsys %>%
  select(AEBODSYS, ord1, ord2) %>%
  mutate(ord3=999)

# Individual term counts
term <- ae_counts(adae, AEBODSYS, AETERM, sort=TRUE, N_counts = header_n) %>%
  mutate(AETERM = paste0('  ', AETERM), ord1=2, ord2=2)

# Bring the data together
combined <- bind_rows(overall, bodsys, pad, term) %>%
  arrange(ord1, AEBODSYS, ord2, desc(ord3), AETERM)

# Build and attach column headers
column_headers <- header_n %>%
  select(-N) %>%
  pivot_wider(names_from = TRT01PN, values_from=labels) %>%
  select(npct_0=`0`, npct_54=`54`, npct_81=`81`) %>%
  mutate(cAEs_0 = '',
         cAEs_54 = '',
         cAEs_81 = '',
         AETERM = '',
         p_low = "Fisher's Exact\\line p-values",
         p_high = '')

# Insert second row of header
column_headers <- bind_rows(column_headers, tibble(
  AETERM = 'System Organ Class/\\line Preferred Term',
  npct_0 = 'n(%)',
  cAEs_0 = '[AEs]',
  npct_54 = 'n(%)',
  cAEs_54 = '[AEs]',
  npct_81 = 'n(%)',
  cAEs_81 = '[AEs]',
  p_low = 'Placebo\\line vs.\\line Low Dose',
  p_high = 'Placebo\\line vs.\\line High Dose'
))

# Attach to final
final <- bind_rows(column_headers, combined) %>%
  select(AETERM, npct_0, cAEs_0, npct_54, cAEs_54, npct_81, cAEs_81, p_low, p_high)


# Make the table ----

ht <- huxtable::as_hux(final) %>%
  huxtable::merge_cells(1, 2:3) %>%
  huxtable::merge_cells(1, 4:5) %>%
  huxtable::merge_cells(1, 6:7) %>%
  huxtable::merge_cells(1, 8:9)
huxtable::bottom_border(ht)[2, ] <- 1
huxtable::valign(ht)[1:2, ] <- 'bottom'
huxtable::bold(ht)[1:2, ] <- TRUE
huxtable::align(ht)[1:2, ] <- 'center'
huxtable::width(ht) <- 1.5
huxtable::escape_contents(ht) <- FALSE
huxtable::col_width(ht) <- c(.3, .1, .07, .1, .07, .1, .07, .09, .1)
huxtable::bottom_padding(ht) <- 0
huxtable::top_padding(ht) <- 0
# ht

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht, header_rows = 2) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-5.01') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-5.01.rtf')
