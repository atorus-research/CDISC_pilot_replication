# t-14-1-01.R
#   CDISC Pilot Table 14-1.01

library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(huxtable)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in the ADLB datasets
adsl <- read_xpt(glue("{adam_lib}/adsl.xpt"))

# Create the total values upfront for quicker summary ----
adsl_ <- adsl %>%
  union(adsl %>%
          mutate(TRT01P = 'Total',
                 TRT01PN = 99)) %>%
  mutate(
    COMPL = ifelse(DCDECOD == "COMPLETED", "Y", "N")
  )

# Calculate the header Ns
header_n <- get_header_n(adsl_)

# Column headers
column_headers <- header_n %>% select(TRT01PN, labels) %>%
  pivot_wider(names_from = TRT01PN, values_from = labels)

# Intent to treat
itt <- sum_subgrp(adsl_, ITTFL, order_var=STUDYID, include.n=FALSE, header_n = header_n) %>%
  mutate(rowlbl1 = "Intent-To-Treat (ITT)")

# Safety
safety <- sum_subgrp(adsl_, SAFFL, order_var=STUDYID, include.n=FALSE, header_n = header_n) %>%
  mutate(rowlbl1 = "Safety")

# Efficacy
efficacy <- sum_subgrp(adsl_, EFFFL, order_var=STUDYID, include.n=FALSE, header_n = header_n) %>%
  mutate(rowlbl1 = "Efficacy")

# Commpleters Week 24
compl_24 <- sum_subgrp(adsl_, COMP24FL, order_var=STUDYID, include.n=FALSE, header_n = header_n) %>%
  mutate(rowlbl1 = "Complete Week 24")

# Study completers
compl <- sum_subgrp(adsl_, COMPL, order_var=STUDYID, include.n=FALSE, header_n = header_n) %>%
  mutate(rowlbl1 = "Complete Study")

# Pull the body together
body <- rbind(itt, safety, efficacy, compl_24, compl) %>%
           filter(rowlbl2 == "Y") %>%
           select(-rowlbl2)

# Cleanup
rm(itt, safety, efficacy, compl_24, compl)

# Attach the header
final <- bind_rows(column_headers, body) %>%
  select(rowlbl1, `0`, `54`, `81`, `99`)

# Make the table
ht <- as_hux(final, add_colnames = FALSE) %>%
  huxtable::set_bold(1, 1:ncol(final), TRUE) %>%
  huxtable::set_align(1, 1:ncol(final), 'center') %>%
  huxtable::set_valign(1, 1:ncol(final), 'bottom') %>%
  huxtable::set_bottom_border(1, 1:ncol(final), 1) %>%
  huxtable::set_width(1.1) %>%
  huxtable::set_escape_contents(FALSE) %>%
  huxtable::set_col_width(c(.4, .15, .15, .15, .15)) %>%
  huxtable::set_wrap(TRUE)

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-1.01') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-1.01.rtf')