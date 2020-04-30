# t-14-3-09.R
#   CDISC Pilot Table 14-3.10

library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(huxtable)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in the ADAS datasets and filter ----
adas <- read_xpt(glue("{adam_lib}/adadas.xpt")) %>%
  filter(EFFFL == "Y" & PARAMCD == "ACTOT" & ITTFL == "Y" & AVISITN %in% c(0, 8, 16, 24) & ANL01FL=="Y") %>%
  mutate(SET = "LOCF") %>%
  select(TRTPN, TRTP, AVISIT, AVISITN, AVAL, BASE, CHG, DTYPE, SET)

# Dataframe to merge display order of visits
visits <- tibble(
  ORD = rep(c(1:8), 3),
  AVISIT = rep(c("Baseline", "Week 8 (Windowed)", "Week 16 (Windowed)", "Week 24 (Windowed)",
             "Week 8 LOCF", "Week 16 LOCF", "Week 24 LOCF", ''), 3),
  TRTPN = c(rep(c(0), 8), rep(54, 8), rep(81, 8))
)

# Join the LOCF and Windowed sets together
step1 <- adas %>%
  bind_rows(adas %>%
    filter(AVISITN != 0 & DTYPE != 'LOCF') %>%
    mutate(SET="WINDOWED")
  ) %>%

  mutate(
    # Format AVISIT for display
    AVISIT =
      case_when(
        SET == 'WINDOWED' & AVISITN != 0 ~ paste(AVISIT, '(Windowed)'),
        SET == 'LOCF' & AVISITN != 0 ~ paste(AVISIT, 'LOCF'),
        TRUE ~ AVISIT
      ),
    # Display of TRTP
    TRTP =
      case_when(
        TRTPN == 0 ~ 'Placebo',
        TRTPN == 54 ~ 'Xan.Low',
        TRTPN == 81 ~ 'Xan.High'
      )
  )

# Get all summaries for aval
aval <- step1 %>%
  group_by(TRTPN, TRTP, AVISITN, AVISIT, SET) %>%
  summarize(
    n = num_fmt(n(), int_len=2, size=2),
    mean = num_fmt(mean(AVAL), digits=1, int_len=2, size=4),
    sd = num_fmt(sd(AVAL), digits=2, int_len=2, size=5),
    md = num_fmt(median(AVAL), digits=1, int_len=2, size=4),
    mn = num_fmt(min(AVAL), int_len=2, size=4),
    mx = num_fmt(max(AVAL), int_len=2, size=4),
  ) %>%
  full_join(visits, by=c("TRTPN", "AVISIT"))

# Get all summaries for chg
chg <- step1 %>%
  group_by(TRTPN, TRTP, AVISITN, AVISIT, SET) %>%
  filter(AVISITN != 0) %>%
  summarize(
    meanc = num_fmt(mean(CHG), digits=1, int_len=1, size=4),
    sdc = num_fmt(sd(CHG), digits=2, int_len=1, size=4),
    mdc = num_fmt(median(CHG), digits=1, int_len=1, size=4),
    mnc = num_fmt(min(CHG), int_len=3, size=4),
    mxc = num_fmt(max(CHG), int_len=2, size=4),
    bmn = num_fmt(mean(BASE), digits=1, int_len=2, size=4),
    bsd = num_fmt(sd(BASE), digits=2, int_len=2, size=5)
  )

# Join to AVAL and CHG results together to create final table
final <- left_join(aval, chg, by=c('TRTPN', 'TRTP', 'AVISITN', 'AVISIT', 'SET')) %>%
  arrange(TRTPN, ORD) %>%
  ungroup() %>%
  mutate(TRTP = ifelse(ORD==1, TRTP, '')) %>%
  select(TRTP, AVISIT, n, mean, sd, md, mn, mx, bmn, bsd, meanc, sdc, mdc, mnc, mxc)

# Create the column headers
header <- tibble(
  TRTP=character(2),
  AVISIT=character(2),
  n=c('', 'nc'),
  mean=c('', 'Mean'),
  sd=c('', 'Std'),
  md=c('', 'Med.'),
  mn=c('', 'Min.'),
  mx=c('', 'Max.'),
  bmn=c('', 'Bsln\\line Mean'),
  bsd=c('', 'Bsln\\line Std'),
  meanc=c('---Change from baseline---', 'Mean'),
  sdc=c('', 'Std'),
  mdc=c('', 'Med.'),
  mnc=c('', 'Min.'),
  mxc=c('', 'Max.')
  )

# Make the table
ht <- as_hux(bind_rows(header, final)) %>%
  huxtable::merge_cells(1, 11:15) %>% # Span header for Change from Baseline
  huxtable::set_bold(1:2, 1:ncol(final), TRUE) %>% # Bold the header
  huxtable::set_align(1:2, 1:ncol(final), 'center') %>% # Align the header
  huxtable::set_align(1:(nrow(final) +2), 3:ncol(final), 'center') %>% # Center all of the numeric cells
  huxtable::set_valign(1:2, 1:ncol(final), 'bottom') %>% # Attach the column headers to the bottom of the cell
  huxtable::set_bottom_border(2, 1:ncol(final), 1) %>% # Bottom border under column header
  huxtable::set_width(1.5) %>% # Take up the whole width of the page
  huxtable::set_escape_contents(FALSE) %>% # Allow RTF strings
  huxtable::set_col_width(c(.09, .19, .05, .06, .06, .06, .05, .05, .06, .06, .06, .05, .05, .05, .06)) # Column widths as a ratio

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht, header_rows=2) %>% pharmaRTF:::titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-3.10') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-3.10.rtf')



