# t-14-2-01.R
#   CDISC Pilot Table 14-2.01

library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Import and explore the data frame ----
adsl <- read_xpt(glue("{adam_lib}/adsl.xpt")) %>%
  filter(ITTFL == "Y") %>%
  mutate(
    RACE_DISPLAY = case_when(
      ETHNIC == 'HISPANIC OR LATINO' ~ 'Hispanic',
      RACE == 'WHITE' ~ 'Caucasian',
      RACE == 'BLACK OR AFRICAN AMERICAN' ~ 'African Descent',
      RACE == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 'Other',
    ),
    RACEN_DISPLAY = case_when(
      ETHNIC == 'HISPANIC OR LATINO' ~ 3,
      RACE == 'WHITE' ~ 1,
      RACE == 'BLACK OR AFRICAN AMERICAN' ~ 2,
      RACE == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 4,
    ),
    SEX =
      case_when(
        SEX == 'M' ~ 'Male',
        SEX == 'F' ~ 'Female'
      ),
    SEXN =
      case_when(
        SEX == 'Male' ~ 1,
        SEX == 'Female' ~ 2
      ),
    DURDSGR1N =
      case_when(
        DURDSGR1 == '<12' ~ 1,
        DURDSGR1 == '>=12' ~ 2
      ),
    DURDSGR1 = paste(DURDSGR1, 'months'),
    BMIBLGR1N =
      case_when(
        BMIBLGR1 == '<25' ~ 1,
        BMIBLGR1 == '25-<30' ~ 2,
        BMIBLGR1 == '>=30' ~ 3
      ),
    AGEGR1 = paste(AGEGR1, 'yrs')
  )

get_meta(adsl)

# Create the total values upfront for quicker summary ----
adsl_ <- adsl %>%
  bind_rows(adsl %>%
              mutate(TRT01P = 'Total',
                     TRT01PN = 99))


# Get the header N's ----
header_n <- get_header_n(adsl_)

## Exploring Age ----

# Descriptive stats
age_1 <- adsl_ %>% desc_stats(AGE)
age_p <- adsl %>% aov_p(AGE ~ TRT01P) # anova

age_1 <- attach_p(age_1, age_p)

# Categorical n counts
age_2 <- adsl_ %>% sum_subgrp(AGEGR1, AGEGR1N, include.n=FALSE, header_n=header_n)

agegrp_p <- adsl %>% chi_p(AGEGR1, TRT01P)
age_2 <- attach_p(age_2, agegrp_p)

age <- rbind(age_1, age_2) %>%
  mutate(rowlbl1 = "Age (y)")

rm(age_1, age_2, age_p, agegrp_p)

## Exploring sex ----
sex = adsl_ %>%
  sum_subgrp(SEX, SEXN, header_n=header_n)

sex_p <- adsl %>% chi_p(SEX, TRT01P)

sex <- attach_p(sex, sex_p) %>%
  mutate(rowlbl1 = 'Sex')

rm(sex_p)

## Exploring race ----
race = adsl_ %>%
  sum_subgrp(RACE_DISPLAY, RACEN_DISPLAY, header_n=header_n) %>%
  rowwise() %>%
  mutate(
    rowlbl1 = "Race (Origin)",
  )

race_p <- adsl %>% chi_p(RACE_DISPLAY, TRT01P)

race <- attach_p(race, race_p)

rm(race_p)

## Exploring MMSE ---
mmse <- adsl_ %>% desc_stats(MMSETOT) %>%
  mutate(
    rowlbl1 = 'MMSE'
  )

mmse_p <- adsl %>% aov_p(MMSETOT ~ TRT01P)

mmse <- attach_p(mmse, mmse_p)

rm(mmse_p)

## Exploring disease duration ----

# Descriptive
durdis_1 <- adsl_ %>% desc_stats(DURDIS)
durdis_1p <- adsl %>% aov_p(DURDIS ~ TRT01P)
durdis_1 <- attach_p(durdis_1, durdis_1p)

# Categorical
durdis_2 <- adsl_ %>% sum_subgrp(DURDSGR1, DURDSGR1N, include.n=FALSE, header_n=header_n)
durdis_2p <- adsl %>% chi_p(DURDSGR1, TRT01P)
durdis_2 <- attach_p(durdis_2, durdis_2p)

durdis <- durdis_1 %>%
  union(durdis_2) %>%
  mutate(
    rowlbl1 = 'Duration of disease '
  ) %>%
  pad_row()

rm(durdis_1, durdis_2, durdis_1p, durdis_2p)

## Years of education ----
educlvl <- adsl_ %>% desc_stats(EDUCLVL) %>%
  mutate(
    rowlbl1 = 'Years of education'
  )
educlvl_p <- adsl %>% aov_p(EDUCLVL ~ TRT01P)
educlvl <- attach_p(educlvl, educlvl_p)

rm(educlvl_p)

## Baseline weight ----
weightbl <- adsl_ %>% desc_stats(WEIGHTBL) %>%
  mutate(
    rowlbl1 = 'Baseline weight(kg)'
  )
weightbl_p <- adsl %>% aov_p(WEIGHTBL ~ TRT01P)
weightbl <- attach_p(weightbl, weightbl_p)

rm(weightbl_p)

## Baseline height ----
heightbl <- adsl_ %>% desc_stats(HEIGHTBL) %>%
  mutate(
    rowlbl1 = 'Baseline height(cm)'
  )
heightbl_p <- adsl %>% aov_p(HEIGHTBL ~ TRT01P)
heightbl <- attach_p(heightbl, heightbl_p)

rm(heightbl_p)

## Baseline BMI ----

# Descriptive
bmi_1 <- adsl_ %>% desc_stats(BMIBL)
bmi_1p <- adsl %>% aov_p(BMIBL ~ TRT01P)
bmi_1 <- attach_p(bmi_1, bmi_1p)

# Categorical
bmi_2 <- adsl_ %>% sum_subgrp(BMIBLGR1, BMIBLGR1N, include.n=FALSE, header_n=header_n)
bmi_2p <- adsl %>% chi_p(BMIBLGR1, TRT01P)
bmi_2 <- attach_p(bmi_2, bmi_2p)

bmi <- rbind(bmi_1, bmi_2) %>%
  mutate(
    rowlbl1 = 'Baseline BMI'
  )

rm(bmi_1, bmi_2, bmi_1p, bmi_2p)

## Stack together final tables ---
final <- rbind(age, sex, race, mmse, durdis, educlvl, weightbl, heightbl, bmi) %>%
  group_by(rowlbl1) %>%
  mutate(ord1 = row_number()) %>%
  ungroup() %>%
  mutate(rowlbl1 = ifelse(ord1 == 1, rowlbl1, ""))

rm(age, sex, race, mmse, durdis, educlvl, weightbl, heightbl, bmi)

# Make and attach column headers
header_n_v <- header_n %>% select(TRT01PN, labels) %>%
  pivot_wider(names_from = TRT01PN, values_from = labels) %>%
  mutate(
    rowlbl1 = '',
    rowlbl2 = '',
    p = 'p-value\\line [1]'
  )

final <- bind_rows(header_n_v, final) %>%
  select(rowlbl1, rowlbl2, `0`, `54`, `81`, `99`, p)

## Table build and RTF Output
ht <- huxtable::as_hux(final, add_colnames = FALSE)

# ht <- as_hux(mtcars, add_colnames = TRUE)
huxtable::bottom_border(ht)[1, ] <- 1
huxtable::valign(ht)[1, ] <- 'bottom'
huxtable::bold(ht)[1, ] <- TRUE
huxtable::align(ht)[1, ] <- 'center'
huxtable::width(ht) <- 1.5
huxtable::escape_contents(ht) <- FALSE
huxtable::col_width(ht) <- c(.2, .2, .12, .12, .12, .12, .12)
huxtable::bottom_padding(ht) <- 0
huxtable::top_padding(ht) <- 0

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-2.01') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top=1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-2.01.rtf')
