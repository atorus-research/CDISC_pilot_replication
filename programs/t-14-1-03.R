# t-14-1-03.R
#   CDISC Pilot Table 14-1.03

library(dplyr)
library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(huxtable)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

# Read in the datasets
adsl <- read_xpt(glue("{adam_lib}/adsl.xpt")) %>%
  filter(ITTFL == "Y")

adsl$SITEGR1 <- ordered(adsl$SITEGR1, c(
  "Pooled\\line Id",
  "701",
  "703",
  "704",
  "705",
  "708",
  "709",
  "710",
  "713",
  "716",
  "718",
  "900",
  "TOTAL"
))
adsl$SITEID <- ordered(adsl$SITEID, c(
  "Site\\line Id",
  "701",
  "703",
  "704",
  "705",
  "708",
  "709",
  "710",
  "713",
  "716",
  "718",
  "702",
  "706",
  "707",
  "711",
  "714",
  "715",
  "717",
  ""
))
adsl$ITTFL <- ordered(adsl$ITTFL, c("Y", "N"))
adsl$EFFFL <- ordered(adsl$EFFFL, c("Y", "N"))
adsl$COMP24FL <- ordered(adsl$COMP24FL, c("Y", "N"))



adsl_grp1 <- adsl %>%
  select(SITEGR1, SITEID, TRT01P, ITTFL) %>%
  group_by(SITEGR1, SITEID, TRT01P, ITTFL) %>%
  filter(ITTFL == "Y") %>%
  summarise(n = n())
adsl_grp1[,4] <- "ITTFL"
names(adsl_grp1)[4] <- "FLFL"
adsl_grp2 <- adsl %>%
  select(SITEGR1, SITEID, TRT01P, EFFFL) %>%
  group_by(SITEGR1, SITEID, TRT01P, EFFFL) %>%
  filter(EFFFL == "Y") %>%
  summarise(n = n())
adsl_grp2[,4] <- "EFFFL"
names(adsl_grp2)[4] <- "FLFL"
adsl_grp3 <- adsl %>%
  select(SITEGR1, SITEID, TRT01P, COMP24FL) %>%
  group_by(SITEGR1, SITEID, TRT01P, COMP24FL) %>%
  filter(COMP24FL == "Y") %>%
  summarise(n = n())
adsl_grp3[,4] <- "COMP24FL"
names(adsl_grp3)[4] <- "FLFL"
adsl_grp4 <- adsl %>%
  select(SITEGR1, SITEID, ITTFL) %>%
  group_by(SITEGR1, SITEID, ITTFL) %>%
  filter(ITTFL == "Y") %>%
  summarise(n = n())
adsl_grp4[,3] <- "ITTFL"
names(adsl_grp4)[3] <- "FLFL"
adsl_grp4$TRT01P <- "Total"
adsl_grp5 <- adsl %>%
  select(SITEGR1, SITEID, EFFFL) %>%
  group_by(SITEGR1, SITEID, EFFFL) %>%
  filter(EFFFL == "Y") %>%
  summarise(n = n())
adsl_grp5[,3] <- "EFFFL"
names(adsl_grp5)[3] <- "FLFL"
adsl_grp5$TRT01P <- "Total"
adsl_grp6 <- adsl %>%
  select(SITEGR1, SITEID, COMP24FL) %>%
  group_by(SITEGR1, SITEID, COMP24FL) %>%
  filter(COMP24FL == "Y") %>%
  summarise(n = n())
adsl_grp6[,3] <- "COMP24FL"
names(adsl_grp6)[3] <- "FLFL"
adsl_grp6$TRT01P <- "Total"

all <- rbind(adsl_grp1, adsl_grp2, adsl_grp3, adsl_grp4, adsl_grp5, adsl_grp6)
all$FLFL <- ordered(all$FLFL, c("ITTFL", "EFFFL", "COMP24FL"))
all$TRT01P <- ordered(all$TRT01P, c(
  "Placebo",
  "Xanomeline Low Dose",
  "Xanomeline High Dose",
  "Total"
))
df <- all %>%
  arrange(SITEGR1, SITEID, TRT01P, FLFL) %>%
  pivot_wider(id_cols = c(SITEGR1, SITEID), names_from = c(TRT01P, FLFL), values_from = c(n), values_fill = list(n = 0)) %>%
  ungroup() %>%
  as.data.frame()

# Stack the total row to the bottom of the data frame
df <-rbind(df, 
          data.frame(
            SITEGR1 = "TOTAL",
            SITEID = "",
            t(apply(df[,3:ncol(df)], 2, sum)), check.names = FALSE
        )
      )

names(df) <- c(
  "Pooled\\line Id",
  "Site\\line Id",
  rep(c("ITT", "Eff", "Com"), 4)
)

df[2:(nrow(df) + 1),] <- df[1:nrow(df),]
df[1,] <- as.list(names(df))
df <- df %>%
  add_row("Pooled\\line Id" = "", .before = 1) %>%
  add_row("Pooled\\line Id" = "", .before = 1)


### Add Headers
headers <- adsl %>%
  group_by(ARM) %>%
  summarise(N = n()) %>%
  mutate(labels = str_replace_all(str_wrap(glue('{ARM} (N={N})'), width=10), "\n", function(x) "\\line "))
headers[4,] <- list(
  ARM = "Total",
  N = nrow(adsl),
  labels = paste0("Total\\line(N=", nrow(adsl), ")")
)

df[1, 3] <- headers[1, "labels"]
df[1, 6] <- headers[3, "labels"]
df[1, 9] <- headers[2, "labels"]
df[1, 12] <- headers[4, "labels"]

ht <- df %>%
  huxtable::as_hux(add_colnames=FALSE) %>%
  merge_cells(1, 3:5) %>%
  set_bottom_border(2, 3:5, 1) %>%
  merge_cells(1, 6:8) %>%
  set_bottom_border(2, 6:8, 1) %>%
  merge_cells(1, 9:11) %>%
  set_bottom_border(2, 9:11, 1) %>%
  merge_cells(1, 12:14) %>%
  set_bottom_border(2, 12:14, 1) %>%
  set_escape_contents(FALSE) %>%
  set_width(1.1) %>%
  set_bottom_border(3, 1:14, 1) %>%
  huxtable::set_align(3, 1:14, "center") %>%
  huxtable::set_valign(3, 1:14, "bottom") %>%
  huxtable::set_align(4:nrow(df), 3:ncol(df), "right") %>%
  huxtable::set_align(1:nrow(df), 1:2, "center") %>%
  huxtable::set_align(1, 1:14, "center") %>%
  huxtable::set_valign(1, 1:14, "bottom") %>%
  huxtable::set_col_width(1:ncol(df), value = c(0.1, 0.1, rep(0.07, 12))) %>%
  huxtable::set_wrap(FALSE)

doc <- rtf_doc(ht, header_rows = 2) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-1.03') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top = 1)

write_rtf(doc, file='./outputs/14-1.03.rtf')

