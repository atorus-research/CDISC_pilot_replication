## Table 14-6.06

library(huxtable)
library(glue)
library(tidyverse)
library(haven)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

pad_row <- function(df, r) {
  #df - dataframe to insert pad
  #r - row number to pad
  for(i in seq(along = r)) {
    if(r[i] + i - 1 < nrow(df)){
      df[seq(r[i] + i, nrow(df) + 1),] <- df[seq(r[i] + (i - 1), nrow(df)),]
      df[r[i] + (i - 1),] <- NA
    } else {
      df[r[i] + (i - 1),] <- NA
    }
  }
  df
}

n_pct <- function(n, pct, n_width=3, pct_width=3) {
  n <- unlist(n)
  pct <- unique(pct)
  # n (%) formatted string. e.g. 50 ( 75%)
  unlist(lapply(n, function(x) {
    if(x == 0) " 0      "
    else {
      as.character(
        # Form the string using glue and format
        glue('{format(x, width=n_width)}({format(round((x/pct) * 100), width=pct_width)}%)')
      )
    }
  }))
}

# Old data used because new data is missing columns
adlbhy <- read_xpt(glue("{adam_lib}/adlbhy.xpt")) %>%
  filter(SAFFL == "Y", PARAMCD %in% c("TRANSHY", "HYLAW"), !is.na(BASE), AVISITN > 0)

adlbhy2 <- adlbhy %>%
  group_by(USUBJID) %>%
  filter(AVAL == max(AVAL)) %>%
  filter(AVISITN == max(AVISITN)) %>%
  ungroup()

adlbhy2$BASE <- ordered(adlbhy2$BASE, c(0, 1))
adlbhy2$AVAL <- ordered(adlbhy2$AVAL, c(0, 1))
adlbhy2$TRTP <- ordered(adlbhy2$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))

total_t <- adlbhy2 %>%
  filter(!is.na(TRTP), !is.na(BASE), PARAMCD == "TRANSHY") %>%
  group_by(TRTP, BASE) %>%
  complete(nesting(TRTP, BASE)) %>%
  summarise(N = n(),
            Nc = num_fmt(n(), size = 2, int_len = 2))
total_tl <- total_t %>%
  mutate("Shift\\line[1]" = ordered("T", c("T", "N", "H"))) %>%
  pivot_wider(id_cols = c("Shift\\line[1]"), names_from = c(TRTP, BASE), values_from = Nc)

adlbhy_t1 <- adlbhy2 %>%
  filter(PARAMCD == "TRANSHY") %>%
  group_by(TRTP, BASE, AVAL) %>%
  complete(nesting(BASE, AVAL)) %>%
  summarise(N = n())
adlbhy_t <- adlbhy_t1 %>%
  mutate(N2 = n_pct(N, total_t[total_t$TRTP == TRTP &
                               total_t$BASE == BASE, "N"], n_width = 2)) %>%
  pivot_wider(id_cols = c("AVAL"), names_from = c("TRTP", "BASE"), values_from = c("N2"))

names(adlbhy_t)[1] <- "Shift\\line[1]"

total_t3 <- total_tl %>%
  rbind(adlbhy_t)

total_t3[, " "] <- "Transaminase 1.5 x ULN"

total_t3[, "p-\\line value\\line[2]"] <- c(
  num_fmt(mantelhaen.test(array(unlist(adlbhy_t1[,"N"]), dim = c(2,3,2)))$p.value, size = 6, int_len = 1, digits = 3)
  , "", "")

adlbhy3 <- adlbhy2 %>%
  filter(PARAMCD %in% c("HYLAW"))

adlbhy3$BASE <- ordered(adlbhy3$BASE, c(0, 1))
adlbhy3$AVAL <- ordered(adlbhy3$AVAL, c(0, 1))
adlbhy3$TRTP <- ordered(adlbhy3$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))

total_b <- adlbhy3 %>%
  filter(!is.na(TRTP), !is.na(BASE)) %>%
  group_by(TRTP, BASE) %>%
  complete(nesting(TRTP, BASE)) %>%
  summarise(N = n(),
            Nc = num_fmt(n(), size = 2, int_len = 2))
total_bl <- total_b %>%
  mutate("Shift\\line[1]" = ordered("T", c("T", "N", "H"))) %>%
  pivot_wider(id_cols = c("Shift\\line[1]"), names_from = c(TRTP, BASE), values_from = Nc)

adlbhy_b1 <- adlbhy3 %>%
  group_by(TRTP, BASE, AVAL) %>%
  complete(nesting(BASE, AVAL)) %>%
  summarise(N = n())
adlbhy_b <- adlbhy_b1 %>%
  mutate(N2 = n_pct(N, total_b[total_b$TRTP == TRTP &
                               total_b$BASE == BASE, "N"], n_width = 2)) %>%
  pivot_wider(id_cols = c("AVAL"), names_from = c("TRTP", "BASE"), values_from = c("N2"))

names(adlbhy_b)[1] <- "Shift\\line[1]"

total_b3 <- total_bl %>%
  rbind(adlbhy_b)
total_b3[, " "] <- "Total Bili 1.5 x ULN and\\line Transaminase 1.5 x ULN"

## FIXME - Different counts???
# total_b3[, "p-\\line value\\line[2]"] <- c(
#   num_fmt(mantelhaen.test(array(unlist(adlbh_b1[,"N"]), dim = c(2,3,2)))$p.value, size = 6, int_len = 1, digits = 3)
#   , "", "")
total_b3[, "p-\\line value\\line[2]"] <- c("", "", "")

## Table construction
# Lots of weird properties for this table so I'm doing it manually
comb <- rbind(total_t3, total_b3)
comb <- comb[, c(8,1,2,3,4,5,6,7,9)]
comb[(comb$`Shift\\line[1]` != "T"), " "] <- ""
comb$`Shift\\line[1]` <- as.character(recode(comb$`Shift\\line[1]`,
                                    "T" = "n",
                                    "N" = "Normal",
                                    "Y" = "High"))
comb2 <- comb %>%
  add_row("Shift\\line[1]" = NA, .before = 1) %>%
  add_row("Shift\\line[1]" = NA, .before = 1)
comb2 <- comb2 %>%
  add_row("Shift\\line[1]" = NA, .before = 6)

names(comb2) <- c(
  " ",
  "Shift\\line[1]",
  "Normal at Baseline",
  "High at Baseline",
  "Normal at Baseline",
  "High at Baseline",
  "Normal at Baseline",
  "High at Baseline",
  "p-\\line value\\line[2]"
)

adsl <- read_xpt(glue("{adam_lib}/adsl.xpt"))
headers <- adsl %>%
  filter(ARM != "Screen Failure") %>%
  group_by(ARM) %>%
  summarise(N = n()) %>%
  mutate(label = paste0(recode(ARM,
                               "Placebo" = "Placebo",
                               "Xanomeline Low Dose" = "Xan. Low",
                               "Xanomeline High Dose" = "Xan. High"), " (N=", N, ")"))

ht <- comb2 %>%
  huxtable::as_hux(add_colnames=TRUE)

ht <- as_hux(pad_row(as.data.frame(ht), c(1,1)), add_colnames = FALSE)
ht[1, 3] <- headers[1, "label"]
ht[1, 5] <- headers[3, "label"]
ht[1, 7] <- headers[2, "label"]

ht2 <- ht %>%
  huxtable::merge_cells(1, 3:4) %>%
  huxtable::merge_cells(1, 5:6) %>%
  huxtable::merge_cells(1, 7:8) %>%
  huxtable::set_bottom_border(2, 3:4, 1) %>%
  huxtable::set_bottom_border(2, 5:6, 1) %>%
  huxtable::set_bottom_border(2, 7:8, 1) %>%
  huxtable::set_bottom_border(3, 1:9, 1)  %>%
  huxtable::set_width(1.5) %>%
  huxtable::set_escape_contents(FALSE) %>%
  huxtable::set_bold(1:3, 1:9, TRUE) %>%
  huxtable::set_valign(1:3, 1:9, "bottom") %>%
  huxtable::set_align(3, 1:9, "center") %>%
  huxtable::set_align(1, 1:9, "center") %>%
  huxtable::set_align(4:12, 9, "right") %>%
  huxtable::set_valign(10, 2:9, "bottom") %>%
  huxtable::set_col_width(1:9, c(0.31, rep(0.09, 7), 0.06))

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht2, header_rows = 3) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-6.06') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top = 1) %>%
  set_header_height(0.75) %>%
  set_footer_height(1)

write_rtf(doc, file='./outputs/14-6.06.rtf')

