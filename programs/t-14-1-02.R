### Table 14-1.02 Summary of End of Study Data

library(dplyr)
library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(pharmaRTF)

source('./programs/config.R')
source('./programs/funcs.R')

#Read in Source and order factors
adsl <- read_xpt(glue("{adam_lib}/adsl.xpt"))
adsl$COMP24FL <- ordered(adsl$COMP24FL, c("Y", "N", NA))
adsl$ARM <- ordered(adsl$ARM, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))
adsl$DCREASCD <- ordered(adsl$DCSREAS, c("Adverse Event",
                                          "Death",
                                          "Lack of Efficacy",
                                          "Lost to Follow-up",
                                          "Withdrew Consent",
                                          "Physician Decision",
                                          "I/E Not Met",
                                          "Protocol Violation",
                                          "Sponsor Decision"))

#### Completion Status Table
comp_stat <- adsl %>%
  group_by(COMP24FL, ARM) %>%
  summarise(n = n())

#Make data.frame for table, unnamed so the cols are named correctly
comp_df <- data.frame(
  "Placebo" = n_pct(unlist(comp_stat[c(1,4), "n"]), sum(unlist(comp_stat[c(1,4), "n"])), mark_lt=FALSE),
  "Xanomeline Low Dose" = n_pct(unlist(comp_stat[c(2,5), "n"]), sum(unlist(comp_stat[c(2,5), "n"])), mark_lt=FALSE),
  "Xanomeline High Dose" = n_pct(unlist(comp_stat[c(3,6), "n"]), sum(unlist(comp_stat[c(3,6), "n"])), mark_lt=FALSE),
  "Total" = c(n_pct(sum(comp_stat[1:3, "n"]), sum(comp_stat[,"n"]), mark_lt=FALSE),
              n_pct(sum(comp_stat[4:6, "n"]), sum(comp_stat[,"n"]), mark_lt=FALSE)),
  row.names = c("\tCompleted Week 24", "\tEarly Termination (prior to Week 24)"),
  #Stop data.frame from adding periods
  check.names = FALSE, stringsAsFactors = FALSE
)
# Add tabs to row.names

# Add missing row.
comp_df["\tMissing", ] <- "  0 (  0%)"

# p-value
comp_p <- fish_p(adsl, adsl$COMP24FL, adsl$ARM)
comp_df <- attach_p(comp_df, comp_p)

#### Reason for Early Termination Table
## By ARM
term_reas <- adsl %>%
  filter(COMP24FL == "N") %>%
  group_by(DCREASCD, ARM) %>%
  complete(nesting(DCREASCD, ARM)) %>%
  summarise(n = n())

## Total
term_reas_tot <- adsl %>%
  filter(COMP24FL == "N", !is.na(DCDECOD)) %>%
  group_by(DCREASCD) %>%
  complete(nesting(DCREASCD, ARM)) %>%
  summarise(n = n())


term_df <- data.frame(
  "Placebo" = n_pct(unlist(term_reas[seq(1, 27, 3), "n"]), sum(adsl %>% filter(ARM == "Placebo") %>% summarise(n = n())), mark_lt=FALSE),
  "Xanomeline Low Dose" = n_pct(unlist(term_reas[seq(2, 27, 3), "n"]), sum(adsl %>% filter(ARM == "Xanomeline Low Dose") %>% summarise(n = n())), mark_lt=FALSE),
  "Xanomeline High Dose" = n_pct(unlist(term_reas[seq(3, 27, 3), "n"]), sum(adsl %>% filter(ARM == "Xanomeline High Dose") %>% summarise(n = n())), mark_lt=FALSE),
  "Total" = n_pct(unlist(term_reas_tot[, "n"]), sum(adsl %>% summarise(n = n())), mark_lt=FALSE),
  row.names = c(
    "\tAdverse Event",
    "\tDeath",
    "\tLack of Efficacy[2]",
    "\tLost to Follow-up",
    "\tSubject decided to withdraw",
    "\tPhysician decided to withdraw subject",
    "\tProtocol criteria not met",
    "\tProtocol violation",
    "\tSponsor decision"
  ),
  #Stop data.frame from adding periods
  check.names = FALSE, stringsAsFactors = FALSE
)
term_df["\tMissing", ] <- "  0 (  0%)"

# p-value
term_p_1 <- adsl %>%
  select(ARM, DCREASCD) %>%
  mutate(loefl = ifelse(DCREASCD %in% "Adverse Event", 1, 0)) %>%
  fish_p(loefl, ARM, width = 6)
term_df <- attach_p(term_df, term_p_1)

term_p_2 <- adsl %>%
  select(ARM, DCREASCD) %>%
  mutate(loefl = ifelse(DCREASCD %in% "Lack of Efficacy", 1, 0)) %>%
  fish_p(ARM ,loefl, width = 6)
term_df["\tLack of Efficacy[2]",] <- attach_p(term_df[3,], term_p_2)


## Add Table lables
comp_df <- add_column(comp_df, " " = row.names(comp_df), .before = 1)
comp_df <- add_row(comp_df, " " = "Completion Status:", .before = 1)
comp_df <- add_row(comp_df, " " = "", .before = 1)

term_df <- add_column(term_df, " " = row.names(term_df), .before = 1)
term_df <- add_row(term_df, " " = "Reason for Early Termination (prior to Week 24):", .before = 1)
term_df <- add_row(term_df, " " = "", .before = 1)

combinedTable <- rbind(comp_df, term_df)
# Rename to get rid of period seperation
names(combinedTable)

headers <- adsl %>%
  group_by(ARM) %>%
  summarise(N = n())
headers_2 <- adsl %>%
  summarise(N = n()) %>%
  mutate(ARM = "Total")
headers_3 <- rbind(headers, headers_2) %>%
  mutate(labels = str_replace_all(str_wrap(glue('{ARM} (N={N})'), width=10), "\n", function(x) "\\line "))
headers_4 <- c(" ", headers_3$labels, "p-value [1]")
names(combinedTable) <- headers_4

ht <- combinedTable %>%
  huxtable::as_hux(add_colnames=TRUE) %>%
  huxtable::set_wrap(FALSE)

huxtable::bottom_border(ht)[1, ] <- 1
huxtable::bold(ht)[1, ] <- TRUE
huxtable::align(ht)[1, ] <- 'center'
huxtable::align(ht)[, 6] <- "center"
huxtable::width(ht) <- 1.5
huxtable::escape_contents(ht) <- FALSE
huxtable::col_width(ht) <- c(.4, .12, .12, .12, .12, .12)
huxtable::bottom_padding(ht) <- 0
huxtable::top_padding(ht) <- 0
huxtable::valign(ht)[1,] <- "bottom"
ht[8,2] <- ""
ht <- huxtable::merge_cells(ht, 8, 1:2)


# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-1.02') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top = 1)

# Write out the RTF
write_rtf(doc, file='./outputs/14-1.02.rtf')

