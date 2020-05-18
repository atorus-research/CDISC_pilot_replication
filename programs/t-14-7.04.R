### Table 14-7.01: Summary of Concomitant Medications (Number of Subjects)

library(dplyr)
library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(pharmaRTF)
library(tibble)

source('./programs/config.R')
source('./programs/funcs.R')

## Modified n_pct function
n_pct <- function(n, pct) {
  # n (%) formatted string. e.g. 50 ( 75%)
  return(
    # Suppress conversion warnings
    as.character(
      # Form the string using glue and format
      glue('{format(n, width=3)} ({format(round((n/pct) * 100))}%)')
    )
  )
}


cm <- read_xpt(glue("{sdtm_lib}/cm.xpt"))
adsl <- read_xpt(glue("{adam_lib}/adsl.xpt"))
adsl$ARM <- ordered(adsl$ARM, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))

## Patients receiving at least one medication
cm_1 <- adsl %>%
  group_by(ARM) %>%
  summarise(n = sum(USUBJID %in% unique(cm$USUBJID)),
            total = n())
cm_res <- n_pct(cm_1$n, cm_1$total)

cm_2 <- data.frame(
  "Therapeutic class, n (%)" = "Patients receiving at least one concomitant medication",
  "Placebo" = cm_res[1],
  "Xanomeline Low Dose" = cm_res[2],
  "Xanomeline High Dose" = cm_res[3],
  stringsAsFactors = FALSE, check.names = FALSE, row.names = FALSE
)

### Table
# Medication classes
cm_class <- sort(unique(cm$CMCLAS))

# By Class
df <- plyr::ldply(cm_class, function(class_i){
  class_by_arm <- as.data.frame(adsl %>%
    group_by(ARM) %>%
    summarise(n = sum(USUBJID %in% unlist(unique(cm[cm$CMCLAS == class_i, "USUBJID"])))))

  df_1 <- data.frame(
    "Therapeutic class, n (%)" = class_i,
    "Placebo" = unname(ifelse(class_by_arm[1, "n"] == 0, "  0", n_pct(class_by_arm[1, "n"], cm_1[1, "total"]))),
    "Xanomeline Low Dose" = unname(ifelse(class_by_arm[2, "n"] == 0, "  0", n_pct(class_by_arm[2, "n"], cm_1[2, "total"]))),
    "Xanomeline High Dose" = unname(ifelse(class_by_arm[3, "n"] == 0, "  0", n_pct(class_by_arm[3, "n"], cm_1[3, "total"]))),
    stringsAsFactors = FALSE, check.names = FALSE, row.names = FALSE
  )

  #Pad Row
  df_1 <- add_row(df_1, "Therapeutic class, n (%)" = "", .before = 1)

  #Coded medication names
  cm_medi <- unlist(unique(cm[cm$CMCLAS == class_i, "CMDECOD"]), use.names = FALSE)

  #By Medication
  df_2 <- plyr::ldply(cm_medi, function(medi_i) {

    medi_by_arm <- as.data.frame(adsl %>%
      group_by(ARM) %>%
      summarise(n = sum(USUBJID %in% unlist(unique(cm[cm$CMDECOD == medi_i, "USUBJID"])))))


    df_3 <- data.frame(
      "Therapeutic class, n (%)" = paste0("\t", unname(medi_i)),
      "Placebo" = unname(ifelse(medi_by_arm[1, "n"] == 0, "  0", n_pct(medi_by_arm[1, "n"], cm_1[1, "total"]))),
      "Xanomeline Low Dose" = unname(ifelse(medi_by_arm[2, "n"] == 0, "  0", n_pct(medi_by_arm[2, "n"], cm_1[2, "total"]))),
      "Xanomeline High Dose" = unname(ifelse(medi_by_arm[3, "n"] == 0, "  0", n_pct(medi_by_arm[3, "n"], cm_1[3, "total"]))),
      stringsAsFactors = FALSE, check.names = FALSE, row.names = FALSE
    )
  })
  ## Order Medications. Order Descending by placebo count and ascending alphabetically
  # radix used because its the only method that supports a decreasing vector.
  df_2 <- df_2[order(df_2$Placebo, df_2[,1], decreasing = c(TRUE, FALSE), method = "radix"),]
  rbind(df_1, df_2)
})


combinedTable <- rbind(cm_2, df)


### Add Headers
headers <- adsl %>%
  group_by(ARM) %>%
  summarise(N = n()) %>%
  mutate(labels = str_replace_all(str_wrap(glue('{ARM} (N={N})'), width=10), "\n", function(x) "\\line "))
names(combinedTable) <- c("Therapeutic class, n (%)", headers$labels)


ht <- combinedTable %>%
  huxtable::as_hux(add_colnames=TRUE)


huxtable::bottom_border(ht)[1, ] <- 1
huxtable::bold(ht)[1, ] <- TRUE
huxtable::align(ht)[1, ] <- 'center'
huxtable::align(ht)[1, 1] <- "left"
huxtable::width(ht) <- 1.2
huxtable::bottom_padding(ht) <- 0
huxtable::top_padding(ht) <- 0
huxtable::col_width(ht) <- c(.6, .15, .15, .15)
huxtable::valign(ht)[1,] <- "bottom"
huxtable::escape_contents(ht) <- FALSE
huxtable::align(ht)[-1,2:4] <- "left"




# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-7.04') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_header_height(1) %>%
  set_column_header_buffer(1,0) %>%
  set_footer_height(1)

write_rtf(doc, file='./outputs/14-7.04.rtf')

