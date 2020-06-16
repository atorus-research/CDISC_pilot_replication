## Table 14-7.02 Summary of Vital Signs Change from Baseline at End of Treatment

library(huxtable)
library(dplyr)
library(glue)
library(tidyverse)
library(haven)
library(pharmaRTF)
library(tibble)

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

adsl <- read_xpt(glue("{adam_lib}/adsl.xpt"))
advs <- read_xpt(glue("{adam_lib}/advs.xpt")) %>%
  filter(SAFFL == "Y" & !is.na(BASE))

advs$EOTFL <- ifelse(advs[,"AVISIT"] == "End of Treatment", "Y", "")
advs$W24FL <- ifelse(advs[, "AVISIT"] == "Week 24", "Y", "")

advs2 <- advs %>%
  filter(EOTFL == "Y" | W24FL == "Y") %>%
  filter(PARAM %in% c("Diastolic Blood Pressure (mmHg)",
                      "Pulse Rate (beats/min)",
                      "Systolic Blood Pressure (mmHg)"))

advs2$PRTFL <- ifelse(advs2[,"EOTFL"] == "Y", "End of Trt.","Week 24")

advs2$TRTP <- ordered(advs2$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))
## Add ordered VISITS to order visits
advs2$AVISIT <- ordered(advs2$AVISIT, c("Baseline", "Week 24", "End of Treatment"))
advs2$PRTFL <- ordered(advs2$PRTFL, c("Week 24", "End of Trt."))
advs2$PARAM <- recode(advs2$PARAM,
                      "Pulse Rate (beats/min)" = "Pulse (bpm)")
advs2$PARAM <- ordered(advs2$PARAM, c("Systolic Blood Pressure (mmHg)",
                                      "Diastolic Blood Pressure (mmHg)",
                                      "Pulse (bpm)"))

advs3 <- advs2 %>%
  group_by(PARAM, ATPT, TRTP, PRTFL) %>%
  summarise(n = n(),
            mean = mean(CHG, na.rm = TRUE),
            sd = sd(CHG, na.rm = TRUE),
            median = median(CHG, na.rm = TRUE),
            min = min(CHG, na.rm = TRUE),
            max = max(CHG, na.rm = TRUE))

advs4 <- add_column(advs3, "N" = apply(advs3,
                                       1,
                                       function(x) {aSum <- sum(adsl[,"ARM"] == x["TRTP"], na.rm = TRUE)
                                       ifelse(aSum == 0, NA, aSum)}),
                    .after = 3)

advs4[!(advs4$PRTFL %in% "Week 24"), "TRTP"] <- NA
advs4[!(advs4$TRTP %in% "Placebo"), "ATPT"] <- NA
advs4[!(advs4$ATPT %in% "AFTER LYING DOWN FOR 5 MINUTES"), "PARAM"] <- NA
advs4[!(advs4$PRTFL %in% "Week 24"), "N"] <- NA



advs4$TRTP <- apply(advs4, 1, function(x) {switch(x["TRTP"],
                                                  "Placebo" = "Placebo",
                                                  "Xanomeline High Dose" = "Xan.High",
                                                  "Xanomeline Low Dose" = "Xan.Low",
                                                  NA)})

advs4$mean <- num_fmt(advs4$mean, digits = 1, size = 3)
advs4$sd <- num_fmt(advs4$sd, digits = 2, size = 4, int_len = 2)
advs4$median <- num_fmt(advs4$median, digits = 1, size = 2, int_len = 2)
advs4$min <- num_fmt(advs4$min, digits = 1, size = 4, int_len = 2)
advs4$max <- num_fmt(advs4$max, digits = 1, size = 4, int_len = 2)


names(advs4) <- c(
  "Measure",
  "Position",
  "Treatment",
  "N",
  "Planned Relative Time",
  "n",
  "Mean",
  "SD",
  "Median",
  "Min.",
  "Max."
)

advs4 <- pad_row(advs4, which(advs4[, "Planned Relative Time"] == "End of Trt.") + 1)

ht <- advs4 %>%
  huxtable::as_hux(add_colnames=TRUE) %>%
  huxtable::set_bold(1, 1:ncol(advs4), TRUE) %>%
  huxtable::set_align(1, 1:ncol(advs4), "center") %>%
  huxtable::set_align(2:nrow(advs4), 3, "center") %>%
  huxtable::set_align(2:nrow(advs4), 4:ncol(advs4), "left") %>%
  huxtable::set_valign(1, 1:ncol(advs4), "bottom") %>%
  huxtable::set_bottom_border(1, 1:ncol(advs4), 1) %>%
  huxtable::set_width(1.45) %>%
  huxtable::set_col_width(1:ncol(advs4), c(0.2, 0.15, 0.19, 0.03, 0.1, 0.03, 0.06, 0.06, 0.06, 0.06, 0.06))
wrap(ht) <- FALSE

doc <- rtf_doc(ht) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  
  reader=example_custom_reader,
  table_number='14-7.02') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_header_height(1) %>%
  set_column_header_buffer(1,0) %>%
  set_footer_height(1.3)

write_rtf(doc, file='./outputs/14-7.02.rtf')

