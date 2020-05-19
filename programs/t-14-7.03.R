### Table 14-7.03 pg 148 Summary of Weight Change from Baseline at End of Treatment


library(dplyr)
library(glue)
library(tidyverse)
library(haven)
library(assertthat)
library(pharmaRTF)
library(tibble)


source('./programs/config.R')
source('./programs/funcs.R')

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

advs <- read_xpt(glue("{adam_lib}/advs.xpt")) %>%
  filter(PARAM == "Weight (kg)")

advs$TRTP <- ordered(advs$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))

advs$EOTFL <- ifelse(advs$AVISITN == 99, "Y", NA)
advs$W24 <- ifelse(advs$AVISITN == 24, "Y", NA)
advs$ABLFL <- ifelse(advs$ABLFL == "Y", "Y", NA)

#Rbinded data.frame
advs2 <- rbind(advs[advs$EOTFL %in% "Y", ], advs[advs$W24 %in% "Y",], advs[advs$ABLFL %in% "Y",])

# Create table for stats
bw_stats <- advs2 %>%
  group_by(TRTP, ABLFL, W24, EOTFL) %>%
  summarise(n = n(),
            Mean = mean(AVAL),
            SD = sd(AVAL),
            Median = median(AVAL),
            Min. = min(AVAL),
            Max. = max(AVAL)) %>%
  arrange()
bw_stats[, 2] <- rep(c("Baseline", "Week 24", "End of Trt."), 3)
bw_stats <- bw_stats[, c(-3, -4)]

bw_stats <- add_column(bw_stats, "Measure" = "Weight (kg)", .before= 1)
bw_stats[unlist(bw_stats[, 3]) != "Baseline", "TRTP"] <- NA
bw_stats[!(bw_stats$TRTP %in% "Placebo"), "Measure"] <- NA

adsl <- read_xpt(glue("{adam_lib}/adsl.xpt"))
bw_stats <- add_column(bw_stats, "N" = apply(bw_stats,
                           1,
                           function(x) {aSum <- sum(adsl[,"ARM"] == x["TRTP"], na.rm = TRUE)
                           ifelse(aSum == 0, NA, aSum)}),
           .before = 3)
# Pad blank rows after End of Trt. rows
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
bw_stats <- pad_row(bw_stats, which(bw_stats$ABLFL == "End of Trt.",) + 1)
names(bw_stats)[4] <- "VISIT"


### Weight Change from Baseline table
# Create table for baseline changes


.blfun = function(x, usubjid = NULL) {
  x <- x[x$USUBJID == unique(usubjid),]
  bl <- as.numeric(x[x$ABLFL %in% "Y", "AVAL"])
  w24 <- as.numeric(x[x$W24 %in% "Y", "AVAL"])
  eot <- as.numeric(x[x$EOTFL %in% "Y", "AVAL"])
  arm <- unique(x$TRTP)
  ## Done this way to make dplyr easier
  c(ifelse(length(w24-bl) == 0, NA, w24-bl),
    ifelse(length(eot-bl) == 0, NA, eot-bl))
}
bw_bl <- advs2 %>%
  select(USUBJID, TRTP, ABLFL, W24, EOTFL, AVAL) %>%
  group_by(USUBJID) %>%
  summarise(`WEEK 24` = .blfun(., USUBJID)[1],
         `End of Trt.` = .blfun(., USUBJID)[2],
         TRTP = unique(TRTP)) %>%
  select(USUBJID, TRTP, `WEEK 24`, `End of Trt.`) %>%
  pivot_longer(c(`WEEK 24`, `End of Trt.`), names_to = "VISIT", values_to = "change")
## Add ordered factor to order arms
bw_bl$TRTP <- ordered(bw_bl$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))
bw_bl$VISIT <- ordered(bw_bl$VISIT,c("WEEK 24", "End of Trt."))

bw_bl_1 <- bw_bl %>%
  group_by(TRTP, VISIT) %>%
  summarise(n = sum(!is.na(change)),
            Mean = mean(change, na.rm = TRUE),
            SD = sd(change, na.rm = TRUE),
            Median = median(change, na.rm = TRUE),
            Min. = min(change, na.rm = TRUE),
            Max. = max(change, na.rm = TRUE)) %>%
  ungroup()
bw_bl_1 <- add_column(bw_bl_1, "Measure" = "Weight Change\\line from Baseline", .before = 1)
bw_bl_1[bw_bl_1$VISIT != "WEEK 24", "TRTP"] <- NA
bw_bl_1[!(bw_bl_1$TRTP %in% "Placebo"), "Measure"] <- NA

bw_bl_1 <- add_column(bw_bl_1, "N" = apply(bw_bl_1,
                                             1,
                                             function(x) {aSum <- sum(adsl[,"ARM"] == x["TRTP"], na.rm = TRUE)
                                             ifelse(aSum == 0, NA, aSum)}),
                       .before = 3)
bw_bl_1 <- pad_row(bw_bl_1, which(bw_bl_1$VISIT == "End of Trt.") + 1)

### Combine Tables and match output
combinedTable <- rbind(bw_stats, bw_bl_1)
names(combinedTable)[2] <- "Treatment"
names(combinedTable)[4] <- "Planned Relative Time"
combinedTable[,"Treatment"] <- apply(combinedTable, 1, function(x){
  switch(x["Treatment"],
         "Placebo" = "Placebo",
         "Xanomeline Low Dose" = "Xan.Low",
         "Xanomeline High Dose" = "Xan.High",
         NA)
})
combinedTable[,"Planned Relative Time"] <- apply(combinedTable, 1, function(x){
  switch(x["Planned Relative Time"],
         "Baseline" = "Baseline",
         "WEEK 24" = "Week 24",
         "Week 24" = "Week 24",
         "End of Trt." = "End of Trt.",
         NA)
})

### Number formatting
class(combinedTable) <- "data.frame"
combinedTable[!is.na(combinedTable$Mean),"Mean"] <- num_fmt(unlist(combinedTable[!is.na(combinedTable$Mean),"Mean"]),
                                                            digits = 1, size = 3, int_len = 2)
combinedTable[!is.na(combinedTable$SD),"SD"] <-  num_fmt(unlist(combinedTable[!is.na(combinedTable$SD),"SD"]),
                                                         digits = 2, size = 3, int_len = 2)
combinedTable[!is.na(combinedTable$Median),"Median"] <-  num_fmt(unlist(combinedTable[!is.na(combinedTable$Median),"Median"]),
                                                                 digits = 1, size = 3, int_len = 2)
combinedTable[!is.na(combinedTable$`Min.`),"Min."] <-  num_fmt(unlist(combinedTable[!is.na(combinedTable$`Min.`),"Min."]),
                                                               digits = 1, size = 3, int_len = 2)
combinedTable[!is.na(combinedTable$`Max.`),"Max."] <-  num_fmt(unlist(combinedTable[!is.na(combinedTable$`Max.`),"Max."]),
                                                               digits = 1, size = 3, int_len = 2)


ht <- combinedTable %>%
  huxtable::as_hux(add_colnames=TRUE)


huxtable::bottom_border(ht)[1, ] <- 1
huxtable::bold(ht)[1, ] <- TRUE
huxtable::align(ht)[1, ] <- 'center'
huxtable::align(ht)[,c(3, 5:10)] <- "center"
huxtable::width(ht) <- 1.5
huxtable::bottom_padding(ht) <- 0
huxtable::top_padding(ht) <- 0
huxtable::col_width(ht) <- c(0.25, 0.15, 0.05, 0.185, 0.04, 0.065, 0.065, 0.065, 0.065, 0.065)
huxtable::valign(ht)[1,] <- "bottom"
huxtable::escape_contents(ht) <- FALSE



# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-7.03') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_header_height(1) %>%
  set_column_header_buffer(1,0)

write_rtf(doc, file='./outputs/14-7.03.rtf')

