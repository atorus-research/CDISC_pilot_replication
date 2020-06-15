# T-14-6.04

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

## Chem
adlbc <- read_xpt(glue("{adam_lib}/adlbc.xpt")) %>%
  filter(SAFFL == "Y", AVISITN != 99)
adlbh <- read_xpt(glue("{adam_lib}/adlbh.xpt")) %>%
  filter(SAFFL == "Y", AVISITN != 99)
comb <- rbind(adlbc, adlbh)

#sort tests
comb$PARAM<- recode(comb$PARAM,
                    "Alanine Aminotransferase (U/L)" = "ALANINE AMINOTRANSFERASE",
                    "Albumin (g/L)" = "ALBUMIN",
                    "Alkaline Phosphatase (U/L)" = "ALKALINE PHOSPHATASE",
                    "Aspartate Aminotransferase (U/L)" = "ASPARTATE AMINOTRANSFERASE",
                    "Bilirubin (umol/L)" = "BILIRUBIN",
                    "Calcium (mmol/L)" = "CALCIUM",
                    "Chloride (mmol/L)" = "CHLORIDE",
                    "Cholesterol (mmol/L)" = "CHOLESTEROL",
                    "Creatine Kinase (U/L)" = "CREATINE KINASE",
                    "Creatinine (umol/L)" = "CREATININE",
                    "Gamma Glutamyl Transferase (U/L)" = "GAMMA GLUTAMYL TRANSFERASE",
                    "Glucose (mmol/L)" = "GLUCOSE",
                    "Phosphate (mmol/L)" = "PHOSPHATE",
                    "Potassium (mmol/L)" = "POTASSIUM",
                    "Protein (g/L)" = "PROTEIN",
                    "Sodium (mmol/L)" = "SODIUM",
                    "Urate (umol/L)" = "URATE",
                    "Blood Urea Nitrogen (mmol/L)" = "UREA NITROGEN",
                    "Basophils (GI/L)" = "BASOPHILS",
                    "Eosinophils (GI/L)" = "EOSINOPHILS",
                    "Ery. Mean Corpuscular HGB Concentration (mmol/L)" = "ERY. MEAN CORPUSCULAR HB CONCENTRATION",
                    "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))" = "ERY. MEAN CORPUSCULAR HEMOGLOBIN",
                    "Ery. Mean Corpuscular Volume (fL)" = "ERY. MEAN CORPUSCULAR VOLUME",
                    "Erythrocytes (TI/L)" = "ERYTHROCYTES",
                    "Hematocrit" = "HEMATOCRIT",
                    "Hemoglobin (mmol/L)" = "HEMOGLOBIN",
                    "Leukocytes (GI/L)" = "LEUKOCYTES",
                    "Lymphocytes (GI/L)" = "LYMPHOCYTES",
                    "Monocytes (GI/L)" = "MONOCYTES",
                    "Platelet (GI/L)" = "PLATELET")
#sort tests
comb$PARAM <-ordered(comb$PARAM, c(
  "ALANINE AMINOTRANSFERASE",
  "ALBUMIN",
  "ALKALINE PHOSPHATASE",
  "ASPARTATE AMINOTRANSFERASE",
  "BILIRUBIN",
  "CALCIUM",
  "CHLORIDE",
  "CHOLESTEROL",
  "CREATINE KINASE",
  "CREATININE",
  "GAMMA GLUTAMYL TRANSFERASE",
  "GLUCOSE",
  "PHOSPHATE",
  "POTASSIUM",
  "PROTEIN",
  "SODIUM",
  "URATE",
  "UREA NITROGEN",
  "BASOPHILS",
  "EOSINOPHILS",
  "ERY. MEAN CORPUSCULAR HB CONCENTRATION",
  "ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "ERY. MEAN CORPUSCULAR VOLUME",
  "ERYTHROCYTES",
  "HEMATOCRIT",
  "HEMOGLOBIN",
  "LEUKOCYTES",
  "LYMPHOCYTES",
  "MONOCYTES",
  "PLATELET"))
comb$TRTP <- ordered(comb$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))
comb$ANRIND <- ordered(comb$ANRIND, c("N", "H"))
comb$BNRIND <- ordered(comb$BNRIND, c("N", "H"))
comb$VISIT <- ordered(comb$VISIT, c(
  "WEEK 2",
  "WEEK 4",
  "WEEK 6",
  "WEEK 8",
  "WEEK 12",
  "WEEK 16",
  "WEEK 20",
  "WEEK 24",
  "WEEK 26"
))

total_ABLFL1 <- comb%>%
  filter(!is.na(VISIT), !is.na(TRTP), !is.na(BNRIND), !is.na(ANRIND), !is.na(PARAM)) %>%
  group_by(PARAM, VISIT, TRTP, BNRIND) %>%
  complete(nesting(TRTP, ABLFL)) %>%
  summarise(N = n())
total_ABLFL <- total_ABLFL1 %>%
  mutate(ANRIND = ordered("T", c("T", "N", "H"))) %>%
  pivot_wider(id_cols = c(PARAM, VISIT, ANRIND), names_from = c(TRTP, BNRIND), values_from = N) %>%
  ungroup()

comb2 <- comb %>%
  filter(!is.na(VISIT), !is.na(TRTP), !is.na(BNRIND), !is.na(ANRIND), !is.na(PARAM)) %>%
  group_by(PARAM, VISIT, TRTP, BNRIND, ANRIND) %>%
  complete(nesting(BNRIND, ANRIND)) %>%
  summarise(N = n()) %>%
  mutate(n2 = n_pct(N, total_ABLFL1[total_ABLFL1$PARAM == PARAM &
                                         total_ABLFL1$VISIT == VISIT   &
                                         total_ABLFL1$TRTP == TRTP     &
                                         total_ABLFL1$BNRIND == BNRIND, "N"], n_width = 2)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(PARAM, VISIT, ANRIND), names_from = c(TRTP, BNRIND), values_from = n2)

comb2$ANRIND <- ordered(comb2$ANRIND, c("T", "N", "H"))

total_ABLFL$Placebo_N <- num_fmt(total_ABLFL$Placebo_N, size = 2, int_len = 2)
total_ABLFL$Placebo_H <- num_fmt(total_ABLFL$Placebo_H, size = 2, int_len = 2)
total_ABLFL$`Xanomeline Low Dose_N` <- num_fmt(total_ABLFL$`Xanomeline Low Dose_N`, size = 2, int_len = 2)
total_ABLFL$`Xanomeline Low Dose_H` <- num_fmt(total_ABLFL$`Xanomeline Low Dose_H`, size = 2, int_len = 2)
total_ABLFL$`Xanomeline High Dose_N` <- num_fmt(total_ABLFL$`Xanomeline High Dose_N`, size = 2, int_len = 2)
total_ABLFL$`Xanomeline High Dose_H` <- num_fmt(total_ABLFL$`Xanomeline High Dose_H`, size = 2, int_len = 2)

comb3 <- comb2 %>%
  rbind(total_ABLFL) %>%
  arrange(PARAM, VISIT, ANRIND)

comb3$VISIT <- as.character(str_extract(comb3$VISIT, "[0-9]+"))
comb3$ANRIND <- as.character(recode(comb3$ANRIND,
                       "T" = "n",
                       "N" = "Normal",
                       "H" = "High"))
names(comb3) <- c(
  "rowlbl",
  "Week",
  "Shift to",
  "Normal at Baseline",
  "High at Baseline",
  "Normal at Baseline",
  "High at Baseline",
  "Normal at Baseline",
  "High at Baseline"
)

comb3 <- comb3[!apply(comb3, 1, function(x) {
  all(x[4:9] ==  " 0      ") & all(x[3] == "High")
}), ]

comb4 <- pad_row(comb3, which(comb3$`Shift to` == "n")) %>%
  add_row('Week' = NA, .before = 1) %>%
  add_row("Week" = NA, .before = 1)
comb4 <- comb4 %>%
  add_row("Week" = NA, .before = 541) %>%
  add_row("Week" = NA, .before = 541)

comb4[,1] <- as.character(comb4$rowlbl)

comb4[!(comb4$`Shift to` %in% "n") , 2] <- NA
comb4[!(comb4$Week %in% "2"), 1] <- NA
comb4[2,1] <- "CHEMISTRY"
comb4[3,1] <- "----------"
comb4[542,1] <- "HEMATOLOGY"
comb4[543,1] <- "----------"


names(comb4) <- c(
  "",
  "Week",
  "Shift to",
  "Normal at Baseline",
  "High at Baseline",
  "Normal at Baseline",
  "High at Baseline",
  "Normal at Baseline",
  "High at Baseline"
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



ht <- comb4 %>%
  as_huxtable(add_colnames = TRUE)

ht <- as_hux(pad_row(as.data.frame(ht), c(1,1)), add_colnames = FALSE)
ht[1, 4] <- headers[1, "label"]
ht[1, 6] <- headers[3, "label"]
ht[1, 8] <- headers[2, "label"]

ht2 <- ht %>%
  huxtable::merge_cells(1, 4:5) %>%
  huxtable::merge_cells(1, 6:7) %>%
  huxtable::merge_cells(1, 8:9) %>%
  huxtable::set_bottom_border(2, 4:5, 1) %>%
  huxtable::set_bottom_border(2, 6:7, 1) %>%
  huxtable::set_bottom_border(2, 8:9, 1) %>%
  huxtable::set_bottom_border(3, 1:9, 1) %>%
  huxtable::set_width(1.4) %>%
  huxtable::set_escape_contents(FALSE) %>%
  huxtable::set_bold(1:3, 1:9, TRUE) %>%
  huxtable::set_valign(1:3, 1:9, "bottom") %>%
  huxtable::set_align(3, 1:9, "center") %>%
  huxtable::set_align(1, 1:9, "center") %>%
  huxtable::set_col_width(1:9, c(0.29, 0.06, 0.07, rep(0.1, 6)))

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht2, header_rows = 3) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-6.04') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top = 1) %>%
  set_header_height(0.75) %>%
  set_footer_height(1)


# Write out the RTF
write_rtf(doc, file='./outputs/14-6.04.rtf')


