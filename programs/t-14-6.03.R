# Table 14-6.03

library(huxtable)
library(dplyr)
library(glue)
library(tidyverse)
library(haven)
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

## Chem
adlbc <- read_xpt(glue("{adam_lib}/adlbcpv.xpt")) %>%
  filter(SAFFL == "Y", ANL01FL == "Y", AVISITN != 99)

adlbc$TRTP <- ordered(adlbc$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))
adlbc$PARAM <- recode(adlbc$PARAM,
  "Albumin (g/L) change from previous visit, relative to normal range" = "ALBUMIN",
  "Alkaline Phosphatase (U/L) change from previous visit, relative to normal range" = "ALKALINE PHOSPHATASE",
  "Alanine Aminotransferase (U/L) change from previous visit, relative to normal range" = "ALANINE AMINOTRANSFERASE",
  "Aspartate Aminotransferase (U/L) change from previous visit, relative to normal range" = "ASPARTATE AMINOTRANSFERASE",
  "Bilirubin (umol/L) change from previous visit, relative to normal range" = "BILIRUBIN",
  "Blood Urea Nitrogen (mmol/L) change from previous visit, relative to normal range" = "UREA NITROGEN",
  "Calcium (mmol/L) change from previous visit, relative to normal range" = "CALCIUM",
  "Cholesterol (mmol/L) change from previous visit, relative to normal range" = "CHOLESTEROL",
  "Creatine Kinase (U/L) change from previous visit, relative to normal range" = "CREATINE KINASE",
  "Chloride (mmol/L) change from previous visit, relative to normal range" = "CHLORIDE",
  "Creatinine (umol/L) change from previous visit, relative to normal range" = "CREATININE",
  "Gamma Glutamyl Transferase (U/L) change from previous visit, relative to normal range" = "GAMMA GLUTAMYL TRANSFERASE",
  "Glucose (mmol/L) change from previous visit, relative to normal range" = "GLUCOSE",
  "Potassium (mmol/L) change from previous visit, relative to normal range" = "POTASSIUM",
  "Sodium (mmol/L) change from previous visit, relative to normal range" = "SODIUM",
  "Phosphate (mmol/L) change from previous visit, relative to normal range" = "PHOSPHATE",
  "Protein (g/L) change from previous visit, relative to normal range" = "PROTEIN",
  "Urate (umol/L) change from previous visit, relative to normal range" = "URATE"
)
adlbc$PARAM <- ordered(adlbc$PARAM, c(
  "ALBUMIN",
  "ALKALINE PHOSPHATASE",
  "ALANINE AMINOTRANSFERASE",
  "ASPARTATE AMINOTRANSFERASE",
  "BILIRUBIN",
  "UREA NITROGEN",
  "CALCIUM",
  "CHOLESTEROL",
  "CREATINE KINASE",
  "CHLORIDE",
  "CREATININE",
  "GAMMA GLUTAMYL TRANSFERASE",
  "GLUCOSE",
  "POTASSIUM",
  "SODIUM",
  "PHOSPHATE",
  "PROTEIN",
  "URATE"
))
adlbc$ANRIND <- ordered(adlbc$ANRIND, c("L", "N", "H"))

adlbc2 <- adlbc %>%
  filter(!is.na(TRTP)) %>%
  group_by(PARAM, TRTP, ANRIND) %>%
  complete(nesting(PARAM, TRTP, ANRIND)) %>%
  summarise(N = n()) %>%
  group_by(PARAM, TRTP) %>%
  mutate(tot = sum(N))

adlbc_pvals <- list()

for(i in seq(nrow(adlbc2)/9)) {
  adlbc_pvals[i] <- round(fisher.test(
    matrix(unlist(adlbc2[((i-1)*9+1):(i*9), "N"]), nrow = 3, ncol = 3, byrow = TRUE)
  )$p.value, 3)
}

adlbc3 <- adlbc2 %>%
  mutate(n_w_pct = n_pct(N, tot, n_width = 2)) %>%
  pivot_wider(id_cols = PARAM,names_from = c(TRTP, ANRIND), values_from = n_w_pct) %>%
  add_column("p-val\\line [1]" = num_fmt(unlist(adlbc_pvals), digits = 3, int_len = 1, size = 5))

### Heme
adlbh <- read_xpt(glue("{adam_lib}/adlbhpv.xpt")) %>%
  filter(SAFFL == "Y", ANL01FL == "Y", AVISITN != 99)

adlbh$TRTP <- ordered(adlbh$TRTP, c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"))
adlbh$PARAM <- recode(adlbh$PARAM,
  "Basophils (GI/L) change from previous visit, relative to normal range" = "BASOPHILS",
  "Eosinophils (GI/L) change from previous visit, relative to normal range" = "EOSINOPHILS",
  "Hematocrit change from previous visit, relative to normal range" = "HEMATOCRIT",
  "Hemoglobin (mmol/L) change from previous visit, relative to normal range" = "HEMOGLOBIN",
  "Lymphocytes (GI/L) change from previous visit, relative to normal range" = "LYMPHOCYTES",
  "Ery. Mean Corpuscular Hemoglobin (fmol(Fe)) change from previous visit, relative to normal range" = "ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "Ery. Mean Corpuscular HGB Concentration (mmol/L) change from previous visit, relative to normal rang" = "ERY. MEAN CORPUSCULAR HB CONCENTRATION",
  "Ery. Mean Corpuscular Volume (fL) change from previous visit, relative to normal range" = "ERY. MEAN CORPUSCULAR VOLUME",
  "Monocytes (GI/L) change from previous visit, relative to normal range" = "MONOCYTES",
  "Platelet (GI/L) change from previous visit, relative to normal range" = "PLATELET",
  "Erythrocytes (TI/L) change from previous visit, relative to normal range" = "ERYTHROCYTES",
  "Leukocytes (GI/L) change from previous visit, relative to normal range" = "LEUKOCYTES"
)
adlbh$PARAM <- ordered(adlbh$PARAM, c(
  "BASOPHILS",
  "EOSINOPHILS",
  "HEMATOCRIT",
  "HEMOGLOBIN",
  "LYMPHOCYTES",
  "ERY. MEAN CORPUSCULAR HEMOGLOBIN",
  "ERY. MEAN CORPUSCULAR HB CONCENTRATION",
  "ERY. MEAN CORPUSCULAR VOLUME",
  "MONOCYTES",
  "PLATELET",
  "ERYTHROCYTES",
  "LEUKOCYTES"
))
adlbh$ANRIND <- ordered(adlbh$ANRIND, c("L", "N", "H"))


adlbh2 <- adlbh %>%
  filter(!is.na(PARAM), !is.na(TRTP)) %>%
  group_by(PARAM, TRTP, ANRIND) %>%
  complete(nesting(PARAM, TRTP, ANRIND)) %>%
  summarise(N = n()) %>%
  group_by(PARAM, TRTP) %>%
  mutate(tot = sum(N))

adlbh_pvals <- list()

for(i in seq(nrow(adlbh2)/9)) {
  adlbh_pvals[i] <- round(fisher.test(
    matrix(unlist(adlbh2[((i-1)*9+1):(i*9), "N"]), nrow = 3, ncol = 3, byrow = TRUE)
  )$p.value, 3)
}

adlbh3 <- adlbh2 %>%
  mutate(n_w_pct = n_pct(N, tot, n_width = 2)) %>%
  pivot_wider(id_cols = PARAM,names_from = c(TRTP, ANRIND), values_from = n_w_pct) %>%
  add_column("p-val\\line [1]" = num_fmt(unlist(adlbh_pvals), digits = 3, int_len = 1, size = 5))

final <- adlbc3 %>%
  ungroup() %>%
  add_row("PARAM" = "----------", .before = 1) %>%
  add_row("PARAM" = "CHEMISTRY", .before = 1) %>%
  add_row("PARAM" = "", .before = 1) %>%
  add_row("PARAM" = "") %>%
  add_row("PARAM" = "") %>%
  add_row("PARAM" = "") %>%
  add_row("PARAM" = "") %>%
  add_row("PARAM" = "") %>%
  add_row("PARAM" = "HEMATOLOGY") %>%
  add_row("PARAM" = "----------") %>%
  rbind(ungroup(adlbh3))


names(final) <- c(
  "",
  "Low",
  "Normal",
  "High",
  "Low",
  "Normal",
  "High",
  "Low",
  "Normal",
  "High",
  "p-val\\line[1]"
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

ht <- final %>%
  as_huxtable(add_colnames = TRUE)

ht <- as_hux(pad_row(as.data.frame(ht), c(1,1)), add_colnames = FALSE)
ht[1, 2] <- headers[1, "label"]
ht[1, 5] <- headers[3, "label"]
ht[1, 8] <- headers[2, "label"]

ht2 <- ht %>%
  huxtable::merge_cells(1, 2:4) %>%
  huxtable::merge_cells(1, 5:7) %>%
  huxtable::merge_cells(1, 8:10) %>%
  huxtable::set_bottom_border(2, 2:4, 1) %>%
  huxtable::set_bottom_border(2, 5:7, 1) %>%
  huxtable::set_bottom_border(2, 8:10, 1) %>%
  huxtable::set_bottom_border(3, 1:11, 1) %>%
  huxtable::set_col_width(1:11, c(0.18, rep(0.082, 10))) %>%
  huxtable::set_width(1.4) %>%
  huxtable::set_escape_contents(FALSE) %>%
  huxtable::set_bold(1:3, 1:11, TRUE) %>%
  huxtable::set_valign(1:3, 1:11, "bottom") %>%
  huxtable::set_align(3, 1:11, "center") %>%
  huxtable::set_align(1, 1:11, "center")

# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht2, header_rows = 3) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-6.03') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top = 1) %>%
  set_footer_height(1.25)


# Write out the RTF
write_rtf(doc, file='./outputs/14-6.03.rtf')

