## 14-6.05

library(huxtable)
library(glue)
library(tidyverse)
library(haven)
library(tibble)
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

adlbc <- read_xpt(glue("{adam_lib}/adlbc.xpt")) %>%
  filter(SAFFL == "Y", ANL01FL == "Y")
adlbh <- read_xpt(glue("{adam_lib}/adlbh.xpt")) %>%
  filter(SAFFL == "Y", ANL01FL == "Y")
comb <- rbind(adlbc, adlbh)


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
comb$BNRIND <- ordered(comb$BNRIND, c("N", "H"))
comb$ANRIND <- ordered(comb$ANRIND, c("N", "H"))

comb <- comb %>%
  filter(!is.na(comb$PARAM), !is.na(comb$TRTP), !is.na(comb$BNRIND), !is.na(comb$ANRIND), AVISITN != 99)

total_bltrfl1 <- comb%>%
  filter(!is.na(TRTP), !is.na(BNRIND), !is.na(ANRIND), !is.na(PARAM)) %>%
  group_by(PARAM, TRTP, BNRIND) %>%
  complete(nesting(TRTP, BNRIND)) %>%
  summarise(N = n())
total_bltrfl <- total_bltrfl1 %>%
  mutate(ANRIND = ordered("T", c("T", "N", "H"))) %>%
  pivot_wider(id_cols = c(PARAM, ANRIND), names_from = c(TRTP, BNRIND), values_from = N) %>%
  ungroup()

comb2 <- comb %>%
  group_by(PARAM, TRTP, BNRIND, ANRIND) %>%
  complete(nesting(BNRIND, ANRIND)) %>%
  summarise(N = n()) %>%
  ungroup()


pvals <- c()
for(i in levels(comb$PARAM)) {
  mat <-  comb[comb$PARAM == i, c("ANRIND", "TRTP", "BNRIND")]

  if(all(mat[, "ANRIND"] == "N")) pvals[i] <- ""
  else {
    pvals[i] <- tryCatch(num_fmt(cmh_p(mat, ANRIND ~ TRTP | BNRIND, alternate = "rmeans"), digits = 3, int_len = 1),
                         error = function(c) ""
             )
  }
}

temp1 <- mat %>%
  group_by(TRTP) %>%
  summarise(n = n())

comb3 <- comb2 %>%
  group_by(PARAM, TRTP, BNRIND) %>%
  mutate(n2 = n_pct(N, total_bltrfl1[total_bltrfl1$PARAM == PARAM &
                                       total_bltrfl1$TRTP == TRTP  &
                                       total_bltrfl1$BNRIND == BNRIND, "N"], n_width = 2)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(PARAM, ANRIND), names_from = c(TRTP, BNRIND), values_from = n2)

comb4 <- comb3[!apply(comb3, 1, function(x) {
  all(x[3:8] ==  " 0      ") & all(x[2] == "H")
}), ]


comb4$ANRIND <- ordered(comb4$ANRIND, c("T", "N", "H"))

total_bltrfl$Placebo_N <- num_fmt(total_bltrfl$Placebo_N, size = 2, int_len = 2)
total_bltrfl$Placebo_H <- num_fmt(total_bltrfl$Placebo_H, size = 2, int_len = 2)
total_bltrfl$`Xanomeline Low Dose_N` <- num_fmt(total_bltrfl$`Xanomeline Low Dose_N`, size = 2, int_len = 2)
total_bltrfl$`Xanomeline Low Dose_H` <- num_fmt(total_bltrfl$`Xanomeline Low Dose_H`, size = 2, int_len = 2)
total_bltrfl$`Xanomeline High Dose_N` <- num_fmt(total_bltrfl$`Xanomeline High Dose_N`, size = 2, int_len = 2)
total_bltrfl$`Xanomeline High Dose_H` <- num_fmt(total_bltrfl$`Xanomeline High Dose_H`, size = 2, int_len = 2)



comb5 <- comb4 %>%
  rbind(total_bltrfl) %>%
  arrange(PARAM, ANRIND)

comb5$ANRIND <- as.character(recode(comb5$ANRIND,
                                    "T" = "n",
                                    "N" = "Normal",
                                    "H" = "High"))

comb5[unlist(comb5[,2] == "n")[,1], 9] <- pvals
comb5 <- pad_row(comb5, which(comb5[,2] == "n")) %>%
  ungroup() %>%
  add_row("PARAM" = NA, .before = 1) %>%
  add_row("PARAM" = NA, .before = 1)
comb5 <- comb5 %>%
  add_row("PARAM" = NA, .before = 65) %>%
  add_row("PARAM" = NA, .before = 65)

comb5[!(unlist(comb5[,2]) %in% "n") , 1] <- NA

comb5 <- comb5[!apply(comb5, 1, function(x) {
  all(x[4:9] ==  " 0      ") & all(x[3] == "High")
}), ]

comb5[,1] <- as.character(comb5$PARAM)
comb5[2,1] <- "CHEMISTRY"
comb5[3,1] <- "----------"
comb5[66,1] <- "HEMATOLOGY"
comb5[67,1] <- "----------"

names(comb5) <- c(
  "",
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

ht <- comb5 %>%
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
  huxtable::set_bottom_border(3, 1:9, 1) %>%
  huxtable::set_width(1.5) %>%
  huxtable::set_escape_contents(FALSE) %>%
  huxtable::set_bold(1:3, 1:9, TRUE) %>%
  huxtable::set_valign(1:3, 1:9, "bottom") %>%
  huxtable::set_align(3, 1:9, "center") %>%
  huxtable::set_align(1, 1:9, "center") %>%
  huxtable::set_align(4:102, 9, "right") %>%
  huxtable::set_col_width(1:9, c(0.31, rep(0.09, 7), 0.06))


# Write into doc object and pull titles/footnotes from excel file
doc <- rtf_doc(ht2, header_rows = 3) %>% titles_and_footnotes_from_df(
  from.file='./data/titles.xlsx',
  reader=example_custom_reader,
  table_number='14-6.05') %>%
  set_font_size(10) %>%
  set_ignore_cell_padding(TRUE) %>%
  set_column_header_buffer(top = 1) %>%
  set_header_height(1) %>%
  set_footer_height(1.3)

# Write out the RTF
write_rtf(doc, file='./outputs/14-6.05.rtf')

