# funcs.R
#   File to store functions separately from data processing

require(tidyverse)
require(glue)

example_custom_reader <- function(..., table_number=NULL) {

  # Make sure that readxl is installed before
  if (suppressWarnings(!require('readxl'))) {
    stop("This reader requires the package `readxl`. Install using `install.packages('readxl')")
  }

  # If a column isn't populated then the type may be guessed wrong so force it
  col_types <- c('text', 'numeric', 'text', 'text', 'text', 'text', 'logical', 'logical', 'text')
  # pass through arguments from ...
  df <- readxl::read_excel(..., col_types=col_types)

  # Subset and return that dataframe
  df[df$table_number==table_number, !names(df) == 'table_number']
}


get_meta <- function(df) {
  # Examines the metadata of a SAS imported tibble
  for (name in names(df)) {
    lab <- attr(df[[name]], 'label')
    name <- str_pad(name, 8)
    cat(name, lab, "\n", sep="\t")
  }
}

get_header_n <- function(.data, trtp = TRT01P, trtpn = TRT01PN) {
  # Extract header N's into a dataframe to be used on merges or for display

  trtp = enquo(trtp)
  trtpn = enquo(trtpn)

  # Get the header N's ----
  .data %>%
    group_by(!!trtp, !!trtpn) %>%
    summarize(N = n()) %>%
    mutate(
      trtp = !!trtp,
      labels = str_replace_all(str_wrap(glue('{trtp} (N={N})'), width=10), "\n", function(x) "\\line ")
      # labels = str_wrap(glue('{TRTP} (N={N})'), width=10)
    ) %>%
    ungroup() %>%
    arrange(!!trtpn) %>%
    select(-!!trtp, -trtp)
}

num_fmt <- function(var, digits=0, size=10, int_len=3) {
  # Formats summary stat strings to align display correctly

  if (is.na(var)) return('')

  # Set nsmall to input digits
  nsmall = digits

  # Incremement digits for to compensate for display
  if (digits > 0) {
    digits = digits + 1
  }

  # Form the string
  return(str_pad(
    format(
      # Round
      round(var, nsmall),
      # Set width of format string
      width=(int_len+digits),
      # Decimals to display
      nsmall=nsmall
    ),
    # Overall width padding
    side='right', size
  ))
}

num_fmt <- Vectorize(num_fmt)

n_pct <- function(n, pct, n_width=3, pct_width=3, digits=0, mark_lt=TRUE) {
  # n (%) formatted string. e.g. 50 ( 75%)
  res <- n / pct

  if (res < .01 & mark_lt) {
    disp <- str_pad('<1', width=pct_width)
  } else {
    disp <- format(round(res * 100, digits=digits), width=pct_width, nsmall=digits)
  }
  return(
    # Suppress conversion warnings
    as.character(
      # Form the string using glue and format
      glue('{format(n, width=n_width)} ({disp}%)')
    )
  )
}

n_pct <- Vectorize(n_pct)

sum_subgrp <- function(.data, subgroup_var, order_var = NULL, include.n=TRUE, pad.row=TRUE, header_n = header_n) {
  # Create n (%) subgroups by TRT01P

  # Pull from adsl with totals
  subgrps <- .data %>%
    # Keep only the gtwo group variables and group byC:\Users\16105\OneDrive - ATorus\Documents\Projects\Explore\test2.rtf
    select(TRT01PN, {{ subgroup_var }}, {{ order_var }}) %>%
    filter(!is.na({{ subgroup_var }})) %>%
    group_by(TRT01PN, {{ subgroup_var }}, {{ order_var }}) %>%
    # Summarize counts
    summarize(
      n = n()
    ) %>%
    arrange(TRT01PN, {{ order_var }}) %>%
    # Merge with big Ns
    left_join(header_n, by = 'TRT01PN') %>%
    rowwise() %>%
    # Create the n (%) string
    mutate(
      res = n_pct(n, N)
    ) %>%
    # Drop unnecessary vars
    select(-n, -N, -labels, -{{ order_var }}) %>%
    # Transpose
    pivot_wider(names_from = TRT01PN, values_from = res) %>%
    # Take care of NA results
    replace(is.na(.), '  0       ') %>%
    # Rename row label column
    rename(rowlbl2 = {{ subgroup_var }})

  if (include.n){
    subgrps <- rbind(desc_stats(.data, {{subgroup_var}}, include='n')[1,], subgrps)
  }

  pad_row(subgrps)

}

desc_stats <- function(.data, var, group = TRT01PN, na.rm=TRUE, int_len=3, size=10, include=c('n', 'Mean', 'SD', 'Median', 'Min', 'Max')) {
  # Provides descriptive statistics of provided variable, by TRT01PN
  # n, Mean, SD, Median, Min, Max

  # Ensure that the include argument was valid
  include = match.arg(include, several.ok=TRUE)

  # This is gonna get wonky - store each summary as an expression
  #TODO: Allow flexibility in significant digits - right now it's hard coded
  summaries <- list(
    n      = rlang::expr(num_fmt(      n()         , digits=0, int_len=int_len, size=size)),
    Mean   = rlang::expr(num_fmt(   mean({{ var }}), digits=1, int_len=int_len, size=size)),
    SD     = rlang::expr(num_fmt(     sd({{ var }}), digits=2, int_len=int_len, size=size)),
    Median = rlang::expr(num_fmt( median({{ var }}), digits=1, int_len=int_len, size=size)),
    Min    = rlang::expr(num_fmt(    min({{ var }}), digits=1, int_len=int_len, size=size)),
    Max    = rlang::expr(num_fmt(    max({{ var }}), digits=1, int_len=int_len, size=size))
  )[include] # this is a named list, so subset based on the input arguments

  # Pull from ADSL with totals
  .data %>%
    # Pick of TRT01PN and the variable of interest
    select({{ group }}, {{ var }}) %>%
    # Filter out missing values
    filter(!is.na({{ var }})) %>%
    # Group by treatment
    group_by({{ group }}) %>%
    # Summarize each statistic and use num_fmt for rounding/formatting
    summarize(!!!summaries) %>% # unpack the expressions into syntax to be evaluated
    # Transpose statistics into one column
    pivot_longer(-{{ group }}, names_to = 'rowlbl2', values_to = 'temp') %>%
    # Transpose treatments into separate columns
    pivot_wider(names_from = {{ group }}, values_from = temp) %>%
    pad_row()
}

invert.list <- function (NL) {
  # Invert a list's items and names (assuming it's key value)
  L <- list()
  for (i in 1:length(NL)) {
    L[[NL[[i]]]] <- names(NL[i])
  }
  L
}

# Count of subjects with an adverse event
ae_counts <- function(.data, ..., N_counts = NULL, sort=FALSE) {

  # Get the grouping
  grouped_data <- .data %>%
    group_by(TRTAN, TRTA, ...) %>%
    select(TRTA, TRTAN, ..., USUBJID)

  # Counts of each subject
  subject_counts <- grouped_data %>%
    distinct() %>%
    summarize(n = n())

  # Count of adverse events
  event_counts <- grouped_data %>%
    summarize(AEs = n())

  # Join the subject and event counts, pivot out by treatment
  counts <- subject_counts %>%
    left_join(event_counts) %>%
    pivot_wider(id_cols=c(...), names_from=TRTAN, values_from=c(n, AEs))

  # If no events for a treatment group, they won't be in the pivoted table, so create
  for (g in unique(N_counts$TRT01PN)) {
    cnames <- c(paste0('n_', g), paste0('AEs_', g))
    if (!all(cnames %in% names(counts))) {
      # If one is missing, they're both missing
      counts[cnames[1]] <- 0
      counts[cnames[2]] <- 0
    }
  }

  # Add in subject counts
  counts['N_0'] <- N_counts[N_counts$TRT01PN == 0, 'N']
  counts['N_54'] <- N_counts[N_counts$TRT01PN == 54, 'N']
  counts['N_81'] <- N_counts[N_counts$TRT01PN == 81, 'N']

  # Fill all NAs with 0
  counts[is.na(counts)] <- 0

  # Find no event counts
  counts['no_event_0'] <- counts$N_0 - counts$n_0
  counts['no_event_54'] <- counts$N_54 - counts$n_54
  counts['no_event_81'] <- counts$N_81 - counts$n_81

  # Calculate p-values
  counts['p_low']  <- apply(counts[, c('n_0', 'n_54', 'no_event_0', 'no_event_54')], MARGIN=1, FUN=fisher_test_ae)
  counts['p_high'] <- apply(counts[, c('n_0', 'n_81', 'no_event_0', 'no_event_81')], MARGIN=1, FUN=fisher_test_ae)

  # Formatting
  counts <- counts %>%
    rowwise() %>%
    mutate(
      npct_0  = ifelse(n_0  > 0, n_pct(n_0,   N_0,   n_width=2, pct_width=4, digits=1), ' 0        '),
      npct_54 = ifelse(n_54 > 0, n_pct(n_54,  N_54,  n_width=2, pct_width=4, digits=1), ' 0        '),
      npct_81 = ifelse(n_81 > 0, n_pct(n_81,  N_81,  n_width=2, pct_width=4, digits=1), ' 0        '),
      cAEs_0  = ifelse(n_0  > 0, paste0('[',AEs_0,  ']'), ''),
      cAEs_54 = ifelse(n_54 > 0, paste0('[',AEs_54, ']'), ''),
      cAEs_81 = ifelse(n_81 > 0, paste0('[',AEs_81, ']'), ''),
      ord3 = n_81 # Use for descending event order
    )

  # Remove numeric columns not used in display
  counts <- counts %>%
    select(-starts_with('n_'), -starts_with('no_'), -starts_with('AEs'))

}

# P-value for anova test
aov_p <- function(.data, forumula) {
  # Run the anova test
  a <- aov(forumula, .data, na.action=na.omit)

  # Extract the P value
  p <- summary(a)[[1]][['Pr(>F)']][1]
  # Format the output

  format(round(p, 4), width=10, nsmall=4)
}

# P-value for chi-squared
chi_p <- function(data, results, categories) {
  # get the arguments as a off of the function call
  arguments <- as.list(match.call())
  # Evaluate the arguments within the dataframe environment
  # This is all just so I can allow acceptance of variables without quotes
  cats <- factor(eval(arguments$categories, data))
  res <- factor(eval(arguments$results, data))

  p <- chisq.test(res, cats)$p.value
  if(round(p, 4) == 0) return("<.0001")
  format(round(p, 4), width=10, nsmall=4)
}

# P-vaule for fisher's test
fish_p <- function(data, results, categories, width = 10) {
  # get the arguments as a off of the function call
  arguments <- as.list(match.call())
  # Evaluate the arguments within the dataframe environment
  # This is all just so I can allow acceptance of variables without quotes
  cats <- factor(eval(arguments$categories, data))
  res <- factor(eval(arguments$results, data))

  p <- fisher.test(res, cats)$p.value
  if(round(p, 4) == 0) return("<.0001")
  format(round(p, 4), width=width, nsmall=4)
}

# Fisher test built for row-wise derivations suited for our AE tables
fisher_test_ae <- function(.data) {

  # If there were no events in either treatment arm then don't compute
  if (sum(.data[1:2]) == 0){
    return('')
  }

  # convert to a 2X2 matrix
  dim(.data) <- c(2, 2)

  # Return the p-value of interest
  p <- fisher.test(.data)$p.value

  # Format the p-values for display
  disp <- num_fmt(p, digits=3, size=5, int_len=1)

  # Post process the display for special representation
  if (p > .99) {
    disp <- '>.99'
  } else if (p < .15) {
    disp <- paste0(disp, '*')
  } else {
    disp <- paste0(disp, ' ')
  }
}

# CMH test p-value with options for alternate hyptheses
# !!! NOTE: To obtain the same p-value used in SAS for this display, a modification had to be made to the vcdExtra library.
#           Please refer to this github issue: https://github.com/friendly/vcdExtra/issues/3
#           And you can access our fork of the library here: https://github.com/mstackhouse/vcdExtra
cmh_p <- function(.data, formula, alternate=c('rmeans', 'cmeans', 'general', 'cor')) {
  # Pull out the hypoth
  alternate <- match.arg(alternate, several.ok=FALSE)

  # Run the test
  res <- vcdExtra::CMHtest(formula, data=.data, overall=TRUE)$ALL

  pvalue <- unlist(res$table[alternate, 'Prob'])
  pvalue
}

# Attach P-value to the first row of a dataframe
attach_p <- function(.data, p_value, digits = 4) {
  # Empty column
  .data[['p']] = character(nrow(.data))
  .data[1, 'p'] = p_value

  .data

}

# Add an empty row
pad_row <- function(.data, n=1) {
  .data[(nrow(.data)+1):(nrow(.data)+n), ] <- ""
  .data
}


## Efficacy functions ----


# Create the summary data portion for 14-3.01 through XXXX ----
summary_data <- function(data, var, week, stub_header) {

  # Get the summary statistics
  df <- data %>%
    # Filter to analsis week
    filter(AVISITN==week) %>%
    # Summary statistics
    group_by(TRTPN) %>%
    summarize(n = n(),
              mean = mean({{var}}),
              sd = sd({{var}}),
              median = median({{var}}),
              min = min({{var}}),
              max = max({{var}})) %>%
    rowwise() %>%
    # Form into display strings
    mutate(
      N = num_fmt(n, size=12, int_len=2),
      mean_sd = as.character(
        glue('{num_fmt(mean, digits=1, int_len=2, size=4)} ({num_fmt(sd, digits=2, int_len=2, size=4)})')
      ),
      med_ran = as.character(
        glue('{num_fmt(median, digits=1, int_len=2, size=4)} ({num_fmt(min, int_len=3, size=3)};{num_fmt(max, int_len=2, size=2)})')
      )
    ) %>%
    # Drop the numeric values
    select(-n, -mean, -sd, -median, -min, -max) %>%
    # Make summary stats vertical
    pivot_longer(c(N, mean_sd, med_ran), names_to = "rowlbl1") %>%
    # Split out treatment groups into separate columns
    pivot_wider(names_from=TRTPN, values_from=value) %>%
    # Fix the row labels
    mutate(rowlbl1 =
      case_when(
        rowlbl1 == 'N' ~ '  n',
        rowlbl1 == 'mean_sd' ~ '  Mean (SD)',
        rowlbl1 == 'med_ran' ~ '  Median (Range)'
      )
    )

  # Add in the stub header
  bind_rows(tibble(rowlbl1=c(stub_header)), df)
}


efficacy_models <- function(data, var=NULL, wk=NULL, model_type='ancova') {

  if (model_type == 'ancova') {
    # Need to set contrasts to work for Type III SS. See analysis results metadata for
    # table 14-3.01. Reference for R here: https://www.r-bloggers.com/anova-%E2%80%93-type-iiiiii-ss-explained/
    op <- options(contrasts = c("contr.sum","contr.poly"))

    # Subset to analyze
    data <- data %>%
      filter(AVISITN == wk)
  }

  data <- data %>%
    mutate(
      TRTPCD = case_when(
        TRTPN == 0 ~ 'Pbo',
        TRTPN == 54 ~ 'Xan_Lo',
        TRTPN == 81 ~ 'Xan_Hi'
      )
    )

  # Create an ordered factor variable for the models
  data['TRTPCD_F'] <- factor(data$TRTPCD, levels=c('Xan_Hi', 'Xan_Lo', 'Pbo'))
  data['AWEEKC'] = factor(data$AVISIT)

  # Set up the models
  if (model_type == 'ancova') {
    if (var == "CHG") {
      model1 <- lm(CHG ~ TRTPN + SITEGR1 + BASE, data=data)
      model2 <- lm(CHG ~ TRTPCD_F + SITEGR1 + BASE, data=data)
    } else {
      model1 <- lm(AVAL ~ TRTPN + SITEGR1, data=data)
      model2 <- lm(AVAL ~ TRTPCD_F + SITEGR1, data=data)
    }
  } else {
    model2 <- lme4::lmer(CHG ~ TRTPCD_F + SITEGR1 + AWEEKC + TRTPCD_F:AWEEKC + BASE + BASE:AWEEKC + (AVISITN | USUBJID),
                         data=data)
  }

  ## Dose Response ---
  # NOTE: For statistics portions, I purposefully did not import the libraries to make it explicitly clear which
  # packages were being used to match P-values.
  # Use `car` package Anova test with type III SS.

  if (model_type == 'ancova') {
    ancova <- car::Anova(model1, type=3)

    # Pull it out into a table
    sect1 <- tibble(rowlbl1=c('p-value(Dose Response) [1][2]'),
                            `81` = c(num_fmt(ancova[2, 'Pr(>F)'], int_len=4, digits=3, size=12))
    ) %>%
      pad_row()
  }

  ## Pairwise Comparisons ----
  # Here's a reference for the emmeans package and how to use it:
  #   https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html
  # Adjustments made are in line with the analysis results metadata in the analysis define
  # and PROC GLM documentation.

  # Linear model but use treatment group as a factor now

  if (model_type == 'ancova') {
    # LS Means and weight proportionately to match OM option on PROC GLM in SAS
    lsm <- emmeans::lsmeans(model2, ~TRTPCD_F, weights='proportional')
  } else { # Mixed model LS means values - Get the dataset from here
    # See analysis results metadata in analysis define (ARM-Leaf0046) for details on implementation in SAS
    # These numbers end up being slightly off, mostly in the P-values. To anyone that finds this and can match the numbers,
    # please feel free to submit a PR and correct our implementation!
    # To the best of my knowledge, I've matched what I could here. It's not explicit, but the default
    # covariance structure in the lme4 package in unstructured
    lsm <- emmeans::lsmeans(model2, ~TRTPCD_F, lmer.df='kenward-roger')
    # Build the section 1 data here instead of above because its from the same model
    sect1 <- as_tibble(lsm) %>%
      rowwise() %>%
      mutate(
        rowlbl1 = 'LS Means (SE)',
        values = as.character(glue('{num_fmt(lsmean, int_len=1, digits=1, size=3)} ({num_fmt(SE, int_len=1, digits=2, size=4)})'))
      ) %>%
      select(rowlbl1, TRTPCD_F, values) %>%
      pivot_wider(id_cols = rowlbl1, names_from=TRTPCD_F, values_from=values) %>%
      rename(`0` = Pbo, `54`=Xan_Lo, `81`=Xan_Hi) %>%
      pad_row()
  }

  # Here on out - it's all the same data manipulation
  # Get pairwise contrast and remove P-values adjustment for multiple groups
  cntrst_p <- emmeans::contrast(lsm, method="pairwise", adjust=NULL)
  # 95% CI
  cntrst_ci <- confint(cntrst_p)

  # merge and convert into dataframe
  pw_data <- as_tibble(summary(cntrst_p)) %>%
    merge(as_tibble(cntrst_ci)) %>%
    rowwise() %>%
    # Create the display strings
    mutate(
      p = num_fmt(p.value, int_len=4, digits=3, size=12),
      diff_se = as.character(
        glue('{num_fmt(estimate, int_len=2, digits=1, size=4)} ({num_fmt(SE, int_len=1, digits=2, size=4)})')
      ),
      ci = as.character(
        glue('({num_fmt(lower.CL, int_len=2, digits=1, size=4)};{num_fmt(upper.CL, int_len=1, digits=1, size=3)})')
      )
    ) %>%
    # Clean out the numeric variables
    select(contrast, p, diff_se, ci) %>%
    # Transpose
    pivot_longer(c('p', 'diff_se', 'ci'), names_to = 'row_label')

  # Subset Xan_Lo - Pbo into table variables
  xan_lo <- pw_data %>%
    filter(contrast == 'Xan_Lo - Pbo') %>%
    # Rename to the table display variable
    select(`54`=value) %>%
    pad_row()

  #Add in rowlbl
  xan_lo['rowlbl1'] <- c('p-value(Xan - Placebo) [1][3]', '  Diff of LS Means (SE)', '  95% CI', '')

  # Subset Xan_hi - Pbo into table variables
  xan_hi <- pw_data %>%
    filter(contrast == 'Xan_Hi - Pbo') %>%
    # Rename to the table display variable
    select(`81`=value) %>%
    pad_row()
  # Add in rowlbl
  xan_hi['rowlbl1'] <- c('p-value(Xan - Placebo) [1][3]', '  Diff of LS Means (SE)', '  95% CI', '')
  xan_hi['ord'] <- c(1,2,3,4) # Order for sorting

  # Subset Xan_Hi - Xan_Lo into table variable
  xan_xan <- pw_data %>%
    filter(contrast == 'Xan_Hi - Xan_Lo') %>%
    # Rename to the table display variable
    select(`81`=value)
  # Add in rowlbl
  xan_xan['rowlbl1'] <- c('p-value(Xan High - Xan Low) [1][3]', '  Diff of LS Means (SE)', '  95% CI')
  xan_xan['ord'] <- c(5,6,7) # Order for sorting

  # Pack it all together
  pw_final <- merge(xan_lo, xan_hi, by='rowlbl1') %>%
    bind_rows(xan_xan) %>%
    arrange(ord)

  # Return the statistics together
  return(bind_rows(sect1, pw_final))

}