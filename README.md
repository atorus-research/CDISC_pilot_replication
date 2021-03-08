# CDISC Pilot Replication in R

## Updates

We've update this repository to be compatible with Huxtable v5.0.0! Huxtable v5.0.0 had some backwards compatibility breaking changes. All updates within this repository are compatible back to Huxtable v4.7.1. See a full list of changes to Huxtable [here](https://hughjonesd.github.io/whats-new-in-huxtable-5.0.0.html). 

The changes most of interest to a user of this repository are:
 - The way indexing was handled has changed. You can no longer index columns like `ht[1:5]`. The new syntax is `ht[1:5, ]`. 
 - Additionally, the `add_columns` argument now changed from a default of `TRUE` instead of `FALSE`. We updated [`config.R`](programs/config.R) to reset the default option to keep the code consistent. You can do this like so:

```
options(huxtable.add_colnames = FALSE)
```

## Introduction
Welcome to the Atorus CDISC Pilot replication repository! 
In 2007, the [original pilot project submission package](https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse) was finalized and released following a review by FDA Staff, where the CDISC data and metadata contained within the package were evaluated for its suitability in meeting the needs and expectations of medical and statistical reviewers. In 2019, the [PHUSE Test Data Factory](https://advance.phuse.global/display/WEL/Test+Dataset+Factory) took on the goal of replicating the SDTM and ADaM data within the CDISC pilot package to match more modern data standards, bringing the ADaM data up to version 1.1. 
Atorus Research has now regenerated the table outputs within the CDISC Pilot Project using the PHUSE Test Data Factory project’s data and the R Programming language. Our motivation behind this project was to:
-	Demonstrate that we were able to obtain matching outputs using R
-	Provide open source code to the public to demonstrate how we were able to do this
-	Demonstrate our first publicly released R package, `pharmaRTF`, in action.
You can find our package pharmaRTF right [here](https://github.com/atorus-research/pharmaRTF)

## Setup Instructions
To obtain the data for this repository, you can download the data from the PHUSE Github Repository, using [this link](https://github.com/phuse-org/phuse-scripts/blob/master/data/adam/TDF_ADaM_v1.0.zip) for ADaM data and [this link](https://github.com/phuse-org/phuse-scripts/blob/master/data/sdtm/TDF_SDTM_v1.0%20.zip)
 for the SDTM.

This repository was programmed using R 3.6. For further system information, see our [session information](SessionInfo.txt).

## Notes on Data
Every effort was made to use best programming practices to recreate the values on the CDISC Pilot displays, however some values on our outputs do not align with the values on the CDISC Pilot displays due to the following reasons:

### General:
-	The ADaMs we used to regenerate the CDISC Pilot displays were the PHUSE CDISC Pilot replication ADaMs following ADaM V1.1.  Since the CDISC Pilot displays were not regenerated using the PHUSE CDISC Pilot replication data there are likely discrepancies between the original CDISC Pilot analysis data and the PHUSE CDISC Pilot replication ADaMs.
-	SAS and R round differently.  While SAS rounds up if the value is 5 or greater, R rounds to the nearest even number.
-	In some circumstances, R packages will not produce a p-value if the the counts within the data are not high enough to make it statistically meaningful. An example of this is BILIRUBIN on Table 14-6.05. The High at Baseline stratum only has a single count, and the `mantelhein.test` function in R requires that each stratum has more than 1 observation. 

### Output Specific Details:
- Table 14-3.07 ADAS Cog (11) - Change from Baseline to Week 24 - Completers at Wk 24-Observed Cases-Windowed
  - Difference in baseline values.  Despite following the analysis results metadata (ARM) in the original CDISC Pilot Define.xml, the baseline counts are off by 1. Following the information available within the original ARM and SAP there is a discrepancy between the available data and the display using both the PHUSE CDISC Pilot replication data and the original CDISC Pilot analysis data.
- Table 14-3.11 ADAS Cog (11) - Repeated Measures Analysis of Change from Baseline to Week 24
  - Differences in p-values.  See ARM in the original CDISC Pilot Define.xml (ARM-Leaf0046) for details on implementation in SAS.  These numbers end up being slightly off.  To anyone that finds this and can match the numbers, please feel free to submit a PR and correct our implementation!  To the best of our knowledge, we've matched what we could. It's not explicit, but the default covariance structure in the lme4 package in unstructured.
- Table 14-3.12 Mean NPI-X Total Score from Week 4 through Week 24 – Windowed
  - Difference in values of Mean of Weeks 4-24. This was programmed using the derived NPTOTMN variable. The ARM in the original CDISC Pilot Define.xml was followed as best as possible to determine the subset.  This means that the counts are a discrepancy with the original CDISC Pilot analysis data, which is no longer available, therefore we were not able to investigate the discrepancy with the current PHUSE CDISC Pilot replication data.  The subsequent statistical summaries therefore also have differences.
- Table 14-6.05 Shifts of Laboratory Values During Treatment, Categorized Based on Threshold Ranges
  -	Difference in the values for BILIRUBIN values for the Xan. Low group and MONOCYTES values for the Placebo group.  The Analysis Reference Range Indicator and Shift variables are used as is from the ADaM which indicates there are likely discrepancies for reference ranges and shifts between the original CDISC Pilot analysis data and the PHUSE CDISC Pilot replication data.  The subsequent statistical summaries therefore also have differences.
- Table 14-6.06 Shifts of Hy's Law Values During Treatment
  -	Difference in the values for Transaminase 1.5 x ULN for the Xan. Low group and Total Bili 1.5 x ULN and Transaminase 1.5 x ULN for all groups.  The Analysis Reference Range Indicator and Shift variables are used as is from the ADaM which indicates there are likely discrepancies for reference ranges and shifts between the original CDISC Pilot analysis data and the PHUSE CDISC Pilot replication data.  The subsequent statistical summaries therefore also have differences.
- Table 14-7.01 Summary of Vital Signs at Baseline and End of Treatment
  -	Difference in values for End of Treatment for all groups.  The End of Treatment flag is used as is from the ADaM which indicates there are likely discrepancies for end of treatment between the original CDISC Pilot analysis data and the PHUSE CDISC Pilot replication data.
- Table 14-7.02 Summary of Vital Signs Change from Baseline at End of Treatment
  - Difference in values throughout table.  The End of Treatment flag is used as is from the ADaM which indicates there are likely discrepancies for end of treatment between the original CDISC Pilot analysis data and the PHUSE CDISC Pilot replication data.
- Table 14-7.03 Summary of Weight Change from Baseline at End of Treatment
  -	Difference in values for End of Treatment for all groups.  The End of Treatment flag is used as is from the ADaM which indicates there are likely discrepancies for end of treatment between the original CDISC Pilot analysis data and the PHUSE CDISC Pilot replication data.
## Notes on R Packages
As many programmers in the R community do, we relied on the [tidyverse](https://www.tidyverse.org/packages/) for much of our data processing. There are a few addition libraries that we used worth mentioning:
-	For the CMH test where testing for the alternate hypothesis that row means differ, the package vcdExtra was used, which is not included in the base distribution of R. We additionally had to make a slight modification to this library, which is available in [this fork]( https://github.com/mstackhouse/vcdExtra) of the package. The update is due to the fact that the `solve` function in R will throw an error when processing large, sparse tables. By replacing `solve` with `MASS::ginv`, the error is bypassed. We as an organization plan to perform further testing of this update and submit the update back to the original author.
-	For Mixed Models, the `lme4` package was used
-	For ANCOVA models, `car` was used and the `emmeans` package  was used to do LSMEANS.
