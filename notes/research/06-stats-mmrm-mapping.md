I now have everything verified. Let me compose the deliverable.

---

# Legacy → Modern Statistical Method Mapping — CDISC Pilot Rebuild

Scope: `/Users/mstackhouse/Documents/repos/CDISC_pilot_replication`. All statistical logic lives in **`programs/funcs.R`** (the engine) and is called by the individual table programs. Below, every legacy signature/argument/default was verified against the actual source; every modern snippet was verified against the installed packages (`mmrm 0.3.15`, `emmeans 2.0.0`, `car 3.1-3`, `tplyr2 0.1.0.9000`, `clinify 0.3.0`) on R 4.5.1, and the mmrm+emmeans pipeline was run end-to-end to confirm the API.

## 1. Inventory of statistics-bearing tables

| Table(s) | Program | Legacy method | Engine function | Modern target |
|---|---|---|---|---|
| 14-3.01, .02, .03, .04, .05, .06, .08, .09, .12 | `t-14-3-0x.R`, `t-14-3-12.R` | **ANCOVA** (`lm`), Type III dose-response via `car::Anova`, LSMEANS + pairwise via `emmeans` | `efficacy_models(..., model_type='ancova')` | `lm` + `emmeans` (keep) |
| 14-3.07 | `t-14-3-07.R` | ANCOVA (same as above; a windowed/subset ADAS analysis) | `efficacy_models(...,'ancova')` | `lm` + `emmeans` |
| **14-3.11** | `t-14-3-11.R` | **Mixed model for repeated measures** via `lme4::lmer` + `emmeans` (Kenward-Roger) | `efficacy_models(..., model_type='repeated')` | **`mmrm` + `emmeans`** |
| 14-3.13 | `t-14-3-13.R` | **CMH row-mean-scores**, stratified, via forked `vcdExtra::CMHtest` | `cmh_p()` | `coin::cmh_test` or CRAN `vcdExtra::CMHtest` |
| 14-6.05 | `t-14-6.05.R` | **CMH row-mean-scores** (shift table), stratified by baseline range | `cmh_p()` | same as above |
| 14-6.06 | `t-14-6.06.R` | **Generalized CMH / general association** via base `stats::mantelhaen.test` | inline | **keep base R** |
| 14-2.01 | `t-14-2-01.R` | one-way **ANOVA** (`aov`) for continuous; `chisq.test`/`fisher.test` for categorical | `aov_p()`, `chi_p()`, `fish_p()` | `lm`/`aov` + `emmeans` (optional); keep chi/fisher |
| 14-5.01, .02 | `t-14-5-0x.R` | **Fisher's exact** per AE term (2×2) | `fisher_test_ae()` | keep `fisher.test` |

Vitals tables 14-7.0x, 14-4.01, and the AE-frequency 14-5.0x carry no model-based inference beyond Fisher's exact.

---

## 2. Exact legacy models

### 2.1 ANCOVA — `efficacy_models(data, var, wk, model_type='ancova')` (`funcs.R:401-546`)

Consumed by e.g. `t-14-3-01.R:35` → `efficacy_models(adas, 'CHG', 24)`. Steps:

1. **Data subset**: `data %>% filter(AVISITN == wk)` — a single analysis visit (e.g. Week 24). For 14-3.12, `wk=98` is the derived "Mean of Weeks 4-24" record.
2. **Type III setup**: `op <- options(contrasts = c("contr.sum","contr.poly"))` (needed so `car::Anova(type=3)` is meaningful for factor terms).
3. **Recode / factor**: `TRTPN {0,54,81}` → `TRTPCD {Pbo, Xan_Lo, Xan_Hi}`; `TRTPCD_F <- factor(TRTPCD, levels=c('Xan_Hi','Xan_Lo','Pbo'))`; `AWEEKC <- factor(AVISIT)`.
4. **Two models** (note the two exist only to separate the numeric dose-response test from the factor-based LSMEANS):
   - `var == "CHG"`:
     - `model1 <- lm(CHG ~ TRTPN + SITEGR1 + BASE, data=data)` — **TRTPN numeric** = dose-response linear term
     - `model2 <- lm(CHG ~ TRTPCD_F + SITEGR1 + BASE, data=data)` — treatment as factor
   - else (`AVAL`, e.g. adcibc; no baseline term):
     - `model1 <- lm(AVAL ~ TRTPN + SITEGR1, data=data)`
     - `model2 <- lm(AVAL ~ TRTPCD_F + SITEGR1, data=data)`
5. **Dose-response p-value**: `ancova <- car::Anova(model1, type=3)`; value taken is `ancova[2, 'Pr(>F)']` (row 2 = the `TRTPN` term). Row label `p-value(Dose Response) [1][2]`, placed in the `81` column.
6. **LSMEANS**: `lsm <- emmeans::lsmeans(model2, ~TRTPCD_F, weights='proportional')` — `weights='proportional'` replicates SAS `PROC GLM ... LSMEANS / OM` (weights marginal means by observed cell frequencies; averages over `SITEGR1`, holds numeric `BASE` at its mean).
7. **Pairwise contrasts**: `cntrst_p <- emmeans::contrast(lsm, method="pairwise", adjust=NULL)` (no multiplicity adjustment); `cntrst_ci <- confint(cntrst_p)` (95% CI).
8. **Placement into display** (`funcs.R:493-543`) — three formatted strings per contrast (`p`, `Diff of LS Means (SE)`, `95% CI`):
   - `Xan_Lo - Pbo` → column **`54`**
   - `Xan_Hi - Pbo` → column **`81`**
   - `Xan_Hi - Xan_Lo` → column **`81`** (rows "p-value(Xan High - Xan Low)")
   - Row labels: `p-value(Xan - Placebo) [1][3]`, `  Diff of LS Means (SE)`, `  95% CI`.

Response = `CHG` (change from baseline in ADAS-Cog `ACTOT`) or `AVAL` (CIBIC). Datasets: `adadas` (`PARAMCD=='ACTOT'`, `EFFFL/ITTFL=='Y'`, `ANL01FL=='Y'`), `adcibc`, `adnpix`.

### 2.2 Repeated-measures MMRM — `efficacy_models(..., model_type='repeated')` (14-3.11)

`t-14-3-11.R:29` → `efficacy_models(adas, 'CHG', 24, model_type='repeated')` on `adadas` filtered `EFFFL=='Y' & PARAMCD=='ACTOT' & ANL01FL=='Y' & DTYPE != 'LOCF' & AVISITN > 0` (all observed post-baseline visits, no LOCF).

- Model (`funcs.R:436-437`):
  ```r
  model2 <- lme4::lmer(
    CHG ~ TRTPCD_F + SITEGR1 + AWEEKC + TRTPCD_F:AWEEKC + BASE + BASE:AWEEKC + (AVISITN | USUBJID),
    data = data)
  ```
- LSMEANS (`funcs.R:472`): `lsm <- emmeans::lsmeans(model2, ~TRTPCD_F, lmer.df='kenward-roger')` — marginalized over weeks and site; produces the **`LS Means (SE)`** row for columns `0`, `54`, `81`.
- Pairwise contrasts + CI: identical downstream code to the ANCOVA path, same column placement.
- **Author's own caveat** (`funcs.R:468-472`): results, especially p-values, do **not** exactly match SAS. The random term `(AVISITN | USUBJID)` with **numeric** `AVISITN` is a random intercept + random slope, **not** an unstructured within-subject covariance (SAS `REPEATED / TYPE=UN`). This is precisely the mismatch that `mmrm` removes.

### 2.3 CMH row-mean-scores — `cmh_p()` (`funcs.R:324-333`)

```r
cmh_p <- function(.data, formula, alternate=c('rmeans','cmeans','general','cor')) {
  alternate <- match.arg(alternate, several.ok=FALSE)
  res <- vcdExtra::CMHtest(formula, data=.data, overall=TRUE)$ALL
  pvalue <- unlist(res$table[alternate, 'Prob'])
  pvalue
}
```
- Requires the **`mstackhouse/vcdExtra` fork** (comment `funcs.R:321-323`) to match SAS; upstream bug = `friendly/vcdExtra#3`.
- **14-3.13** (`t-14-3-13.R:106-119`): per visit (8/16/24), `cmh_p(AVAL ~ TRTP | SITEGR1)` with default `alternate='rmeans'`. `AVAL` = ordinal CIBIC score (1-7, rows), `TRTP` = treatment (columns), `SITEGR1` = stratum. This is SAS PROC FREQ CMH **"Row Mean Scores Differ"** — do mean CIBIC scores differ across treatments, controlling for site. One p-value placed on the "n" row per visit block.
- **14-6.05** (`t-14-6.05.R:142`): per lab parameter, `cmh_p(mat, ANRIND ~ TRTP | BNRIND, alternate="rmeans")`, where `mat = comb[comb$PARAM==i, c("ANRIND","TRTP","BNRIND")]`. Shift-from-baseline (`ANRIND` Normal/High) vs treatment, stratified by baseline range indicator `BNRIND`. Skipped (`""`) when all `ANRIND=="N"`; wrapped in `tryCatch(..., error=function(c) "")`.

### 2.4 Generalized CMH — `mantelhaen.test` (14-6.06)

`t-14-6.06.R:83`:
```r
num_fmt(mantelhaen.test(array(unlist(adlbhy_t1[,"N"]), dim = c(2,3,2)))$p.value,
        size = 6, int_len = 1, digits = 3)
```
- `adlbhy` `PARAMCD %in% c("TRANSHY","HYLAW")`, last-worst per subject. Array is **2 (post AVAL 0/1) × 3 (TRTP) × 2 (BASE strata)** → base R **generalized CMH general-association** test. (The second block, `HYLAW`/total bilirubin, is stubbed out at `t-14-6.06.R:118-122` — "Different counts", left blank.)
- No continuity correction is applied for I×J×K with dims > 2×2×K (base R `correct=` only affects 2×2×K).

### 2.5 Demographics inference — 14-2.01

- Continuous (`aov_p`, `funcs.R:254-263`): `a <- aov(formula, .data, na.action=na.omit)`; p = `summary(a)[[1]][['Pr(>F)']][1]`. Applied as one-way ANOVA `~ TRT01P` for `AGE`, `MMSETOT`, `DURDIS`, `EDUCLVL`, `WEIGHTBL`, `HEIGHTBL`, `BMIBL`.
- Categorical: `chi_p` (`chisq.test(res, cats)$p.value`, `<.0001` floor) and `fish_p` (`fisher.test`) against `TRT01P`.

### 2.6 AE Fisher — `fisher_test_ae()` (`funcs.R:294-318`), `fish_p`

2×2 `fisher.test` per AE term, treatment vs placebo; formats `>.99` and flags `*` when `p < .15`.

---

## 3. Modern equivalents

### 3.1 MMRM (replaces `lme4::lmer` for 14-3.11) — **required path**

**Verified** `mmrm 0.3.15` API (from `mmrm.Rd`, `mmrm_control.Rd`, `covariance_types.Rd`, `emmeans_support.Rd`, and a live run):

- `mmrm(formula, data, weights = NULL, covariance = NULL, reml = TRUE, control = mmrm_control(...), ...)`.
- Formula carries **exactly one** special covariance term. Covariance types (`cov_types()`): `us` (unstructured = SAS `TYPE=UN`), `cs/csh`, `ar1/ar1h`, `toep/toeph`, `ad/adh`, spatial `sp_exp`. **For all non-spatial structures the time variable must be a factor.**
- `mmrm_control(method = c("Satterthwaite","Kenward-Roger","Residual","Between-Within"), vcov = NULL, ...)`. `method` = d.f. adjustment; `vcov` = coefficient-covariance adjustment. Defaults: Satterthwaite→Asymptotic, Kenward-Roger→Kenward-Roger. **The help explicitly states: for an Unstructured covariance, `"Kenward-Roger"` differs from SAS; use `vcov = "Kenward-Roger-Linear"` to better match SAS.** `method`/`vcov`/`reml` can be passed straight to `mmrm(...)` via `...`.

**Direct translation of the 14-3.11 model** (SAS `PROC MIXED; MODEL CHG = TRT SITEGR1 AVISIT TRT*AVISIT BASE BASE*AVISIT; REPEATED AVISIT / SUBJECT=USUBJID TYPE=UN;`):

```r
library(mmrm); library(emmeans); library(dplyr); library(broom)

adas <- adas %>%
  filter(EFFFL == "Y", PARAMCD == "ACTOT", ANL01FL == "Y",
         DTYPE != "LOCF", AVISITN > 0) %>%
  mutate(
    TRTPCD_F = factor(
      case_when(TRTPN == 0 ~ "Pbo", TRTPN == 54 ~ "Xan_Lo", TRTPN == 81 ~ "Xan_Hi"),
      levels = c("Pbo", "Xan_Lo", "Xan_Hi")),   # Pbo first => cleaner Xan - Pbo contrasts
    AWEEKC   = factor(AVISIT),                    # visit MUST be a factor for us()
    SITEGR1  = factor(SITEGR1),
    USUBJID  = factor(USUBJID))

fit <- mmrm(
  formula = CHG ~ SITEGR1 + TRTPCD_F * AWEEKC + BASE * AWEEKC + us(AWEEKC | USUBJID),
  data    = adas,
  reml    = TRUE,
  method  = "Kenward-Roger",        # d.f. method
  vcov    = "Kenward-Roger-Linear"  # matches SAS for TYPE=UN (per mmrm docs)
)
```

`TRTPCD_F * AWEEKC` expands to the treatment main effect, week main effect, and treatment×week interaction; `BASE * AWEEKC` gives baseline + baseline×week — matching the legacy `lmer` fixed-effects exactly, while `us(AWEEKC | USUBJID)` supplies the true unstructured covariance the `lmer` term lacked.

**LSMEANS at Week 24** (emmeans support is built into mmrm — confirmed):

```r
emm <- emmeans(fit, ~ TRTPCD_F | AWEEKC, weights = "proportional")  # SAS OM analogue
lsm_wk24 <- as.data.frame(emm) %>% filter(AWEEKC == "Week 24")
# columns: TRTPCD_F, AWEEKC, emmean (LS mean), SE, df, lower.CL, upper.CL
```

**Treatment contrasts at Week 24** (estimate / SE / CI / p, unadjusted):

```r
emm24 <- emmeans(fit, ~ TRTPCD_F, at = list(AWEEKC = "Week 24"), weights = "proportional")

pw <- contrast(emm24, method = list(
  "Xan_Lo - Pbo"    = c(Pbo = -1, Xan_Lo = 1, Xan_Hi = 0),
  "Xan_Hi - Pbo"    = c(Pbo = -1, Xan_Lo = 0, Xan_Hi = 1),
  "Xan_Hi - Xan_Lo" = c(Pbo = 0,  Xan_Lo = -1, Xan_Hi = 1)
), adjust = NULL)

pw_df <- as.data.frame(summary(pw, infer = c(TRUE, TRUE)))
# -> contrast, estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value  (all needed cells)
```

Notes verified in the live run: `summary(pw, infer = c(TRUE, TRUE))` returns CI **and** p-value together; the explicit named contrast list removes any ambiguity about direction/order (safer than `pairs(..., reverse=TRUE)`).

**Programmatic single-contrast alternative** (no emmeans) — `df_1d(object, contrast)` returns `list(est, se, df, t_stat, p_val)`; `df_md` for multi-row F-tests. Useful if you prefer to build contrast vectors against `fit$beta_est` directly.

**Diagnostics/glue**: `tidy(fit, conf.int = TRUE, conf.level = 0.95)` → coefficient tibble (`term, estimate, std.error, df, statistic, p.value, conf.low, conf.high`); `glance(fit)` for fit stats; `VarCorr(fit)` for the estimated covariance.

### 3.2 ANCOVA (14-3.01 … 14-3.12) — `lm` + `emmeans` (unchanged engine, modernized)

`emmeans` operates on `lm` directly; no special support needed.

```r
library(emmeans); library(car); library(dplyr)

d24 <- adas %>% filter(AVISITN == 24) %>%
  mutate(TRTPCD_F = factor(case_when(TRTPN==0~"Pbo", TRTPN==54~"Xan_Lo", TRTPN==81~"Xan_Hi"),
                           levels = c("Pbo","Xan_Lo","Xan_Hi")),
         SITEGR1  = factor(SITEGR1))

## LS means + pairwise (treatment as factor)
m_lsm <- lm(CHG ~ TRTPCD_F + SITEGR1 + BASE, data = d24)   # drop BASE for AVAL/CIBIC tables
emm   <- emmeans(m_lsm, ~ TRTPCD_F, weights = "proportional")
pw    <- contrast(emm, method = "pairwise", adjust = NULL)
ci    <- confint(pw)                    # 95% CI on the differences

## Dose-response linear trend (treatment as numeric dose)
op    <- options(contrasts = c("contr.sum", "contr.poly"))   # for Type III on factor terms
m_dr  <- lm(CHG ~ TRTPN + SITEGR1 + BASE, data = d24)
p_dose <- car::Anova(m_dr, type = 3)["TRTPN", "Pr(>F)"]
options(op)
```

Faithful to the legacy: keep `car::Anova(type=3)` for the dose-response term. (For the single numeric `TRTPN` term the Type III F equals `summary(m_dr)`'s t² and yields the identical p-value, so `summary(m_dr)$coefficients["TRTPN","Pr(>|t|)"]` is an equivalent, dependency-free substitute — the `contr.sum` option does not change a numeric main-effect test, it only matters if you Type-III-test the `SITEGR1` factor.)

### 3.3 CMH row-mean-scores (14-3.13, 14-6.05) — maintained replacement for the vcdExtra fork

Neither `coin` nor `vcdExtra` is currently installed here, so exact SAS-parity should be **validated on install** — flagging that explicitly. Two maintained options:

**Option A — `coin::cmh_test()`** (actively maintained on CRAN). SAS "Row Mean Scores Differ" = assign ordinal scores to the response, nominal grouping to treatment, stratify with `| block`:

```r
library(coin)
# 14-3.13, per visit: ordinal CIBIC (1-7) response, treatment groups, site strata
d8 <- cbic %>% filter(AVISITN == 8) %>%
  mutate(AVAL = ordered(AVAL), TRTP = factor(TRTP), SITEGR1 = factor(SITEGR1))
res <- cmh_test(AVAL ~ TRTP | SITEGR1, data = d8,
                scores = list(AVAL = 1:nlevels(d8$AVAL)))  # integer row scores => row-mean-scores stat
pvalue(res)
```
With an `ordered` response + integer `scores`, `coin::cmh_test` computes the mean-score (nonzero-location) CMH statistic — the analogue of SAS's row-mean-scores-differ. Validate that `pvalue(res)` reproduces the legacy value before locking it in.

**Option B — CRAN `vcdExtra::CMHtest()`** (still maintained by M. Friendly). Keeps the exact `cmh_p()` interface and the `"rmeans"/"cmeans"/"cor"/"general"` statistic selection. Before relying on it, confirm whether `friendly/vcdExtra#3` (the reason for the fork) is fixed in the current CRAN release; if so, the fork is no longer needed:
```r
vcdExtra::CMHtest(AVAL ~ TRTP | SITEGR1, data = d8, overall = TRUE)$ALL$table["rmeans", "Prob"]
```

Recommendation: prefer **coin** (clean CRAN dependency, no fork), but keep `vcdExtra::CMHtest` as the drop-in if numeric parity with the fork is required and the upstream fix has landed.

### 3.4 Generalized CMH (14-6.06) — keep base R

`stats::mantelhaen.test()` is base R and maintained; no replacement needed. Modern-tidy wrapper only:
```r
arr <- array(unlist(adlbhy_t1[, "N"]), dim = c(2, 3, 2))  # response x treatment x stratum
p   <- mantelhaen.test(arr, correct = FALSE)$p.value      # general association for I x J x K
```

### 3.5 Demographics / AE — keep base R

One-way ANOVA `aov`/`lm` (optionally surface via `emmeans`/`broom::tidy`), `chisq.test`, `fisher.test` are all base/maintained. No migration required; only wrap outputs into the tplyr2/clinify assembly (below).

---

## 4. Assembling model output into a display + tplyr2 / clinify

**tplyr2 build output schema** (verified by running `tplyr_build`): columns are `rowlabel1`, `rowlabel2`, … (row labels), `res1`, `res2`, `res3`, … (one result column **per column-variable value, in tplyr2's column order** — not `var1_<name>` as in Tplyr v1), plus `ord_layer_1`, …, `ord_layer_index`.

**Key architectural point.** `group_analyze()` (signature `group_analyze(target_var, by=NULL, where=NULL, analyze_fn, settings=layer_settings())`) calls `analyze_fn(.data, .target_var)` **once per column × by group**, i.e. it only ever sees one treatment column's slice at a time. That is perfect for self-contained per-group statistics (geometric mean, custom ratios), but it is **not** suitable for model-based LSMEANS/contrasts, which must be fit **once on the full dataset** (they borrow strength across arms via the covariance and covariates, and the contrasts are cross-column differences). Therefore:

- **Descriptive block** (the legacy `summary_data`: N, Mean (SD), Median (Range)) → native tplyr2 `group_desc()`.
- **Model block** (LS Means (SE), Diff of LS Means (SE), 95% CI, p-values, dose-response, CMH p) → fit once, format into a data.frame that **matches the tplyr2 schema**, then **row-bind** below the descriptive build. This mirrors exactly the legacy `bind_rows(column_headers, summary_portion, model_portion)`.

### 4.1 Descriptive block via tplyr2

```r
library(tplyr2)
spec <- tplyr_spec(
  cols   = "TRTP",
  layers = tplyr_layers(
    group_desc("AVAL", by = "Week 24",
      settings = layer_settings(
        format_strings = list(
          "n"              = f_str("xx",              "n"),
          "Mean (SD)"      = f_str("xx.x (xx.xx)",    "mean", "sd"),
          "Median (Range)" = f_str("xx.x (xx;xx)",    "median", "min", "max")
        )))
  ))
desc <- tplyr_build(spec, adas %>% filter(AVISITN == 24))
# -> rowlabel1, rowlabel2, res1, res2, res3, ord_layer_1, ord_layer_index
```

### 4.2 Model block: emmeans → schema-matched rows

Turn the mmrm/emmeans objects (Section 3.1) into a `res1/res2/res3` frame. Confirm which treatment each `resN` maps to (tplyr2 orders columns by the sorted/factor order of `TRTP`); build the model frame in the **same** order.

```r
fmt <- function(x, d = 1) formatC(round(x, d), format = "f", digits = d)

# LS Means (SE) row -> one cell per arm
lsm_row <- lsm_wk24 %>%
  transmute(TRTPCD_F,
            cell = sprintf("%s (%s)", fmt(emmean,1), fmt(SE,2))) %>%
  tidyr::pivot_wider(names_from = TRTPCD_F, values_from = cell)

model_block <- tibble::tibble(
  rowlabel1 = "Model",
  rowlabel2 = c("LS Means (SE)",
                "p-value(Xan - Placebo)", "  Diff of LS Means (SE)", "  95% CI",
                "p-value(Xan High - Xan Low)", "  Diff of LS Means (SE)", "  95% CI"),
  # res1 = Placebo, res2 = Xan Low, res3 = Xan High  (verify against tplyr2 column order!)
  res1 = c(lsm_row$Pbo,  "", "", "",  "", "", ""),
  res2 = c(lsm_row$Xan_Lo,
           fmt(pw_df$p.value[pw_df$contrast=="Xan_Lo - Pbo"],3),
           sprintf("%s (%s)", fmt(pw_df$estimate[pw_df$contrast=="Xan_Lo - Pbo"],1),
                              fmt(pw_df$SE[pw_df$contrast=="Xan_Lo - Pbo"],2)),
           sprintf("(%s;%s)", fmt(pw_df$lower.CL[pw_df$contrast=="Xan_Lo - Pbo"],1),
                              fmt(pw_df$upper.CL[pw_df$contrast=="Xan_Lo - Pbo"],1)),
           "", "", ""),
  res3 = c(lsm_row$Xan_Hi,
           fmt(pw_df$p.value[pw_df$contrast=="Xan_Hi - Pbo"],3),
           sprintf("%s (%s)", fmt(pw_df$estimate[pw_df$contrast=="Xan_Hi - Pbo"],1),
                              fmt(pw_df$SE[pw_df$contrast=="Xan_Hi - Pbo"],2)),
           sprintf("(%s;%s)", fmt(pw_df$lower.CL[pw_df$contrast=="Xan_Hi - Pbo"],1),
                              fmt(pw_df$upper.CL[pw_df$contrast=="Xan_Hi - Pbo"],1)),
           fmt(pw_df$p.value[pw_df$contrast=="Xan_Hi - Xan_Lo"],3),
           sprintf("%s (%s)", fmt(pw_df$estimate[pw_df$contrast=="Xan_Hi - Xan_Lo"],1),
                              fmt(pw_df$SE[pw_df$contrast=="Xan_Hi - Xan_Lo"],2)),
           sprintf("(%s;%s)", fmt(pw_df$lower.CL[pw_df$contrast=="Xan_Hi - Xan_Lo"],1),
                              fmt(pw_df$upper.CL[pw_df$contrast=="Xan_Hi - Xan_Lo"],1)))
) %>%
  mutate(ord_layer_index = max(desc$ord_layer_index) + 1,
         ord_layer_1 = row_number())

final <- dplyr::bind_rows(desc, model_block) %>%
  dplyr::arrange(ord_layer_index, ord_layer_1) %>%
  tplyr2::apply_row_masks()
```

`apply_row_masks()`, `apply_conditional_format()`, `collapse_row_labels()` are exported tplyr2 post-processors usable on the combined frame.

**Alternative for self-contained stats** (e.g. a CMH p-value that only needs one column's slice, or a geometric mean) — a `group_analyze()` layer keeps it inside the tplyr2 build. Two modes (from `vignettes/analyze.Rmd`): *format-strings mode* (return a one-row numeric data.frame; format via `layer_settings(format_strings=...)`) and *pre-formatted mode* (return `data.frame(row_label=, formatted=)`). Use this only when the statistic is genuinely computable from the per-group `.data` — not for cross-arm model contrasts.

### 4.3 clinify rendering

`clinify` (verified exports) is flextable-backed and accepts any presentation-ready data.frame; you set headers/widths by column name (`res1/res2/res3` here). From `clinify/vignettes/end_to_end.Rmd`:

```r
library(clinify)
header_n <- setNames(desc_header_ns$n, desc_header_ns$TRTP)  # N per arm

ct <- clintable(final) |>
  clin_auto_page("ord_layer_index", drop = TRUE) |>
  clin_group_pad("rowlabel1", when = "notempty") |>
  clin_column_headers(
    rowlabel1 = "", rowlabel2 = "",
    res1 = sprintf("Placebo\n(N=%s)",   header_n["Placebo"]),
    res2 = c("Xanomeline", sprintf("Low Dose\n(N=%s)",  header_n["Xanomeline Low Dose"])),
    res3 = c("Xanomeline", sprintf("High Dose\n(N=%s)", header_n["Xanomeline High Dose"]))
  ) |>
  flextable::align(j = c("res1","res2","res3"), align = "center") |>
  clin_col_widths(rowlabel1 = .17, rowlabel2 = .30, res1 = .176, res2 = .176, res3 = .176) |>
  clin_add_titles(list(
    c("Protocol: CDISCPILOT01", "Page {PAGE} of {NUMPAGES}"),
    c("Table 14-3.11"),
    c("ADAS-Cog (11) - Change from Baseline (MMRM)"))) |>
  clin_add_footnotes(list(c("Source: programs/t-14-3-11.R",
                            format(Sys.time(), "%H:%M %A, %B %d, %Y"))))

write_clindoc(ct, "outputs/14-3.11.docx")
```

`clintable()`/`as_clintable()` build the object, `clin_*` helpers add paging/headers/widths/titles/footnotes, flextable functions can be piped in directly, and `write_clindoc()` (or `as_clindoc()`) emits the `.docx`. `{PAGE}`/`{NUMPAGES}` tokens are auto-detected for page numbering.

---

## 5. Migration checklist (per method)

- **MMRM (14-3.11)**: `lme4::lmer(... + (AVISITN|USUBJID))` → `mmrm(CHG ~ SITEGR1 + TRTPCD_F*AWEEKC + BASE*AWEEKC + us(AWEEKC|USUBJID), method="Kenward-Roger", vcov="Kenward-Roger-Linear")`; `emmeans(fit, ~TRTPCD_F | AWEEKC)` for LS means; `contrast(..., adjust=NULL)` + `summary(..., infer=c(TRUE,TRUE))` for diffs/CI/p. This also fixes the known SAS mismatch (true UN covariance instead of random slope).
- **ANCOVA (14-3.0x/.12)**: keep `lm` + `emmeans(weights="proportional")` + `contrast(adjust=NULL)`; keep `car::Anova(type=3)` for the numeric-dose response term.
- **CMH row-mean-scores (14-3.13, 14-6.05)**: replace forked `vcdExtra::CMHtest` with `coin::cmh_test(ordered(y) ~ trt | strata, scores=...)` (preferred, CRAN) or CRAN `vcdExtra::CMHtest(..., )["rmeans","Prob"]` if the upstream fork-bug fix is confirmed. **Numeric parity must be validated** — neither package is installed in this environment.
- **Generalized CMH (14-6.06)** and **ANOVA/chi/fisher (14-2.01, 14-5.0x)**: no change; already base R.
- **Assembly**: descriptive rows via tplyr2 `group_desc`/`group_count`; model rows fit once and row-bound in the `rowlabel*/res1..resN/ord_*` schema; render via clinify `clintable() |> clin_* |> write_clindoc()`. Reserve `group_analyze()` for statistics computable from a single per-group `.data` slice, never for cross-arm model contrasts.

Key source files: engine `programs/funcs.R` (lines cited above); ANCOVA callers `programs/t-14-3-01.R`…`t-14-3-09.R`, `t-14-3-12.R`; MMRM `programs/t-14-3-11.R`; CMH `programs/t-14-3-13.R`, `programs/t-14-6.05.R`; generalized CMH `programs/t-14-6.06.R`; demographics `programs/t-14-2-01.R`. tplyr2/clinify references: `/Users/mstackhouse/Documents/repos/tplyr2/vignettes/analyze.Rmd`, `/Users/mstackhouse/Documents/repos/clinify/vignettes/end_to_end.Rmd`.