# medication-mortality-biobank

## Introduction
Code for Perez, Mak et al., "A multinational pharmacoepidemiology atlas of prescription medication associations with
all-cause mortality across six population cohorts" (2026).

Individual-level data cannot be shared and stay in each cohort's secure environment (UK Biobank RAP, All of Us
Researcher Workbench, FinnGen Sandbox, CDARS, JMDC, NHANES). The scripts below were run inside those environments;
only the TDC Cox function ships with synthetic data and runs anywhere.

## Contents

| Path | Manuscript | Runs in |
|---|---|---|
| `TDC_Cox.R` | Time-dependent covariate (TDC) Cox screen, Models 1–3 and strata (Figure 2) | anywhere (synthetic data); each cohort |
| `active_comparator/run_ac_fast.R` | New-user active-comparator Cox models (Figure S2) | each cohort |
| `meta_analysis/01_meta_analysis.R` | Random-effects (REML) meta-analysis of AOU, UKBB, FinnGen and CDARS; Bonferroni within group | anywhere (cohort-level HRs) |
| `figures/02_figure2.R` | Figure 2, panels A–E | anywhere (cohort-level HRs) |
| `negative_controls/nco_cause_specific_finngen.R` | Cause-specific mortality as negative-control outcomes (Figure 3C) | FinnGen Sandbox |
| `figures/plot_cause_specific.R` | Figure 3C forest plot | anywhere (summary table) |
| `target_trial_emulation/aou/vzv_tte_aou.ipynb` | VZV vaccine target trial emulation vs matched unvaccinated (Figure S3, S4A) | All of Us Workbench |
| `target_trial_emulation/aou/vzv_tte_aou_ac.ipynb` + `vzv_tte_aou_ac_models.R` | VZV vaccine vs pneumococcal or Td/Tdap vaccine recipients (Figure S3, S4B–C) | All of Us Workbench |
| `target_trial_emulation/aou/e2_tte_aou.ipynb` + `e2_tte_aou_models.R` | Estradiol emulation (women aged 45–59; too few deaths to report) | All of Us Workbench |

## TDC Cox model

```
run_subset_fast(subset_df, drug_list, group_name, progressor_fn, formula)
```
Input: one row per person × drug with `ID, Gender, YearEnrolled, AgeEnroll, Time, Status, Normalized, TimeToDrug`
(time in years from enrollment). Follow-up before the first prescription counts as unexposed and follow-up after it
as exposed (`survival::tmerge`). Set `MIN_DEATHS = 15` for real data.

## Meta-analysis and Figure 2

```
Rscript meta_analysis/01_meta_analysis.R ResultsByCohort.csv MetaAnalysis_recomputed.csv
Rscript figures/02_figure2.R ResultsByCohort.csv MetaAnalysis_recomputed.csv Figure2.pdf
```
`ResultsByCohort.csv` holds the cohort-level output of `TDC_Cox.R` (`drug, HR, conf.low, conf.high, p.value, p.adj,
Group, Bio`).

## Target trial emulation (All of Us)
The Python notebooks build the cohort in BigQuery (`WORKSPACE_CDR`) and write a model file inside the Workbench; the R
scripts fit weighted Cox models (`survival`), weighted Kaplan–Meier risks and bootstrap CIs. Set `COMPARATOR` at the
top of each R script. Row-level model files must not leave the Workbench, and counts of 1–20 must not be exported
(All of Us Data and Statistics Dissemination Policy).

## Requirements
R ≥ 4.2 with `survival`, `dplyr`, `purrr`, `broom`, `progressr`, `data.table`, `ggplot2`, `patchwork`, and `meta` or
`metafor` (meta-analysis); Python 3 with `pandas`, `numpy`, `statsmodels`, `lifelines`, `matplotlib` and
`google-cloud-bigquery` (All of Us).

## Licence and citation
Code is released under the MIT licence (`LICENSE`). Citation metadata for the archived release is in
`CITATION.cff`.

## Correspondence

Kevin Perez\
Institute of Oncology Research (IOR), Bellinzona, Switzerland\
email: kevin.perez@ior.usi.ch

## Summary results
`results/ResultsByCohort.csv` – cohort-level TDC Cox HRs (95% CI) for every drug, model and stratum.
`results/MetaAnalysis.csv` – pooled random-effects estimates (k cohorts, I², τ, Q).
`results/aou/` – target trial emulation HRs, 3-year risks and covariate balance (standardized mean differences).
No participant counts are included.
