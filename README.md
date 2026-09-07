# medication-mortality-biobank

## Introduction
In this repository we provide the code used to derive medication-mortality associations using a Time-Depedent Covariate (TDC) Cox model. 

Manuscript.
Perez et al., 2026 "A multinational pharmacoepidemiology atlas of prescription medication associations with all-cause mortality across six population cohorts".

## TDC Cox model

```
run_subset_fast(subset_df, drug_list, group_name, progressor_fn, formula)
```
The main function run_subset_fast (TDC_Cox.R) is provided with synthetic data.
