library(dplyr)
library(purrr)
library(survival)
library(broom)

# New-user active-comparator design (Figure S2).
# Patients initiating drug A or comparator B are followed from first initiation (index) until death or end of
# follow-up. Same-day starters are dropped; patients who later start the other drug stay in their first arm
# ("itt") or are censored at the switch ("censor"). Cox model for A vs B adjusted for age at index, time from
# cohort entry to index and sex. Requires >= min_events deaths in each arm.
#
# Input: same data.frame as TDC_Cox.R (one row per person x drug):
#   ID, Gender, YearEnrolled, AgeEnroll, Time, Status, Normalized, TimeToDrug (years from enrollment)
# Output includes per-arm counts (n_A, n_B, ev_A, ev_B): check cohort small-cell rules before exporting them
# (All of Us: no counts 1-20; FinnGen: no counts < 5).

run_ac_fast <- function(subset_df,
                        drug_list,
                        comparator,
                        group_name,
                        progressor_fn,
                        formula    = NULL,
                        on_switch  = c("itt", "censor"),
                        min_events = 15,
                        time_unit  = 365.25) {   # scale of AgeIndex/LagEnroll terms; the HR for arm is unaffected

  on_switch <- match.arg(on_switch)

  baseline <- subset_df %>%
    distinct(ID, Gender, YearEnrolled, AgeEnroll, Time, Status, .keep_all = TRUE)

  # first initiation of each drug, per patient
  first_use <- subset_df %>%
    group_by(ID, Normalized) %>%
    summarise(t_first = min(TimeToDrug, na.rm = TRUE), .groups = "drop")

  comp <- first_use %>%
    filter(Normalized == comparator) %>%
    select(ID, t_comp = t_first)

  map_dfr(setdiff(drug_list, comparator), function(drug) {
    progressor_fn()

    trt <- first_use %>%
      filter(Normalized == drug) %>%
      select(ID, t_trt = t_first)

    dat <- baseline %>%
      left_join(trt,  by = "ID") %>%
      left_join(comp, by = "ID") %>%
      filter(!is.na(t_trt) | !is.na(t_comp)) %>%          # initiators of A or B only
      filter(!(!is.na(t_trt) & !is.na(t_comp) &
                 t_trt == t_comp)) %>%                      # drop same-day starters
      mutate(
        first_is_trt = coalesce(t_trt, Inf) < coalesce(t_comp, Inf),
        index        = pmin(t_trt, t_comp, na.rm = TRUE),
        t_other      = if_else(first_is_trt, t_comp, t_trt),   # switch time, if any
        arm          = factor(if_else(first_is_trt, "A", "B"), levels = c("B", "A"))
      ) %>%
      mutate(
        stop_time = if (on_switch == "censor")
          pmin(Time, coalesce(t_other, Inf)) else Time,
        status    = if_else(stop_time >= Time, as.integer(Status), 0L),
        fu        = stop_time - index,
        AgeIndex  = AgeEnroll + index / time_unit,
        LagEnroll = index / time_unit
      ) %>%
      filter(fu > 0)

    # require events in BOTH arms, not just overall
    n_ev <- c(A = sum(dat$status[dat$arm == "A"]),
              B = sum(dat$status[dat$arm == "B"]))
    if (any(n_ev < min_events)) return(NULL)

    f <- if (!is.null(formula) && n_distinct(dat$Gender) > 1) {
      as.formula(formula)
    } else if (n_distinct(dat$Gender) > 1) {
      as.formula("Surv(fu, status) ~ arm + AgeIndex + LagEnroll + Gender")
    } else {
      as.formula("Surv(fu, status) ~ arm + AgeIndex + LagEnroll")
    }

    model <- coxph(f, data = dat)
    res   <- tidy(model, conf.int = TRUE) %>% filter(term == "armA")

    tibble(
      drug       = drug,
      comparator = comparator,
      n_A        = sum(dat$arm == "A"),
      n_B        = sum(dat$arm == "B"),
      ev_A       = n_ev[["A"]],
      ev_B       = n_ev[["B"]],
      HR         = exp(res$estimate),
      conf.low   = exp(res$conf.low),
      conf.high  = exp(res$conf.high),
      p.value    = res$p.value,
      Group      = group_name
    )
  })
}

# Comparisons in the manuscript (Figure S2): target drug vs each comparator
comparisons <- list(
  "vzv vaccine" = c("pneumococcal vaccine", "influenza vaccine", "hepatitis A vaccine", "hepatitis B vaccine",
                    "DTP vaccine"),
  "DTP vaccine" = c("pneumococcal vaccine", "influenza vaccine", "hepatitis A vaccine", "hepatitis B vaccine",
                    "vzv vaccine"),
  "estradiol"   = c("estriol", "tibolone", "conjugated estrogens")
)

# res <- imap_dfr(comparisons, function(comps, target) {
#   map_dfr(comps, function(comp) {
#     run_ac_fast(df, drug_list = target, comparator = comp, group_name = "All",
#                 progressor_fn = function() invisible(NULL), on_switch = "itt")
#   })
# })
# Pooled estimates: random-effects meta-analysis across cohorts, as in meta_analysis/01_meta_analysis.R.
