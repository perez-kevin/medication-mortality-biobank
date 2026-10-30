# nco_cause_specific_finngen.R
# Cause-specific death as negative-control outcomes, using your existing run_subset_fast().
#
# Trick: no new model code is needed. For each cause we redefine Status so that
#   Status = 1 only if the person died of that cause,
# while Time stays the same. Deaths from other causes are then censored at their death
# date - this is the standard cause-specific hazard model.
#
# Requires in memory (from TDC_Cox.R):
#   df               your FinnGen analysis data (ID, Gender, YearEnrolled, AgeEnroll, Time, Status,
#                    Normalized, TimeToDrug)
#   run_subset_fast  your function, unchanged
# and a death table with FINNGENID and the ICD-10 cause of death.

library(dplyr); library(purrr); library(progressr)

# ---- settings --------------------------------------------------------------------
death <- read.csv("death_table.csv")          # <- your death table
ID_COL  <- "FINNGENID"                         # id column in the death table
COD_COL <- "COD"                               # ICD-10 underlying cause column
drugs   <- c("vzv vaccine", "atovaquone / proguanil", "typhoid vaccine", "estradiol")
formula <- "Surv(tstart, tstop, status) ~ using + AgeEnroll + Gender + YearEnrolled"  # Model 1

# ICD-10 regexes (codes with or without dot, e.g. W19, W19.9, W199)
causes <- c(
  "All-cause (reference)"     = ".*",
  "External causes (V01-Y98)" = "^[VWXY]",
  "Falls (W00-W19)"           = "^W[01][0-9]",
  "Transport (V01-V99)"       = "^V[0-9]",
  "Suicide (X60-X84)"         = "^X[67][0-9]|^X8[0-4]",
  "Cancer (C00-D48)"          = "^C|^D[0-4]",
  "Circulatory (I00-I99)"     = "^I"
)
# ----------------------------------------------------------------------------------

death <- death %>%
  transmute(ID = .data[[ID_COL]], COD = toupper(trimws(.data[[COD_COL]]))) %>%
  distinct(ID, .keep_all = TRUE)               # one underlying cause per person

d <- df %>% left_join(death, by = "ID")

# sanity checks
pp <- d %>% distinct(ID, Status, COD)
cat("Deaths (Status = 1):", sum(pp$Status == 1),
    "| with a cause of death:", sum(pp$Status == 1 & !is.na(pp$COD) & pp$COD != ""),
    "| cause recorded but Status = 0:", sum(pp$Status == 0 & !is.na(pp$COD)), "\n")
for (k in names(causes)[-1])
  cat(sprintf("  %-28s %6d deaths\n", k, sum(pp$Status == 1 & grepl(causes[[k]], pp$COD))))

# run the usual TDC Cox once per cause
res <- imap_dfr(causes, function(rx, cause) {
  dc <- d %>% mutate(Status = as.integer(Status == 1 & !is.na(COD) & grepl(rx, COD)))
  # all-cause row keeps every death, also those without a recorded cause
  if (cause == "All-cause (reference)") dc$Status <- d$Status
  out <- with_progress({
    p <- progressor(steps = length(drugs))
    run_subset_fast(dc, drugs, cause, p, formula = formula)
  })
  # events among exposed, to judge precision
  ev <- map_int(drugs, function(dr) {
    x <- dc %>% filter(Normalized == dr, Status == 1, TimeToDrug < Time)
    n_distinct(x$ID)
  })
  out %>% left_join(tibble(drug = drugs, deaths_exposed = ev), by = "drug")
})

res <- res %>% rename(outcome = Group) %>%
  select(drug, outcome, deaths_exposed, HR, conf.low, conf.high, p.value) %>%
  arrange(drug, match(outcome, names(causes)))
print(res, n = Inf)
write.csv(res, "nco_cause_specific_finngen.csv", row.names = FALSE)

# Reading the result
#  - Falls / transport / suicide / external causes: the drugs cannot plausibly prevent these.
#    HR ~ 1  -> no healthy-user signal on this probe.
#    HR < 1, close to the all-cause HR -> the mortality association is largely selection.
#  - Typhoid vaccine is a built-in check: it should look like Malarone if both reflect
#    the healthy-traveller profile.
#  - Few deaths_exposed (< ~10) -> very wide CI; report but do not over-interpret.
#  - MIN_DEATHS in run_subset_fast: set to 0 or 5 for these rarer outcomes,
#    otherwise drugs with < 15 cause-specific deaths are silently dropped.
