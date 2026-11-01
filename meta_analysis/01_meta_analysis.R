# 01_meta_analysis.R
# Random-effects meta-analysis of cohort-level Cox HRs (ResultsByCohort.csv).
#  - Inverse-variance pooling of log(HR); SE derived from the 95% CI
#  - Between-cohort variance by REML (as in Methods)
#  - NHANES and JMDC excluded from pooling (kept in the cohort file / plots)
#  - Pooled only if the drug was analysed in >= 2 of the remaining cohorts
#  - Bonferroni correction within each analysis group (Group), over the
#    number of drugs meta-analysed in that group
#
# Usage: Rscript 01_meta_analysis.R ResultsByCohort.csv MetaAnalysis_recomputed.csv

suppressPackageStartupMessages({
  library(dplyr)
})

args     <- commandArgs(trailingOnly = TRUE)
in_file  <- if (length(args) >= 1) args[1] else "ResultsByCohort.csv"
out_file <- if (length(args) >= 2) args[2] else "MetaAnalysis_recomputed.csv"

EXCLUDE_FROM_META <- c("NHANES", "JMDC")
MIN_COHORTS       <- 2
Z                 <- qnorm(0.975)

res <- read.csv(in_file, stringsAsFactors = FALSE) %>%
  mutate(
    logHR = log(HR),
    se    = (log(conf.high) - log(conf.low)) / (2 * Z)
  )

# sanity checks on input
stopifnot(all(res$HR > 0), all(res$conf.low <= res$HR), all(res$HR <= res$conf.high))
dups <- res %>% count(drug, Group, Bio) %>% filter(n > 1)
if (nrow(dups) > 0) stop("Duplicate drug x Group x cohort rows in input")

use_meta <- requireNamespace("meta", quietly = TRUE)
if (!use_meta) {
  if (!requireNamespace("metafor", quietly = TRUE))
    stop("Install either the 'meta' or the 'metafor' package")
  message("Package 'meta' not found: using metafor::rma (same REML random-effects model)")
}

pool_one <- function(y, se) {
  if (use_meta) {
    m <- meta::metagen(TE = y, seTE = se, sm = "HR", common = FALSE,
                       random = TRUE, method.tau = "REML")
    data.frame(logHR = m$TE.random, se = m$seTE.random,
               lo = m$lower.random, hi = m$upper.random,
               p.value = m$pval.random, I2 = m$I2, tau = m$tau,
               Q = m$Q, p.Q = m$pval.Q)
  } else {
    m <- metafor::rma(yi = y, sei = se, method = "REML",
                      control = list(stepadj = 0.5, maxiter = 1000))
    data.frame(logHR = as.numeric(m$b), se = m$se,
               lo = m$ci.lb, hi = m$ci.ub,
               p.value = m$pval, I2 = m$I2 / 100, tau = sqrt(m$tau2),
               Q = m$QE, p.Q = m$QEp)
  }
}

to_pool <- res %>%
  filter(!Bio %in% EXCLUDE_FROM_META) %>%
  group_by(drug, Group) %>%
  filter(n() >= MIN_COHORTS) %>%
  ungroup()

meta_res <- to_pool %>%
  group_by(drug, Group) %>%
  group_modify(~ cbind(
    pool_one(.x$logHR, .x$se),
    k        = nrow(.x),
    cohorts  = paste(sort(.x$Bio), collapse = ";")
  )) %>%
  ungroup() %>%
  mutate(HR = exp(logHR), conf.low = exp(lo), conf.high = exp(hi)) %>%
  group_by(Group) %>%
  mutate(n_tests = n(),
         p.adj   = pmin(1, p.value * n_tests)) %>%
  ungroup() %>%
  select(Drug = drug, Group, k, cohorts, HR, conf.low, conf.high,
         p.value, p.adj, n_tests, I2, tau, Q, p.Q) %>%
  arrange(Group, HR)

# Pooled estimates must lie within the range of the cohort estimates
chk <- meta_res %>%
  inner_join(to_pool %>% group_by(drug, Group) %>%
               summarise(mn = min(HR), mx = max(HR), .groups = "drop"),
             by = c("Drug" = "drug", "Group")) %>%
  filter(HR < mn * 0.999 | HR > mx * 1.001)
if (nrow(chk) > 0) warning(nrow(chk), " pooled HRs outside cohort range - check input")

write.csv(meta_res, out_file, row.names = FALSE)

cat("Cohorts pooled: ", paste(setdiff(unique(res$Bio), EXCLUDE_FROM_META), collapse = ", "), "\n")
cat("Rows written:   ", nrow(meta_res), "\n")
print(meta_res %>% count(Group, name = "drugs_meta_analysed"))
