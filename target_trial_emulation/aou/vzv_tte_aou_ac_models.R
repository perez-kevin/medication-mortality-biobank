## Active-comparator TTE (All of Us): models, risks, E-value, plot
## Input: vzv_aou_ac_<COMPARATOR>_model.csv from the Python notebook (cell 5)
library(survival); library(data.table)
COMPARATOR <- 'pneumo'          # 'pneumo' or 'tdap'
B_BOOT     <- 50                # as in the manuscript; use 20 for a quick test run
RISK_YEARS <- 3

d <- fread(sprintf('vzv_aou_ac_%s_model.csv', COMPARATOR))
d[, one := 1]
hr <- function(f, label) {
  s <- summary(f)$coefficients; r <- grep('treated', rownames(s))
  data.table(analysis = label, term = rownames(s)[r], HR = exp(s[r, 'coef']),
             lo = exp(s[r, 'coef'] - 1.96 * s[r, 'robust se']),
             hi = exp(s[r, 'coef'] + 1.96 * s[r, 'robust se']))
}
ev_n <- function(x, ev = 'death') x[, .(events = sum(get(ev))), by = treated][order(-treated), paste(events, collapse = ' vs ')]

## Hazard ratios
res <- list(
  hr(coxph(Surv(time, death) ~ treated, d, weights = one, cluster = person_id), 'Unweighted'),
  hr(coxph(Surv(time, death) ~ treated, d, weights = w0,  cluster = person_id), 'ITT, PS no SRH'),
  hr(coxph(Surv(time, death) ~ treated, d, weights = w,   cluster = person_id), 'ITT, PS + SRH (primary)'),
  hr(coxph(Surv(time, death) ~ treated + age + female + cci + lvis + ldrug + flu + travel, d,
           weights = w, cluster = person_id), 'Weighted + adjusted'),
  hr(coxph(Surv(pp_time, pp_death) ~ treated, d, weights = w, cluster = person_id), 'Per-protocol'))
sp <- survSplit(Surv(time, death) ~ treated + w + person_id, data = d, cut = c(90, 365), episode = 'period')
sp$period <- factor(sp$period, labels = c('0-90d', '91-365d', '>1y'))
res$win <- hr(coxph(Surv(tstart, time, death) ~ treated:period, sp, weights = w, cluster = person_id), 'By window')
for (vt in c('recombinant', 'live')) {
  x <- d[treated == 0 | vacc_type == vt]
  if (x[treated == 1, sum(death)] >= 20)
    res[[vt]] <- hr(coxph(Surv(time, death) ~ treated, x, weights = w, cluster = person_id), paste('ITT, zoster', vt))
}
## zoster given alone (no comparator vaccine the same day)
x <- d[treated == 0 | same_day_cmp == 0]
res$alone <- hr(coxph(Surv(time, death) ~ treated, x, weights = w, cluster = person_id), 'ITT, zoster not co-administered')
results <- rbindlist(res)
print(results, digits = 3)
cat('Deaths zoster vs comparator:', ev_n(d), '\n')
print(as.data.table(table(sp$period, sp$treated, sp$death)[, , '1']))   # deaths per window and arm

## Absolute risk (weighted KM), bootstrap by person
T <- RISK_YEARS * 365.25
risk_at <- function(x) {
  s <- summary(survfit(Surv(time, death) ~ treated, data = x, weights = w), times = T, extend = TRUE)
  r <- setNames(1 - s$surv, sub('treated=', '', s$strata))
  c(risk_cmp = r[['0']], risk_zos = r[['1']], rd = r[['1']] - r[['0']])
}
pt <- risk_at(d)
set.seed(2026)
ids <- unique(d$person_id); setkey(d, person_id)
boot <- t(replicate(B_BOOT, risk_at(d[.(sample(ids, replace = TRUE)), allow.cartesian = TRUE])))
risks <- data.table(measure = names(pt), estimate = pt,
                    lo = apply(boot, 2, quantile, .025), hi = apply(boot, 2, quantile, .975))
print(risks, digits = 4)

## E-value (deaths rare: RR ~ HR)
ev <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
p <- results[analysis == 'ITT, PS + SRH (primary)']
lim <- if (p$lo <= 1 & p$hi >= 1) 1 else ev(if (p$HR < 1) p$hi else p$lo)
cat('E-value point:', round(ev(p$HR), 2), ' CI limit:', round(lim, 2), '\n')

## Plot (5 years)
km <- survfit(Surv(time / 365.25, death) ~ treated, data = d, weights = w)
cmp_lab <- c(pneumo = 'Pneumococcal vaccine', tdap = 'Td/Tdap vaccine')[COMPARATOR]
png(sprintf('vzv_aou_ac_%s_cumulative_mortality.png', COMPARATOR), width = 1600, height = 1200, res = 200)
plot(km, fun = 'event', col = c('grey40', '#2a7f62'), lwd = 2, conf.int = FALSE, xlim = c(0, 5),
     xlab = 'Years since vaccination', ylab = 'Cumulative mortality')
legend('topleft', c(cmp_lab, 'Zoster vaccine'), col = c('grey40', '#2a7f62'), lwd = 2, bty = 'n')
dev.off()

fwrite(results, sprintf('vzv_aou_ac_%s_hazard_ratios.csv', COMPARATOR))
fwrite(risks,   sprintf('vzv_aou_ac_%s_risk.csv', COMPARATOR))
## AoU dissemination policy: no counts between 1 and 20 in anything you export.
