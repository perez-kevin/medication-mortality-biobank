## Estradiol TTE (All of Us): models, positive controls, risks, E-value, plot
## Input: e2_aou_<TAG>_model.csv from the Python notebook (cell 5)
library(survival); library(data.table)
TAG        <- 'vaginal_45_60'   # must match TAG in the notebook, e.g. 'nonuser_45_60', 'vaginal_60_80'
B_BOOT     <- 200               # 20 for a quick test
RISK_YEARS <- 3

d <- fread(sprintf('e2_aou_%s_model.csv', TAG))
d[, one := 1]
hr <- function(f, label, n_ev = NA_character_) {
  s <- summary(f)$coefficients; r <- grep('treated', rownames(s))
  data.table(analysis = label, term = rownames(s)[r], HR = exp(s[r, 'coef']),
             lo = exp(s[r, 'coef'] - 1.96 * s[r, 'robust se']),
             hi = exp(s[r, 'coef'] + 1.96 * s[r, 'robust se']), events_E2_vs_cmp = n_ev)
}
evs <- function(x, ev) x[, .(n = sum(get(ev))), by = treated][order(-treated), paste(n, collapse = ' vs ')]

## All-cause mortality
res <- list(
  hr(coxph(Surv(time, death) ~ treated, d, weights = one, cluster = person_id), 'Death: unweighted', evs(d, 'death')),
  hr(coxph(Surv(time, death) ~ treated, d, weights = w0,  cluster = person_id), 'Death: PS no SRH'),
  hr(coxph(Surv(time, death) ~ treated, d, weights = w,   cluster = person_id), 'Death: PS + SRH (primary)'),
  hr(coxph(Surv(time, death) ~ treated + age + cci + lvis + ldrug + flu + mammo + hx_hyster, d,
           weights = w, cluster = person_id), 'Death: weighted + adjusted'),
  hr(coxph(Surv(pp_time, pp_death) ~ treated, d, weights = w, cluster = person_id), 'Death: per-protocol', evs(d, 'pp_death')))
sp <- survSplit(Surv(time, death) ~ treated + w + person_id, data = d, cut = c(90, 365), episode = 'period')
sp$period <- factor(sp$period, labels = c('0-90d', '91-365d', '>1y'))
res$win <- tryCatch(hr(coxph(Surv(tstart, time, death) ~ treated:period, sp, weights = w, cluster = person_id),
                       'Death: by window'), error = function(e) NULL)

## Secondary outcomes (prevalent cases excluded per outcome; death censors)
## Expected: VTE and breast cancer HR > 1 (positive controls); CHD, stroke, hip fracture HR <= 1
for (k in c('vte', 'breast_ca', 'chd', 'stroke', 'hip_fx')) {
  x <- d[get(paste0('prev_', k)) == 0]
  setnames(x, c(paste0('t_', k), paste0('e_', k)), c('tk', 'ek'))
  if (x[, sum(ek)] >= 10)
    res[[k]] <- hr(coxph(Surv(tk, ek) ~ treated, x, weights = w, cluster = person_id), paste('Outcome:', k), evs(x, 'ek'))
}
results <- rbindlist(res, fill = TRUE)
print(results, digits = 3)

## Absolute mortality risk (weighted KM), bootstrap by person
T <- RISK_YEARS * 365.25
risk_at <- function(x) {
  s <- summary(survfit(Surv(time, death) ~ treated, data = x, weights = w), times = T, extend = TRUE)
  r <- setNames(1 - s$surv, sub('treated=', '', s$strata))
  c(risk_cmp = r[['0']], risk_e2 = r[['1']], rd = r[['1']] - r[['0']])
}
pt <- risk_at(d)
set.seed(2026)
ids <- unique(d$person_id); setkey(d, person_id)
boot <- t(replicate(B_BOOT, risk_at(d[.(sample(ids, replace = TRUE)), allow.cartesian = TRUE])))
risks <- data.table(measure = names(pt), estimate = pt,
                    lo = apply(boot, 2, quantile, .025), hi = apply(boot, 2, quantile, .975))
print(risks, digits = 4)

## E-value for the primary mortality HR
ev <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }
p <- results[analysis == 'Death: PS + SRH (primary)']
lim <- if (p$lo <= 1 & p$hi >= 1) 1 else ev(if (p$HR < 1) p$hi else p$lo)
cat('E-value point:', round(ev(p$HR), 2), ' CI limit:', round(lim, 2), '\n')

## Plot
km <- survfit(Surv(time / 365.25, death) ~ treated, data = d, weights = w)
cmp_lab <- if (grepl('^vaginal', TAG)) 'Vaginal estrogen' else 'No systemic estrogen'
png(sprintf('e2_aou_%s_cumulative_mortality.png', TAG), width = 1600, height = 1200, res = 200)
plot(km, fun = 'event', col = c('grey40', '#8a3b8f'), lwd = 2, conf.int = FALSE, xlim = c(0, 5),
     xlab = 'Years since initiation', ylab = 'Cumulative mortality')
legend('topleft', c(cmp_lab, 'Systemic estradiol'), col = c('grey40', '#8a3b8f'), lwd = 2, bty = 'n')
dev.off()

fwrite(results, sprintf('e2_aou_%s_hazard_ratios.csv', TAG))
fwrite(risks,   sprintf('e2_aou_%s_risk.csv', TAG))
## AoU dissemination policy: no counts between 1 and 20 in anything you export.
