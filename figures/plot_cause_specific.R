# plot_cause_specific.R
# Forest plot of cause-specific mortality (FinnGen) for selected drugs.
# Usage: Rscript plot_cause_specific.R nco_cause_specific_finngen.csv cause_specific_finngen.pdf
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "nco_cause_specific_finngen.csv"
out_pdf <- if (length(args) >= 2) args[2] else "cause_specific_finngen.pdf"

d <- read.csv(in_csv, check.names = FALSE, stringsAsFactors = FALSE)
names(d) <- sub("^deaths.exposed$|^deaths exposed$", "deaths_exposed", names(d))

d <- d %>%
  filter(!is.na(HR)) %>%                                   # drop suppressed rows (N < 5)
  mutate(
    cause = sub(" \\(.*", "", outcome),
    cause = recode(cause, "All-cause" = "All-cause"),
    type  = case_when(
      cause == "All-cause"                                  ~ "All-cause",
      cause %in% c("Cancer", "Circulatory")                 ~ "Disease-specific",
      TRUE                                                  ~ "Negative control\n(not plausibly affected)"
    ),
    label = sprintf("%.2f (%.2f, %.2f)  n=%s", HR, conf.low, conf.high, deaths_exposed)
  )

order <- c("All-cause", "Circulatory", "Cancer", "External causes", "Falls", "Transport", "Suicide")
d$cause <- factor(d$cause, levels = rev(order[order %in% d$cause]))
d$type  <- factor(d$type, levels = c("All-cause", "Disease-specific",
                                     "Negative control\n(not plausibly affected)"))

p <- ggplot(d, aes(y = cause, x = log(HR), colour = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = log(conf.low), xmax = log(conf.high)), height = 0, linewidth = 0.6) +
  geom_point(size = 2.4) +
  geom_text(aes(x = 0.95, label = label), hjust = 0, size = 2.5, colour = "grey20") +
  facet_wrap(~ drug, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("All-cause" = "black", "Disease-specific" = "#D55E00",
                                 "Negative control\n(not plausibly affected)" = "#0072B2"),
                      name = NULL) +
  coord_cartesian(xlim = c(-1.8, 0.75), clip = "off") +
  labs(x = "log(HR), cause-specific mortality", y = NULL,
       title = "Cause-specific mortality (FinnGen)",
       caption = "HR (95% CI), n = deaths among exposed. Model 1; deaths from other causes censored.\nCauses with < 5 deaths among exposed not shown.") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom", legend.text = element_text(size = 7.5),
        plot.margin = margin(5, 105, 5, 5),
        plot.caption = element_text(size = 6.5, hjust = 0))

ggsave(out_pdf, p, width = 4.8, height = 4.6, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", out_pdf), p, width = 4.8, height = 4.6, dpi = 250)
message("Saved ", out_pdf)
