# 02_figure2.R
# Recreates Figure 2 (panels A-E) from the cohort results and the recomputed meta-analysis.
#   A  all meta-analysed drugs (Total, Model 1)
#   B  selected drugs: every cohort (incl. NHANES and JMDC, shown as open circles = not pooled) + meta
#   C  Total / Male / Female, Model 1 (meta)
#   D  Models 1-3 (meta)
#   E  pre-existing conditions (meta)
# Selected drugs = lowest pooled HR in "Total, Model 1", drug combinations removed,
# list stops at STOP_AT (estradiol).
#
# Usage: Rscript 02_figure2.R ResultsByCohort.csv MetaAnalysis_recomputed.csv Figure2.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args       <- commandArgs(trailingOnly = TRUE)
cohort_csv <- if (length(args) >= 1) args[1] else "ResultsByCohort.csv"
meta_csv   <- if (length(args) >= 2) args[2] else "MetaAnalysis_recomputed.csv"
out_pdf    <- if (length(args) >= 3) args[3] else "Figure2.pdf"

# ---- settings -----------------------------------------------------------------
REF_GROUP      <- "Total, Model 1"
STOP_AT        <- "estradiol"                  # last drug shown on panels B-E
NOT_POOLED     <- c("NHANES", "JMDC")           # shown in B, not in the meta-analysis
COMBO_PATTERN  <- " / | and |\\+|,"             # how combinations are recognised
COMBO_KEEP     <- c("atovaquone / proguanil")  # single products kept despite matching
MAX_P_NOMINAL  <- 1                            # e.g. 0.05 to also require nominal significance
XLIM_MAIN      <- c(-2, 0.5)                   # log(HR) range for panels B-E (wide = CIs look smaller)
XLIM_SMALL     <- c(-1, 1)                     # log(HR) range for panels F-G
STATINS        <- c("atorvastatin", "rosuvastatin")
FOCUS_HRT      <- c("estradiol", "tibolone")

display_name <- c(
  "vzv vaccine"            = "Shingles vaccine",
  "atovaquone / proguanil" = "Malarone",
  "dtp vaccine"            = "DTP vaccine",
  "hepatitis a vaccine"    = "hepatitis A vaccine",
  "hepatitis b vaccine"    = "hepatitis B vaccine"
)
lab <- function(x) ifelse(x %in% names(display_name), display_name[x], x)

cohort_cols <- c(Meta = "black", JMDC = "#00BFC4", FINNGEN = "#7CAE00", CDARS = "#C49A00",
                 NHANES = "#C77CFF", AOU = "#F8766D", UKBB = "#FF61CC")
# ------------------------------------------------------------------------------

coh  <- read.csv(cohort_csv, stringsAsFactors = FALSE)
meta <- read.csv(meta_csv,   stringsAsFactors = FALSE)

is_combo <- function(d) grepl(COMBO_PATTERN, d) & !d %in% COMBO_KEEP

# ---- drug selection -----------------------------------------------------------
ref <- meta %>%
  filter(Group == REF_GROUP, !is_combo(Drug), p.value <= MAX_P_NOMINAL) %>%
  arrange(HR)
if (!STOP_AT %in% ref$Drug) stop(STOP_AT, " not in the reference meta-analysis")
sel <- ref$Drug[seq_len(which(ref$Drug == STOP_AT))]
message("Drugs shown: ", paste(sel, collapse = ", "))

# ---- shared forest-plot helper -------------------------------------------------
forest <- function(df, colour_var, colours, title, dodge = 0.13,
                   shape_var = NULL, legend_rows = 2, drugs = sel,
                   key_labels = function(x) sub(", Model 1$", "", x)) {
  lv   <- names(colours)
  ypos <- setNames(seq_along(drugs), rev(drugs))   # first drug at the top
  df   <- df %>%
    filter(Drug %in% drugs, .data[[colour_var]] %in% lv) %>%
    mutate(key = factor(.data[[colour_var]], levels = lv),
           off = (as.integer(key) - (length(lv) + 1) / 2) * dodge,
           y   = ypos[Drug] + off)
  p <- ggplot(df, aes(y = y, colour = key)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = log(conf.low), xmax = log(conf.high)), height = 0, linewidth = 0.5)
  p <- if (is.null(shape_var)) {
    p + geom_point(aes(x = log(HR)), size = 1.8)
  } else {
    p + geom_point(aes(x = log(HR), shape = .data[[shape_var]], size = .data[[shape_var]]),
                   stroke = 0.8) +
      scale_shape_manual(values = c(Meta = 16, Pooled = 16, `Not pooled` = 1), guide = "none") +
      scale_size_manual(values = c(Meta = 2.8, Pooled = 1.8, `Not pooled` = 1.8), guide = "none")
  }
  p +
    scale_colour_manual(values = colours, drop = TRUE, name = NULL, labels = key_labels,
                        guide = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    scale_y_continuous(breaks = ypos, labels = lab(names(ypos)),
                       expand = expansion(add = 0.6)) +
    labs(title = title, x = "log(HR), mortality", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7))
}

# ---- A: all drugs ---------------------------------------------------------------
a_df <- meta %>% filter(Group == REF_GROUP) %>% arrange(HR) %>%
  mutate(rank = rev(row_number()), shown = Drug %in% sel)
pA <- ggplot(a_df, aes(y = rank)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = log(conf.low), xmax = log(conf.high)),
                 height = 0, colour = "grey55", linewidth = 0.25) +
  geom_point(aes(x = log(HR)), size = 0.7) +
  geom_point(data = filter(a_df, shown), aes(x = log(HR)), colour = "#D55E00", size = 1.4) +
  scale_y_continuous(expand = expansion(add = 2)) +
  labs(title = "All drugs", x = "log(HR), mortality",
       caption = "orange: drugs shown in B-E",
       y = paste0("drugs (n = ", nrow(a_df), ")")) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5, size = 10),
        plot.caption = element_text(size = 6.5, hjust = 0))

# ---- B: by cohort + meta --------------------------------------------------------
b_df <- bind_rows(
  coh %>% filter(Group == REF_GROUP) %>%
    transmute(Drug = drug, HR, conf.low, conf.high, src = Bio,
              pooled = ifelse(Bio %in% NOT_POOLED, "Not pooled", "Pooled")),
  meta %>% filter(Group == REF_GROUP) %>%
    transmute(Drug, HR, conf.low, conf.high, src = "Meta", pooled = "Meta")
)
b_cols <- cohort_cols[c("Meta", sort(unique(b_df$src[b_df$Drug %in% sel & b_df$src != "Meta"])))]
pB <- forest(b_df, "src", b_cols, "All biobanks", dodge = 0.11, shape_var = "pooled",
             key_labels = identity) +
  labs(caption = "open circles: not included in meta-analysis") +
  theme(plot.caption = element_text(size = 6.5, hjust = 0))

# ---- C: sex ---------------------------------------------------------------------
c_cols <- c(`Total, Model 1` = "black", `Male, Model 1` = "#00BA38", `Female, Model 1` = "#619CFF")
pC <- forest(meta, "Group", c_cols, "Sex-specific", legend_rows = 1)

# ---- D: adjustment --------------------------------------------------------------
d_cols <- c(`Total, Model 1` = "black", `Total, Model 2` = "#00BA38", `Total, Model 3` = "#619CFF")
pD <- forest(meta, "Group", d_cols, "Covariate adjustment", legend_rows = 1,
             key_labels = function(x) sub("^Total, ", "", x)) +
  labs(caption = "Model 1: sex + age; Model 2: + CCI;\nModel 3: + smoking, BMI, SES") +
  theme(plot.caption = element_text(size = 6.5, hjust = 0))

# ---- E: pre-existing conditions -------------------------------------------------
e_cols <- c(`Total, Model 1` = "black", `Cancer, Model 1` = "#B79F00", `CVD, Model 1` = "#00BA38",
            `Diabetes, Model 1` = "#00BFC4", `Smoker, Model 1` = "#619CFF",
            `Healthy, Model 1` = "#F564E3")
pE <- forest(meta, "Group", e_cols, "Pre-existing conditions", dodge = 0.12)

# ---- F-G: statins and estradiol across adjustment models (meta) ----------------
model_lab <- function(x) sub("^Total, ", "", x)
pF <- forest(meta, "Group", d_cols, "Statins", drugs = STATINS, dodge = 0.2,
             legend_rows = 1, key_labels = model_lab)
pG <- forest(meta, "Group", d_cols, "Hormone therapy", drugs = FOCUS_HRT, dodge = 0.2,
             legend_rows = 1, key_labels = model_lab)

xl_main  <- coord_cartesian(xlim = XLIM_MAIN)
xl_small <- coord_cartesian(xlim = XLIM_SMALL)

# warn if any CI is cut by the fixed axis
cut <- bind_rows(b_df, meta) %>% filter(Drug %in% sel) %>%
  filter(log(conf.low) < XLIM_MAIN[1] | log(conf.high) > XLIM_MAIN[2])
if (nrow(cut) > 0) message(nrow(cut), " CI(s) extend beyond XLIM_MAIN and are truncated")

# panel A spans all rows; F-G are small panels at the bottom
design <- "ABC\nADE\nAFG"
fig <- wrap_plots(A = pA, B = pB + xl_main, C = pC + xl_main,
                  D = pD + xl_main, E = pE + xl_main,
                  F = pF + xl_small, G = pG + xl_small, design = design,
                  heights = c(1, 1, 0.42)) +
  plot_annotation(tag_levels = "A")

ggsave(out_pdf, fig, width = 12, height = 10.5, device = cairo_pdf)
ggsave(sub("\\.pdf$", ".png", out_pdf), fig, width = 12, height = 10.5, dpi = 200)
ggsave(sub("\\.pdf$", "_panelE.pdf", out_pdf), pE + xl_main, width = 4.2, height = 4.6, device = cairo_pdf)  # Figure2_panelE
ggsave(sub("\\.pdf$", "_panelE.png", out_pdf), pE + xl_main, width = 4.2, height = 4.6, dpi = 250)
message("Saved ", out_pdf)
