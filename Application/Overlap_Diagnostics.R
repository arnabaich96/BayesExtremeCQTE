#!/usr/bin/env Rscript
# Overlap_Diagnostics.R
# Propensity-score overlap and covariate balance for the diabetes application.
# Reads Application/output/PS_main.rds (written by MixGpd_Diabetes.R): posterior-mean PS from the
# Bayesian logistic model, treatment (obesity, BMI >= 30) and covariates of the main analysis data set.
# Outputs: Results/Application/Diabetes/Figures/fig_ps_overlap_bw.png, fig_love_plot_bw.png
#          Results/Application/Diabetes/Tables/Overlap_Summary.csv, Balance_SMD.csv, table_balance.tex

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

fig_dir   <- here("Results/Application/Diabetes/Figures")
table_dir <- here("Results/Application/Diabetes/Tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ps_obj <- readRDS(here("Application/output/PS_main.rds"))
ps <- ps_obj$ps; Tr <- ps_obj$T; X <- ps_obj$X

# ---- Overlap summary ----
overlap <- data.frame(
  arm    = c("Non-obese", "Obese"),
  n      = c(sum(Tr == 0), sum(Tr == 1)),
  min_ps = c(min(ps[Tr == 0]), min(ps[Tr == 1])),
  q05_ps = c(quantile(ps[Tr == 0], 0.05), quantile(ps[Tr == 1], 0.05)),
  median = c(median(ps[Tr == 0]), median(ps[Tr == 1])),
  q95_ps = c(quantile(ps[Tr == 0], 0.95), quantile(ps[Tr == 1], 0.95)),
  max_ps = c(max(ps[Tr == 0]), max(ps[Tr == 1])),
  below_010 = c(sum(ps[Tr == 0] < 0.1), sum(ps[Tr == 1] < 0.1)),
  above_090 = c(sum(ps[Tr == 0] > 0.9), sum(ps[Tr == 1] > 0.9))
)
write.csv(overlap, file.path(table_dir, "Overlap_Summary.csv"), row.names = FALSE)
print(overlap, digits = 3)

# ---- Mirrored histogram of the PS by arm ----
breaks <- seq(0, 1, by = 0.025)
h <- bind_rows(lapply(c(0, 1), function(a) {
  cnt <- hist(ps[Tr == a], breaks = breaks, plot = FALSE)$counts
  data.frame(mid = head(breaks, -1) + 0.0125, count = if (a == 1) cnt else -cnt,
             arm = if (a == 1) "Obese (BMI ≥ 30)" else "Non-obese (BMI < 30)")
}))
p_overlap <- ggplot(h, aes(x = mid, y = count, fill = arm)) +
  geom_col(width = 0.025, colour = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = c(0.1, 0.9), linetype = "dashed", linewidth = 0.3) +
  scale_fill_grey(start = 0.35, end = 0.8, name = NULL) +
  scale_y_continuous(labels = abs) +
  labs(x = "Estimated propensity score", y = "Number of subjects") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "fig_ps_overlap_bw.png"), p_overlap, width = 6.5, height = 3.6, dpi = 300, bg = "white")

# ---- Standardized mean differences, unweighted and inverse-probability weighted ----
w <- ifelse(Tr == 1, 1 / ps, 1 / (1 - ps))
wmean <- function(x, w) sum(w * x) / sum(w)
wvar  <- function(x, w) sum(w * (x - wmean(x, w))^2) / sum(w)
smd <- function(x, w) {
  m1 <- wmean(x[Tr == 1], w[Tr == 1]); m0 <- wmean(x[Tr == 0], w[Tr == 0])
  s  <- sqrt((var(x[Tr == 1]) + var(x[Tr == 0])) / 2)    # pooled unweighted SD in both cases
  (m1 - m0) / s
}
labels <- c(stab.glu = "Stabilized glucose", ratio = "Cholesterol/HDL ratio", location = "Location (Louisa)",
            age = "Age", gender = "Gender (male)", frame = "Body frame", bp.1s = "Systolic BP",
            bp.1d = "Diastolic BP", time.ppn = "Time since last meal")
bal <- data.frame(
  covariate  = labels[names(X)],
  smd_raw    = sapply(X, smd, w = rep(1, length(Tr))),
  smd_ipw    = sapply(X, smd, w = w),
  vr_raw     = sapply(X, function(x) var(x[Tr == 1]) / var(x[Tr == 0])),
  vr_ipw     = sapply(X, function(x) wvar(x[Tr == 1], w[Tr == 1]) / wvar(x[Tr == 0], w[Tr == 0])),
  row.names  = NULL
)
write.csv(bal, file.path(table_dir, "Balance_SMD.csv"), row.names = FALSE)
print(bal, digits = 3)
ess <- function(w) sum(w)^2 / sum(w^2)
cat("Effective sample size of IPW weights: obese", round(ess(w[Tr == 1]), 1),
    "| non-obese", round(ess(w[Tr == 0]), 1), "\n")
writeLines(sprintf("%s & %.3f & %.3f & %.2f & %.2f \\\\", bal$covariate, bal$smd_raw, bal$smd_ipw, bal$vr_raw, bal$vr_ipw),
           file.path(table_dir, "table_balance.tex"))

# ---- Love plot ----
love <- bal %>%
  select(covariate, Unweighted = smd_raw, `IPW weighted` = smd_ipw) %>%
  pivot_longer(-covariate, names_to = "sample", values_to = "smd") %>%
  mutate(covariate = factor(covariate, levels = bal$covariate[order(abs(bal$smd_raw))]))
p_love <- ggplot(love, aes(x = abs(smd), y = covariate, shape = sample)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 2.3) +
  scale_shape_manual(values = c(Unweighted = 1, `IPW weighted` = 16), name = NULL) +
  labs(x = "Absolute standardized mean difference", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "fig_love_plot_bw.png"), p_love, width = 6, height = 3.8, dpi = 300, bg = "white")
cat("Overlap diagnostics written to", fig_dir, "and", table_dir, "\n")
