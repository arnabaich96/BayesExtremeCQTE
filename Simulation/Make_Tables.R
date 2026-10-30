# LaTeX tables of the simulation study (Table 1, Supplementary Tables S2 and S3). Bias = estimate - truth;
# "Overall" averages the 9 tau x 9 profile cells, "Tail" the 5 x 9 cells with tau > 0.95.
library(here)
library(dplyr)

tab_dir <- here("Results/Simulation/Tables")
methods <- c(IPW = "IPW", ExtremeQR = "Extreme QR", SkewT = "Skew-$t$",
             MixGPD_gamma = "MixGammaGPD", MixGPD_normal = "MixInvNormGPD", MixGPD_lognormal = "MixLogNormGPD")

read_method <- function(method, tag) {
  f <- if (startsWith(method, "MixGPD_")) {
    file.path(tab_dir, sprintf("MixGPD_%s_%s.csv", tag, sub("MixGPD_", "", method)))
  } else {
    file.path(tab_dir, sprintf("%s_%s.csv", method, tag))
  }
  if (!file.exists(f)) return(NULL)
  read.csv(f)
}

summarize_cells <- function(df) {
  tail <- df$y_index > 0.95
  c(cov_all = mean(df$coverage), cov_tail = mean(df$coverage[tail]),
    bias_all = mean(df$bias),    bias_tail = mean(df$bias[tail]),
    sd_all = mean(df$sd),        sd_tail = mean(df$sd[tail]),
    len_all = mean(df$length),   len_tail = mean(df$length[tail]))
}

method_rows <- function(tag) {
  rows <- lapply(names(methods), function(m) {
    df <- read_method(m, tag)
    if (is.null(df)) return(NULL)
    data.frame(Method = methods[[m]], t(summarize_cells(df)), check.names = FALSE)
  })
  do.call(rbind, rows)
}

# averages of 10^4 or more in absolute value (Skew-t at extreme levels) are printed as bounds
fmt_num <- function(x) {
  ifelse(x >= 1e4, "$>10^4$", ifelse(x <= -1e4, "$<-10^4$", sprintf("%.3f", x)))
}
fmt_rows <- function(tab, first_cols = "Method") {
  num <- setdiff(names(tab), first_cols)
  apply(tab, 1, function(r) {
    paste0(paste(r[first_cols], collapse = " & "), " & ",
           paste(fmt_num(as.numeric(r[num])), collapse = " & "), " \\\\")
  })
}

write_tex <- function(lines, file) {
  writeLines(lines, file.path(tab_dir, file))
  cat("Wrote", file, "\n")
}

# ---------- Table 1: S1, n = 500 ----------
t1 <- method_rows("spliced_500")
if (!is.null(t1)) write_tex(fmt_rows(t1), "table1_spliced_500.tex")

# ---------- Misspecified scenarios S2 and S3 ----------
mis <- lapply(c(S2 = "burr_500", S3 = "lnmix_500"), method_rows)
mis_lines <- unlist(lapply(names(mis), function(s) {
  if (is.null(mis[[s]])) return(NULL)
  lab <- if (s == "S2") "S2 (Burr)" else "S3 (lognormal mixture)"
  c(sprintf("\\multicolumn{9}{l}{\\textit{%s}} \\\\", lab), fmt_rows(mis[[s]]), "\\midrule")
}))
if (length(mis_lines)) write_tex(head(mis_lines, -1), "table_misspec_500.tex")

# ---------- Metrics recomputed from per-replicate files (common subsets of replicates) ----------
rep_metrics <- function(tag, kern, reps, truth) {
  cells <- lapply(reps, function(i) {
    res <- readRDS(here(sprintf("Simulation/output/qte_%s_%s_rep%03d.rds", tag, kern, i)))$qte
    do.call(rbind, lapply(seq_along(res), function(q) {
      data.frame(rep = i, q = q, x = seq_len(nrow(res[[q]])), res[[q]], truth = truth[, q])
    }))
  })
  d <- do.call(rbind, cells)
  d %>% group_by(q, x) %>%
    summarise(coverage = mean(lower <= truth & truth <= upper), bias = mean(estimate - truth),
              sd = sd(estimate), length = mean(length.ci), .groups = "drop") %>%
    mutate(y_index = as.numeric(colnames(truth))[q])
}
n_reps_done <- function(tag, kern) {
  sum(file.exists(here(sprintf("Simulation/output/qte_%s_%s_rep%03d.rds", tag, kern, 1:100))))
}

# ---------- Threshold sensitivity (MixGammaGPD, S1, first 50 replicates) ----------
if (n_reps_done("spliced_500_lu50", "gamma") > 0) {
  truth <- readRDS(here("Synthetic_Data/Data_spliced_500.rds"))$True_CQTE
  thr_tags <- c("50th" = "spliced_500_lu50", "60th (main)" = "spliced_500",
                "70th" = "spliced_500_lu70", "80th" = "spliced_500_lu80")
  thr_tags <- thr_tags[sapply(thr_tags, n_reps_done, kern = "gamma") > 0]
  reps <- seq_len(min(sapply(thr_tags, n_reps_done, kern = "gamma")))
  thr_tab <- do.call(rbind, lapply(names(thr_tags), function(l) {
    data.frame(Level = l, t(summarize_cells(rep_metrics(thr_tags[[l]], "gamma", reps, truth))), check.names = FALSE)
  }))
  write_tex(fmt_rows(thr_tab, "Level"), "table_threshold_sim.tex")
  write.csv(thr_tab, file.path(tab_dir, "threshold_sensitivity_sim.csv"), row.names = FALSE)
  cat("Threshold sensitivity over", length(reps), "replicates\n"); print(thr_tab, digits = 3)
}

# ---------- Plug-in vs propagated PS on the same replicates (MixGammaGPD, S1) ----------
ps_file <- here("Simulation/output/qte_spliced_500_ps5_gamma_rep001.rds")
if (file.exists(ps_file)) {
  truth <- readRDS(here("Synthetic_Data/Data_spliced_500.rds"))$True_CQTE
  reps  <- seq_len(sum(file.exists(here(sprintf("Simulation/output/qte_spliced_500_ps5_gamma_rep%03d.rds", 1:100)))))
  ps_tab <- rbind(
    data.frame(PS = "Plug-in", t(summarize_cells(rep_metrics("spliced_500", "gamma", reps, truth)))),
    data.frame(PS = "Propagated ($K=5$)", t(summarize_cells(rep_metrics("spliced_500_ps5", "gamma", reps, truth)))))
  write_tex(fmt_rows(ps_tab, "PS"), "table_ps_propagation_sim.tex")
  write.csv(ps_tab, file.path(tab_dir, "ps_propagation_sim.csv"), row.names = FALSE)
  cat("PS comparison over", length(reps), "replicates\n"); print(ps_tab)
}

for (x in list(t1 = t1, S2 = mis$S2, S3 = mis$S3)) if (!is.null(x)) print(x, digits = 3)
