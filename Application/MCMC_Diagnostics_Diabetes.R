# Traceplots, effective sample sizes and R-hat of xi, sigma and the threshold u at the arm-wise median
# profile, for the main fits of the application (Supplementary Section S6).
suppressPackageStartupMessages({
  library(here)
  library(coda)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

fig_dir <- here("Results/Application/Diabetes/Figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

KEY_PARAMS <- c("xi", "sigma", "u_rep")
PARAM_LABELS <- c(xi = "\u03be", sigma = "\u03c3", u_rep = "u")

source(here("Functions/Prior_Distributions.R"))   # median_lnormTrunc()

# Arm-wise median design row (PS first, then the covariates), from the stage-1 output of MixGpd_Diabetes.R
ps_obj <- readRDS(here("Application/output/PS_main.rds"))
y_obs  <- readRDS(here("Application/data/Diabetes_Data.rds"))$glyhb
build_median_X <- function(arm) {
  X_design <- cbind(PS = ps_obj$ps, as.matrix(ps_obj$X))
  idx      <- if (arm == "Treatment") which(ps_obj$T == 1) else which(ps_obj$T == 0)
  apply(X_design[idx, , drop = FALSE], 2, median)
}
arm_l_u <- function(arm, u_lower_q = 0.60) {
  unname(quantile(y_obs[ps_obj$T == as.integer(arm == "Treatment")], u_lower_q))
}

# u_rep: conditional median of the l_u-truncated lognormal threshold at the median profile, as in Sim()
add_u_rep <- function(mcmc_obj, X_median, l_u) {
  lapply(mcmc_obj, function(ch) {
    mat  <- as.matrix(ch)
    d_cols <- grep("^delta1\\[", colnames(mat), value = TRUE)
    d_mat  <- mat[, d_cols, drop = FALSE]
    u_vals <- median_lnormTrunc(as.numeric(d_mat %*% X_median), 1 / sqrt(mat[, "var2"]), l_u)
    cbind(mat, u_rep = u_vals)
  })
}

# ---- helpers ----

mcmc_to_long <- function(mcmc_obj, params) {
  params_present <- intersect(params, colnames(mcmc_obj[[1]]))
  bind_rows(lapply(seq_along(mcmc_obj), function(ch) {
    mat <- mcmc_obj[[ch]][, params_present, drop = FALSE]
    as.data.frame(mat) |>
      mutate(iteration = seq_len(nrow(mat)), chain = factor(ch)) |>
      pivot_longer(
        cols = -c(iteration, chain),
        names_to = "parameter",
        values_to = "value"
      )
  }))
}

make_traceplot <- function(long_df, title_str) {
  long_df |>
    mutate(parameter = factor(parameter, levels = KEY_PARAMS)) |>
    ggplot(aes(x = iteration, y = value,
               color = chain, group = chain)) +
    geom_line(alpha = 0.7, linewidth = 0.4) +
    facet_wrap(
      ~parameter, scales = "free_y", ncol = 3,
      labeller = as_labeller(PARAM_LABELS)
    ) +
    scale_color_grey(start = 0.1, end = 0.7) +
    labs(title = title_str,
         x = "Iteration", y = "Value", color = "Chain") +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold", hjust = 0.5)
    )
}

compute_ess_df <- function(trt_mat, con_mat, params) {
  params_t <- intersect(params, colnames(trt_mat[[1]]))
  params_c <- intersect(params, colnames(con_mat[[1]]))
  ess_t <- effectiveSize(as.mcmc.list(lapply(trt_mat, function(m) as.mcmc(m[, params_t]))))
  ess_c <- effectiveSize(as.mcmc.list(lapply(con_mat, function(m) as.mcmc(m[, params_c]))))
  bind_rows(
    tibble(parameter = names(ess_t), ESS = as.numeric(ess_t), arm = "Treatment"),
    tibble(parameter = names(ess_c), ESS = as.numeric(ess_c), arm = "Control")
  ) |>
    mutate(parameter = factor(parameter, levels = params))
}

make_ess_plot <- function(ess_df, title_str) {
  x_labels <- setNames(
    sapply(levels(ess_df$parameter), function(p) PARAM_LABELS[[p]] %||% p),
    levels(ess_df$parameter)
  )
  ess_df |>
    ggplot(aes(x = parameter, y = ESS, fill = arm)) +
    geom_col(position = position_dodge(0.7), width = 0.6,
             color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
    scale_fill_grey(start = 0.4, end = 0.85) +
    scale_x_discrete(labels = unname(x_labels)) +
    labs(title = title_str,
         x = "Parameter", y = "Effective Sample Size", fill = "Arm") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(face = "bold", hjust = 0.5),
      axis.text.x     = element_text(angle = 30, hjust = 1)
    )
}

# ---- Precompute median covariate vectors (same PS model for both arms) ----
message("Building representative covariate profiles...")
X_trt_median <- build_median_X("Treatment")
X_con_median <- build_median_X("Control")
message("  Treatment median u profile: ", paste(round(X_trt_median, 3), collapse = ", "))
message("  Control   median u profile: ", paste(round(X_con_median, 3), collapse = ", "))

# ---- main loop over kernels ----

for (kern in c("normal", "gamma", "lognormal")) {
  message("Processing kernel: ", kern)

  trt_file <- here("Application/output", paste0(kern, "_Treatment_Samples.rds"))
  con_file <- here("Application/output", paste0(kern, "_Control_Samples.rds"))

  if (!file.exists(trt_file) || !file.exists(con_file)) {
    message("  Missing RDS for kernel ", kern, " — skipping.")
    next
  }

  trt_mcmc <- readRDS(trt_file)
  con_mcmc <- readRDS(con_file)

  # Augment chains with u_rep
  trt_aug <- add_u_rep(trt_mcmc, X_trt_median, arm_l_u("Treatment"))
  con_aug <- add_u_rep(con_mcmc, X_con_median, arm_l_u("Control"))

  # Gelman-Rubin statistics for the key parameters
  rhat <- sapply(list(Treatment = trt_aug, Control = con_aug), function(aug) {
    gelman.diag(as.mcmc.list(lapply(aug, function(m) as.mcmc(m[, KEY_PARAMS]))), autoburnin = FALSE)$psrf[, 1]
  })
  message("  R-hat:\n", paste(capture.output(print(round(rhat, 3))), collapse = "\n"))
  write.csv(data.frame(parameter = rownames(rhat), rhat, kernel = kern),
            file.path(here("Results/Application/Diabetes/Tables"), paste0("Rhat_", kern, ".csv")), row.names = FALSE)

  # --- Traceplot: treatment arm ---
  p_trace_trt <- make_traceplot(
    mcmc_to_long(trt_aug, KEY_PARAMS),
    paste0("Traceplots — Treatment arm (", kern, " kernel)")
  )
  out_trt <- file.path(fig_dir, paste0("Traceplot_Treatment_", kern, "_BW.png"))
  ggsave(out_trt, p_trace_trt, width = 9, height = 3, dpi = 300, bg = "white")
  message("  Saved: ", out_trt)

  # --- Traceplot: control arm ---
  p_trace_con <- make_traceplot(
    mcmc_to_long(con_aug, KEY_PARAMS),
    paste0("Traceplots — Control arm (", kern, " kernel)")
  )
  out_con <- file.path(fig_dir, paste0("Traceplot_Control_", kern, "_BW.png"))
  ggsave(out_con, p_trace_con, width = 9, height = 3, dpi = 300, bg = "white")
  message("  Saved: ", out_con)

  # --- ESS bar chart ---
  ess_df  <- compute_ess_df(trt_aug, con_aug, KEY_PARAMS)
  p_ess   <- make_ess_plot(ess_df, paste0("Effective Sample Size — ", kern, " kernel"))
  out_ess <- file.path(fig_dir, paste0("ESS_", kern, "_BW.png"))
  ggsave(out_ess, p_ess, width = 6, height = 4, dpi = 300, bg = "white")
  message("  Saved: ", out_ess)
}

message("--- MCMC diagnostics complete ---")
