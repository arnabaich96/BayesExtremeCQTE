# MixGpd.R
# Simulation study, MixGPD estimators: fits the two-stage MixGPD model to every replicate of a
# simulated data set (Functions/Report.R) and writes the summary over replicates to
# Results/Simulation/Tables/MixGPD_<scenario>_<N><suffix>_<kern>.csv.
# Per-replicate results are kept in Simulation/output, so an interrupted run resumes.
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(brms)
})
source(here("Functions/Report.R"))

# scenario:   data-generating scenario ("spliced" = S1, "burr" = S2, "lnmix" = S3)
# u_lower_q:  quantile level of the arm-specific threshold lower bound l_u
# n_ps_draws: 0 = plug-in PS; K > 0 = propagate PS uncertainty with K posterior PS draws
# suffix:     appended to the output tag (e.g. "_lu70", "_ps5")
# MCMC defaults are those of the article: 3 chains of 5000 iterations, 1000 burn-in, thinning 3.
run_simulation <- function(N, kern, scenario = "spliced", u_lower_q = 0.60, n_ps_draws = 0,
                           n_rep = 100, suffix = "", n_cores = 24,
                           n_iter = 5000, n_burnin = 1000, n_thin = 3, n_chains = 3, seed = 20) {

  CP = 0.90 ; n_mix = 5 ; epsilon = 0.01
  p <- c(0.75, 0.5, 0.25, 0.05, 10/N ,5/N, 1/(N^1.25), 1/(N * log(N)), 1/(N^1.5))
  p_QTE <- 1 - p
  tag <- paste0(scenario, "_", N, suffix)
  # Load data
  DATA <- readRDS(here(sprintf("Synthetic_Data/Data_%s_%d.rds", scenario, N)))
  DATA$Y <- DATA$Y[, 1:n_rep, drop = FALSE]
  DATA$T <- DATA$T[, 1:n_rep, drop = FALSE]

  assign("DATA", DATA, envir = .GlobalEnv)
  assign("p_QTE", p_QTE, envir = .GlobalEnv)

  Final <- Report(DATA, n_mix, p_QTE, epsilon, CP, mcmc.plot = FALSE, kern,
                  n_iter, n_burnin, n_chains, n_thin, n_cores, n_cores, seed,
                  tag = tag, u_lower_q = u_lower_q, n_ps_draws = n_ps_draws)

  dir.create(here("Results/Simulation/Tables"), showWarnings = FALSE, recursive = TRUE)
  write.csv(Final$QTE, file = here(sprintf("Results/Simulation/Tables/MixGPD_%s_%s.csv", tag, kern)), row.names = FALSE)
  invisible(Final)
}
