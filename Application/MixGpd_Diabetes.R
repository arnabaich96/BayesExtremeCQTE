# MixGpd_Diabetes.R
# Purpose: Run the MixGPD analysis of the diabetes data (effect of obesity, BMI >= 30, on HbA1c)
# and all sensitivity analyses reported in the article:
#   main   : imputation 1, posterior-mean PS, l_u = 60th percentile, three kernels
#   lu50/lu70/lu80 : threshold lower bound at the 50th/70th/80th percentile
#   prop   : cut posterior stacking 10 fits, fit k = imputed data set k + one posterior PS draw
#   trim   : subjects with PS in [0.1, 0.9] (Crump et al., 2009)
#   noglu  : stab.glu removed from the PS and outcome models (possible mediator)
#   xisd10 : prior sd of xi increased from 0.05 to 0.1 (prior sensitivity)
# The confounders are listed in `covariates` below; BMI, which defines the treatment, is not among them.
# Every fit has its own seed, so the results do not depend on the number of cores (CQTE_CORES).

# ------------------------- 1. Packages and Setup -------------------------
packages <- c("nimble", "here", "dplyr", "brms", "foreach", "doParallel")
invisible(lapply(packages, require, character.only = TRUE))
source(here("Functions/Report_RealData.R"))

# ------------------------- 2. Model and MCMC Definitions -------------------------
# longer chains than in the simulation: the GPD scale mixes more slowly on these data
CP = 0.90; n_mix = 5; epsilon = 0.01; n_iter = 15000; n_burnin = 5000
n_thin = 5; n_chains = 3; seed = 2024; n_imp = 10
n_cores <- as.integer(Sys.getenv("CQTE_CORES", min(24, max(1, parallel::detectCores() - 1))))

out_dir   <- here("Application/output")
draw_dir  <- file.path(out_dir, "draws")
table_dir <- here("Results/Application/Diabetes/Tables")
dir.create(draw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

covariates <- c("stab.glu", "ratio", "location", "age", "gender", "frame", "bp.1s", "bp.1d", "time.ppn")

# ------------------------- 3. Data and prediction grid -------------------------
MI <- readRDS(here("Application/data/Diabetes_Data_MI.rds"))
training_data <- MI[[1]]
N_train <- nrow(training_data)

get_mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}
make_grid <- function(dat, covs, grid_cont = c("age", "stab.glu", "ratio"), grid_cat = "gender") {
  grid_cont <- intersect(grid_cont, covs)
  grid_values <- lapply(grid_cont, function(v) unname(round(quantile(dat[[v]], c(0.25, 0.5, 0.75)), 2)))
  names(grid_values) <- grid_cont
  grid_values[[grid_cat]] <- sort(unique(dat[[grid_cat]]))
  grid <- expand.grid(grid_values)
  for (v in setdiff(covs, names(grid))) {
    grid[[v]] <- if (v %in% c("location", "frame")) get_mode(dat[[v]]) else median(dat[[v]])
  }
  grid[, covs]
}
grid_main  <- make_grid(training_data, covariates)
covs_noglu <- setdiff(covariates, "stab.glu")
grid_noglu <- make_grid(training_data, covs_noglu)
cat("Prediction grid:", nrow(grid_main), "profiles (", nrow(grid_noglu), "without stab.glu )\n")

p_QTE_diabetes <- 1 - c(0.75, 0.5, 0.25, 0.05, 10/N_train, 5/N_train, 1/(N_train^1.25),
                        1/(N_train * log(N_train)), 1/(N_train^1.5))

# ------------------------- 4. Stage 1: propensity scores -------------------------
print("--- Stage 1: Bayesian propensity-score models ---")
ps <- vector("list", n_imp)
ps[[1]] <- fit_ps_app(training_data$obese, training_data[, covariates], grid_main, n_draws = 1, seed = seed)
for (k in 2:n_imp) {
  ps[[k]] <- fit_ps_app(MI[[k]]$obese, MI[[k]][, covariates], grid_main, n_draws = 1,
                        seed = seed + k, ps_model = ps[[1]]$model)
}
ps_noglu <- fit_ps_app(training_data$obese, training_data[, covs_noglu], grid_noglu, seed = seed)
saveRDS(list(ps = ps[[1]]$train, T = training_data$obese, X = training_data[, covariates],
             ps_model = ps[[1]]$model),
        file.path(out_dir, "PS_main.rds"))

# ------------------------- 5. Stage 2 jobs -------------------------
design <- function(ps_train, ps_new, dat, grid, covs) {
  list(X = cbind(PS = ps_train, dat[, covs]), new = cbind(PS = ps_new, grid[, covs]))
}
jobs <- list()
for (kern in c("gamma", "normal", "lognormal")) {
  jobs[[length(jobs) + 1]] <- list(name = "main", kern = kern, imp = 1, u_q = 0.60,
                                   d = design(ps[[1]]$train, ps[[1]]$new, training_data, grid_main, covariates),
                                   keep = rep(TRUE, N_train), grid = grid_main)
}
for (kern in c("gamma", "normal")) {
  for (u_q in c(0.50, 0.70, 0.80)) {
    jobs[[length(jobs) + 1]] <- list(name = sprintf("lu%02d", 100 * u_q), kern = kern, imp = 1, u_q = u_q,
                                     d = design(ps[[1]]$train, ps[[1]]$new, training_data, grid_main, covariates),
                                     keep = rep(TRUE, N_train), grid = grid_main)
  }
  for (k in seq_len(n_imp)) {
    jobs[[length(jobs) + 1]] <- list(name = sprintf("prop%02d", k), kern = kern, imp = k, u_q = 0.60,
                                     d = design(ps[[k]]$train_draws[1, ], ps[[k]]$new_draws[1, ], MI[[k]], grid_main, covariates),
                                     keep = rep(TRUE, N_train), grid = grid_main)
  }
  jobs[[length(jobs) + 1]] <- list(name = "noglu", kern = kern, imp = 1, u_q = 0.60,
                                   d = design(ps_noglu$train, ps_noglu$new, training_data, grid_noglu, covs_noglu),
                                   keep = rep(TRUE, N_train), grid = grid_noglu)
  jobs[[length(jobs) + 1]] <- list(name = "xisd10", kern = kern, imp = 1, u_q = 0.60,
                                   d = design(ps[[1]]$train, ps[[1]]$new, training_data, grid_main, covariates),
                                   keep = rep(TRUE, N_train), grid = grid_main, hyper = list(xi_sd = 0.1))
}
jobs[[length(jobs) + 1]] <- list(name = "trim", kern = "gamma", imp = 1, u_q = 0.60,
                                 d = design(ps[[1]]$train, ps[[1]]$new, training_data, grid_main, covariates),
                                 keep = ps[[1]]$train >= 0.1 & ps[[1]]$train <= 0.9, grid = grid_main)
cat("Trimming keeps", sum(jobs[[length(jobs)]]$keep), "of", N_train, "subjects\n")
cat("Running", length(jobs), "stage-2 fits on", n_cores, "cores\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)
job_info <- foreach(job = jobs, .packages = c("here", "nimble", "evd", "coda"),
                    .export = c("fit_mixgpd_app", "summarize_cqte_app"), .combine = rbind) %dopar% {
  dat  <- MI[[job$imp]]
  keep <- job$keep
  fit  <- fit_mixgpd_app(dat$glyhb[keep], dat$obese[keep], job$d$X[keep, , drop = FALSE], job$d$new,
                         job$kern, n_mix, p_QTE_diabetes, epsilon, n_iter, n_burnin, n_chains, n_thin,
                         u_lower_q = job$u_q, seed = seed,
                         hyper_override = if (is.null(job$hyper)) list() else job$hyper)
  tag <- sprintf("%s_%s", job$name, job$kern)
  saveRDS(list(Q_trt = fit$Q_trt, Q_con = fit$Q_con), file.path(draw_dir, paste0("draws_", tag, ".rds")))
  if (job$name == "main") {
    saveRDS(fit$samples_trt, file.path(out_dir, paste0(job$kern, "_Treatment_Samples.rds")))
    saveRDS(fit$samples_con, file.path(out_dir, paste0(job$kern, "_Control_Samples.rds")))
  }
  if (!grepl("^prop", job$name)) {
    s <- summarize_cqte_app(fit$Q_trt, fit$Q_con, job$grid, p_QTE_diabetes, CP)
    write.csv(s$profiles, file.path(table_dir, sprintf("MixGPD_%s_Profiles.csv", tag)), row.names = FALSE)
    write.csv(s$strata,   file.path(table_dir, sprintf("MixGPD_%s_Strata.csv", tag)), row.names = FALSE)
  }
  data.frame(job = job$name, kern = job$kern, n = sum(keep),
             l_u_trt = fit$l_u[["treatment"]], l_u_con = fit$l_u[["control"]],
             exceed_trt = fit$exceed[["treatment"]], exceed_con = fit$exceed[["control"]],
             ncomp_trt = fit$n_components[["treatment"]], ncomp_con = fit$n_components[["control"]])
}
stopCluster(cl)
write.csv(job_info, file.path(table_dir, "MixGPD_job_info.csv"), row.names = FALSE)

# ------------------------- 6. Cut posterior: stack the propagated fits -------------------------
for (kern in c("gamma", "normal")) {
  Q_trt <- Q_con <- vector("list", length(p_QTE_diabetes))
  for (k in seq_len(n_imp)) {
    dr <- readRDS(file.path(draw_dir, sprintf("draws_prop%02d_%s.rds", k, kern)))
    for (j in seq_along(p_QTE_diabetes)) {
      Q_trt[[j]] <- rbind(Q_trt[[j]], dr$Q_trt[[j]])
      Q_con[[j]] <- rbind(Q_con[[j]], dr$Q_con[[j]])
    }
  }
  s <- summarize_cqte_app(Q_trt, Q_con, grid_main, p_QTE_diabetes, CP)
  write.csv(s$profiles, file.path(table_dir, sprintf("MixGPD_prop_%s_Profiles.csv", kern)), row.names = FALSE)
  write.csv(s$strata,   file.path(table_dir, sprintf("MixGPD_prop_%s_Strata.csv", kern)), row.names = FALSE)
}
print("--- MixGpd analysis of the diabetes data complete ---")
