# run_simulation.R
# Simulation study of the article: n = 500, 100 replicates per scenario.
# Usage (from the repository root):  Rscript run_simulation.R [step ...]
#   data       simulated data of S1-S3 (Synthetic_Data/) and the true quantile surfaces (Suppl.)
#   S1, S2, S3 MixGPD (Gamma, inverse-Gaussian and lognormal kernels), IPW and Extreme QR on
#              S1 (spliced), S2 (Burr XII) and S3 (lognormal mixture)
#   threshold  MixGammaGPD on S1 with l_u at the 50th, 70th and 80th percentiles (first 50 replicates)
#   ps         MixGammaGPD on S1 with K = 5 posterior PS draws (first 50 replicates)
#   skewt      Skew-t (bootstrap) on S1-S3
#   tables     LaTeX tables, MCMC diagnostics summary and heatmaps
# Without arguments all steps run in this order. MixGPD fits are saved per replicate in
# Simulation/output, so an interrupted step resumes where it stopped.
# The environment variable CQTE_CORES sets the number of parallel workers (default min(24, cores - 1)).
library(here)

steps <- commandArgs(trailingOnly = TRUE)
all_steps <- c("data", "S1", "S2", "S3", "threshold", "ps", "skewt", "tables")
if (!length(steps)) steps <- all_steps
if (length(bad <- setdiff(steps, all_steps))) stop("Unknown step(s): ", paste(bad, collapse = ", "))

n_cores <- as.integer(Sys.getenv("CQTE_CORES", min(24, max(1, parallel::detectCores() - 1))))
N <- 500
p_QTE <- 1 - c(0.75, 0.5, 0.25, 0.05, 10/N, 5/N, 1/(N^1.25), 1/(N * log(N)), 1/(N^1.5))
dir.create(here("Results/Simulation/Tables"), showWarnings = FALSE, recursive = TRUE)

# Each script runs in its own R process
run_script <- function(file) {
  cat(sprintf("\n[%s] Rscript %s\n", format(Sys.time(), "%H:%M:%S"), file))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(here(file)))
  if (status != 0) stop(file, " failed with exit status ", status)
}

timing_file <- here("Results/Simulation/Methods_timings.csv")
timed <- function(label, expr) {
  cat(sprintf("\n[%s] Running: %s\n", format(Sys.time(), "%H:%M:%S"), label))
  t0 <- Sys.time()
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) { message("Error in ", label, ": ", e$message); FALSE })
  row <- data.frame(Method = label, StartTime = format(t0, "%Y-%m-%d %H:%M:%S"),
                    Duration_minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), OK = ok)
  write.table(row, timing_file, sep = ",", row.names = FALSE, append = file.exists(timing_file),
              col.names = !file.exists(timing_file))
  cat(sprintf("Finished %s in %.1f minutes\n", label, row$Duration_minutes))
}

if ("data" %in% steps) {
  run_script("Simulation/Generate_Scenarios.R")
  run_script("Simulation/Plot_True_Quantile_Surfaces.R")
}

fit_steps <- intersect(steps, c("S1", "S2", "S3", "threshold", "ps", "skewt"))
if (length(fit_steps)) {
  source(here("Simulation/MixGpd.R"))
  MixGPD_run_simulation <- run_simulation
  source(here("Simulation/Competing.R"))
  Competing_run_simulation <- run_simulation

  run_competitors <- function(scenario, methods) {
    for (method in methods) {
      timed(paste0(method, "_", scenario),
            Competing_run_simulation(N, methods = method, n_cores = n_cores, p_QTE = p_QTE, alpha = 0.1,
                                     scenario = scenario))
    }
  }
  run_scenario <- function(scenario) {
    for (kern in c("gamma", "normal", "lognormal")) {
      timed(paste0("MixGPD_", kern, "_", scenario),
            MixGPD_run_simulation(N, kern, scenario = scenario, n_cores = n_cores))
    }
    run_competitors(scenario, c("IPW", "ExtremeQR"))
  }

  if ("S1" %in% steps) run_scenario("spliced")
  if ("S2" %in% steps) run_scenario("burr")
  if ("S3" %in% steps) run_scenario("lnmix")
  if ("threshold" %in% steps) {
    for (u_q in c(0.50, 0.70, 0.80)) {
      sfx <- sprintf("_lu%02d", 100 * u_q)
      timed(paste0("MixGPD_gamma_spliced", sfx),
            MixGPD_run_simulation(N, "gamma", scenario = "spliced", u_lower_q = u_q, n_rep = 50,
                                  suffix = sfx, n_cores = n_cores))
    }
  }
  if ("ps" %in% steps) {
    timed("MixGPD_gamma_spliced_ps5",
          MixGPD_run_simulation(N, "gamma", scenario = "spliced", n_ps_draws = 5, n_rep = 50,
                                suffix = "_ps5", n_cores = n_cores))
  }
  if ("skewt" %in% steps) {
    for (scenario in c("spliced", "burr", "lnmix")) run_competitors(scenario, "SkewT")
  }
}

if ("tables" %in% steps) {
  run_script("Simulation/Make_Tables.R")
  run_script("Simulation/Summarize_Diagnostics.R")
  run_script("Simulation/Analysis.R")
}
cat("\nAll requested simulation steps are complete.\n")
