# Simulated data sets of scenarios S1-S3 (n = 500, 100 replicates, seed 321) and their true CQTEs;
# output: Synthetic_Data/Data_<scenario>_500.rds
library(here)
source(here("Simulation/my_data_sim.R"))
dir.create(here("Synthetic_Data"), showWarnings = FALSE)
n_cores <- as.integer(Sys.getenv("CQTE_CORES", min(10, max(1, parallel::detectCores() - 1))))

N        <- 500
r        <- 100
beta     <- c(-0.03, -0.4, 0.12, 0.55)
x_index  <- seq(0.1, 0.9, length.out = 9)
p_QTE    <- 1 - c(0.75, 0.5, 0.25, 0.05, 10/N, 5/N, 1/(N^1.25), 1/(N * log(N)), 1/(N^1.5))

for (scenario in c("spliced", "burr", "lnmix")) {
  message("Generating scenario ", scenario, " (N = ", N, ")")
  Data <- sim_data_scenario(scenario, r = r, N = N, n_cores = n_cores, beta = beta,
                            p_QTE = p_QTE, x_index = x_index, seed = 321)
  saveRDS(Data, file = here(sprintf("Synthetic_Data/Data_%s_%d.rds", scenario, N)))
  cat("\n==== ", scenario, ": TRUE CQTE (rows: covariate profiles; cols: tau) ====\n")
  print(Data$True_CQTE)
  cat("treated fraction:", round(mean(as.matrix(Data$T)), 3),
      " | outcome quantiles (50/95/99%):", round(quantile(as.matrix(Data$Y), c(.5, .95, .99)), 2), "\n")
}

# check that S1 has no point mass at the threshold u(x)
S1 <- readRDS(here(sprintf("Synthetic_Data/Data_spliced_%d.rds", N)))
u_build <- function(par, PS, X) exp(as.numeric(par$del_ps * PS + X %*% par$del_x))
u1 <- u_build(treatment_params, S1$PS, S1$X); u0 <- u_build(control_params, S1$PS, S1$X)
Ym <- as.matrix(S1$Y); Tm <- as.matrix(S1$T)
u_obs <- ifelse(Tm == 1, u1, u0)
cat("\nS1 share of y within 1e-4 of u(x):", mean(abs(Ym - u_obs) < 1e-4),
    " | share above u(x):", round(mean(Ym > u_obs), 3), "\n")
print(S1$True_tail)
