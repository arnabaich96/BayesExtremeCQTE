# MCMC and mixture diagnostics of the MixGPD simulation fits (Supplementary Table S4), from the
# diagnostics files written by Functions/Report.R.
library(here)
library(dplyr)

diag_dir <- here("Simulation/output/Diagnostics")
tab_dir  <- here("Results/Simulation/Tables")
kern_lab <- c(gamma = "MixGammaGPD", normal = "MixInvNormGPD", lognormal = "MixLogNormGPD")
scen_lab <- c(spliced = "S1", burr = "S2", lnmix = "S3")

files <- list.files(diag_dir, pattern = "^diag_(spliced|burr|lnmix)_500_(gamma|normal|lognormal)\\.csv$", full.names = TRUE)
if (!length(files)) stop("No diagnostics files found in ", diag_dir)

d <- bind_rows(lapply(files, function(f) {
  m <- regmatches(basename(f), regexec("^diag_(spliced|burr|lnmix)_500_(gamma|normal|lognormal)\\.csv$", basename(f)))[[1]]
  read.csv(f) %>% mutate(scenario = m[2], kern = m[3])
}))

summ <- d %>%
  group_by(scenario, kern) %>%
  summarise(
    fits            = n(),
    rhat_sigma_med  = median(rhat_sigma, na.rm = TRUE),
    rhat_xi_med     = median(rhat_xi, na.rm = TRUE),
    share_rhat_1.1  = mean(pmax(rhat_sigma, rhat_xi) < 1.1, na.rm = TRUE),
    ess_sigma_med   = median(ess_sigma, na.rm = TRUE),
    ess_xi_med      = median(ess_xi, na.rm = TRUE),
    ncomp_mean      = mean(n_components, na.rm = TRUE),
    ncomp_max       = max(n_components, na.rm = TRUE),
    xi_med_trt      = median(xi[arm == "treatment"], na.rm = TRUE),
    xi_med_con      = median(xi[arm == "control"], na.rm = TRUE),
    .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = names(scen_lab)),
         kern = factor(kern, levels = names(kern_lab))) %>%
  arrange(scenario, kern)

write.csv(summ, file.path(tab_dir, "Diagnostics_Summary.csv"), row.names = FALSE)
lines <- with(summ, sprintf("%s & %s & %.2f & %.2f & %.2f & %.0f & %.0f & %.2f & %.2f \\\\",
                            scen_lab[as.character(scenario)], kern_lab[as.character(kern)],
                            rhat_sigma_med, rhat_xi_med, share_rhat_1.1, ess_sigma_med, ess_xi_med,
                            ncomp_mean, ncomp_max))
writeLines(lines, file.path(tab_dir, "table_sim_diagnostics.tex"))
print(as.data.frame(summ), digits = 3)
