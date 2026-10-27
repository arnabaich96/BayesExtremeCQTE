# run_application.R
# Application of the article: effect of obesity (BMI >= 30) on HbA1c in the Vanderbilt diabetes data.
# Usage (from the repository root):  Rscript run_application.R
#   1. Application/Prepare_Data.R            download (once), clean and multiply impute the data
#   2. Application/MixGpd_Diabetes.R         PS model and all MixGPD fits (main and sensitivity analyses)
#   3. Application/Overlap_Diagnostics.R     PS overlap figure, balance table and Love plot
#   4. Application/MCMC_Diagnostics_Diabetes.R  traceplots, effective sample sizes and R-hat
#   5. Application/Export_Interval_Tables.R  LaTeX tables of the stratified CQTEs and of the sensitivity analyses
# The environment variable CQTE_CORES sets the number of parallel workers (default min(24, cores - 1)).
library(here)

run_script <- function(file) {
  cat(sprintf("\n[%s] Rscript %s\n", format(Sys.time(), "%H:%M:%S"), file))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(here(file)))
  if (status != 0) stop(file, " failed with exit status ", status)
}

for (f in c("Application/Prepare_Data.R", "Application/MixGpd_Diabetes.R", "Application/Overlap_Diagnostics.R",
            "Application/MCMC_Diagnostics_Diabetes.R", "Application/Export_Interval_Tables.R")) {
  run_script(f)
}
cat("\nApplication complete. Tables: Results/Application/Diabetes/Tables; figures: Results/Application/Diabetes/Figures\n")
