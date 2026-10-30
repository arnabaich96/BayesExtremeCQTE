# Application of the article: effect of obesity (BMI >= 30) on HbA1c in the Vanderbilt diabetes data.
# Usage: Rscript run_application.R (runs the five scripts below in order).
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
