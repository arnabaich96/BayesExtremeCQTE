# Installs the R packages of the pipeline; the article used R 4.5.3 and the versions below.
# Usage: Rscript install_packages.R [--exact]; --exact installs the listed versions.
versions <- c(
  "brms"         = "2.22.0",
  "rstan"        = "2.32.7",
  "nimble"       = "1.3.0",
  "coda"         = "0.19-4.1",
  "doParallel"   = "1.0.17",
  "dplyr"        = "1.1.4",
  "evd"          = "2.3-7.1",
  "foreach"      = "1.5.2",
  "future"       = "1.40.0",
  "future.apply" = "1.11.3",
  "ggplot2"      = "3.5.2",
  "here"         = "1.0.1",
  "mice"         = "3.18.0",
  "progressr"    = "0.15.1",
  "purrr"        = "1.0.4",
  "quantreg"     = "6.1",
  "readr"        = "2.1.5",
  "scales"       = "1.4.0",
  "sn"           = "2.1.3",
  "stringr"      = "1.5.1",
  "tibble"       = "3.2.1",
  "tidyr"        = "1.3.1",
  "tidyverse"    = "2.0.0"
)
repos <- "https://cloud.r-project.org"
exact <- "--exact" %in% commandArgs(trailingOnly = TRUE)

installed <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
if (exact) {
  if (is.na(installed("remotes"))) install.packages("remotes", repos = repos)
  for (p in names(versions)) {
    if (!identical(package_version(installed(p), strict = FALSE), package_version(versions[[p]]))) {
      remotes::install_version(p, version = versions[[p]], repos = repos, upgrade = "never")
    }
  }
} else {
  missing <- names(versions)[is.na(vapply(names(versions), installed, ""))]
  if (length(missing)) install.packages(missing, repos = repos)
}

now <- vapply(names(versions), installed, "")
differs <- is.na(now) | package_version(ifelse(is.na(now), "0", now)) != package_version(versions)
cat("R", as.character(getRversion()), "(article: 4.5.3)\n")
print(data.frame(article = versions, installed = now, differs = differs), right = FALSE)
