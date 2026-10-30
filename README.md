# BayesExtremeCQTE

R code for the article *A semiparametric Bayesian approach for estimating extreme conditional quantile
treatment effect for heavy-tailed data*. It reproduces the simulation study and the application to the
Vanderbilt diabetes data.

For each treatment arm, the conditional outcome distribution is a Dirichlet process mixture (Gamma,
inverse-Gaussian or lognormal kernel) for the central part, joined to a generalized Pareto distribution
above a covariate-dependent latent threshold. The covariates and an estimated propensity score (PS) enter
as predictors. The PS is estimated first, by a Bayesian logistic regression (`brms`), and the outcome
model is then fitted by MCMC with `nimble`. Conditional quantile treatment effects (CQTEs) at all
quantile levels, including extreme levels, come from the same posterior.

## Requirements

- R 4.5.3 was used for the article. Install the packages with `Rscript install_packages.R`.
  Add `--exact` to install the versions used for the article.
- A C++ toolchain, because `nimble` and `rstan` (the back end of `brms`) compile models: Rtools on Windows,
  the Xcode command line tools on macOS.
- Internet access the first time the application runs, to download the public data file.
- Several cores. The article used 24 parallel workers (AMD Ryzen 9 9950X, 128 GB RAM). Set the number
  of workers with the environment variable `CQTE_CORES`. The default is the smaller of 24 and the number
  of cores minus one.

## Running the pipeline

Run all commands from the repository root, so that `here::here()` resolves paths correctly.

```sh
Rscript install_packages.R
Rscript run_application.R          # about 70 minutes with 24 workers
Rscript run_simulation.R           # about 1.5 days with 24 workers
```

`run_simulation.R` also runs individual steps, for example `Rscript run_simulation.R data S1 tables`:

| Step | What it does |
|---|---|
| `data` | Simulates scenarios S1 (spliced), S2 (Burr XII) and S3 (lognormal mixture), with n = 500 and 100 replicates, and plots the true quantile surfaces |
| `S1`, `S2`, `S3` | MixGPD with the three kernels, IPW and Extreme QR |
| `threshold` | MixGammaGPD on S1 with the threshold bound at the 50th, 70th and 80th percentiles (first 50 replicates) |
| `ps` | MixGammaGPD on S1 with the PS uncertainty propagated through K = 5 posterior PS draws (first 50 replicates) |
| `skewt` | Skew-t with bootstrap intervals on S1 to S3 |
| `tables` | LaTeX tables, the MCMC diagnostics summary and the heatmaps |

Each replicate of a MixGPD fit is saved in `Simulation/output/`, so a step that is interrupted resumes
where it stopped. Delete that folder to start a step from scratch.

## Repository structure

```
run_simulation.R, run_application.R   entry points
install_packages.R                    package installation and the versions used for the article
Functions/
  Model_Code.R                        NIMBLE model code and the prior hyperparameters
  Mix_Gpd_Gamma.R, Mix_Gpd_Normal.R,  density and quantile function of the spliced mixture for the
  Mix_Gpd_LN.R                        Gamma, inverse-Gaussian and lognormal kernels
  Prior_Distributions.R               truncated lognormal distribution of the latent threshold
  Sim.R                               MCMC for one arm, and the posterior conditional quantiles
  Report.R                            two-stage fit of all replicates of a simulated data set
  Report_RealData.R                   two-stage fit and posterior summaries for the application
  other_method.R, extremal_qr_inference.R, skewt_method_helpers.R
                                      competing methods: IPW, Extreme QR and Skew-t
Simulation/
  my_data_sim.R                       data-generating mechanisms and true CQTEs
  Generate_Scenarios.R                simulated data sets
  Plot_True_Quantile_Surfaces.R       true quantile surfaces in S1
  MixGpd.R, Competing.R               fits of MixGPD and of the competing methods
  Make_Tables.R, Summarize_Diagnostics.R, Analysis.R
                                      tables, MCMC diagnostics and heatmaps
Application/
  Prepare_Data.R                      download, treatment definition and imputation
  MixGpd_Diabetes.R                   PS model and all MixGPD fits (main and sensitivity analyses)
  Overlap_Diagnostics.R               PS overlap and covariate balance
  MCMC_Diagnostics_Diabetes.R         traceplots, effective sample sizes and R-hat
  Export_Interval_Tables.R            stratified CQTE tables and the sensitivity table
```

The repository tracks only code. Data, fitted models, tables and figures are written to `Synthetic_Data/`,
`Simulation/output/`, `Application/data/`, `Application/output/` and `Results/`, which git ignores.

## Outputs and the article

| Output (under `Results/`) | In the article |
|---|---|
| `Simulation/Tables/table1_spliced_500.tex` | Simulation summary for S1 |
| `Simulation/Tables/table_misspec_500.tex` | Misspecified scenarios S2 and S3 |
| `Simulation/Tables/table_threshold_sim.tex`, `table_ps_propagation_sim.tex` | Threshold bound and PS propagation in the simulation |
| `Simulation/Tables/table_sim_diagnostics.tex` | MCMC diagnostics of the simulation (Supplementary Material) |
| `Simulation/Figures/Competing_*_N500_*` and `Mixgpd_*_N500_*` | Heatmaps of bias, coverage, SD and interval width (Supplementary Material) |
| `Simulation/Figures/Quantile_Surfaces_*` | True quantile surfaces (Supplementary Material) |
| `Application/Diabetes/Figures/fig_ps_overlap_bw.png` | PS overlap |
| `Application/Diabetes/Tables/table_balance.tex` | Covariate balance |
| `Application/Diabetes/Tables/table_app_gender.tex`, `table_app_age.tex`, `table_app_glu.tex` | CQTE by gender, age and stabilized glucose |
| `Application/Diabetes/Tables/table_app_sensitivity.tex` | Sensitivity analyses of the application |
| `Application/Diabetes/Figures/fig_love_plot_bw.png`, `Traceplot_*`, `ESS_*` | Love plot and MCMC diagnostics (Supplementary Material) |

## Application data

`Application/Prepare_Data.R` downloads `diabetes.csv` from the
[Vanderbilt Biostatistics Datasets](https://hbiostat.org/data/) (Willems et al., 1997). The treatment is
obesity, defined as BMI >= 30 kg/m^2, with BMI computed from weight and height. The outcome is HbA1c
(`glyhb`). The confounders are age, gender, body frame, location, stabilized glucose, the
cholesterol/HDL ratio, systolic and diastolic blood pressure and the time since the last meal.

BMI is used only to define the treatment. It is not a confounder and is not kept in the analysis data,
because the treatment is a deterministic function of BMI and adjusting for it would violate the overlap
assumption. Waist and hip circumference, which are alternative measurements of adiposity, are not used
either. The 19 subjects with missing HbA1c, height or weight are excluded (n = 384). The remaining missing
covariate values are imputed by chained equations with predictive mean matching, and the analysis uses the
first completed data set.

## Reproducibility

- The simulated data sets and every application fit have fixed seeds. With R 4.5.3 and the package
  versions in `install_packages.R`, a rerun of the application reproduces its tables and figures in the
  article exactly, and `Rscript run_simulation.R data` reproduces the simulated data sets exactly.
- Each simulation replicate is seeded by its index, and the competing methods use parallel-safe random
  number streams. The simulation runs reported in the article did not fix the seeds of the parallel
  workers, so a rerun reproduces the reported simulation summaries up to Monte Carlo error.
