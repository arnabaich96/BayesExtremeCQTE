# Competing.R
# Simulation study, competing methods (IPW, Extreme QR, Skew-t; Functions/other_method.R).
# Writes the summary over replicates to Results/Simulation/Tables/<method>_<scenario>_<N>.csv.
# ------------------------- Packages -------------------------
suppressPackageStartupMessages({
  library(quantreg)
  library(future)
  library(future.apply)
  library(here)
  library(progressr)
  library(sn)
})
source(here("Functions/other_method.R"))
handlers(global = TRUE)
handlers("txtprogressbar")

# ------------------------- Simulation -------------------------
run_method_loop <- function(method, Y, D, X, new_X, p_QTE, qte_true, N , alpha, n_cores = 1, tag = N) {
  n_rep <- ncol(Y)
  method_results <- list()
  colnames(new_X) <- paste0("X", seq_len(ncol(new_X)))
  rownames(new_X) <- as.character(rownames(new_X))

  old_warn <- getOption("warn")
  options(warn = -1)

  # SkewT fits do not depend on tau: fit once per replicate and predict all tau levels
  if (method == "SkewT") {
    if (n_cores > 1) { plan(multisession, workers = n_cores) } else { plan(sequential) }
    res_all <- future_lapply(1:n_rep, function(i) {
      qte_skewt_multi(Y[, i], D[, i], X, new_X, taus = p_QTE, true_qte = as.matrix(qte_true), alpha = alpha)
    }, future.seed = TRUE)
  }

  for (p_index in seq_along(p_QTE)) {
    q <- p_QTE[p_index]
    true <- as.numeric(qte_true[, p_index])

    if (n_cores > 1) { plan(multisession, workers = n_cores) } else { plan(sequential) }

    if (method == "SkewT") res <- lapply(res_all, `[[`, p_index) else with_progress({
      p <- progressor(along = 1:n_rep)
      res <- future_lapply(1:n_rep, function(i) {
        p(sprintf("Replicate %d", i))
        Yi <- Y[, i]
        Di <- D[, i]
        switch(method,
               ExtremeQR = qte_extremeqr(Yi, Di, X, new_X, tau = q, true_qte = true, alpha = alpha),
               IPW = qte_ipw(Yi, Di, X, new_X, tau = q, true_qte = true, alpha = alpha))
      }, future.seed = TRUE)
    })

    # failed replicates (NA) are dropped; n_valid records how many remain
    out_df <- data.frame(
      x_index = rownames(new_X),
      estimate = rowMeans(sapply(res, function(x) x$estimate), na.rm = TRUE),
      bias = rowMeans(sapply(res, function(x) x$bias), na.rm = TRUE),
      sd = apply(sapply(res, function(x) x$estimate), 1, sd, na.rm = TRUE),
      ise = sqrt(rowMeans(sapply(res, function(x) x$bias^2), na.rm = TRUE)),
      coverage = rowMeans(sapply(res, function(x) x$coverage), na.rm = TRUE),
      length = rowMeans(sapply(res, function(x) x$length_ci), na.rm = TRUE),
      n_valid = rowSums(!is.na(sapply(res, function(x) x$estimate)))
    )
    method_results[[p_index]] <- out_df
    gc()
  }
  plan(sequential)

  final_table <- do.call(rbind, Map(function(df, p) {
    df$y_index <- round(p, 4)
    df
  }, method_results, p_QTE))

  final_table <- final_table[, c("y_index", "x_index", "estimate", "bias", "ise", "sd", "coverage", "length", "n_valid")]

  dir.create(here("Results/Simulation/Tables"), showWarnings = FALSE, recursive = TRUE)
  invisible(suppressWarnings(
    write.csv(final_table, file = paste0(here("Results/Simulation/Tables"), "/", method, "_", tag, ".csv"), row.names = FALSE)
  ))
  options(warn = old_warn)
}

# methods: any of "IPW", "ExtremeQR", "SkewT"; n_rep: number of replicates used (default all)
run_simulation <- function(N, methods, n_cores, p_QTE, alpha, scenario = "spliced", n_rep = NULL, seed = 2024) {

  DATA <- readRDS(here(sprintf("Synthetic_Data/Data_%s_%d.rds", scenario, N)))
  if (is.null(n_rep)) n_rep <- ncol(DATA$Y)
  Y <- DATA$Y[, 1:n_rep, drop = FALSE]
  D <- DATA$T[, 1:n_rep, drop = FALSE]
  X <- as.data.frame(DATA$X)
  new_X <- as.data.frame(DATA$new_X)
  qte_true <- DATA$True_CQTE

  colnames(X) <- paste0("X", seq_len(ncol(X)))
  colnames(new_X) <- colnames(X)

  for (method in methods) {
    set.seed(seed)
    run_method_loop(method, Y, D, X, new_X, p_QTE, qte_true, N, alpha, n_cores,
                    tag = paste0(scenario, "_", N))
  }
}
