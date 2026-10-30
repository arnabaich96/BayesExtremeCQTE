# Two-stage MixGPD fit of every replicate of a simulated data set. u_lower_q: quantile level of the bound l_u;
# n_ps_draws: 0 = posterior-mean PS, K > 0 = one stage-2 fit per posterior PS draw, with the CQTE draws pooled.
Report <- function(DATA, n_mix, p_QTE, epsilon, CP = NULL, mcmc.plot = FALSE, kern =c("gamma","normal","lognormal"),
                   n_iter, n_burnin, n_chains, n_thin, n_cores,
                   block_size = 10, seed, tag = "spliced_500", u_lower_q = 0.60, n_ps_draws = 0) {
  print('Importing Functions...')
  source(here("Functions/Sim.R"))
  source(here("Functions/Model_Code.R"))
  t1 <- Sys.time()
  set.seed(seed)
  # ------------------------- Define NIMBLE Models -------------------------
  code_trt  <- mixgpd_code(kern)
  code_con  <- mixgpd_code(kern)
  hyper_trt <- mixgpd_hyper("sim", kern, "treatment")
  hyper_con <- mixgpd_hyper("sim", kern, "control")

  # ------------------------- Parallel Execution -------------------------
  library(foreach)
  library(doParallel)

  n_rep <- ncol(DATA$Y)

  DATA$X <- as.data.frame(DATA$X)
  colnames(DATA$X) <- paste0("X", seq_len(ncol(DATA$X)))
  new_X <- as.data.frame(DATA$new_X)
  colnames(new_X) <- colnames(DATA$X)
  qte_true <- DATA$True_CQTE

  out_dir <- here("Simulation/output")
  dir.create(file.path(out_dir, "Diagnostics"), showWarnings = FALSE, recursive = TRUE)

  # PS model compiled once (Bayesian logistic regression, N(0, 2.5^2) prior on the slopes)
  D0 <- data.frame(y = rbinom(nrow(DATA$X), 1, 0.5), DATA$X)
  dummy_model <- brm(
    y ~ .,
    data = D0,
    family = bernoulli(link = "logit"),
    prior = set_prior("normal(0, 2.5)", class = "b"),
    refresh = 0
  )

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)

  # replicates already on disk are skipped; a failed batch is retried on a fresh cluster up to three times
  rep_file <- function(i) sprintf(file.path(out_dir, "qte_%s_%s_rep%03d.rds"), tag, kern, i)
  todo <- which(!file.exists(vapply(1:n_rep, rep_file, "")))
  for (attempt in 1:3) {
  if (!length(todo)) break
  batch_ids <- split(todo, ceiling(seq_along(todo) / block_size))

  for (b in seq_along(batch_ids)) {
    cat(sprintf("Running batch %d of %d (attempt %d)...\n", b, length(batch_ids), attempt))
    res <- tryCatch(foreach(i = batch_ids[[b]],
                   .packages = c("brms", "dplyr","evd","here","coda"),
                   .export = c("Sim","code_trt","code_con","hyper_trt","hyper_con")) %dopar% {
                     set.seed(seed + i)   # replicate-specific seed, independent of the worker
                     # --- Fit PS on this replication ---
                     D <- data.frame(y = as.numeric(DATA$T[, i]), DATA$X)
                     model <- suppressMessages(suppressWarnings(
                       update(dummy_model, newdata = D, recompile = FALSE, refresh = 0)
                     ))
                     if (n_ps_draws > 0) {
                       draw_ids <- sample(ndraws(model), n_ps_draws)
                       ps_train <- posterior_epred(model, draw_ids = draw_ids)                   # K x n
                       ps_new   <- posterior_epred(model, newdata = new_X, draw_ids = draw_ids)  # K x G
                     } else {
                       ps_train <- matrix(fitted(model)[, "Estimate"], nrow = 1)
                       ps_new   <- matrix(fitted(model, newdata = new_X)[, "Estimate"], nrow = 1)
                     }

                     Data <- data.frame(y = DATA$Y[, i], T = DATA$T[, i])
                     trt  <- Data$T == 1
                     l_u_trt <- unname(quantile(Data$y[trt],  u_lower_q))
                     l_u_con <- unname(quantile(Data$y[!trt], u_lower_q))

                     Q_trt <- Q_con <- vector("list", length(p_QTE))
                     diag_rows <- list()
                     for (k in seq_len(nrow(ps_train))) {
                       # --- Build predictor matrices: PS first, then original X ---
                       X_design     <- cbind(PS = ps_train[k, ], as.matrix(DATA$X))
                       new_X_design <- cbind(PS = ps_new[k, ],   as.matrix(new_X))

                       treatment_data <- list(y = Data$y[trt],  X = X_design[trt, , drop = FALSE],
                                              full_X = new_X_design, l_u = l_u_trt)
                       control_data   <- list(y = Data$y[!trt], X = X_design[!trt, , drop = FALSE],
                                              full_X = new_X_design, l_u = l_u_con)

                       out_treat <- Sim(treatment_data, code_trt, epsilon, n_mix, p_QTE,
                                        mcmc.plot, kern, n_iter, n_burnin, n_chains, n_thin, TRUE,
                                        hyper = hyper_trt)
                       out_con   <- Sim(control_data, code_con, epsilon, n_mix, p_QTE,
                                        mcmc.plot, kern, n_iter, n_burnin, n_chains, n_thin, TRUE,
                                        hyper = hyper_con)

                       # pool the stage-2 draws over the PS draws
                       for (j in seq_along(p_QTE)) {
                         Q_trt[[j]] <- rbind(Q_trt[[j]], out_treat$estimated_quantiles[[j]])
                         Q_con[[j]] <- rbind(Q_con[[j]], out_con$estimated_quantiles[[j]])
                       }

                       # compact convergence / tail diagnostics (full samples are not stored)
                       for (arm in c("treatment", "control")) {
                         out <- if (arm == "treatment") out_treat else out_con
                         mc  <- out$MCMC_samples[, c("sigma", "xi"), drop = FALSE]
                         ess <- coda::effectiveSize(mc)
                         rh  <- tryCatch(coda::gelman.diag(mc, autoburnin = FALSE)$psrf[, 1],
                                         error = function(e) c(sigma = NA, xi = NA))
                         diag_rows[[length(diag_rows) + 1]] <- data.frame(
                           rep = i, ps_draw = k, arm = arm,
                           sigma = out$GPD_parameters$sigma, xi = out$GPD_parameters$xi,
                           l_u = out$GPD_parameters$l_u,
                           ess_sigma = ess[["sigma"]], ess_xi = ess[["xi"]],
                           rhat_sigma = rh[1], rhat_xi = rh[2],
                           n_components = mean(out$n_components),
                           row.names = NULL)
                       }
                     }

                     # --- CQTE per tau on the evaluation profiles: posterior mean and equal-tailed CP interval ---
                     qte_list <- lapply(seq_along(p_QTE), function(j) {
                       qte <- Q_trt[[j]] - Q_con[[j]]
                       estimate <- apply(qte, 2, mean, na.rm = TRUE)
                       post_median <- apply(qte, 2, median, na.rm = TRUE)
                       lower <- apply(qte, 2, quantile, probs = 0.5 * (1 - CP), na.rm = TRUE)
                       upper <- apply(qte, 2, quantile, probs = 1 - 0.5 * (1 - CP), na.rm = TRUE)
                       length.ci <- upper - lower
                       data.frame(
                         estimate = estimate,
                         lower = lower,
                         upper = upper,
                         length.ci = length.ci,
                         post_median = post_median
                       )
                     })
                     list(qte = qte_list, diag = do.call(rbind, diag_rows))
                   }, error = function(e) {
                     message("Batch failed: ", conditionMessage(e))
                     NULL
                   })

    if (is.null(res)) {
      # a dead worker leaves the cluster unusable: replace it
      try(stopCluster(cl), silent = TRUE)
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      next
    }
    for (j in seq_along(batch_ids[[b]])) {
      saveRDS(res[[j]], file = rep_file(batch_ids[[b]][j]))
      gc()
    }
  }
  todo <- which(!file.exists(vapply(1:n_rep, rep_file, "")))
  }
  if (length(todo)) stop("Replicates still missing after 3 attempts: ", paste(todo, collapse = ", "))

  # ---- Post-processing ----
  reps <- lapply(1:n_rep, function(i) readRDS(rep_file(i)))
  qte <- lapply(reps, `[[`, "qte")
  diag_df <- do.call(rbind, lapply(reps, `[[`, "diag"))
  write.csv(diag_df, file.path(out_dir, "Diagnostics", sprintf("diag_%s_%s.csv", tag, kern)), row.names = FALSE)

  x_vals <- as.numeric(rownames(qte_true))
  n_x <- length(x_vals)
  n_q <- length(p_QTE)

  coverage_array <- array(0, dim = c(n_rep, n_q, n_x))
  length_array   <- array(0, dim = c(n_rep, n_q, n_x))
  estimate_array <- array(0, dim = c(n_rep, n_q, n_x))
  bias_array     <- array(0, dim = c(n_rep, n_q, n_x))

  for (r in seq_len(n_rep)) {
    for (q in seq_along(p_QTE)) {
      for (x in seq_along(x_vals)) {
        result <- qte[[r]][[q]][x, ]
        true_val <- qte_true[x, q]
        coverage_array[r, q, x] <- as.integer(result["lower"] <= true_val & true_val <= result["upper"])
        length_array[r, q, x]   <- as.numeric(result["length.ci"])
        estimate_array[r, q, x] <- as.numeric(result["estimate"])
        bias_array[r, q, x] <- as.numeric(result["estimate"]) - true_val
      }
    }
  }

  coverage_prob <- apply(coverage_array, c(2,3), mean)
  avg_length_ci <- apply(length_array, c(2,3), mean)
  avg_estimate  <- apply(estimate_array, c(2,3), mean)
  sd_estimate   <- apply(estimate_array, c(2,3), sd)
  bias <- avg_estimate - t(qte_true)                     # estimate - truth
  ise <- apply(bias_array, c(2,3), function(x) sqrt(mean(x^2, na.rm = TRUE)))   # RMSE

  result_df <- expand.grid(
    y_index = round(p_QTE, 5),
    x_index = x_vals
  ) %>%
    mutate(
      coverage = as.vector(coverage_prob),
      length = as.vector(avg_length_ci),
      estimate = as.vector(avg_estimate),
      bias = as.vector(bias),
      sd = as.vector(sd_estimate),
      ise = as.vector(ise),
    ) %>%
    round(4)

  list(QTE = result_df, Diagnostics = diag_df, Total_execution_time = Sys.time() - t1)
}
