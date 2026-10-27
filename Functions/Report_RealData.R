# Functions/Report_RealData.R
# Two-stage MixGPD analysis of a real data set:
#   fit_ps_app()         stage 1, Bayesian logistic propensity-score model (brms)
#   fit_mixgpd_app()     stage 2, MixGPD outcome model per arm given one PS vector
#   summarize_cqte_app() posterior summaries of the CQTE draws per profile and per stratum
# Stacking the stage-2 draws over several PS draws (and imputed data sets) gives the cut posterior
# that propagates stage-1 uncertainty; see Application/MixGpd_Diabetes.R.

# Stage 1. Returns posterior-mean PS for the data and the grid, plus n_draws posterior PS draws.
fit_ps_app <- function(T, X, newdata, n_draws = 0, seed = 2024, ps_model = NULL) {
  D_ps <- data.frame(T = as.numeric(T), X)
  if (is.null(ps_model)) {
    ps_model <- brm(T ~ ., data = D_ps, family = bernoulli(link = "logit"),
                    prior = set_prior("normal(0, 2.5)", class = "b"),
                    refresh = 0, seed = seed)
  } else {
    ps_model <- update(ps_model, newdata = D_ps, recompile = FALSE, refresh = 0, seed = seed)
  }
  newdata_X <- newdata[, names(X), drop = FALSE]
  out <- list(model     = ps_model,
              train     = fitted(ps_model)[, "Estimate"],
              new       = fitted(ps_model, newdata = newdata_X)[, "Estimate"])
  if (n_draws > 0) {
    set.seed(seed)
    draw_ids <- sample(ndraws(ps_model), n_draws)
    out$train_draws <- posterior_epred(ps_model, draw_ids = draw_ids)                    # K x n
    out$new_draws   <- posterior_epred(ps_model, newdata = newdata_X, draw_ids = draw_ids) # K x G
  }
  out
}

# Stage 2. y, T: outcome and treatment; X_design, new_X_design: PS first, then the covariates.
# hyper_override: hyperparameters replacing those of mixgpd_hyper("app", ...), e.g. list(xi_sd = 0.1)
fit_mixgpd_app <- function(y, T, X_design, new_X_design, kern = c("gamma","normal","lognormal"),
                           n_mix, p_QTE, epsilon, n_iter, n_burnin, n_chains, n_thin,
                           u_lower_q = 0.60, seed = 2024, hyper_override = list()) {
  kern <- match.arg(kern)
  set.seed(seed)
  source(here("Functions/Sim.R"))
  source(here("Functions/Model_Code.R"))

  trt <- T == 1
  X_design <- as.matrix(X_design); new_X_design <- as.matrix(new_X_design)
  l_u <- c(treatment = unname(quantile(y[trt], u_lower_q)), control = unname(quantile(y[!trt], u_lower_q)))

  out_treat <- Sim(list(y = y[trt], X = X_design[trt, , drop = FALSE], full_X = new_X_design, l_u = l_u[["treatment"]]),
                   mixgpd_code(kern), epsilon, n_mix, p_QTE, mcmc.plot = FALSE, kern,
                   n_iter, n_burnin, n_chains, n_thin, TRUE,
                   hyper = modifyList(mixgpd_hyper("app", kern, "treatment"), hyper_override))
  out_con   <- Sim(list(y = y[!trt], X = X_design[!trt, , drop = FALSE], full_X = new_X_design, l_u = l_u[["control"]]),
                   mixgpd_code(kern), epsilon, n_mix, p_QTE, mcmc.plot = FALSE, kern,
                   n_iter, n_burnin, n_chains, n_thin, TRUE,
                   hyper = modifyList(mixgpd_hyper("app", kern, "control"), hyper_override))

  list(Q_trt = out_treat$estimated_quantiles, Q_con = out_con$estimated_quantiles,
       samples_trt = out_treat$MCMC_samples, samples_con = out_con$MCMC_samples,
       l_u = l_u,
       n_components = c(treatment = mean(out_treat$n_components), control = mean(out_con$n_components)),
       exceed = c(treatment = mean(y[trt] > l_u[["treatment"]]), control = mean(y[!trt] > l_u[["control"]])))
}

# Posterior summaries of CQTE = Q_trt - Q_con. Q_trt, Q_con: lists over tau of (draws x profiles).
# Stratum summaries average the CQTE over the profiles sharing a level, draw by draw.
# The point estimate is the posterior mean; the posterior median is reported as well.
summarize_cqte_app <- function(Q_trt, Q_con, newdata, p_QTE, CP = 0.90, strata = c("age", "stab.glu", "gender")) {
  lo <- 0.5 * (1 - CP); hi <- 1 - lo
  prof <- strat <- list()
  for (j in seq_along(p_QTE)) {
    cqte_draws <- Q_trt[[j]] - Q_con[[j]]
    prof[[j]] <- cbind(newdata, data.frame(
      estimate       = apply(cqte_draws, 2, mean, na.rm = TRUE),
      post_median    = apply(cqte_draws, 2, median, na.rm = TRUE),
      sd             = apply(cqte_draws, 2, sd, na.rm = TRUE),
      lower_cred_int = apply(cqte_draws, 2, quantile, probs = lo, na.rm = TRUE),
      upper_cred_int = apply(cqte_draws, 2, quantile, probs = hi, na.rm = TRUE),
      prob_positive  = apply(cqte_draws, 2, function(d) mean(d > 0, na.rm = TRUE)),
      quantile       = p_QTE[j]))
    for (v in intersect(strata, names(newdata))) {
      for (lev in sort(unique(newdata[[v]]))) {
        avg <- rowMeans(cqte_draws[, newdata[[v]] == lev, drop = FALSE], na.rm = TRUE)
        strat[[length(strat) + 1]] <- data.frame(
          variable = v, level = lev, quantile = p_QTE[j],
          estimate = mean(avg, na.rm = TRUE),
          post_median = median(avg, na.rm = TRUE),
          lower_cred_int = unname(quantile(avg, lo, na.rm = TRUE)),
          upper_cred_int = unname(quantile(avg, hi, na.rm = TRUE)),
          prob_positive = mean(avg > 0, na.rm = TRUE))
      }
    }
  }
  list(profiles = do.call(rbind, prof), strata = do.call(rbind, strat))
}
