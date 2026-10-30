# MCMC fit of the MixGPD model for one arm and posterior conditional quantiles on the profiles full_X.
# data: list(y, X, full_X, l_u), predictors PS first then X (no intercept); hyper: from mixgpd_hyper().
Sim <- function(data, code, epsilon, n_mix, p_QTE, mcmc.plot = FALSE, kern =c("gamma","normal","lognormal"),
                n_iter = 7500, n_burnin = 2500, n_chains = 1, n_thin = 1,
                silent = FALSE, parallel = FALSE, hyper = list())
{
  
  if (silent) {
    sink(tempfile())
    on.exit({ sink() }, add = TRUE)
  }
  
  suppressMessages({
    suppressWarnings({
      if (kern=="normal")
        source(here("Functions/Mix_Gpd_Normal.R"), echo = FALSE, verbose = FALSE)
      else if (kern=="lognormal")
        source(here("Functions/Mix_Gpd_LN.R"), echo = FALSE, verbose = FALSE)
      else
        source(here("Functions/Mix_Gpd_Gamma.R"), echo = FALSE, verbose = FALSE)
      source(here("Functions/Prior_Distributions.R"), echo = FALSE, verbose = FALSE)
    })
  })
  if (parallel) require(future.apply)

  my.list <- list()
  full_X <- as.matrix(data$full_X)   # rows = grid points; cols = predictors (PS + X)
  X_mat  <- as.matrix(data$X)        # rows = obs; cols = predictors (PS + X)
  l_u    <- data$l_u

  # --- Model constants include p, the threshold bound and the hyperparameters ---
  consts <- c(list(n = length(data$y), m = n_mix, p = ncol(X_mat), l_u = l_u), hyper)

  # --- Initial values: latent locations on the scale of the data, u above l_u ---
  y_med <- median(data$y)
  make_inits <- function() {
    ini <- list(
      alpha     = 1,
      eta       = runif(consts$m, 2, 4),
      theta1    = matrix(0, nrow = consts$p, ncol = consts$m),
      mu        = matrix(y_med * runif(consts$n * consts$m, 0.5, 1.5), nrow = consts$n),
      u         = rep(max(quantile(data$y, 0.90, na.rm = TRUE), 1.05 * l_u), consts$n),
      sigma     = max(sd(data$y), 0.1),
      xi        = runif(1, 0.05, 0.3),
      delta1    = rep(0, consts$p),
      v         = rbeta(consts$m - 1, 1, 1),
      var1      = 1,
      var2      = 1
    )
    if (kern == "lognormal") {
      ini$mu  <- matrix(rnorm(consts$n * consts$m, log(y_med), 0.3), nrow = consts$n)
      ini$eta <- runif(consts$m, 0.3, 0.8)
    }
    ini
  }

  # --- Run MCMC ---
  # Chain-specific inits must give a finite log density (1 - F(u) can underflow for poor starts)
  nimble_out <- suppressMessages(suppressWarnings({
    model  <- nimbleModel(code, constants = consts, data = list(y = data$y, X = X_mat),
                          inits = make_inits(), check = FALSE, calculate = FALSE)
    cmodel <- compileNimble(model)
    chain_inits <- lapply(seq_len(n_chains), function(ch) {
      for (attempt in 1:50) {
        ini <- make_inits()
        for (nm in names(ini)) cmodel[[nm]] <- ini[[nm]]
        if (is.finite(cmodel$calculate())) return(ini)
      }
      stop("No initial values with a finite log density were found.")
    })
    conf <- configureMCMC(model, monitors = c("alpha","weights","eta","theta1","delta1",
                                              "sigma","xi","var1","var2"), print = FALSE)
    # GPD scale and shape are strongly correlated a posteriori: add a joint random-walk update
    conf$addSampler(target = c("sigma", "xi"), type = "RW_block")
    cmcmc <- compileNimble(buildMCMC(conf), project = model)
    runMCMC(cmcmc, niter = n_iter, nburnin = n_burnin, nchains = n_chains, thin = n_thin,
            inits = chain_inits, progressBar = FALSE, samplesAsCodaMCMC = TRUE)
  }))
  
  samples <- as.data.frame(as.matrix(nimble_out))
  my.list$MCMC_samples <- nimble_out
  # normalize column names
  colnames(samples) <- gsub("\\s+", "", colnames(samples))
  
  # --- Column names of the monitored parameters ---
  .param_indexers <- function(m, p) list(
    W  = paste0("weights[", 1:m, "]"),
    E  = paste0("eta[",     1:m, "]"),
    T1 = as.vector(outer(1:p, 1:m, function(k, j) paste0("theta1[", k, ",", j, "]"))),
    D1 = paste0("delta1[",  1:p, "]")
  )
  
  # --- Per-tau conditional quantiles: latent locations and thresholds at their conditional medians ---
  # (exp(lp) for Gamma/IG means, lp for the lognormal meanlog, the l_u-truncated lognormal median for u)
  .compute_Q_for_tau <- function(p_tau, samples, full_X, idx, m, p) {
    Xmat <- as.matrix(full_X); K <- nrow(Xmat); nD <- nrow(samples)
    needed <- c(idx$W, idx$E, idx$T1, idx$D1, "sigma", "xi", "var2")
    miss <- setdiff(needed, colnames(samples))
    if (length(miss)) stop("Missing MCMC columns: ", paste(miss, collapse=", "))

    Wmat <- as.matrix(samples[, idx$W,  drop = FALSE])    # nD × m
    Emat <- as.matrix(samples[, idx$E,  drop = FALSE])    # nD × m
    T1bg <- as.matrix(samples[, idx$T1, drop = FALSE])    # nD × (p*m)
    D1m  <- as.matrix(samples[, idx$D1, drop = FALSE])    # nD × p
    sig  <- samples[, "sigma"]
    xi   <- samples[, "xi"]
    sd_u <- 1 / sqrt(samples[, "var2"])

    Q_p  <- matrix(NA_real_, nrow = nD, ncol = K)

    for (k in 1:K) {
      xk <- Xmat[k, ]                                   # length p

      # precompute mu and u for all draws
      mu_all <- matrix(NA_real_, nrow = nD, ncol = m)
      for (j in 1:m) {
        cols <- ((j - 1) * p + 1):(j * p)
        lp_mu <- as.numeric(T1bg[, cols, drop = FALSE] %*% xk)
        mu_all[, j] <- if (kern == "lognormal") lp_mu else exp(lp_mu)
      }
      lp_u_all <- as.numeric(D1m %*% xk)
      u_all    <- median_lnormTrunc(lp_u_all, sd_u, l_u)
      
      for (d in 1:nD) {
        Q_p[d, k] <- tryCatch(
          as.numeric(qTotal(p_tau,
                            Wmat[d, ],
                            mu_all[d, ],
                            Emat[d, ],
                            u_all[d],
                            sig[d],
                            xi[d])),
          error = function(e) NA_real_
        )
      }
    }
    Q_p
  }
  
  idx <- .param_indexers(consts$m, consts$p)
  if (parallel) {
    plan(multisession, workers = length(p_QTE))
    Q <- future.apply::future_lapply(
      p_QTE,
      function(pt) .compute_Q_for_tau(pt, samples, full_X, idx, consts$m, consts$p),
      future.seed = TRUE
    )
    plan(sequential)
  } else {
    Q <- lapply(
      p_QTE,
      function(pt) .compute_Q_for_tau(pt, samples, full_X, idx, consts$m, consts$p)
    )
  }
  my.list$estimated_quantiles <- Q

  # draw-wise number of components with weight above epsilon (all m components enter the CQTE)
  my.list$n_components <- rowSums(as.matrix(samples[, idx$W, drop = FALSE]) > epsilon)

  # posterior means of the threshold coefficients and tail parameters
  D1 <- idx$D1
  delta1_mean <- colMeans(as.matrix(samples[, D1, drop = FALSE]))
  names(delta1_mean) <- paste0("delta1_", seq_len(consts$p))
  my.list$GPD_parameters <- data.frame(
    t(delta1_mean),
    sigma  = mean(samples[, "sigma"]),
    xi     = mean(samples[, "xi"]),
    l_u    = l_u,
    row.names = NULL,
    check.names = FALSE
  )
  
  return(my.list)
}
