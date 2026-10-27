# my_data_sim.R
# Data-generating mechanisms of the simulation study and the true CQTEs:
#   S1 "spliced": Gamma-mixture (treatment) / lognormal-mixture (control) bulk with a GPD tail above
#                 a covariate-dependent threshold (my_data_sim)
#   S2 "burr", S3 "lnmix": misspecified scenarios without the spliced structure (sim_data_scenario)
# Sourced by Simulation/Generate_Scenarios.R and Simulation/Plot_True_Quantile_Surfaces.R.
suppressPackageStartupMessages({
  library(here)
  library(parallel)
  library(doParallel)
  library(evd)
  library(ggplot2)
})

## ============================================================
## Helpers for Gamma-bulk (treatment) and Lognormal-bulk (control)
## ============================================================

# Gamma-mixture quantile for bulk (single draw)
mix_q_gamma <- function(mu_vec, p, eta, wt) {
  id <- sample.int(length(eta), 1, prob = wt)
  qgamma(p, shape = eta[id], rate = eta[id] / mu_vec[id])
}

# Lognormal-mixture quantile for bulk (single draw)
# mu_vec: component means; eta: sdlog
mix_q_lnorm <- function(mu_vec, p, eta, wt) {
  id <- sample.int(length(eta), 1, prob = wt)
  sdlog   <- eta[id]
  meanlog <- log(mu_vec[id]) - 0.5 * sdlog^2
  qlnorm(p, meanlog = meanlog, sdlog = sdlog)
}

# Gamma-mixture CDF at thresholds u (for treatment arm)
mix_cdf_gamma_at_u <- function(wt, eta, mu, u) {
  K <- length(wt)
  cdf_mat <- sapply(1:K, function(k) pgamma(u, shape = eta[k], rate = eta[k] / mu[, k]))
  as.numeric(cdf_mat %*% wt)
}

# Lognormal-mixture CDF at thresholds u (for control arm)
mix_cdf_lnorm_at_u <- function(wt, eta, mu, u) {
  K <- length(wt)
  cdf_mat <- sapply(1:K, function(k) {
    sdlog   <- eta[k]
    meanlog <- log(mu[, k]) - 0.5 * sdlog^2
    plnorm(u, meanlog = meanlog, sdlog = sdlog)
  })
  as.numeric(cdf_mat %*% wt)
}

# Gamma-bulk + GPD tail quantile (treatment)
qTotal_gamma <- function(p, weights, mu, eta, u, sigma, xi) {
  P <- sum(pgamma(u, shape = eta, rate = eta / mu) * weights)
  
  if (p < P) {
    f <- function(x) p - sum(weights * pgamma(x, shape = eta, rate = eta / mu))
    r <- uniroot(f, c(0, max(u * 10, 1)), maxiter = 1e4)
    x <- r$root
  } else {
    a <- (p - P) / (1 - P)
    x <- qgpd(a, loc = u, scale = sigma, shape = xi)
  }
  x
}

# Lognormal-bulk + GPD tail quantile (control)
qTotal_lnorm <- function(p, weights, mu, eta, u, sigma, xi) {
  # mu: component means; eta: sdlog
  F_mix <- function(x) {
    s <- 0
    for (k in seq_along(weights)) {
      sdlog   <- eta[k]
      meanlog <- log(mu[k]) - 0.5 * sdlog^2
      s <- s + weights[k] * plnorm(x, meanlog = meanlog, sdlog = sdlog)
    }
    s
  }
  
  P <- F_mix(u)
  
  if (p < P) {
    f_root <- function(x) F_mix(x) - p
    upper <- max(u * 10, 1)
    r <- uniroot(f_root, c(1e-8, upper), maxiter = 1e4)
    r$root
  } else {
    a <- (p - P) / (1 - P)
    qgpd(a, loc = u, scale = sigma, shape = xi)
  }
}

## ============================================================
## Data generator (NO intercepts; exp links; PS+X in bulk & tail)
## ============================================================
# treatment_params:
#   wt, eta, theta_ps, theta_X, del_ps, del_x, sigma, xi
#   (bulk is Gamma mixture)
# control_params:
#   same structure, but bulk is Lognormal mixture (eta = sdlog)
#
# Bulk mean:   mu_ij = exp(theta_ps_j * PS_i + X_i^T theta_X_{·j})
# Tail thresh: u_i   = exp(del_ps * PS_i + X_i^T del_x)
#
my_data_sim <- function(r, N, n_cores, beta,
                        treatment_params,
                        control_params,
                        p_QTE, x_index, seed = 10) {
  
  set.seed(seed)
  out <- list()
  
  ## --- covariates & PS ---
  p_x <- 3 # number of continuous covariates (excluding categorical)
  X   <- matrix(round(runif(N * p_x, 1, 5), 3), nrow = N, ncol = p_x)
  X   <- cbind(sample.int(4, size = N, replace = TRUE), X)  # categorical covariate with 4 levels
  colnames(X) <- paste0("X", 1:{p_x+1})
  out$X  <- X
  PS     <- plogis(X %*% beta)
  out$PS <- as.numeric(PS)
  
  ## --- unpack params ---
  unpack <- function(par) {
    list(
      wt       = par$wt,
      eta      = par$eta,
      theta_ps = par$theta_ps,
      theta_X  = par$theta_X,      # p_x+1 x K
      del_ps   = par$del_ps,
      del_x    = par$del_x,        # length p_x+1
      sigma    = par$sigma,
      xi       = par$xi
    )
  }
  tr <- unpack(treatment_params)  # Gamma bulk
  co <- unpack(control_params)    # Lognormal bulk
  K  <- length(tr$wt)
  stopifnot(
    length(tr$wt)==K, length(tr$eta)==K, length(tr$theta_ps)==K, all(dim(tr$theta_X)==c(ncol(X), K)),
    length(co$wt)==K, length(co$eta)==K, length(co$theta_ps)==K, all(dim(co$theta_X)==c(ncol(X), K)),
    length(tr$del_x)==ncol(X), length(co$del_x)==ncol(X)
  )
  
  ## --- bulk means & tail thresholds (NO intercepts) ---
  mu_mat <- function(theta_ps, theta_X, PS, X) {
    log_mu <- outer(as.numeric(PS), theta_ps) + X %*% theta_X  # N x K
    exp(log_mu)
  }
  mu1 <- mu_mat(tr$theta_ps, tr$theta_X, PS, X)  # treatment: Gamma means
  mu0 <- mu_mat(co$theta_ps, co$theta_X, PS, X)  # control: Lognormal means
  
  u_build <- function(del_ps, del_x, PS, X) {
    exp(as.numeric(del_ps * PS + X %*% del_x))   # length N
  }
  u1 <- u_build(tr$del_ps, tr$del_x, PS, X)
  u0 <- u_build(co$del_ps, co$del_x, PS, X)
  
  ## --- mixture CDF at thresholds ---
  p1 <- mix_cdf_gamma_at_u(tr$wt, tr$eta, mu1, u1)     # treatment: Gamma bulk
  p0 <- mix_cdf_lnorm_at_u(co$wt, co$eta, mu0, u0)     # control: Lognormal bulk
  
  ## --- sample potential outcomes in parallel ---
  out$Y <- data.frame(matrix(nrow = N, ncol = r))
  out$T <- data.frame(matrix(nrow = N, ncol = r))
  print_id <- sample.int(r, size = 1)
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)
  
  clusterEvalQ(cl, { library(evd) })
  clusterExport(cl, c("qTotal_gamma","qTotal_lnorm","mu1","mu0","u1","u0","tr","co","N"),
                envir = environment())

  for (rep_id in 1:r) {
    tA <- rbinom(N, 1, PS)
    out$T[, rep_id] <- tA

    # potential outcomes by inverse CDF of the spliced bulk + GPD law, i.e. the same
    # distribution whose quantiles define True_CQTE below
    U0 <- runif(N); U1 <- runif(N)

    # y(0): control, Lognormal bulk
    y0 <- parSapply(cl, seq_len(N), function(i) {
      qTotal_lnorm(U0[i], co$wt, mu0[i, ], co$eta, u0[i], co$sigma, co$xi)
    })

    # y(1): treatment, Gamma bulk
    y1 <- parSapply(cl, seq_len(N), function(i) {
      qTotal_gamma(U1[i], tr$wt, mu1[i, ], tr$eta, u1[i], tr$sigma, tr$xi)
    })
    
    y <- numeric(N)
    y[tA==1] <- y1[tA==1]
    y[tA==0] <- y0[tA==0]
    out$Y[, rep_id] <- y
    
    if (rep_id == print_id && r > 1) {
      D <- data.frame(y0 = y0, y1 = y1, y = y)
      out$Densities_wgpd <- ggplot(D) +
        geom_density(aes(x = y0, color = "Control")) +
        geom_density(aes(x = y1, color = "Treatment")) +
        geom_density(aes(x = y,  color = "Observed")) +
        geom_vline(aes(xintercept = mean(u0), color = "Control"), linetype = "dashed") +
        geom_vline(aes(xintercept = mean(u1), color = "Treatment"), linetype = "dashed") +
        theme_classic() + xlab("Value") + ggtitle("Density with bulk/tail split")
    }
  }
  names(out$Y) <- NULL; names(out$T) <- NULL
  
  ## --- true QTE on PS+X grid + threshold/prob table ---
  if (!is.null(p_QTE)) {
    new_X <- apply(X, 2, function(x) quantile(x, x_index, na.rm = TRUE))  # rows = x_index
    rownames(new_X) <- x_index
    new_PS <- plogis(new_X %*% beta)
    
    grid_mu <- function(theta_ps, theta_X, new_PS, new_X) {
      log_mu <- outer(as.numeric(new_PS), theta_ps) + new_X %*% theta_X
      exp(log_mu)  # G x K
    }
    mu1_grid <- grid_mu(tr$theta_ps, tr$theta_X, new_PS, new_X) # Gamma means
    mu0_grid <- grid_mu(co$theta_ps, co$theta_X, new_PS, new_X) # Lognormal means
    
    u_grid <- function(del_ps, del_x, new_PS, new_X) {
      as.numeric(exp(del_ps * new_PS + new_X %*% del_x))    # length G
    }
    u1_grid <- u_grid(tr$del_ps, tr$del_x, new_PS, new_X)
    u0_grid <- u_grid(co$del_ps, co$del_x, new_PS, new_X)
    
    # TRUE QTE: treatment uses Gamma qTotal, control uses Lognormal qTotal
    true_Treat <- sapply(as.numeric(p_QTE), function(p) {
      sapply(seq_along(new_PS), function(g) {
        tryCatch(
          qTotal_gamma(p,
                       tr$wt,
                       mu1_grid[g, ],
                       tr$eta,
                       u1_grid[g],
                       tr$sigma,
                       tr$xi),
          error = function(e) NA_real_
        )
      })
    })
    true_Control <- sapply(as.numeric(p_QTE), function(p) {
      sapply(seq_along(new_PS), function(g) {
        tryCatch(
          qTotal_lnorm(p,
                       co$wt,
                       mu0_grid[g, ],
                       co$eta,
                       u0_grid[g],
                       co$sigma,
                       co$xi),
          error = function(e) NA_real_
        )
      })
    })
    
    true_QTE <- round(true_Treat - true_Control, 3)
    rownames(true_QTE) <- x_index
    colnames(true_QTE) <- round(p_QTE, 5)
    out$new_X <- new_X
    out$True_CQTE <- true_QTE
    
    # 9 x 4 table: u1, p1, u0, p0 (per tau_x)
    G <- length(new_PS)
    p_below_trt <- mix_cdf_gamma_at_u(tr$wt, tr$eta, mu1_grid, u1_grid)
    p_below_con <- mix_cdf_lnorm_at_u(co$wt, co$eta, mu0_grid, u0_grid)
    
    tail_table <- data.frame(
      u1 = u1_grid,
      p1 = round(p_below_trt, 4),
      u0 = u0_grid,
      p0 = round(p_below_con, 4)
    )
    rownames(tail_table) <- x_index
    out$True_tail <- tail_table
  }
  
  out
}

assign("my_data_sim", my_data_sim, .GlobalEnv)

## ============================================================
## Example parameterization (NO intercepts)
## ============================================================
K   <- 3        # mixture components
p_x <- 4        # number of covariates (including categorical X1)

# Treatment: ALL Gamma components
treatment_params <- list(
  wt       = c(0.25, 0.35, 0.40),
  eta      = c(3, 2, 3),                     # Gamma shapes
  theta_ps = c(0.8, 0.5, 0.4)/2,             # PS slopes (bulk)
  theta_X  = matrix(c(                       # X slopes (p_x x K)
    1.25,  0.8,  1.1,
    -0.25, -0.2,  0.3,
    0.4,   0.2, -0.15,
    0.3,   0.1, -0.2)/2,
    nrow = p_x, ncol = K, byrow = TRUE),
  del_ps   =  2,
  del_x    = c(0.6, -0.26, 0.30, 0.14)/2,
  sigma    = 3.5,
  xi       = 0.30
)

# Control: ALL Lognormal components (eta = sdlog)
control_params <- list(
  wt       = c(0.55, 0.20, 0.25),
  eta      = c(0.4, 0.5, 0.35),              # sdlog (example values)
  theta_ps = c(0.4, 0.3, 0.5)/2,
  theta_X  = matrix(c(
    0.3,  0.4,  0.2,
    -0.4, -0.2,  0.2,
    0.3,  0.15, -0.1,
    -0.1,  0.2,  0.05),
    nrow = p_x, ncol = K, byrow = TRUE),
  del_ps   = 1.5,
  del_x    = c(0.15, -0.10, 0.25, 0.10)/2,
  sigma    = 2.0,
  xi       = 0.25
)

## ============================================================
## Misspecified scenarios (no spliced bulk/GPD structure)
## ============================================================
# S2 "burr":  Y(a) | x ~ Burr XII(c_a, k_a) with scale lambda_a(x) = exp(gam_ps * PS + x' gam_x).
#             Regularly varying tail with index xi_a = 1/(c_a k_a), but no threshold above which the
#             excesses are exactly GPD, and the tail scale changes with x.
# S3 "lnmix": Y(a) | x ~ three-component lognormal mixture with no GPD tail (xi = 0 domain,
#             slow convergence of excesses to the GPD limit).
burr_params <- list(
  treatment = list(c = 3, k = 1.1, gam_ps = 0.3, gam_x = c(0.30, -0.10, 0.15, 0.05)),  # xi = 0.30
  control   = list(c = 4, k = 1.0, gam_ps = 0.2, gam_x = c(0.20, -0.15, 0.10, 0.10))   # xi = 0.25
)
lnmix_params <- list(
  treatment = list(wt = treatment_params$wt, sdlog = c(0.60, 0.45, 0.70),
                   theta_ps = treatment_params$theta_ps, theta_X = treatment_params$theta_X),
  control   = list(wt = control_params$wt, sdlog = c(0.50, 0.60, 0.45),
                   theta_ps = control_params$theta_ps, theta_X = control_params$theta_X)
)

qburr <- function(p, lambda, c, k) lambda * ((1 - p)^(-1 / k) - 1)^(1 / c)

# Quantile of a lognormal mixture with component means mu (mean parameterisation)
q_lnmix <- function(p, wt, mu, sdlog) {
  meanlog <- log(mu) - 0.5 * sdlog^2
  qs <- qlnorm(p, meanlog, sdlog)
  if (diff(range(qs)) < 1e-10) return(qs[1])
  uniroot(function(x) sum(wt * plnorm(x, meanlog, sdlog)) - p, range(qs), tol = 1e-10)$root
}

# Same covariates and propensity score in every scenario (my_data_sim draws them first with set.seed)
sim_covariates <- function(N, beta, seed) {
  set.seed(seed)
  p_x <- 3
  X   <- matrix(round(runif(N * p_x, 1, 5), 3), nrow = N, ncol = p_x)
  X   <- cbind(sample.int(4, size = N, replace = TRUE), X)
  colnames(X) <- paste0("X", 1:(p_x + 1))
  list(X = X, PS = as.numeric(plogis(X %*% beta)))
}

sim_data_scenario <- function(scenario = c("spliced", "burr", "lnmix"), r, N, n_cores, beta,
                              p_QTE, x_index, seed = 321) {
  scenario <- match.arg(scenario)
  if (scenario == "spliced") {
    out <- my_data_sim(r, N, n_cores, beta, treatment_params, control_params, p_QTE, x_index, seed)
    out$scenario <- scenario
    return(out)
  }

  cv <- sim_covariates(N, beta, seed)
  X <- cv$X; PS <- cv$PS
  new_X  <- apply(X, 2, function(x) quantile(x, x_index, na.rm = TRUE))
  rownames(new_X) <- x_index
  new_PS <- as.numeric(plogis(new_X %*% beta))

  if (scenario == "burr") {
    lam <- function(par, ps, x) as.numeric(exp(par$gam_ps * ps + x %*% par$gam_x))
    tr <- burr_params$treatment; co <- burr_params$control
    lam1 <- lam(tr, PS, X); lam0 <- lam(co, PS, X)
    draw <- function() {
      list(y1 = qburr(runif(N), lam1, tr$c, tr$k), y0 = qburr(runif(N), lam0, co$c, co$k))
    }
    q_trt <- function(p, g) qburr(p, lam(tr, new_PS[g], new_X[g, , drop = FALSE]), tr$c, tr$k)
    q_con <- function(p, g) qburr(p, lam(co, new_PS[g], new_X[g, , drop = FALSE]), co$c, co$k)
  } else {
    mus <- function(par, ps, x) exp(outer(as.numeric(ps), par$theta_ps) + x %*% par$theta_X)  # N x K
    tr <- lnmix_params$treatment; co <- lnmix_params$control
    mu1 <- mus(tr, PS, X); mu0 <- mus(co, PS, X)
    rmix <- function(mu, par) {
      id <- sample.int(length(par$wt), N, replace = TRUE, prob = par$wt)
      m  <- mu[cbind(seq_len(N), id)]
      rlnorm(N, log(m) - 0.5 * par$sdlog[id]^2, par$sdlog[id])
    }
    draw <- function() list(y1 = rmix(mu1, tr), y0 = rmix(mu0, co))
    q_trt <- function(p, g) q_lnmix(p, tr$wt, mus(tr, new_PS[g], new_X[g, , drop = FALSE]), tr$sdlog)
    q_con <- function(p, g) q_lnmix(p, co$wt, mus(co, new_PS[g], new_X[g, , drop = FALSE]), co$sdlog)
  }

  out <- list(X = X, PS = PS, scenario = scenario)
  out$Y <- data.frame(matrix(nrow = N, ncol = r))
  out$T <- data.frame(matrix(nrow = N, ncol = r))
  for (rep_id in 1:r) {
    tA <- rbinom(N, 1, PS)
    yy <- draw()
    out$T[, rep_id] <- tA
    out$Y[, rep_id] <- ifelse(tA == 1, yy$y1, yy$y0)
  }
  names(out$Y) <- NULL; names(out$T) <- NULL

  true_Treat   <- sapply(as.numeric(p_QTE), function(p) sapply(seq_along(new_PS), function(g) q_trt(p, g)))
  true_Control <- sapply(as.numeric(p_QTE), function(p) sapply(seq_along(new_PS), function(g) q_con(p, g)))
  true_QTE <- round(true_Treat - true_Control, 3)
  rownames(true_QTE) <- x_index
  colnames(true_QTE) <- round(p_QTE, 5)
  out$new_X <- new_X
  out$True_CQTE <- true_QTE
  out
}

