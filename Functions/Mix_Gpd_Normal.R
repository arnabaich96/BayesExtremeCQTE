# Inverse-Gaussian mixture spliced with a GPD tail: density dTotal, sampler rTotal, quantile qTotal.

library(nimble)
library(evd)  # for qgpd in qTotal (R-level only)

## ---- Density ----

dTotal <- nimbleFunction(
  run = function(x      = double(0),
                 weights = double(1),
                 mu      = double(1),   # IG means > 0
                 eta     = double(1),   # IG shapes λ_j > 0
                 u       = double(0),
                 sigma   = double(0),
                 xi      = double(0),
                 log     = integer(0, default = 0)) {
    
    returnType(double(0))
    
    K <- length(weights)
    prob <- 0.0
    
    if (x < u) {
      ## ---------------- BULK: Inverse-Gaussian mixture ----------------
      if (x <= 0.0) {
        prob <- 0.0
      } else {
        for (j in 1:K) {
          mu_j  <- mu[j]
          lam_j <- eta[j]
          
          if (mu_j > 0.0) {
            if (lam_j > 0.0) {
              ## inverse-Gaussian density
              pi_val <- 3.141592653589793
              cst  <- sqrt(lam_j / (2.0 * pi_val * x^3))
              expo <- -lam_j * (x - mu_j) * (x - mu_j) /
                (2.0 * mu_j * mu_j * x)
              prob <- prob + weights[j] * cst * exp(expo)
            }
          }
        }
      }
      
    } else {
      ## ---------------- TAIL: GPD with mixing weight 1 - P, P = Pr(Y < u) ----------------
      P <- 0.0
      if (u > 0.0) {
        for (j in 1:K) {
          mu_j  <- mu[j]
          lam_j <- eta[j]
          
          if (mu_j > 0.0) {
            if (lam_j > 0.0) {
              ## inverse-Gaussian CDF at u
              s   <- sqrt(lam_j / u)
              t1  <- s * (u / mu_j - 1.0)
              t2  <- -s * (u / mu_j + 1.0)
              F1  <- pnorm(t1, 0.0, 1.0)
              F2  <- pnorm(t2, 0.0, 1.0)
              P   <- P + weights[j] * (F1 + exp(2.0 * lam_j / mu_j) * F2)
            }
          }
        }
      }
      
      ## GPD tail density
      z <- (x - u) / sigma
      
      ## avoid && and composite logicals
      cond <- 0
      if (xi < 0.0) {
        if ((1.0 + xi * z) <= 0.0) {
          cond <- 1
        }
      }
      
      if (cond == 1) {
        prob <- 0.0   ## outside finite endpoint
      } else {
        if (abs(xi) > 1e-06) {
          prob <- (1.0 - P) * (1.0 / sigma) * (1.0 + xi * z)^(-(1.0 + 1.0 / xi))
        } else {
          prob <- (1.0 - P) * (1.0 / sigma) * exp(-z)
        }
      }
    }
    
    if (log == 0) return(prob)
    else          return(log(prob))
  }
)

c_dTotal <- compileNimble(dTotal)

## ---- Sampler ----

rTotal <- nimbleFunction(
  run = function(n      = integer(0),
                 weights = double(1),
                 mu      = double(1),   # IG means
                 eta     = double(1),   # IG shapes λ_j
                 u       = double(0),
                 sigma   = double(0),
                 xi      = double(0)) {
    
    returnType(double(0))
    
    K <- length(weights)
    
    ## 1) Mixture mass below threshold: P = Pr(Y < u) under IG mixture
    P <- 0.0
    if (u > 0.0) {
      for (j in 1:K) {
        mu_j  <- mu[j]
        lam_j <- eta[j]
        
        if (mu_j > 0.0) {
          if (lam_j > 0.0) {
            s   <- sqrt(lam_j / u)
            t1  <- s * (u / mu_j - 1.0)
            t2  <- -s * (u / mu_j + 1.0)
            F1  <- pnorm(t1, 0.0, 1.0)
            F2  <- pnorm(t2, 0.0, 1.0)
            P   <- P + weights[j] * (F1 + exp(2.0 * lam_j / mu_j) * F2)
          }
        }
      }
    }
    
    y <- runif(1)
    x <- 0.0
    tailFlag <- 0L
    
    ## ---------- BULK BRANCH ----------
    if (y < P) {
      if (u > 0.0) {
        ## component masses under truncation: w_j * F_j(u)
        comp_mass  <- numeric(K)
        total_mass <- 0.0
        for (j in 1:K) {
          mu_j  <- mu[j]
          lam_j <- eta[j]
          
          if (mu_j > 0.0) {
            if (lam_j > 0.0) {
              s   <- sqrt(lam_j / u)
              t1  <- s * (u / mu_j - 1.0)
              t2  <- -s * (u / mu_j + 1.0)
              F1  <- pnorm(t1, 0.0, 1.0)
              F2  <- pnorm(t2, 0.0, 1.0)
              comp_mass[j] <- weights[j] * (F1 + exp(2.0 * lam_j / mu_j) * F2)
            } else {
              comp_mass[j] <- 0.0
            }
          } else {
            comp_mass[j] <- 0.0
          }
          total_mass <- total_mass + comp_mass[j]
        }
        
        if (total_mass <= 0.0) {
          tailFlag <- 1L
        } else {
          ## normalize conditional mixture weights given Y < u
          for (j in 1:K) {
            comp_mass[j] <- comp_mass[j] / total_mass
          }
          
          ## sample component index
          u_comp <- runif(1)
          cum    <- 0.0
          chosen <- 0L
          for (j in 1:K) {
            cum <- cum + comp_mass[j]
            if (chosen == 0L) {
              if (u_comp <= cum) {
                chosen <- j
              }
            }
          }
          if (chosen == 0L) {
            chosen <- K
          }
          
          ## rejection sampling from IG(mu_j, lam_j) truncated at u (finite loop, no break)
          accepted <- 0L
          x_candidate <- 0.0
          for (iter in 1:10000) {
            if (accepted == 0L) {
              mu_j  <- mu[chosen]
              lam_j <- eta[chosen]
              
              ## Michael–Schucany–Haas sampler for IG(μ, λ)
              z   <- rnorm(1, 0.0, 1.0)
              y2  <- z * z
              mu2 <- mu_j * mu_j
              t   <- mu_j +
                (mu2 * y2) / (2.0 * lam_j) -
                (mu_j / (2.0 * lam_j)) *
                sqrt(4.0 * mu_j * lam_j * y2 + mu2 * y2 * y2)
              u_r <- runif(1)
              
              if (u_r <= mu_j / (mu_j + t)) {
                x_candidate <- t
              } else {
                x_candidate <- mu2 / t
              }
              
              if (x_candidate < u) {
                if (x_candidate > 0.0) {
                  x        <- x_candidate
                  accepted <- 1L
                }
              }
            }
          }
          
          ## If somehow never accepted, fall back to mid-bulk
          if (accepted == 0L) {
            x <- u / 2.0
          }
        }
      } else {
        tailFlag <- 1L
      }
    } else {
      ## ---------- TAIL BRANCH ----------
      tailFlag <- 1L
    }
    
    if (tailFlag == 1L) {
      ## GPD tail by inversion of the splice, a = (1 - y) / (1 - P)
      denom <- 1.0 - P
      if (denom <= 0.0) denom <- 1e-12
      a <- (1.0 - y) / denom
      if (a < 0.0) a <- 0.0
      if (a > 1.0) a <- 1.0
      
      if (abs(xi) > 1e-06) {
        x <- u + (sigma / xi) * (a^(-xi) - 1.0)
      } else {
        x <- u - sigma * log(a)
      }
    }
    
    return(x)
  }
)

c_rTotal <- compileNimble(rTotal)

## ---- Quantile function (R level) ----

## R-level IG CDF: scalar q, vector (mu, lambda)
pinvnorm <- function(q, mu, lambda) {
  mu     <- as.numeric(mu)
  lambda <- as.numeric(lambda)
  K      <- length(mu)
  out    <- numeric(K)
  
  if (q <= 0) {
    ## all mass at 0 for IG
    return(out)  ## 0 for all components
  }
  
  x <- q
  s   <- sqrt(lambda / x)
  t1  <- s * (x / mu - 1)
  t2  <- -s * (x / mu + 1)
  F1  <- pnorm(t1)
  F2  <- pnorm(t2)
  out <- F1 + exp(2 * lambda / mu) * F2
  out
}

## Quantile Function: IG mixture + GPD tail
qTotal <- function(p, weights, mu, eta, u, sigma, xi) {
  ## mu: IG means; eta: IG shapes (lambda_j)
  
  mu      <- as.numeric(mu)
  eta     <- as.numeric(eta)
  weights <- as.numeric(weights)
  weights <- weights / sum(weights)
  
  ## Mixture mass below threshold
  Fj_u <- pinvnorm(u, mu = mu, lambda = eta)
  P    <- sum(weights * Fj_u)
  
  if (p < P) {
    ## -------- bulk quantile: solve sum w_j F_j(x) = p on (0, u) --------
    f_root <- function(x) {
      sum(weights * pinvnorm(x, mu = mu, lambda = eta)) - p
    }
    ## bracket quantile in (0, u); assumes bulk mass is below u
    r <- uniroot(f_root, interval = c(1e-10, u), maxiter = 1e4)
    x <- r$root
  } else {
    ## -------- tail quantile: GPD on (u, ∞) --------
    denom <- 1 - P
    if (denom <= 0) denom <- 1e-12
    a <- (p - P) / denom
    if (a < 0) a <- 0
    if (a > 1) a <- 1
    x <- qgpd(a, loc = u, scale = sigma, shape = xi)
  }
  
  x
}
