# DP mixture and Generalized Pareto Distribution
library(nimble)
library(evd)
# Density Function

dTotal <- nimbleFunction(
  run = function(x      = double(0),
                 weights = double(1),
                 mu      = double(1),
                 eta     = double(1),   # shape
                 u       = double(0),
                 sigma   = double(0),
                 xi      = double(0),
                 log     = integer(0, default = 0)) {
    
    returnType(double(0))
    
    K <- length(weights)
    prob <- 0.0
    
    if (x < u) {
      ## Bulk: Gamma mixture, mean mu[j], shape eta[j], so scale = mu[j]/eta[j]
      for (j in 1:K) {
        prob <- prob + dgamma(x,
                              shape = eta[j],
                              scale = mu[j] / eta[j]) * weights[j]
      }
      
    } else {
      ## Tail: GPD with mixing weight = 1 - P(Y < u)
      P <- 0.0
      for (j in 1:K) {
        P <- P + weights[j] * pgamma(u,
                                     shape = eta[j],
                                     scale = mu[j] / eta[j])
      }
      
      z <- (x - u) / sigma
      
      ## Optional: guard finite endpoint for xi < 0
      cond <- (xi < 0) & ((1 + xi * z) <= 0)
      if (cond) {
        prob <- 0.0
      } else {
        if (abs(xi) > 1e-06) {
          prob <- (1 - P) * (1 / sigma) * (1 + xi * z)^(-(1 + 1 / xi))
        } else {
          prob <- (1 - P) * (1 / sigma) * exp(-z)
        }
      }
    }
    
    if (log == 0) return(prob)
    else          return(log(prob))
  }
)

 compileNimble(dTotal)

# Generating Function
rTotal <- nimbleFunction(
  run = function(n      = integer(0),
                 weights = double(1),
                 mu      = double(1),
                 eta     = double(1),
                 u       = double(0),
                 sigma   = double(0),
                 xi      = double(0)) {
    
    returnType(double(0))
    
    if (n != 1) {
      print("rTotal only allows n = 1;")
    }
    
    K <- length(weights)
    
    ## 1) Mixture mass below threshold: P = Pr(Y < u)
    P <- 0.0
    for (j in 1:K) {
      P <- P + weights[j] * pgamma(u,
                                   shape = eta[j],
                                   scale = mu[j] / eta[j])
    }
    
    y <- runif(1)
    x <- 0.0
    tailFlag <- 0L
    
    ## ---------- BULK BRANCH ----------
    if (y < P) {
      ## component masses under truncation: w_j * F_j(u)
      comp_mass <- numeric(K)
      total_mass <- 0.0
      for (j in 1:K) {
        comp_mass[j] <- weights[j] * pgamma(u,
                                            shape = eta[j],
                                            scale = mu[j] / eta[j])
        total_mass <- total_mass + comp_mass[j]
      }
      
      if (total_mass == 0.0) {
        ## numerically no bulk -> send to tail
        tailFlag <- 1L
      } else {
        ## normalize to conditional mixture weights given Y < u
        for (j in 1:K) {
          comp_mass[j] <- comp_mass[j] / total_mass
        }
        
        ## sample component index using these weights
        u_comp <- runif(1)
        cum <- 0.0
        chosen <- 0L
        for (j in 1:K) {
          cum <- cum + comp_mass[j]
          if ((chosen == 0L) & (u_comp <= cum)) {
            chosen <- j
          }
        }
        
        ## truncated Gamma within chosen component: F_j(X) ~ U(0, F_j(u))
        Fj_u <- pgamma(u,
                       shape = eta[chosen],
                       scale = mu[chosen] / eta[chosen])
        
        if (Fj_u == 0.0) {
          tailFlag <- 1L
        } else {
          u_trunc <- runif(1) * Fj_u
          x <- qgamma(u_trunc,
                      shape = eta[chosen],
                      scale = mu[chosen] / eta[chosen])
        }
      }
      
    } else {
      tailFlag <- 1L
    }
    
    ## ---------- TAIL BRANCH ----------
    if (tailFlag == 1L) {
      ## a = (1 - y) / (1 - P) solves F(x) = y for the splice
      a <- (1 - y) / (1 - P)
      if (abs(xi) > 1e-06) {
        x <- u + (sigma / xi) * (a^(-xi) - 1)
      } else {
        x <- u - sigma * log(a)
      }
    }
    
    return(x)
  }
)

compileNimble(rTotal)


# Quantile Function
qTotal <- function(p,weights,mu,eta,u,sigma,xi)
{
  P = sum(pgamma(u,shape = eta , scale = mu/eta)*weights)
  if(p<P){
    r = uniroot(function(x) p - sum(weights*pgamma(x, shape = eta, scale = mu/eta))
              ,c(0,1e10),maxiter = 1e4)
    x = r$root
    }else{
    a = (p-P)/(1-P)
    x = qgpd(a, loc = u,scale = sigma, shape =  xi)
  }
  return(x)
}



