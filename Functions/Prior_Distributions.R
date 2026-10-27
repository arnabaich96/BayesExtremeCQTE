# Custom NIMBLE distributions used by the MixGPD model code (Model_Code.R).
library(nimble)

# ---------------------------------------------------------------------------
# Lognormal truncated below at `lower`, for the covariate-dependent threshold u(x) >= l_u.
# The normalising constant 1 - F(lower) is computed on the log scale; NIMBLE's generic
# T() truncation underflows when exp(meanlog) is far below `lower`.
# ---------------------------------------------------------------------------
dlnormTrunc <- nimbleFunction(
  run = function(x       = double(0),
                 meanlog = double(0),
                 sdlog   = double(0),
                 lower   = double(0),
                 log     = integer(0, default = 0)) {
    returnType(double(0))
    if (x < lower) {
      if (log == 1) return(-Inf)
      else          return(0.0)
    }
    log_tail <- pnorm(log(lower), mean = meanlog, sd = sdlog, lower.tail = FALSE, log.p = TRUE)
    logdens  <- dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE) - log_tail
    if (log == 1) return(logdens)
    else          return(exp(logdens))
  }
)

rlnormTrunc <- nimbleFunction(
  run = function(n       = integer(0),
                 meanlog = double(0),
                 sdlog   = double(0),
                 lower   = double(0)) {
    returnType(double(0))
    if (n != 1) print("rlnormTrunc only allows n = 1;")
    p_low <- plnorm(lower, meanlog = meanlog, sdlog = sdlog)
    x <- qlnorm(p_low + runif(1) * (1 - p_low), meanlog = meanlog, sdlog = sdlog)
    if (x < lower) x <- lower
    return(x)
  }
)

# R-level conditional median of the truncated lognormal, used for prediction in Sim()
median_lnormTrunc <- function(meanlog, sdlog, lower) {
  log_tail <- pnorm(log(lower), meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)
  # median solves P(U > m) = (1 - F(lower)) / 2, evaluated on the log scale
  pmax(lower, exp(qnorm(log_tail - log(2), meanlog, sdlog, lower.tail = FALSE, log.p = TRUE)))
}
