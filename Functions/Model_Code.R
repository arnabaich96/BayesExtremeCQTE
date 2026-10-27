# NIMBLE model code shared by the simulation (Report.R) and the application (Report_RealData.R).
# Gamma and inverse-Gaussian kernels share one template; dTotal comes from the kernel file
# sourced in Sim(). The lognormal kernel places the component location on the log scale.
# Hyperparameters and the threshold lower bound l_u enter as constants (see mixgpd_hyper());
# dlnormTrunc is defined in Prior_Distributions.R.
library(nimble)

mixgpd_code <- function(kern = c("gamma", "normal", "lognormal")) {
  kern <- match.arg(kern)
  if (kern == "lognormal") {
    nimbleCode({
      # dimensions: n (obs), m (mixture comps), p (predictors = PS + X)
      alpha ~ dgamma(shape = 1, scale = alpha_scale)
      for (i in 1:(m-1)) {
        v[i] ~ dbeta(1, alpha)
      }
      weights[1:m] <- stick_breaking(v[1:(m-1)])

      # mixture component parameters
      for (j in 1:m) {
        eta[j] ~ dgamma(shape = eta_shape, scale = 1)
        theta1[1, j] ~ dnorm(mean = 0, sd = theta_sd1)   # PS coefficient
        for (k in 2:p) {
          theta1[k, j] ~ dnorm(mean = 0, sd = theta_sd)  # X coefficients
        }
      }

      # threshold regression coefficients (PS+X, no intercept)
      delta1[1] ~ dnorm(mean = 0, sd = delta_sd1)
      for (k in 2:p) {
        delta1[k] ~ dnorm(mean = 0, sd = delta_sd)
      }

      # tail parameters (informative priors of the original analysis; see mixgpd_hyper())
      xi ~ T(dnorm(mean = 0.25, sd = xi_sd), -0.5, Inf)
      sigma ~ dgamma(shape = 2, scale = sigma_scale)

      # dispersions of the latent locations and thresholds
      var1 ~ dinvgamma(shape = 2, scale = 1)
      var2 ~ dinvgamma(shape = 2, scale = 1)

      for (i in 1:n) {
        for (j in 1:m) {
          lp_mu[i, j] <- inprod(theta1[1:p, j], X[i, 1:p])  # linear predictor
          mu[i, j]    ~ dnorm(mean = lp_mu[i, j], sd = 1/sqrt(var1))  # meanlog
        }

        # threshold: lognormal around linear predictor, truncated below at l_u
        lp_u[i] <- inprod(delta1[1:p], X[i, 1:p])
        u[i] ~ dlnormTrunc(meanlog = lp_u[i], sdlog = 1/sqrt(var2), lower = l_u)

        y[i] ~ dTotal(weights[1:m], mu[i, 1:m], eta[1:m], u[i], sigma, xi)
      }
    })
  } else {
    nimbleCode({
      alpha ~ dgamma(shape = 1, scale = alpha_scale)
      for (i in 1:(m-1)) {
        v[i] ~ dbeta(1, alpha)
      }
      weights[1:m] <- stick_breaking(v[1:(m-1)])

      for (j in 1:m) {
        eta[j] ~ dgamma(shape = eta_shape, scale = 1)
        theta1[1, j] ~ dnorm(mean = 0, sd = theta_sd1)
        for (k in 2:p) {
          theta1[k, j] ~ dnorm(mean = 0, sd = theta_sd)
        }
      }

      delta1[1] ~ dnorm(mean = 0, sd = delta_sd1)
      for (k in 2:p) {
        delta1[k] ~ dnorm(mean = 0, sd = delta_sd)
      }

      xi ~ T(dnorm(mean = 0.25, sd = xi_sd), -0.5, Inf)
      sigma ~ dgamma(shape = 2, scale = sigma_scale)

      var1 ~ dinvgamma(shape = 2, scale = 1)
      var2 ~ dinvgamma(shape = 2, scale = 1)

      for (i in 1:n) {
        for (j in 1:m) {
          lp_mu[i, j] <- inprod(theta1[1:p, j], X[i, 1:p])
          mu[i, j]    ~ dlnorm(meanlog = lp_mu[i, j], sdlog = 1/sqrt(var1))  # component mean
        }

        lp_u[i] <- inprod(delta1[1:p], X[i, 1:p])
        u[i] ~ dlnormTrunc(meanlog = lp_u[i], sdlog = 1/sqrt(var2), lower = l_u)

        y[i] ~ dTotal(weights[1:m], mu[i, 1:m], eta[1:m], u[i], sigma, xi)
      }
    })
  }
}

# Prior hyperparameters actually used. "sim" = simulation study, "app" = diabetes application.
mixgpd_hyper <- function(setting = c("sim", "app"), kern = c("gamma", "normal", "lognormal"),
                         arm = c("treatment", "control")) {
  setting <- match.arg(setting); kern <- match.arg(kern); arm <- match.arg(arm)
  trt <- arm == "treatment"
  # GPD priors: xi ~ N(0.25, xi_sd^2) truncated to (-0.5, Inf); sigma ~ Gamma(2, scale = sigma_scale)
  if (setting == "sim") {
    list(alpha_scale = 1, eta_shape = if (trt) 2 else 1,
         theta_sd1 = 1, theta_sd = 1, delta_sd1 = 1, delta_sd = 1,
         xi_sd = 0.1, sigma_scale = if (trt) 2 else 1)
  } else {
    ln_trt <- trt && kern == "lognormal"
    list(alpha_scale = if (trt) 0.75 else 0.5, eta_shape = if (trt) 2 else 1,
         theta_sd1 = if (ln_trt) 2 else 1,
         theta_sd  = if (ln_trt) 1 else 0.5,
         delta_sd1 = if (ln_trt) 2 else 1,
         delta_sd  = if (trt && kern != "gamma") 1 else 0.5,
         xi_sd = 0.05, sigma_scale = if (trt) 2 else 1)
  }
}
