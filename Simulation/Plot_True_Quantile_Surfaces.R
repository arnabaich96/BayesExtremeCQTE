# True quantile surfaces of the two arms and of the CQTE in scenario S1 at tau = 0.5 and 0.95
# (Supplementary Figures S1 and S2), over covariate ranges and medians of an example sample (N = 1000).
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
source(here("Simulation/my_data_sim.R"))

beta_ex <- c(-0.03, -0.4, 0.12, 0.55)
X_ex    <- sim_covariates(N = 1e3, beta = beta_ex, seed = 123)$X

plot_true_quantile_surfaces <- function(
    X_mat,
    beta,
    tr_par,
    co_par,
    p,                     # quantile level
    out_dir  = here::here("Results/Simulation/Figures/"),
    n_grid   = 100
){
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  colnames(X_mat) <- paste0("X", 1:ncol(X_mat))
  x_meds <- apply(X_mat, 2, median)

  ## true Q_tau at a covariate vector
  compute_q_true <- function(x_vec, beta, par, p, family) {
    ps <- plogis(sum(x_vec * beta))
    K  <- length(par$wt)

    mu_vec <- sapply(1:K, function(k){
      lp <- par$theta_ps[k] * ps + sum(x_vec * par$theta_X[, k])
      exp(lp)
    })

    u_val <- exp(par$del_ps * ps + sum(x_vec * par$del_x))

    if (family == "gamma") {
      qTotal_gamma(p, par$wt, mu_vec, par$eta, u_val, par$sigma, par$xi)
    } else {
      qTotal_lnorm(p, par$wt, mu_vec, par$eta, u_val, par$sigma, par$xi)
    }
  }

  compute_surface_for_grid <- function(grid){
    n_g <- nrow(grid)
    q_trt <- q_con <- q_qte <- numeric(n_g)

    for(i in seq_len(n_g)){
      x_vec <- as.numeric(grid[i, ])
      names(x_vec) <- colnames(grid)

      qt <- compute_q_true(x_vec, beta, tr_par, p, "gamma")
      qc <- compute_q_true(x_vec, beta, co_par, p, "lnorm")

      q_trt[i] <- qt
      q_con[i] <- qc
      q_qte[i] <- qt - qc
    }

    tibble::tibble(grid, q_trt, q_con, q_qte)
  }

  legend_label <- bquote(Q[.(p)])
  fill_scale   <- scale_fill_gradient(low = "white", high = "black")

  plot_pair <- function(v1, v2, fixed, xlab, ylab, file_tag) {
    rng1 <- quantile(X_mat[, v1], c(0.05, 0.95), na.rm = TRUE)
    rng2 <- quantile(X_mat[, v2], c(0.05, 0.95), na.rm = TRUE)
    grid <- expand.grid(seq(rng1[1], rng1[2], length.out = n_grid),
                        seq(rng2[1], rng2[2], length.out = n_grid))
    names(grid) <- c(v1, v2)
    for (v in fixed) grid[[v]] <- x_meds[v]
    grid <- grid[, c("X1","X2","X3","X4")]

    surf_long <- compute_surface_for_grid(grid) %>%
      tidyr::pivot_longer(
        cols = c(q_trt, q_con, q_qte),
        names_to = "scenario",
        values_to = "Qv"
      ) %>%
      dplyr::mutate(
        scenario = dplyr::recode(
          scenario,
          q_trt = "Treatment",
          q_con = "Control",
          q_qte = "QTE"
        ),
        scenario = factor(scenario, levels = c("Treatment","Control","QTE")),
        x_plot = .data[[v1]],
        y_plot = .data[[v2]]
      )

    plt <- ggplot(surf_long, aes(x = x_plot, y = y_plot)) +
      geom_tile(aes(fill = Qv)) +
      geom_contour(
        aes(z = Qv),
        colour   = "black",
        linewidth= 1.0,
        alpha    = 0.9,
        breaks   = pretty(surf_long$Qv, n = 10)
      ) +
      facet_grid(. ~ scenario) +
      fill_scale +
      labs(
        x    = xlab,
        y    = ylab,
        fill = legend_label
      ) +
      theme_classic() +
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_blank()
      )
    ggsave(
      file.path(out_dir, sprintf("Quantile_Surfaces_%s_p%.2f_BW.png", file_tag, p)),
      plt, width = 10, height = 4, dpi = 300
    )
    plt
  }

  invisible(list(
    panel12 = plot_pair("X1", "X2", c("X3", "X4"), expression(x[1]), expression(x[2]), "X1X2"),
    panel34 = plot_pair("X3", "X4", c("X1", "X2"), expression(x[3]), expression(x[4]), "X3X4")
  ))
}

for (p in c(0.5, 0.95)) {
  message("Plotting the true quantile surfaces for p = ", p)
  plot_true_quantile_surfaces(X_ex, beta_ex, treatment_params, control_params, p)
}
