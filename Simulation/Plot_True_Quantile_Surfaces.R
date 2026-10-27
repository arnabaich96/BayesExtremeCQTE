# Plot_True_Quantile_Surfaces.R
# True tau-quantile surfaces of the treatment and control outcomes and of the CQTE in scenario S1
# (Supplementary Material), from one example data set with N = 1000 and 5 replicates.
# Output: Results/Simulation/Figures/Quantile_Surfaces_X1X2_p<tau>_<BW|Color>.png and ..._X3X4_...
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
})
source(here("Simulation/my_data_sim.R"))

## ============================================================
## Example data set: print the true CQTE and the threshold table
## ============================================================
n_rep_ex   <- 5
N_ex       <- 1e3
beta_ex    <- c(-0.03, -0.4, 0.12, 0.55)
x_index_ex <- seq(0.1, 0.9, length.out = 9)
p_ex       <- c(0.75, 0.5, 0.25, 0.05, 10/N_ex, 5/N_ex,
                1/(N_ex^1.25), 1/(N_ex * log(N_ex)), 1/(N_ex^1.5))
p_QTE_ex   <- 1 - p_ex

message("Running single example to print TRUE QTE and threshold table...")
Data_ex <- my_data_sim(n_rep_ex, N_ex,
                       n_cores = max(1, parallel::detectCores()-1),
                       beta = beta_ex,
                       treatment_params = treatment_params,
                       control_params   = control_params,
                       p_QTE = p_QTE_ex, x_index = x_index_ex, seed = 123)

cat("\n==== TRUE QTE (rows: tau_x quantiles of X; cols: tau_y) ====\n")
print(Data_ex$True_CQTE)

cat("\n==== TRUE THRESHOLD TABLE (9 x 4: u1, p1, u0, p0) ====\n")
print(Data_ex$True_tail)


## ============================================================
## TRUE tau-quantile surfaces (Treatment / Control / QTE)
## ============================================================

plot_true_quantile_surfaces <- function(
    DATA,
    beta,
    tr_par,
    co_par,
    p,                     # quantile level
    out_dir  = here::here("Results/Simulation/Figures/"),
    n_grid   = 100,
    x1_range = NULL,
    x2_range = NULL,
    x3_range = NULL,
    x4_range = NULL,
    palette  = c("color","bw")
){
  palette <- match.arg(palette)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  X_mat <- DATA$X
  colnames(X_mat) <- paste0("X", 1:ncol(X_mat))
  x_meds <- apply(X_mat, 2, median)
  
  ## helper: true Q_tau for given x-vector
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
  suffix <- ifelse(palette == "color", "Color", "BW")
  
  fill_scale <- if (palette == "color") {
    scale_fill_viridis_c(option = "C")
  } else {
    scale_fill_gradient(low = "white", high = "black")
  }
  contour_col  <- if (palette == "color") "white" else "black"
  contour_size <- if (palette == "color") 0.9 else 1.0
  
  ## PANEL 1: X1, X2 vary
  if (is.null(x1_range)) {
    rng1 <- quantile(X_mat[, "X1"], c(0.05, 0.95), na.rm = TRUE)
  } else rng1 <- x1_range
  
  if (is.null(x2_range)) {
    rng2 <- quantile(X_mat[, "X2"], c(0.05, 0.95), na.rm = TRUE)
  } else rng2 <- x2_range
  
  grid_12 <- expand.grid(
    X1 = seq(rng1[1], rng1[2], length.out = n_grid),
    X2 = seq(rng2[1], rng2[2], length.out = n_grid)
  )
  grid_12$X3 <- x_meds["X3"]
  grid_12$X4 <- x_meds["X4"]
  grid_12 <- grid_12[, c("X1","X2","X3","X4")]
  
  surf_12 <- compute_surface_for_grid(grid_12)
  
  surf_12_long <- surf_12 %>%
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
      x_plot = X1,
      y_plot = X2
    )
  
  contour_breaks_12 <- pretty(surf_12_long$Qv, n = 10)
  
  p_12 <- ggplot(surf_12_long, aes(x = x_plot, y = y_plot)) +
    geom_tile(aes(fill = Qv)) +
    geom_contour(
      aes(z = Qv),
      colour   = contour_col,
      linewidth= contour_size,
      alpha    = 0.9,
      breaks   = contour_breaks_12
    ) +
    facet_grid(. ~ scenario) +
    fill_scale +
    labs(
      x    = expression(x[1]),
      y    = expression(x[2]),
      fill = legend_label
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_blank()
    )
  ggsave(
    file.path(out_dir, sprintf("Quantile_Surfaces_X1X2_p%.2f_%s.png", p, suffix)),
    p_12, width = 10, height = 4, dpi = 300
  )
  
  ## PANEL 2: X3, X4 vary
  if (is.null(x3_range)) {
    rng3 <- quantile(X_mat[, "X3"], c(0.05, 0.95), na.rm = TRUE)
  } else rng3 <- x3_range
  
  if (is.null(x4_range)) {
    rng4 <- quantile(X_mat[, "X4"], c(0.05, 0.95), na.rm = TRUE)
  } else rng4 <- x4_range
  
  grid_34 <- expand.grid(
    X3 = seq(rng3[1], rng3[2], length.out = n_grid),
    X4 = seq(rng4[1], rng4[2], length.out = n_grid)
  )
  grid_34$X1 <- x_meds["X1"]
  grid_34$X2 <- x_meds["X2"]
  grid_34 <- grid_34[, c("X1","X2","X3","X4")]
  
  surf_34 <- compute_surface_for_grid(grid_34)
  
  surf_34_long <- surf_34 %>%
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
      x_plot = X3,
      y_plot = X4
    )
  
  contour_breaks_34 <- pretty(surf_34_long$Qv, n = 10)
  
  p_34 <- ggplot(surf_34_long, aes(x = x_plot, y = y_plot)) +
    geom_tile(aes(fill = Qv)) +
    geom_contour(
      aes(z = Qv),
      colour   = contour_col,
      linewidth= contour_size,
      alpha    = 0.9,
      breaks   = contour_breaks_34
    ) +
    facet_grid(. ~ scenario) +
    fill_scale +
    labs(
      x    = expression(x[3]),
      y    = expression(x[4]),
      fill = legend_label
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_blank()
    )
  ggsave(
    file.path(out_dir, sprintf("Quantile_Surfaces_X3X4_p%.2f_%s.png", p, suffix)),
    p_34, width = 10, height = 4, dpi = 300
  )
  
  invisible(list(panel12 = p_12, panel34 = p_34))
}

# --- plot TRUE quantile surfaces for example data --- 
for(p in c(0.25,0.5,0.75,0.95)){
  for(type in c("color","bw")){
    message("Plotting TRUE quantile surfaces for p = ", p, ", palette = ", type, " ...")
    plot_true_quantile_surfaces(
      DATA   = Data_ex,
      beta   = beta_ex,
      tr_par = treatment_params,
      co_par = control_params,
      p      = p,
      palette = type
    )
  }
}

