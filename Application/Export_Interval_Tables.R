# LaTeX tables of the application (Tables 3-5 and Supplementary Table S7) from the saved CQTE draws.
# Each cell: posterior mean (90% credible interval) [P(CQTE > 0)] of the CQTE averaged draw by draw over profiles.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

tab_dir  <- here("Results/Application/Diabetes/Tables")
draw_dir <- here("Application/output/draws")
CP <- 0.90
lo_p <- (1 - CP) / 2; hi_p <- 1 - lo_p

load_cqte <- function(tags, kern) {
  Q <- NULL
  for (tag in tags) {
    d <- readRDS(file.path(draw_dir, sprintf("draws_%s_%s.rds", tag, kern)))
    cq <- Map(`-`, d$Q_trt, d$Q_con)
    Q <- if (is.null(Q)) cq else Map(rbind, Q, cq)
  }
  Q
}
grid_of <- function(job, kern) {
  prof <- read.csv(file.path(tab_dir, sprintf("MixGPD_%s_%s_Profiles.csv", job, kern)))
  taus <- sort(unique(prof$quantile))
  list(grid = prof[prof$quantile == taus[1], ], taus = taus)
}
summ <- function(draws, idx) {
  avg <- rowMeans(draws[, idx, drop = FALSE], na.rm = TRUE)
  c(est = mean(avg, na.rm = TRUE), lo = unname(quantile(avg, lo_p, na.rm = TRUE)),
    hi = unname(quantile(avg, hi_p, na.rm = TRUE)), p = mean(avg > 0, na.rm = TRUE))
}
fmt  <- function(x) formatC(x, format = "f", digits = 2)
cell <- function(s, prob = TRUE) {
  out <- paste0(fmt(s[["est"]]), " (", fmt(s[["lo"]]), ", ", fmt(s[["hi"]]), ")")
  if (prob) out <- paste0(out, " [", fmt(s[["p"]]), "]")
  out
}
tau_label <- function(t) ifelse(t < 0.99, sub("0+$", "", formatC(t, format = "f", digits = 3)),
                                formatC(t, format = "f", digits = 4))

jobs_by <- list(main = "main", lu50 = "lu50", lu70 = "lu70", lu80 = "lu80", trim = "trim",
                noglu = "noglu", xisd10 = "xisd10")

# ---------- Stratified tables (main analysis, MixInvNormGPD and MixGammaGPD) ----------
strat_table <- function(var, labels, file) {
  g <- grid_of("main", "normal"); grid <- g$grid; taus <- g$taus
  Q <- list(normal = load_cqte("main", "normal"), gamma = load_cqte("main", "gamma"))
  levels <- sort(unique(grid[[var]]))
  rows <- sapply(seq_along(taus), function(j) {
    cells <- unlist(lapply(c("normal", "gamma"), function(k) {
      sapply(levels, function(lev) cell(summ(Q[[k]][[j]], which(grid[[var]] == lev))))
    }))
    paste0(tau_label(taus[j]), " & ", paste(cells, collapse = " & "), " \\\\")
  })
  writeLines(rows, file.path(tab_dir, file))
  cat("Wrote", file, "| levels of", var, ":", paste(labels[as.character(levels)], collapse = ", "), "\n")
}
g0 <- grid_of("main", "gamma")$grid
strat_table("gender", c("1" = "Female", "2" = "Male"), "table_app_gender.tex")
strat_table("age", setNames(as.character(sort(unique(g0$age))), sort(unique(g0$age))), "table_app_age.tex")
strat_table("stab.glu", setNames(as.character(sort(unique(g0$stab.glu))), sort(unique(g0$stab.glu))), "table_app_glu.tex")
cat("Grid levels: age", sort(unique(g0$age)), "| stab.glu", sort(unique(g0$stab.glu)),
    "| ratio", sort(unique(g0$ratio)), "\n")

# ---------- Sensitivity table: CQTE averaged over all profiles ----------
sens_taus_idx <- c(2, 4, 6, 7, 9)   # tau = 0.5, 0.95, 1 - 5/n, 1 - n^-1.25, 1 - n^-1.5
sens_rows <- list(
  c("Main analysis (plug-in PS, $l_u$ at 60th percentile)", "main"),
  c("$l_u$ at 50th percentile", "lu50"),
  c("$l_u$ at 70th percentile", "lu70"),
  c("$l_u$ at 80th percentile", "lu80"),
  c("Trimmed to PS in $[0.1, 0.9]$", "trim"),
  c("Without stabilized glucose", "noglu"),
  c("Prior standard deviation of $\\xi$ equal to 0.1", "xisd10"))
sens_csv <- list()
sens_lines <- unlist(lapply(c(gamma = "MixGammaGPD", normal = "MixInvNormGPD"), function(kname) {
  k <- names(which(c(gamma = "MixGammaGPD", normal = "MixInvNormGPD") == kname))
  out <- c(sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", length(sens_taus_idx) + 1, kname))
  for (r in sens_rows) {
    tags <- jobs_by[[r[2]]]
    if (!all(file.exists(file.path(draw_dir, sprintf("draws_%s_%s.rds", tags, k))))) next
    Q <- load_cqte(tags, k)
    s <- lapply(sens_taus_idx, function(j) summ(Q[[j]], seq_len(ncol(Q[[j]]))))
    out <- c(out, paste0(r[1], " & ", paste(sapply(s, cell, prob = FALSE), collapse = " & "), " \\\\"))
    sens_csv[[length(sens_csv) + 1]] <<- data.frame(kernel = kname, analysis = r[2],
                                                    tau = grid_of("main", k)$taus[sens_taus_idx],
                                                    do.call(rbind, s))
  }
  c(out, "\\midrule")
}))
writeLines(head(sens_lines, -1), file.path(tab_dir, "table_app_sensitivity.tex"))
write.csv(do.call(rbind, sens_csv), file.path(tab_dir, "App_Sensitivity_Summary.csv"), row.names = FALSE)
cat("Wrote table_app_sensitivity.tex\n")
print(do.call(rbind, sens_csv), digits = 3)
