# Heatmaps of the simulation metrics over the 9 quantile levels x 9 covariate profiles (Supplementary
# Section S3), from Results/Simulation/Tables/*.csv; output in Results/Simulation/Figures.
library(here)
library(tidyverse)   # dplyr, readr, purrr, stringr, ggplot2

# ---------- Helpers ----------

.get_method_raw <- function(basename_str) sub("_.*$", "", basename_str)

.pretty_method_label <- function(method_raw) {
  mlow <- tolower(method_raw)
  dplyr::case_when(
    mlow == "mixgpd"       ~ "Mix-GPD",
    mlow == "extremeqr"    ~ "Extreme QR",
    mlow == "skewt"        ~ "Skew t",
    mlow == "ipw"          ~ "IPW",
    TRUE                   ~ method_raw
  )
}

.resolve_metric <- function(metric, available_names) {
  m <- tolower(metric)
  candidates <- switch(
    m,
    "bias"     = c("bias"),
    "coverage" = c("coverage"),
    "length"   = c("length_ci", "length"),
    "sd"       = c("SD", "sd"),
    c(metric)
  )
  hit <- candidates[candidates %in% available_names]
  if (!length(hit)) {
    stop(
      "Metric column for '", metric, "' not found. Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable: ", paste(available_names, collapse = ", ")
    )
  }
  hit[[1]]
}

# ---------- Load all tables once ----------
load_all_sim_tables <- function(table_dir) {
  if (!dir.exists(table_dir)) stop("Directory not found: ", table_dir)
  all_csv <- list.files(table_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(all_csv)) stop("No CSV files in: ", table_dir)

  data_all <- all_csv %>%
    purrr::map_dfr(function(path) {
      base <- basename(path)
      readr::read_csv(path, show_col_types = FALSE) %>%
        mutate(
          filename    = path,
          file_base   = base,
          method_raw  = .get_method_raw(base),
          scenario    = stringr::str_match(base, "_(spliced|burr|lnmix)_")[, 2],
          sample_size = suppressWarnings(
            as.numeric(stringr::str_extract(base, "\\d+"))
          )
        )
    })

  message("Loaded ", nrow(data_all), " rows from ", length(all_csv), " files.")
  message("Methods: ",
          paste(sort(unique(data_all$method_raw)), collapse = ", "))
  invisible(data_all)
}

# ---------- Broadcast singletons across rows/cols ----------
.broadcast_singletons <- function(heat_df) {
  # rows
  heat_df <- heat_df %>%
    dplyr::group_by(method, tau_y) %>%
    dplyr::mutate(
      non_na_row = sum(!is.na(value)),
      row_val    = ifelse(non_na_row == 1,
                          value[which(!is.na(value))[1]],
                          NA_real_)
    ) %>%
    dplyr::mutate(
      value = ifelse(non_na_row == 1, row_val, value)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-non_na_row, -row_val)

  # columns
  heat_df <- heat_df %>%
    dplyr::group_by(method, tau_x) %>%
    dplyr::mutate(
      non_na_col = sum(!is.na(value)),
      col_val    = ifelse(non_na_col == 1,
                          value[which(!is.na(value))[1]],
                          NA_real_)
    ) %>%
    dplyr::mutate(
      value = ifelse(non_na_col == 1, col_val, value)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-non_na_col, -col_val)

  heat_df
}

.format_quantile_labels <- function(values) {
  vapply(values, function(value) {
    if (value < 0.99) {
      format(value, trim = TRUE, nsmall = 2)
    } else {
      format(round(value, 4), trim = TRUE, nsmall = 4)
    }
  }, character(1))
}

# ---------- Heatmap data: one row of panels per method family ----------

.build_heat_df <- function(df, metric) {
  df <- df %>%
    filter(!is.na(x_num), !is.na(y_num))

  # common tau_y across "method"
  y_by_method <- df %>%
    group_by(method) %>%
    summarise(y_vals = list(sort(unique(y_num))), .groups = "drop") %>%
    pull(y_vals)
  common_y <- Reduce(intersect, y_by_method)
  common_y <- sort(common_y)

  # common tau_x
  x_by_method <- df %>%
    group_by(method) %>%
    summarise(x_vals = list(sort(unique(x_num))), .groups = "drop") %>%
    pull(x_vals)
  common_x <- Reduce(intersect, x_by_method)
  common_x <- sort(common_x)

  x_ref <- x_by_method[[1]]
  y_ref <- y_by_method[[1]]

  if (!length(common_x)) {
    keep_x <- seq_len(min(vapply(x_by_method, length, integer(1))))
    x_labels <- .format_quantile_labels(x_ref[keep_x])
    df <- df %>%
      group_by(method) %>%
      mutate(x_key = match(x_num, sort(unique(x_num)))) %>%
      ungroup() %>%
      filter(x_key %in% keep_x)
  } else {
    keep_x <- common_x
    x_labels <- .format_quantile_labels(keep_x)
    df <- df %>%
      mutate(x_key = x_num) %>%
      filter(x_key %in% keep_x)
  }

  if (!length(common_y)) {
    keep_y <- seq_len(min(vapply(y_by_method, length, integer(1))))
    y_labels <- .format_quantile_labels(y_ref[keep_y])
    df <- df %>%
      group_by(method) %>%
      mutate(y_key = match(y_num, sort(unique(y_num)))) %>%
      ungroup() %>%
      filter(y_key %in% keep_y)
  } else {
    keep_y <- common_y
    y_labels <- .format_quantile_labels(keep_y)
    df <- df %>%
      mutate(y_key = y_num) %>%
      filter(y_key %in% keep_y)
  }

  df <- df %>%
    mutate(
      tau_x = factor(x_key, levels = keep_x, labels = x_labels),
      tau_y = factor(y_key, levels = keep_y, labels = y_labels)
    )

  mcol <- .resolve_metric(metric, names(df))

  heat_df <- df %>%
    group_by(method, tau_x, tau_y) %>%
    summarise(value = mean(.data[[mcol]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(method,
                    tau_y,
                    tau_x,
                    fill = list(value = NA_real_))

  .broadcast_singletons(heat_df)
}

.resolve_y_positions <- function(data_all, sample_size, y_positions = NULL) {
  y_ref <- data_all %>%
    filter(!is.na(sample_size), sample_size == !!sample_size) %>%
    mutate(y_num = as.numeric(as.character(y_index))) %>%
    pull(y_num) %>%
    unique() %>%
    sort()

  if (is.null(y_positions)) {
    return(y_ref)
  }

  y_ref[y_positions]
}

.plot_competing_family <- function(
    data_all,
    sample_size,
    metric,
    methods,
    file_prefix,
    y_positions = NULL,
    save_dir,
    scenario = "spliced"
) {
  data_all <- data_all %>% filter(!is.na(.data$scenario), .data$scenario == !!scenario,
                                  !grepl("_(lu|ps)\\d+", file_base))
  keep_y <- .resolve_y_positions(data_all, sample_size, y_positions = y_positions)
  methods_lc <- tolower(methods)
  method_levels <- .pretty_method_label(methods)

  df <- data_all %>%
    filter(
      !is.na(sample_size),
      sample_size == !!sample_size,
      tolower(method_raw) %in% methods_lc
    ) %>%
    mutate(
      method = factor(.pretty_method_label(method_raw), levels = method_levels),
      x_num  = as.numeric(as.character(x_index)),
      y_num  = as.numeric(as.character(y_index))
    ) %>%
    filter(y_num %in% keep_y)

  if (!nrow(df)) {
    stop(
      "No rows for methods {", paste(methods, collapse = ", "),
      "} at N = ", sample_size
    )
  }

  heat_df <- .build_heat_df(df, metric)
  # Skew-t averages at extreme tau reach 1e20 and more: colour limits from the other cells, larger ones saturated
  limits  <- range(heat_df$value[abs(heat_df$value) < 1e4], na.rm = TRUE)
  mcol    <- .resolve_metric(metric, names(df))

  fill_scale <- scale_fill_gradient(
    low = "white", high = "black",
    limits = limits,
    name = stringr::str_to_title(mcol),
    oob = scales::squish
  )
  suffix <- "bw"
  base_theme <- theme_bw(base_size = 14)

  label_fill <- if (tolower(metric) == "sd" || metric == "SD") "SD" else stringr::str_to_title(mcol)
  n_methods <- dplyr::n_distinct(heat_df$method)
  plot_width <- max(12, 3.5 * n_methods)
  p <- ggplot(heat_df, aes(x = tau_x, y = tau_y, fill = value)) +
    geom_raster() +
    fill_scale +
    facet_wrap(~ method, ncol = n_methods, nrow = 1) +
    labs(
      x = expression("Covariate Quantile Index"),
      y = expression("Outcome Quantile Index:"~tau),
      fill = label_fill
    ) +
    base_theme +
    scale_x_discrete(drop = FALSE) +
    theme(
      strip.text  = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title  = element_text(size = 18, face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_line(colour = "grey85"),
      panel.grid.minor.x = element_blank()
    ) +
    if (tolower(metric) == "sd" || metric == "SD") guides(fill = guide_legend(title = "SD")) else NULL

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  metric_tag <- gsub("\\s+", "", tolower(mcol))
  base_name <- file.path(
    save_dir,
    paste0(file_prefix, "_", scenario, "_N", sample_size, "_", metric_tag, "_", suffix)
  )

  ggsave(paste0(base_name, ".png"),
         plot = p, width = plot_width, height = 5, dpi = 300, bg = "white")
  message("Saved competing-family figure: ", base_name, ".png")

  p
}

# ---- 1) Competing methods: IPW, ExtremeQR, SkewT (1x3 row) ----
plot_competing_row <- function(
    data_all,
    sample_size,
    metric    = "coverage",
    save_dir  = here("Results/Simulation/Figures"),
    scenario  = "spliced"
) {
  .plot_competing_family(
    data_all = data_all,
    sample_size = sample_size,
    metric = metric,
    methods = c("IPW", "ExtremeQR", "SkewT"),
    file_prefix = "Competing",
    save_dir = save_dir,
    scenario = scenario
  )
}

# ---- 2) Three MixGPD models: gamma vs normal vs lognormal (1x3 row) ----
plot_mixgpd_row <- function(
    table_dir = here("Results/Simulation/Tables"),
    sample_size,
    metric   = "coverage",
    save_dir = here("Results/Simulation/Figures"),
    scenario = "spliced"
) {
  f_gamma    <- file.path(table_dir, paste0("MixGPD_", scenario, "_", sample_size, "_gamma.csv"))
  f_normal   <- file.path(table_dir, paste0("MixGPD_", scenario, "_", sample_size, "_normal.csv"))
  f_lognorm  <- file.path(table_dir, paste0("MixGPD_", scenario, "_", sample_size, "_lognormal.csv"))

  if (!file.exists(f_gamma) || !file.exists(f_normal) || !file.exists(f_lognorm)) {
    stop("Missing one or more MixGPD files (gamma/normal/lognormal) for N = ", sample_size)
  }

  gamma_df <- readr::read_csv(f_gamma, show_col_types = FALSE) %>%
    mutate(
      method = "MixGammaGPD",
      x_num  = as.numeric(as.character(x_index)),
      y_num  = as.numeric(as.character(y_index))
    )

  normal_df <- readr::read_csv(f_normal, show_col_types = FALSE) %>%
    mutate(
      method = "MixInvNormGPD",
      x_num  = as.numeric(as.character(x_index)),
      y_num  = as.numeric(as.character(y_index))
    )

  lognorm_df <- readr::read_csv(f_lognorm, show_col_types = FALSE) %>%
    mutate(
      method = "MixLogNormGPD",
      x_num  = as.numeric(as.character(x_index)),
      y_num  = as.numeric(as.character(y_index))
    )

  df <- dplyr::bind_rows(gamma_df, normal_df, lognorm_df)

  heat_df <- .build_heat_df(df, metric)
  limits  <- range(heat_df$value[abs(heat_df$value) < 1e4], na.rm = TRUE)
  mcol    <- .resolve_metric(metric, names(df))

  fill_scale <- scale_fill_gradient(
    low = "white", high = "black",
    limits = limits,
    name   = stringr::str_to_title(mcol),
    oob    = scales::squish
  )
  suffix <- "bw"
  base_theme <- theme_bw(base_size = 14)

  label_fill <- if (tolower(metric) == "sd" || metric == "SD") "SD" else stringr::str_to_title(mcol)
  p <- ggplot(heat_df, aes(x = tau_x, y = tau_y, fill = value)) +
    geom_raster() +
    fill_scale +
    facet_wrap(~ method, ncol = 3, nrow = 1) +
    labs(
       x = expression("Covariate Quantile Index"),
      y = expression("Outcome Quantile Index:"~tau),
      fill = label_fill
    ) +
    base_theme +
    scale_x_discrete(labels = function(v) sprintf("%.1f", as.numeric(v))) +
    theme(
      strip.text  = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",size = 16),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title  = element_text(size = 18, face = "bold"),
      legend.position = "right",
      panel.grid.major.x = element_line(colour = "grey85"),
      panel.grid.minor.x = element_blank()
    ) +
    if (tolower(metric) == "sd" || metric == "SD") guides(fill = guide_legend(title = "SD")) else NULL

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  metric_tag <- gsub("\\s+", "", tolower(mcol))
  base_name <- file.path(
    save_dir,
    paste0("Mixgpd_", scenario, "_N", sample_size, "_", metric_tag, "_", suffix)
  )

  ggsave(paste0(base_name, ".png"),
         plot = p, width = 12, height = 5, dpi = 300, bg = "white")
  message("Saved MixGPD figure: ", base_name, ".png")

  p
}

# ---------- Figures of the Supplementary Material: S1 all four metrics, S2 and S3 coverage and bias ----------
table_dir <- here("Results/Simulation/Tables")
save_dir  <- here("Results/Simulation/Figures")
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
all_tables <- load_all_sim_tables(table_dir)
figures <- list(spliced = c("coverage", "bias", "SD", "length"),
                burr    = c("coverage", "bias"),
                lnmix   = c("coverage", "bias"))
for (scen in names(figures)) {
  for (met in figures[[scen]]) {
    message("=== scenario = ", scen, " | metric = ", met, " ===")
    tryCatch({
      plot_competing_row(data_all = all_tables, sample_size = 500, metric = met, save_dir = save_dir, scenario = scen)
      plot_mixgpd_row(table_dir = table_dir, sample_size = 500, metric = met, save_dir = save_dir, scenario = scen)
    }, error = function(e) message("Skipped (", scen, ", ", met, "): ", e$message))
  }
}
