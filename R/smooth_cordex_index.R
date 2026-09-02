################################################################################
# Internal single-index EuroCORDEX implementation.
# Run R/smooth_cordex.R rather than executing this file directly.
# EBM-inspired two-way ensemble smoother + posterior predictive sampling
################################################################################

if (!exists("project_dir", inherits = FALSE)) {
  stop("project_dir must be supplied by R/smooth_cordex.R.")
}
if (!exists("var_name", inherits = FALSE) || !(var_name %in% c("fwd", "wsmax"))) {
  stop("var_name must be supplied as 'fwd' or 'wsmax' by R/smooth_cordex.R.")
}

library(TimSPEC)
data(SSP585data)  # should load ERF585
source(file.path(project_dir, "R", "smooth_helpers.R"), local = TRUE)

################################################################################
# 0. User settings
################################################################################

data_dir <- file.path(project_dir, "PrecipData")

region_to_use <- "London"
year_min <- 1950
year_max <- 2080
init_year_min <- 1950
init_year_max <- 1979

# For final dissertation runs, use at least 1000. For debugging, use 250.
N_pps <- 1000
random_seed <- 2000

out_root <- file.path(project_dir, "outputs", "CORDEX_smooth")
out_dir <- file.path(out_root, paste0("CORDEX_", var_name))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fig_width <- 12.5
fig_height <- 7.5
fig_res <- 300

show_outer_title <- TRUE

################################################################################
# 1. Input files
################################################################################

haduk_candidates <- c(
  file.path(data_dir, "HadUK_1950_2021.csv"),
  file.path(data_dir, "HadUK.csv")
)
haduk_file <- haduk_candidates[file.exists(haduk_candidates)][1]
if (is.na(haduk_file)) {
  stop("Cannot find HadUK_1950_2021.csv or HadUK.csv in data_dir.")
}

haduk <- read.csv(haduk_file, na.strings = c("  NA", "NA", "", "NaN"))

csv_cordex <- list.files(
  data_dir,
  pattern = "^CORDEX_.*\\.csv$",
  full.names = TRUE
)

cordex_subdir <- file.path(data_dir, "CORDEX-data")
if (length(csv_cordex) == 0 && dir.exists(cordex_subdir)) {
  csv_cordex <- list.files(
    cordex_subdir,
    pattern = "^CORDEX_.*\\.csv$",
    full.names = TRUE
  )
}

csv_cordex <- sort(csv_cordex)
if (length(csv_cordex) == 0) {
  stop("No CORDEX_*.csv files found in data_dir or data_dir/CORDEX-data.")
}

parse_cordex_file <- function(file_name) {
  stem <- tools::file_path_sans_ext(basename(file_name))
  m <- regexec("^CORDEX_(.+)_([^_]+)_(r[0-9]+i[0-9]+p[0-9]+)$", stem)
  parsed <- regmatches(stem, m)[[1]]

  if (length(parsed) != 4) {
    stop(
      "Unexpected CORDEX file name: ", basename(file_name),
      "\nExpected: CORDEX_<GCM>_<RCM>_<run>.csv"
    )
  }

  gcm <- parsed[2]
  rcm <- parsed[3]
  run <- parsed[4]

  data.frame(
    file_name = file_name,
    gcm = gcm,
    rcm = rcm,
    run = run,
    member_name = paste("CORDEX", gcm, rcm, run, sep = "_"),
    stringsAsFactors = FALSE
  )
}

cordex_info <- do.call(rbind, lapply(csv_cordex, parse_cordex_file))

if (any(duplicated(cordex_info$member_name))) {
  stop("Duplicated CORDEX member names found after parsing file names.")
}

GCMIDs <- cordex_info$gcm
RCMIDs <- cordex_info$rcm
GCMgroup <- as.numeric(as.factor(GCMIDs))
RCMgroup <- as.numeric(as.factor(RCMIDs))

Groups <- matrix(c(RCMgroup, GCMgroup), ncol = 2)
colnames(Groups) <- c("RCM", "GCM")

cat("HadUK file:", haduk_file, "\n")
cat("Number of CORDEX files found:", length(csv_cordex), "\n")
cat("GCM counts:\n")
print(table(cordex_info$gcm))
cat("RCM counts:\n")
print(table(cordex_info$rcm))

################################################################################
# 2. Seasons, priors and initialisation
################################################################################

seasons <- c("Annual", "DJF", "MAM", "JJA", "SON")

season_label <- c(
  Annual = "Annual",
  DJF    = "winter",
  MAM    = "spring",
  JJA    = "summer",
  SON    = "autumn"
)

panel_title <- c(
  Annual = "Annual",
  DJF    = "Winter",
  MAM    = "Spring",
  JJA    = "Summer",
  SON    = "Autumn"
)

season_keys <- unname(season_label[seasons])
make_season_list <- function() {
  setNames(vector("list", length(season_keys)), season_keys)
}

# TimSPEC EnsEBM2waytrend.modeldef() order is:
#   alpha,
#   log(sigsq[0]),
#   log(tausq[0]),
#   log(tausq[w]),
#   log(sigsq[1]),
#   log(tausq[1]),
#   logit(phi[0]),
#   logit(phi[1]).
priors <- rbind(
  c( 1.00, 1.00),
  c(-3.70, 3.45),
  c(-8.30, 3.45),
  c(-8.30, 3.45),
  c(-3.70, 3.45),
  c(-8.30, 3.45),
  c( 0.00, 5.00),
  c( 0.00, 5.00)
)

# fwd keeps the original diffuse value used in the draft CORDEX fwd analysis.
# wsmax is fitted on cube-root scale, so use the smaller initial variance described
# in Chapter 5: initial SD 0.5 -> variance 0.25.
kappa <- if (var_name == "wsmax") 0.25 else 25
plot_units <- if (var_name == "wsmax") "days" else "proportion"

################################################################################
# 3. Utility functions
################################################################################

fill_partial_ensemble_na <- function(df, ens_cols) {
  if (length(ens_cols) == 0) {
    return(df)
  }

  for (col in ens_cols) {
    df[[col]] <- to_num(df[[col]])
  }

  ens_mat <- as.matrix(df[, ens_cols, drop = FALSE])
  for (i in seq_len(nrow(ens_mat))) {
    na_idx <- is.na(ens_mat[i, ])
    if (any(na_idx) && !all(na_idx)) {
      ens_mat[i, na_idx] <- mean(ens_mat[i, !na_idx], na.rm = TRUE)
    }
  }

  df[, ens_cols] <- as.data.frame(ens_mat)
  df
}

################################################################################
# 4. Data preparation
################################################################################

build_merged_season_data <- function(region, season) {
  obs <- haduk[
    haduk$Region == region & haduk$Season == season,
    c("Year", var_name)
  ]

  if (nrow(obs) == 0) {
    stop("No HadUK rows for region=", region, ", season=", season, ", variable=", var_name)
  }

  obs$Year <- to_num(obs$Year)
  obs[[2]] <- to_num(obs[[2]])

  if (var_name == "wsmax") {
    obs[[2]] <- pmax(obs[[2]], 0)^(1 / 3)
  }

  colnames(obs) <- c("Year", "HadUK")
  obs <- obs[obs$Year >= year_min & obs$Year <= year_max, ]
  obs <- obs[order(obs$Year), ]

  merged <- obs

  for (j in seq_len(nrow(cordex_info))) {
    file_name <- cordex_info$file_name[j]
    data_ensemble <- read.csv(file_name, na.strings = c("  NA", "NA", "", "NaN"))

    required_cols <- c("Year", "Season", "Region", var_name)
    missing_cols <- setdiff(required_cols, names(data_ensemble))
    if (length(missing_cols) > 0) {
      stop(
        "Missing required columns in ", basename(file_name), ": ",
        paste(missing_cols, collapse = ", ")
      )
    }

    sub <- data_ensemble[
      data_ensemble$Region == region & data_ensemble$Season == season,
      c("Year", var_name)
    ]

    sub$Year <- to_num(sub$Year)
    sub[[2]] <- to_num(sub[[2]])

    if (var_name == "wsmax") {
      sub[[2]] <- pmax(sub[[2]], 0)^(1 / 3)
    }

    colnames(sub) <- c("Year", cordex_info$member_name[j])
    sub <- sub[sub$Year >= year_min & sub$Year <= year_max, ]
    sub <- sub[order(sub$Year), ]

    merged <- merge(merged, sub, by = "Year", all = TRUE)
  }

  full_years <- data.frame(Year = year_min:year_max)
  merged <- merge(full_years, merged, by = "Year", all.x = TRUE)
  merged <- merged[order(merged$Year), ]

  for (j in 2:ncol(merged)) {
    merged[[j]] <- to_num(merged[[j]])
  }

  ens_cols <- if (ncol(merged) >= 3) 3:ncol(merged) else integer(0)
  merged <- fill_partial_ensemble_na(merged, ens_cols)
  merged
}

make_m0 <- function(mu0_init) {
  n_state <- 3 * (length(unique(GCMgroup)) + length(unique(RCMgroup)))
  m0 <- rep(0, n_state)
  m0[1] <- mu0_init
  m0
}

################################################################################
# 5. Plotting functions
################################################################################

open_panel_png <- function(fig_path) {
  png(
    filename = fig_path,
    width = fig_width * fig_res,
    height = fig_height * fig_res,
    res = fig_res,
    bg = "white"
  )
  par(
    mfrow = c(2, 3),
    mar = c(4.2, 4.2, 2.2, 1.1),
    oma = c(0, 0, 1.8, 0),
    bg = "white",
    cex.main = 1.0
  )
}

close_panel_png <- function(outer_title) {
  plot.new()
  if (show_outer_title) {
    mtext(outer_title, outer = TRUE, cex = 1.05, font = 2)
  }
  dev.off()
}

plot_model_panel <- function(Data, Smooth = NULL, Samples = NULL, season, use_wsmax_days) {
  this_title <- panel_title[[season]]

  if (use_wsmax_days) {
    plot_wsmax_day_panel(
      Data = Data,
      Smooth = Smooth,
      Samples = Samples,
      Groups = Groups,
      DatColours = c("black", "coral3"),
      PredColours = c("darkblue", "darkgoldenrod"),
      PlotConsensus = TRUE,
      SmoothPIs = !is.null(Samples),
      EnsTransp = 0.2,
      Units = "days",
      plot.title = this_title,
      LegPos = "topleft"
    )
  } else if (is.null(Samples)) {
    SmoothPlot(
      Data,
      Smooth,
      DatColours = c("black", "coral3"),
      Groups = Groups,
      PredColours = c("darkblue", "darkgoldenrod"),
      PlotConsensus = TRUE,
      EnsTransp = 0.2,
      Units = plot_units,
      plot.title = this_title,
      LegPos = "topleft"
    )
  } else {
    SmoothPlot(
      Data,
      Samples = Samples,
      DatColours = c("black", "coral3"),
      Groups = Groups,
      PredColours = c("darkblue", "darkgoldenrod"),
      PlotConsensus = TRUE,
      SmoothPIs = TRUE,
      EnsTransp = 0.2,
      Units = plot_units,
      plot.title = this_title,
      LegPos = "topleft"
    )
  }
}

plot_grid <- function(fig_path, merged_data_list, fit_list = NULL, samples_list = NULL, outer_title) {
  open_panel_png(fig_path)
  for (season in seasons) {
    show_sea <- season_label[[season]]
    plot_model_panel(
      Data = merged_data_list[[show_sea]],
      Smooth = if (is.null(fit_list)) NULL else fit_list[[show_sea]],
      Samples = if (is.null(samples_list)) NULL else samples_list[[show_sea]],
      season = season,
      use_wsmax_days = (var_name == "wsmax")
    )
  }
  close_panel_png(outer_title)
  invisible(fig_path)
}

plot_cumulative_weights <- function(w, main) {
  w <- as.numeric(w)
  w <- w[is.finite(w)]

  if (length(w) == 0 || sum(w) <= 0) {
    plot.new()
    title(main = main)
    text(0.5, 0.5, "No valid weights")
    return(invisible(NULL))
  }

  w <- sort(w, decreasing = TRUE)
  x <- 100 * seq_along(w) / length(w)
  y <- 100 * cumsum(w) / sum(w)

  plot(
    x, y,
    type = "l",
    lwd = 2,
    col = "darkblue",
    xlim = c(0, 100),
    ylim = c(0, 100),
    xlab = "Cumulative % of samples",
    ylab = "Cumulative % of total weight",
    main = main
  )
  grid()
  invisible(NULL)
}

################################################################################
# 6. Fitting, PPS and outputs
################################################################################

make_2080_summary <- function(dat, pps, season_code, season_name, mu0_init) {
  row_id <- which(dat$Year == 2080)
  if (length(row_id) != 1) {
    stop("Year 2080 is not uniquely present in the data.")
  }

  obs_samples <- extract_obs_matrix(pps)[row_id, ]
  w <- get_weights(pps)
  obs_model_scale <- obs_samples

  if (var_name == "wsmax") {
    obs_report_scale <- wsmax_to_days(obs_samples)
    invalid_fraction <- mean(obs_model_scale < 0, na.rm = TRUE)
    report_units <- "days"
  } else {
    obs_report_scale <- obs_samples
    invalid_fraction <- mean(obs_model_scale < 0 | obs_model_scale > 1, na.rm = TRUE)
    report_units <- "fraction"
  }

  ok <- is.finite(obs_report_scale) & is.finite(w)
  x <- obs_report_scale[ok]
  ww <- w[ok]
  ww <- ww / sum(ww)

  qs <- weighted_quantile(x, ww, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  mu <- sum(ww * x)
  sdev <- sqrt(sum(ww * (x - mu)^2))
  ess <- 1 / sum(ww^2)

  ens_2080 <- as.numeric(dat[row_id, -(1:2)])
  if (var_name == "wsmax") ens_2080 <- wsmax_to_days(ens_2080)

  data.frame(
    index = var_name,
    season_code = season_code,
    season_name = season_name,
    units = report_units,
    n_ens = nrow(Groups),
    mu0_init = mu0_init,
    kappa = kappa,
    raw_ensemble_mean_2080 = mean(ens_2080, na.rm = TRUE),
    raw_ensemble_sd_2080 = sd(ens_2080, na.rm = TRUE),
    posterior_mean_2080 = mu,
    posterior_sd_2080 = sdev,
    posterior_p025_2080 = qs[1],
    posterior_q1_2080 = qs[2],
    posterior_median_2080 = qs[3],
    posterior_q3_2080 = qs[4],
    posterior_p975_2080 = qs[5],
    importance_ess = ess,
    invalid_fraction_model_scale = invalid_fraction
  )
}

merged_data_list <- make_season_list()
EnsEBMFit_list <- make_season_list()
PPSApprox_list <- make_season_list()
PPS_list <- make_season_list()
m0_list <- make_season_list()
Xt_list <- make_season_list()
mu0_table <- data.frame()
theta_table <- data.frame()
summary_table <- data.frame()

cat("\n============================================================\n")
cat("Running EuroCORDEX postprocessing for index:", var_name, "\n")
cat("Region:", region_to_use, "\n")
cat("Output folder:", out_dir, "\n")
cat("============================================================\n")

for (season in seasons) {
  show_sea <- season_label[[season]]
  cat("\n--- Fitting season:", season, "(", show_sea, ") ---\n")

  dat <- build_merged_season_data(region_to_use, season)
  Xt_use <- get_erf_vector(dat$Year)

  mu0_init <- mean(
    dat$HadUK[dat$Year >= init_year_min & dat$Year <= init_year_max],
    na.rm = TRUE
  )

  if (!is.finite(mu0_init)) {
    stop("mu0_init is not finite for ", var_name, " / ", season)
  }

  m0 <- make_m0(mu0_init)

  fit <- EnsEBM2waytrendSmooth(
    as.matrix(dat[, -1]),
    Xt = Xt_use,
    Groups = Groups,
    m0 = m0,
    kappa = kappa,
    prior.pars = priors,
    interactions = "none",
    messages = TRUE
  )

  merged_data_list[[show_sea]] <- dat
  EnsEBMFit_list[[show_sea]] <- fit
  m0_list[[show_sea]] <- m0
  Xt_list[[show_sea]] <- Xt_use

  mu0_table <- rbind(
    mu0_table,
    data.frame(
      index = var_name,
      season_code = season,
      season_name = show_sea,
      n_ens = nrow(Groups),
      mu0_init = mu0_init,
      kappa = kappa
    )
  )

  theta_est <- fit$Theta$par
  theta_table <- rbind(
    theta_table,
    data.frame(
      index = var_name,
      season_code = season,
      season_name = show_sea,
      parameter = names(theta_est),
      MAP = as.numeric(theta_est)
    )
  )
}

# MAP smoother figure.
if (var_name == "wsmax") {
  map_name <- paste0("CORDEX_", gsub(" ", "", region_to_use), "_", var_name, "_Obs_MAP_smooth_days.png")
} else {
  map_name <- paste0("CORDEX_", gsub(" ", "", region_to_use), "_", var_name, "_Obs_MAP_smooth.png")
}
map_fig <- file.path(out_dir, map_name)
plot_grid(
  fig_path = map_fig,
  merged_data_list = merged_data_list,
  fit_list = EnsEBMFit_list,
  outer_title = paste0("CORDEX ", var_name, ": MAP smoother")
)
cat("Saved MAP smoother figure:", map_fig, "\n")

# Posterior predictive sampling.
set.seed(random_seed)

for (season in seasons) {
  show_sea <- season_label[[season]]
  cat("\n--- PPS for season:", season, "(", show_sea, ") ---\n")

  PPSApprox_list[[show_sea]] <- PostPredSample(
    ModelBundle = EnsEBMFit_list[[show_sea]],
    Build = EnsEBM2waytrend.modeldef,
    N = N_pps,
    m0 = m0_list[[show_sea]],
    kappa = kappa,
    Xt = Xt_list[[show_sea]],
    WhichEls = 1,
    ReplaceAll = FALSE,
    Importance = FALSE,
    NonNeg = (var_name == "wsmax"),
    messages = FALSE
  )

  PPS_list[[show_sea]] <- PostPredSample(
    ModelBundle = EnsEBMFit_list[[show_sea]],
    Importance = TRUE,
    CheckMax = TRUE,
    ReplaceOnFail = TRUE,
    Build = EnsEBM2waytrend.modeldef,
    N = N_pps,
    m0 = m0_list[[show_sea]],
    kappa = kappa,
    Xt = Xt_list[[show_sea]],
    WhichEls = 1,
    ReplaceAll = FALSE,
    NonNeg = (var_name == "wsmax"),
    messages = FALSE
  )

  summary_table <- rbind(
    summary_table,
    make_2080_summary(
      dat = merged_data_list[[show_sea]],
      pps = PPS_list[[show_sea]],
      season_code = season,
      season_name = show_sea,
      mu0_init = m0_list[[show_sea]][1]
    )
  )
}

# PPS Laplace approximation figure.
if (var_name == "wsmax") {
  pps_approx_name <- paste0("CORDEX_", gsub(" ", "", region_to_use), "_", var_name, "_smooth_PPSapprox_days.png")
} else {
  pps_approx_name <- paste0("CORDEX_", gsub(" ", "", region_to_use), "_", var_name, "_smooth_PPSapprox.png")
}
pps_approx_fig <- file.path(out_dir, pps_approx_name)
plot_grid(
  fig_path = pps_approx_fig,
  merged_data_list = merged_data_list,
  samples_list = PPSApprox_list,
  outer_title = paste0("CORDEX ", var_name, ": Laplace approximation")
)
cat("Saved PPS approximation figure:", pps_approx_fig, "\n")

# PPS importance sampling figure.
if (var_name == "wsmax") {
  pps_name <- paste0("CORDEX_", gsub(" ", "", region_to_use), "_", var_name, "_smooth_PPS_days.png")
} else {
  pps_name <- paste0("CORDEX_", gsub(" ", "", region_to_use), "_", var_name, "_smooth_PPS.png")
}
pps_fig <- file.path(out_dir, pps_name)
plot_grid(
  fig_path = pps_fig,
  merged_data_list = merged_data_list,
  samples_list = PPS_list,
  outer_title = paste0("CORDEX ", var_name, ": importance sampling")
)
cat("Saved PPS importance figure:", pps_fig, "\n")

# Pairwise importance-weight plots.
weight_figs <- list()
for (season in seasons) {
  show_sea <- season_label[[season]]
  fig_name <- paste0("CORDEX_", show_sea, gsub(" ", "", region_to_use), "_", var_name, "_weight.png")
  fig_path <- file.path(out_dir, fig_name)

  png(filename = fig_path, width = 2000, height = 2000, res = 300, bg = "white")
  par(bg = "white")
  if (!is.null(PPS_list[[show_sea]]$Thetas) && !is.null(PPS_list[[show_sea]]$Weights)) {
    PlotImportanceWts(PPS_list[[show_sea]]$Thetas, PPS_list[[show_sea]]$Weights)
  } else {
    plot.new()
    title(main = paste0(panel_title[[season]], ": no weights"))
  }
  dev.off()

  weight_figs[[show_sea]] <- fig_path
  cat("Saved weight plot:", fig_path, "\n")
}

# Cumulative importance-weight figure.
cumulative_fig <- file.path(
  out_dir,
  paste0("CORDEX_", gsub(" ", "", region_to_use), "_", var_name, "_cumulative.png")
)
open_panel_png(cumulative_fig)
for (season in seasons) {
  show_sea <- season_label[[season]]
  if (!is.null(PPS_list[[show_sea]]$Weights) && !is.null(PPS_list[[show_sea]]$Weights$w)) {
    plot_cumulative_weights(PPS_list[[show_sea]]$Weights$w, main = panel_title[[season]])
  } else {
    plot.new()
    title(main = paste0(panel_title[[season]], ": no weights"))
  }
}
close_panel_png(paste0("CORDEX ", var_name, ": cumulative weights"))
cat("Saved cumulative weight figure:", cumulative_fig, "\n")

# Tables and RDS.
write.csv(
  mu0_table,
  file.path(out_dir, paste0("CORDEX_London_", var_name, "_mu0_initial_values.csv")),
  row.names = FALSE
)

write.csv(
  theta_table,
  file.path(out_dir, paste0("CORDEX_London_", var_name, "_MAP_theta.csv")),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(out_dir, paste0("CORDEX_London_", var_name, "_posterior_2080_summary.csv")),
  row.names = FALSE
)

result <- list(
  index = var_name,
  region = region_to_use,
  merged = merged_data_list,
  fit = EnsEBMFit_list,
  pps_approx = PPSApprox_list,
  pps = PPS_list,
  m0 = m0_list,
  Xt = Xt_list,
  priors = priors,
  kappa = kappa,
  groups = Groups,
  cordex_info = cordex_info,
  mu0_table = mu0_table,
  theta_table = theta_table,
  summary_2080 = summary_table,
  figures = list(
    map = map_fig,
    pps_approx = pps_approx_fig,
    pps_importance = pps_fig,
    weights = weight_figs,
    cumulative_weights = cumulative_fig
  )
)

saveRDS(
  result,
  file.path(out_dir, paste0("CORDEX_London_", var_name, "_results.rds"))
)

cat("\nFinished EuroCORDEX index:", var_name, "\n")
print(summary_table)
cat("Main output folder:", out_dir, "\n")
