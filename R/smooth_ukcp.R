################################################################################
# UKCP18 postprocessing for London fwd and wsmax.
# EBM-inspired ensemble smoother + posterior predictive sampling
#
# Main corrections relative to the draft scripts:
#   1. UKCP18 is treated as an unstructured/exchangeable ensemble: Groups = NULL.
#   2. No artificial UKCP18_A / UKCP18_B grouping is used.
#   3. fwd and cube-root(wsmax) use different priors and kappa values.
#   4. m0[1] is set from the early HadUK mean, season by season.
#   5. wsmax is modelled on cube-root scale, so m0[1] is also on cube-root scale.
#   6. ERF is matched explicitly to the Year column of the merged data.
#   7. PPS uses the season-specific m0, kappa, Xt and actual number of UKCP members.
################################################################################

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else if (basename(getwd()) == "R") {
  normalizePath(getwd())
} else {
  normalizePath(file.path(getwd(), "R"))
}
project_dir <- normalizePath(file.path(script_dir, ".."))

library(TimSPEC)
data(SSP585data)  # should load ERF585
source(file.path(script_dir, "smooth_helpers.R"), local = TRUE)

################################################################################
# 0. User settings
################################################################################

# All inputs and outputs are kept inside the project.
data_dir <- file.path(project_dir, "PrecipData")

# The script will use HadUK.csv if present; otherwise HadUK_1950_2021.csv.
haduk_candidates <- c(
  file.path(data_dir, "HadUK.csv"),
  file.path(data_dir, "HadUK_1950_2021.csv")
)
haduk_file <- haduk_candidates[file.exists(haduk_candidates)][1]
if (is.na(haduk_file)) {
  stop("Cannot find HadUK.csv or HadUK_1950_2021.csv in data_dir.")
}

# fwd stays on proportion scale; wsmax is fitted on cube-root scale and
# converted back to days in plots and summaries.
indices_to_run <- c("fwd", "wsmax")

region_to_use <- "London"
year_min <- 1950
year_max <- 2080
init_year_min <- 1950
init_year_max <- 1979

# For final dissertation runs, use at least 1000. For debugging, use 250.
N_pps <- 1000
random_seed <- 2000

out_root <- file.path(project_dir, "outputs", "UKCP_smooth")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

fig_width <- 10
fig_height <- 7
fig_res <- 300

################################################################################
# 1. Fixed season mapping
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

show_outer_title <- TRUE

################################################################################
# 2. Priors and initialisation variances
################################################################################

# TimSPEC EnsEBMtrend.modeldef() order is:
#   alpha,
#   log(sigsq[0]),
#   log(tausq[0]),
#   log(tausq[w]),
#   log(sigsq[1]),
#   log(tausq[1]),
#   logit(phi[0]),
#   logit(phi[1]).

priors_fwd <- rbind(
  c(  1.00, 0.25),  # alpha
  c( -5.46, 1.18),  # log(sigsq[0])
  c(-12.61, 1.53),  # log(tausq[0])
  c(-12.61, 1.53),  # log(tausq[w])
  c( -5.46, 1.18),  # log(sigsq[1])
  c(-12.61, 1.53),  # log(tausq[1])
  c(  0.00, 5.00),  # logit(phi[0])
  c(  0.00, 5.00)   # logit(phi[1])
)

priors_wsmax <- rbind(
  c(  1.00, 0.25),  # alpha
  c( -3.69, 1.18),  # log(sigsq[0])
  c( -9.69, 1.18),  # log(tausq[0])
  c( -9.69, 1.18),  # log(tausq[w])
  c( -3.69, 1.18),  # log(sigsq[1])
  c( -9.69, 1.18),  # log(tausq[1])
  c(  0.00, 5.00),  # logit(phi[0])
  c(  0.00, 5.00)   # logit(phi[1])
)

prior_by_index <- list(
  fwd = priors_fwd,
  wsmax = priors_wsmax
)

# kappa is the diffuse initialisation variance used by TimSPEC.
# fwd: initial SD 0.20 -> variance 0.04.
# cube-root(wsmax): initial SD 0.50 -> variance 0.25.
kappa_by_index <- list(
  fwd = 0.04,
  wsmax = 0.25
)

units_by_index <- list(
  fwd = "proportion",
  wsmax = "days"
)

################################################################################
# 3. Utility functions
################################################################################

get_member_id <- function(path) {
  sub("^UKCP18_([0-9]+)\\.csv$", "\\1", basename(path))
}

get_ukcp_files <- function(data_dir) {
  files <- list.files(
    data_dir,
    pattern = "^UKCP18_[0-9]+\\.csv$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    stop("No UKCP18_*.csv files found in data_dir.")
  }
  ids <- get_member_id(files)
  files[order(as.integer(ids))]
}

check_partial_ensemble_missing <- function(Ymat, years) {
  if (ncol(Ymat) <= 2) return(invisible(TRUE))

  ens_na <- is.na(Ymat[, -1, drop = FALSE])
  partial <- apply(ens_na, 1, function(z) any(z) && !all(z))

  if (any(partial)) {
    bad <- years[partial]
    stop(
      "Some rows have partially missing UKCP18 ensemble members. ",
      "Do not fill them manually with row means. Inspect these years first: ",
      paste(bad, collapse = ", ")
    )
  }

  invisible(TRUE)
}


################################################################################
# 4. Data preparation
################################################################################

prepare_merged_data <- function(haduk, ukcp_files, region, season, var_name,
                                transform_wsmax = FALSE,
                                year_min = 1950, year_max = 2080) {
  obs <- haduk[
    haduk$Region == region & haduk$Season == season,
    c("Year", var_name)
  ]

  if (nrow(obs) == 0) {
    stop("No HadUK rows for region=", region, ", season=", season, ", variable=", var_name)
  }

  obs$Year <- to_num(obs$Year)
  obs[[var_name]] <- to_num(obs[[var_name]])

  if (transform_wsmax) {
    obs[[var_name]] <- obs[[var_name]]^(1 / 3)
  }

  colnames(obs) <- c("Year", "HadUK")
  obs <- obs[obs$Year >= year_min & obs$Year <= year_max, ]
  obs <- obs[order(obs$Year), ]

  merged <- obs

  for (path in ukcp_files) {
    id <- get_member_id(path)
    ens <- read.csv(path, na.strings = c("  NA", "NA", "", "NaN"))

    sub <- ens[
      ens$Region == region & ens$Season == season,
      c("Year", var_name)
    ]

    if (nrow(sub) == 0) {
      stop("No UKCP rows in ", basename(path), " for region=", region, ", season=", season)
    }

    sub$Year <- to_num(sub$Year)
    sub[[var_name]] <- to_num(sub[[var_name]])

    if (transform_wsmax) {
      sub[[var_name]] <- sub[[var_name]]^(1 / 3)
    }

    colnames(sub) <- c("Year", paste0("UKCP18_", id))
    sub <- sub[sub$Year >= year_min & sub$Year <= year_max, ]
    sub <- sub[order(sub$Year), ]

    merged <- merge(merged, sub, by = "Year", all = TRUE)
  }

  full_years <- data.frame(Year = year_min:year_max)
  merged <- merge(full_years, merged, by = "Year", all.x = TRUE)
  merged <- merged[order(merged$Year), ]

  # Force all data columns to numeric after merging.
  for (j in 2:ncol(merged)) {
    merged[[j]] <- to_num(merged[[j]])
  }

  merged
}

################################################################################
# 5. Fitting, PPS, plotting and summaries
################################################################################

make_m0 <- function(n_ens, mu0_init) {
  m0 <- rep(0, 3 * (n_ens + 1))
  m0[1] <- mu0_init
  m0
}

make_2080_summary <- function(dat, pps, var_name, season_code, season_name,
                              mu0_init, n_ens, kappa) {
  row_id <- which(dat$Year == 2080)
  if (length(row_id) != 1) {
    stop("Year 2080 is not uniquely present in the data.")
  }

  obs_samples <- extract_obs_matrix(pps)[row_id, ]
  w <- get_weights(pps)

  # Keep a copy on modelling scale for diagnostics.
  obs_model_scale <- obs_samples

  if (var_name == "wsmax") {
    # Convert cube-root scale posterior predictive samples back to days.
    # Use pmax() only for the final reporting scale, not for model fitting.
    obs_report_scale <- pmax(obs_samples, 0)^3
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
  if (var_name == "wsmax") ens_2080 <- pmax(ens_2080, 0)^3

  data.frame(
    index = var_name,
    season_code = season_code,
    season_name = season_name,
    units = report_units,
    n_ens = n_ens,
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

plot_map_grid <- function(merged_list, fit_list, var_name, out_dir) {
  suffix <- if (var_name == "wsmax") "MAP_smooth_days" else "MAP_smooth"
  fig_path <- file.path(out_dir, paste0("UKCP18_London_", var_name, "_", suffix, ".png"))

  png(fig_path, width = fig_width * fig_res, height = fig_height * fig_res, res = fig_res)
  oldpar <- par(no.readonly = TRUE)
  on.exit({ par(oldpar); dev.off() }, add = TRUE)

  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 1.6, 0))

  for (season in seasons) {
    show_sea <- season_label[[season]]
    this_title <- panel_title[[season]]

    if (var_name == "wsmax") {
      plot_wsmax_day_panel(
        Data = merged_list[[show_sea]],
        Smooth = fit_list[[show_sea]],
        Groups = NULL,
        DatColours = c("black", "coral3"),
        PredColours = c("darkblue", "darkgoldenrod"),
        PlotConsensus = TRUE,
        EnsTransp = 0.20,
        Units = "days",
        plot.title = this_title,
        LegPos = "topleft"
      )
    } else {
      SmoothPlot(
        Data = merged_list[[show_sea]],
        Smooth = fit_list[[show_sea]],
        DatColours = c("black", "coral3"),
        PredColours = c("darkblue", "darkgoldenrod"),
        PlotConsensus = TRUE,
        EnsTransp = 0.20,
        Units = units_by_index[[var_name]],
        plot.title = this_title,
        LegPos = "topleft"
      )
    }
  }

  plot.new()
  if (show_outer_title) {
    mtext(
      paste0("UKCP18 ", var_name, ": MAP smoother"),
      outer = TRUE,
      cex = 1.05,
      font = 2
    )
  }

  invisible(fig_path)
}

plot_pps_grid <- function(merged_list, samples_list, var_name, out_dir, suffix, title_suffix) {
  suffix_use <- if (var_name == "wsmax") paste0(suffix, "_days") else suffix
  fig_path <- file.path(out_dir, paste0("UKCP18_London_", var_name, "_", suffix_use, ".png"))

  png(fig_path, width = fig_width * fig_res, height = fig_height * fig_res, res = fig_res)
  oldpar <- par(no.readonly = TRUE)
  on.exit({ par(oldpar); dev.off() }, add = TRUE)

  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 1.6, 0))

  for (season in seasons) {
    show_sea <- season_label[[season]]
    this_title <- panel_title[[season]]

    if (var_name == "wsmax") {
      plot_wsmax_day_panel(
        Data = merged_list[[show_sea]],
        Samples = samples_list[[show_sea]],
        Groups = NULL,
        DatColours = c("black", "coral3"),
        PredColours = c("darkblue", "darkgoldenrod"),
        PlotConsensus = TRUE,
        SmoothPIs = TRUE,
        EnsTransp = 0.20,
        Units = "days",
        plot.title = this_title,
        LegPos = "topleft"
      )
    } else {
      SmoothPlot(
        Data = merged_list[[show_sea]],
        Samples = samples_list[[show_sea]],
        DatColours = c("black", "coral3"),
        PredColours = c("darkblue", "darkgoldenrod"),
        PlotConsensus = TRUE,
        SmoothPIs = TRUE,
        EnsTransp = 0.20,
        Units = units_by_index[[var_name]],
        plot.title = this_title,
        LegPos = "topleft"
      )
    }
  }

  plot.new()
  if (show_outer_title) {
    mtext(
      paste0("UKCP18 ", var_name, ": ", title_suffix),
      outer = TRUE,
      cex = 1.05,
      font = 2
    )
  }

  invisible(fig_path)
}

save_importance_weight_plots <- function(samples_list, var_name, out_dir, region = "London") {
  fig_paths <- list()

  for (season in seasons) {
    show_sea <- season_label[[season]]
    pps <- samples_list[[show_sea]]

    fig_path <- file.path(
      out_dir,
      paste0("UKCP18_", show_sea, region, "_", var_name, "_weight.png")
    )

    if (!is.null(pps$Thetas) && !is.null(pps$Weights)) {
      png(fig_path, width = 2000, height = 2000, res = 300)
      PlotImportanceWts(pps$Thetas, pps$Weights)
      dev.off()
    } else {
      png(fig_path, width = 2000, height = 2000, res = 300)
      plot.new()
      title(main = paste0(panel_title[[season]], ": no weights"))
      dev.off()
    }

    fig_paths[[show_sea]] <- fig_path
  }

  invisible(fig_paths)
}

plot_cumulative_weight_grid <- function(samples_list, var_name, out_dir) {
  fig_path <- file.path(out_dir, paste0("UKCP18_London_", var_name, "_cumulative_weights.png"))

  png(fig_path, width = fig_width * fig_res, height = fig_height * fig_res, res = fig_res)
  oldpar <- par(no.readonly = TRUE)
  on.exit({ par(oldpar); dev.off() }, add = TRUE)

  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  for (season in seasons) {
    show_sea <- season_label[[season]]
    pps <- samples_list[[show_sea]]

    if (!is.null(pps$Weights) && !is.null(pps$Weights$w)) {
      CumulativeWeightPlot(pps$Weights$w, main = panel_title[[season]])
      # title(main = paste0(var_name, " - ", show_sea))
    } else {
      plot.new()
      title(main = paste0(panel_title[[season]], ": no weights"))
    }
  }

  plot.new()
  mtext(
    paste0("UKCP18 ", var_name, ": cumulative weights"),
    outer = TRUE,
    cex = 1.1
  )

  invisible(fig_path)
}

run_ukcp_index <- function(var_name, haduk, ukcp_files, region = "London") {
  if (!(var_name %in% c("fwd", "wsmax"))) {
    stop("var_name must be either 'fwd' or 'wsmax'.")
  }

  transform_wsmax <- (var_name == "wsmax")
  priors <- prior_by_index[[var_name]]
  kappa <- kappa_by_index[[var_name]]

  out_dir <- file.path(out_root, paste0("UKCP18_", var_name))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  merged_list <- list()
  fit_list <- list()
  pps_approx_list <- list()
  pps_list <- list()
  m0_list <- list()
  Xt_list <- list()
  mu0_table <- data.frame()
  theta_table <- data.frame()
  summary_table <- data.frame()

  cat("\n============================================================\n")
  cat("Running UKCP18 postprocessing for index:", var_name, "\n")
  cat("Region:", region, "\n")
  cat("Output folder:", out_dir, "\n")
  cat("============================================================\n")

  for (season in seasons) {
    show_sea <- season_label[[season]]
    cat("\n--- Fitting season:", season, "(", show_sea, ") ---\n")

    dat <- prepare_merged_data(
      haduk = haduk,
      ukcp_files = ukcp_files,
      region = region,
      season = season,
      var_name = var_name,
      transform_wsmax = transform_wsmax,
      year_min = year_min,
      year_max = year_max
    )

    Ymat <- as.matrix(dat[, -1])
    check_partial_ensemble_missing(Ymat, dat$Year)

    n_ens <- ncol(Ymat) - 1
    Xt_use <- get_erf_vector(dat$Year)

    mu0_init <- mean(
      dat$HadUK[dat$Year >= init_year_min & dat$Year <= init_year_max],
      na.rm = TRUE
    )

    if (!is.finite(mu0_init)) {
      stop("mu0_init is not finite for ", var_name, " / ", season)
    }

    m0 <- make_m0(n_ens = n_ens, mu0_init = mu0_init)

    fit <- EnsEBMtrendSmooth(
      Y = Ymat,
      Xt = Xt_use,
      Groups = NULL,
      compact = FALSE,
      m0 = m0,
      kappa = kappa,
      prior.pars = priors,
      messages = TRUE
    )

    merged_list[[show_sea]] <- dat
    fit_list[[show_sea]] <- fit
    m0_list[[show_sea]] <- m0
    Xt_list[[show_sea]] <- Xt_use

    mu0_table <- rbind(
      mu0_table,
      data.frame(
        index = var_name,
        season_code = season,
        season_name = show_sea,
        n_ens = n_ens,
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

  # MAP smoother plots.
  map_fig <- plot_map_grid(merged_list, fit_list, var_name, out_dir)
  cat("Saved MAP smoother figure:", map_fig, "\n")

  # Posterior predictive sampling.
  set.seed(random_seed)

  for (season in seasons) {
    show_sea <- season_label[[season]]
    cat("\n--- PPS for season:", season, "(", show_sea, ") ---\n")

    dat <- merged_list[[show_sea]]
    n_ens <- ncol(dat) - 2

    pps_approx <- PostPredSample(
      ModelBundle = fit_list[[show_sea]],
      Build = EnsEBMtrend.modeldef,
      N = N_pps,
      m0 = m0_list[[show_sea]],
      kappa = kappa,
      Xt = Xt_list[[show_sea]],
      NEnsTS = n_ens,
      WhichEls = 1,
      ReplaceAll = FALSE,
      Importance = FALSE,
      NonNeg = TRUE,
      messages = FALSE
    )

    pps <- PostPredSample(
      ModelBundle = fit_list[[show_sea]],
      Build = EnsEBMtrend.modeldef,
      N = N_pps,
      m0 = m0_list[[show_sea]],
      kappa = kappa,
      Xt = Xt_list[[show_sea]],
      NEnsTS = n_ens,
      WhichEls = 1,
      ReplaceAll = FALSE,
      Importance = TRUE,
      CheckMax = TRUE,
      ReplaceOnFail = TRUE,
      NonNeg = TRUE,
      messages = FALSE
    )

    pps_approx_list[[show_sea]] <- pps_approx
    pps_list[[show_sea]] <- pps

    summary_table <- rbind(
      summary_table,
      make_2080_summary(
        dat = dat,
        pps = pps,
        var_name = var_name,
        season_code = season,
        season_name = show_sea,
        mu0_init = m0_list[[show_sea]][1],
        n_ens = n_ens,
        kappa = kappa
      )
    )
  }

  approx_fig <- plot_pps_grid(
    merged_list = merged_list,
    samples_list = pps_approx_list,
    var_name = var_name,
    out_dir = out_dir,
    suffix = "PPS_Laplace_approx",
    title_suffix = "Laplace approximation"
  )
  cat("Saved PPS approximation figure:", approx_fig, "\n")

  pps_fig <- plot_pps_grid(
    merged_list = merged_list,
    samples_list = pps_list,
    var_name = var_name,
    out_dir = out_dir,
    suffix = "PPS_importance",
    title_suffix = "posterior with importance sampling"
  )
  cat("Saved PPS importance figure:", pps_fig, "\n")

  weight_figs <- save_importance_weight_plots(
    samples_list = pps_list,
    var_name = var_name,
    out_dir = out_dir,
    region = region
  )
  cat("Saved importance weight figures:\n")
  print(unlist(weight_figs))

  cumulative_fig <- plot_cumulative_weight_grid(
    samples_list = pps_list,
    var_name = var_name,
    out_dir = out_dir
  )
  cat("Saved cumulative weight figure:", cumulative_fig, "\n")

  write.csv(
    mu0_table,
    file.path(out_dir, paste0("UKCP18_London_", var_name, "_mu0_initial_values.csv")),
    row.names = FALSE
  )

  write.csv(
    theta_table,
    file.path(out_dir, paste0("UKCP18_London_", var_name, "_MAP_theta.csv")),
    row.names = FALSE
  )

  write.csv(
    summary_table,
    file.path(out_dir, paste0("UKCP18_London_", var_name, "_posterior_2080_summary.csv")),
    row.names = FALSE
  )

  result <- list(
    index = var_name,
    region = region,
    merged = merged_list,
    fit = fit_list,
    pps_approx = pps_approx_list,
    pps = pps_list,
    m0 = m0_list,
    Xt = Xt_list,
    priors = priors,
    kappa = kappa,
    mu0_table = mu0_table,
    theta_table = theta_table,
    summary_2080 = summary_table,
    figures = list(
      map = map_fig,
      pps_approx = approx_fig,
      pps_importance = pps_fig,
      weights = weight_figs,
      cumulative_weights = cumulative_fig
    )
  )

  saveRDS(
    result,
    file.path(out_dir, paste0("UKCP18_London_", var_name, "_results.rds"))
  )

  cat("\nFinished index:", var_name, "\n")
  print(summary_table)

  result
}

################################################################################
# 6. Run analysis
################################################################################

run_ukcp_smoothing <- function() {
  haduk <- read.csv(haduk_file, na.strings = c("  NA", "NA", "", "NaN"))
  ukcp_files <- get_ukcp_files(data_dir)

  cat("HadUK file:", haduk_file, "\n")
  cat("Number of UKCP18 files found:", length(ukcp_files), "\n")
  cat("UKCP18 members:", paste(basename(ukcp_files), collapse = ", "), "\n")

  all_results <- setNames(vector("list", length(indices_to_run)), indices_to_run)
  for (idx in indices_to_run) {
    all_results[[idx]] <- run_ukcp_index(
      var_name = idx,
      haduk = haduk,
      ukcp_files = ukcp_files,
      region = region_to_use
    )
  }

  saveRDS(
    all_results,
    file.path(out_root, "UKCP18_London_fwd_wsmax_all_results.rds")
  )

  cat("\nAll requested UKCP18 postprocessing runs finished.\n")
  cat("Main output folder:", out_root, "\n")
  invisible(all_results)
}

if (sys.nframe() == 0L) run_ukcp_smoothing()
