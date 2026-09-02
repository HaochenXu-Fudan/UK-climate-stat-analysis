# Observation-only EBM trend smoothing for HadUK.
# fwd stays on proportion scale; wsmax is modelled on cube-root scale and plotted in days.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else if (file.exists(file.path(getwd(), "climate_common.R"))) {
  normalizePath(getwd())
} else {
  normalizePath(file.path(getwd(), "R"))
}
project_dir <- normalizePath(file.path(script_dir, ".."))
source(file.path(script_dir, "climate_common.R"), local = TRUE)

library(TimSPEC)
data(SSP585data)

indices_to_run <- c("fwd", "wsmax")
region_to_use <- "London"
year_min <- 1950
year_max <- 2021
init_year_min <- 1950
init_year_max <- 1979

fig_width <- 10
fig_height <- 7
fig_res <- 300

prior_by_index <- list(
  fwd = rbind(
    c(-5.46, 1.18),
    c(-12.61, 1.53),
    c(0.00, 5.00)
  ),
  wsmax = rbind(
    c(-3.69, 1.18),
    c(-9.69, 1.18),
    c(0.00, 5.00)
  )
)
kappa_by_index <- c(fwd = 0.04, wsmax = 0.25)

plot_wsmax_observation_panel <- function(dat, fit, title_text) {
  pred <- dlm.ObsPred(fit)
  mu_model <- fit$Smooth$s[-1, 1]
  se_model <- pred$SE[, 1]
  mu_days <- to_report_scale(mu_model, "wsmax")
  lower <- to_report_scale(mu_model - 1.96 * se_model, "wsmax")
  upper <- to_report_scale(mu_model + 1.96 * se_model, "wsmax")
  obs_days <- to_report_scale(dat$value, "wsmax")

  ylim <- range(c(obs_days, lower, upper), finite = TRUE)
  plot(
    dat$Year,
    obs_days,
    type = "l",
    col = "black",
    xlab = "Year",
    ylab = "days",
    ylim = ylim,
    main = title_text
  )
  polygon(
    c(dat$Year, rev(dat$Year)),
    c(lower, rev(upper)),
    border = NA,
    col = adjustcolor("darkblue", alpha.f = 0.12)
  )
  lines(dat$Year, obs_days, col = "black")
  lines(dat$Year, mu_days, col = "darkblue", lwd = 3)
}

smooth_observation_index <- function(haduk, index) {
  validate_index(index)
  out_dir <- ensure_dir(file.path(output_root, "observation_smooth", index))
  fig_path <- file.path(
    out_dir,
    paste0("HadUK_", safe_file_label(region_to_use), "_", index, "_smooth_5seasons.png")
  )

  png(
    fig_path,
    width = fig_width * fig_res,
    height = fig_height * fig_res,
    res = fig_res,
    bg = "white"
  )
  oldpar <- par(no.readonly = TRUE)
  on.exit({ par(oldpar); dev.off() }, add = TRUE)
  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 1.8, 0))

  fits <- list()
  data_used <- list()

  for (season in climate_seasons) {
    dat <- read_index_series(
      haduk,
      region_to_use,
      season,
      index,
      model_scale = TRUE
    )
    dat <- dat[
      dat$Year >= year_min & dat$Year <= year_max &
        is.finite(dat$Year) & is.finite(dat$value),
      ,
      drop = FALSE
    ]
    if (nrow(dat) == 0L) stop("No observation data for ", index, " / ", season)

    Xt <- align_erf(dat$Year, ERF585)
    mu0 <- mean(
      dat$value[dat$Year >= init_year_min & dat$Year <= init_year_max],
      na.rm = TRUE
    )
    m0 <- c(mu0, 0, 1)
    fit <- EBMtrendSmooth(
      dat$value,
      Xt = Xt,
      m0 = m0,
      kappa = unname(kappa_by_index[[index]]),
      UsePhi = TRUE,
      prior.pars = prior_by_index[[index]]
    )

    fits[[season]] <- fit
    data_used[[season]] <- dat
    if (index == "wsmax") {
      plot_wsmax_observation_panel(dat, fit, season_labels[[season]])
    } else {
      plot_data <- data.frame(Year = dat$Year, HadUK = dat$value)
      SmoothPlot(
        plot_data,
        fit,
        DatColours = c("black", NA),
        PredColours = c("darkblue", NA),
        Units = "proportion",
        plot.title = season_labels[[season]]
      )
    }
  }

  plot.new()
  mtext(
    paste0("HadUK observation-only EBM trend - ", region_to_use, " (", index, ")"),
    outer = TRUE,
    cex = 1.05,
    font = 2
  )

  result <- list(
    index = index,
    region = region_to_use,
    data = data_used,
    fit = fits,
    priors = prior_by_index[[index]],
    kappa = unname(kappa_by_index[[index]]),
    figure = fig_path
  )
  saveRDS(result, file.path(out_dir, paste0("HadUK_", region_to_use, "_", index, "_smooth.rds")))
  invisible(result)
}

run_observation_smoothing <- function() {
  haduk <- read_haduk(prefer_filtered = TRUE)
  results <- setNames(vector("list", length(indices_to_run)), indices_to_run)
  for (index in indices_to_run) {
    results[[index]] <- smooth_observation_index(haduk, index)
  }
  invisible(results)
}

if (sys.nframe() == 0L) run_observation_smoothing()
