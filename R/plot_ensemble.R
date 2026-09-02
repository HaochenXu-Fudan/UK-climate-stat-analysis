# Plot HadUK observations, ensemble mean and all ensemble members.
# Change the settings below to select sources, indices or regions.

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

sources_to_run <- c("UKCP18", "CORDEX")
indices_to_run <- c("fwd", "wsmax")
regions_to_run <- c("London")
year_min <- 1950
year_max <- 2080

fig_width <- 10
fig_height <- 7
fig_res <- 300

plot_ensemble_grid <- function(haduk, source, index, region) {
  files <- find_ensemble_files(source)
  out_dir <- ensure_dir(file.path(output_root, "ensemble_raw", source, index))
  fig_path <- file.path(
    out_dir,
    paste0(source, "_", safe_file_label(region), "_", index, "_raw_5seasons.png")
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

  for (season in climate_seasons) {
    dat <- build_ensemble_matrix(
      haduk = haduk,
      files = files,
      source = source,
      region = region,
      season = season,
      index = index,
      model_scale = FALSE,
      year_min = year_min,
      year_max = year_max
    )

    if (nrow(dat) == 0L || all(is.na(dat$HadUK))) {
      plot.new()
      title(main = paste0(season_labels[[season]], ": no data"))
      next
    }

    n_ens <- max(ncol(dat) - 2L, 0L)
    PlotEnsTS(
      dat,
      Colours = c("black", "red", rep("coral3", n_ens)),
      EnsTransp = 0.15,
      Units = index_units(index),
      main = season_labels[[season]]
    )
  }

  plot.new()
  mtext(
    paste0(source, " ensemble vs HadUK - ", region, " (", index, ")"),
    outer = TRUE,
    cex = 1.05,
    font = 2
  )
  cat("Saved:", fig_path, "\n")
  invisible(fig_path)
}

run_ensemble_plots <- function() {
  haduk <- read_haduk(prefer_filtered = TRUE)
  outputs <- list()
  for (source in sources_to_run) {
    for (index in indices_to_run) {
      validate_index(index)
      for (region in regions_to_run) {
        key <- paste(source, index, region, sep = "/")
        outputs[[key]] <- plot_ensemble_grid(haduk, source, index, region)
      }
    }
  }
  invisible(outputs)
}

if (sys.nframe() == 0L) run_ensemble_plots()
