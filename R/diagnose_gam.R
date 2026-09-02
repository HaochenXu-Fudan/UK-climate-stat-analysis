# GAM residual and distribution diagnostics for HadUK, UKCP18 and CORDEX.
# fwd is fitted on its original proportion scale; wsmax is fitted on cube-root scale.

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

library(mgcv)

sources_to_run <- c("HadUK", "UKCP18", "CORDEX")
indices_to_run <- c("fwd", "wsmax")
regions_to_run <- c("London")

fig_width <- 10
fig_height <- 7
fig_res <- 300

save_qq_grid <- function(dat, source, member, region, index, out_dir) {
  member_clean <- safe_file_label(member)
  fig_path <- file.path(
    out_dir,
    paste0(source, "_", member_clean, "_", safe_file_label(region), "_", index, "_GAM_residual_QQ.png")
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
    sub <- read_index_series(
      dat = dat,
      region = region,
      season = season,
      index = index,
      model_scale = TRUE
    )
    sub <- sub[is.finite(sub$Year) & is.finite(sub$value), , drop = FALSE]

    if (nrow(sub) < 8L || length(unique(sub$value)) < 3L) {
      plot.new()
      title(main = paste0(season_labels[[season]], ": insufficient data"))
      next
    }

    spline_k <- min(10L, max(3L, floor(nrow(sub) / 4L)))
    fit <- gam(value ~ s(Year, k = spline_k), data = sub, method = "REML")
    residual <- residuals(fit)
    qqnorm(residual, main = season_labels[[season]])
    qqline(residual, col = "red", lwd = 1.5)
  }

  plot.new()
  scale_note <- if (index == "wsmax") "cube-root scale" else "proportion scale"
  mtext(
    paste0(source, " ", member, " - GAM residual QQ (", index, ", ", scale_note, ")"),
    outer = TRUE,
    cex = 0.95,
    font = 2
  )
  invisible(fig_path)
}

save_haduk_distribution <- function(haduk, region, index, model_scale) {
  scale_name <- if (model_scale && index == "wsmax") "cube_root" else "raw"
  out_dir <- ensure_dir(file.path(output_root, "diagnostics", "HadUK", index))
  fig_path <- file.path(
    out_dir,
    paste0("HadUK_", safe_file_label(region), "_", index, "_", scale_name, "_distribution.png")
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
    sub <- read_index_series(haduk, region, season, index, model_scale)
    x <- sub$value[is.finite(sub$value)]
    if (length(x) < 2L) {
      plot.new()
      title(main = paste0(season_labels[[season]], ": no data"))
      next
    }

    hist(
      x,
      probability = TRUE,
      breaks = 20,
      col = "lightblue",
      border = "white",
      xlab = if (model_scale && index == "wsmax") "cube-root(days)" else index,
      main = season_labels[[season]]
    )
    mu <- mean(x)
    sigma <- sd(x)
    if (is.finite(sigma) && sigma > 0) {
      curve(dnorm(x, mean = mu, sd = sigma), add = TRUE, col = "red", lwd = 2)
    }
  }

  plot.new()
  mtext(
    paste0("HadUK distribution - ", region, " (", index, ", ", scale_name, ")"),
    outer = TRUE,
    cex = 1.0,
    font = 2
  )
  invisible(fig_path)
}

run_gam_diagnostics <- function() {
  haduk <- read_haduk(prefer_filtered = TRUE)
  outputs <- list()

  for (index in indices_to_run) {
    validate_index(index)
    for (region in regions_to_run) {
      if ("HadUK" %in% sources_to_run) {
        out_dir <- ensure_dir(file.path(output_root, "diagnostics", "HadUK", index))
        key <- paste("HadUK", index, region, "qq", sep = "/")
        outputs[[key]] <- save_qq_grid(haduk, "HadUK", "observations", region, index, out_dir)
        outputs[[paste("HadUK", index, region, "raw", sep = "/")]] <-
          save_haduk_distribution(haduk, region, index, model_scale = FALSE)
        if (index == "wsmax") {
          outputs[[paste("HadUK", index, region, "cube_root", sep = "/")]] <-
            save_haduk_distribution(haduk, region, index, model_scale = TRUE)
        }
      }

      for (source in intersect(sources_to_run, c("UKCP18", "CORDEX"))) {
        files <- find_ensemble_files(source)
        out_dir <- ensure_dir(file.path(output_root, "diagnostics", source, index))
        for (path in files) {
          member <- read.csv(path, na.strings = c("  NA", "NA", "", "NaN"))
          member_name <- ensemble_member_name(path, source)
          key <- paste(source, member_name, index, region, sep = "/")
          outputs[[key]] <- save_qq_grid(
            member,
            source,
            member_name,
            region,
            index,
            out_dir
          )
        }
      }
    }
  }

  cat("Saved", length(outputs), "diagnostic figures under", file.path(output_root, "diagnostics"), "\n")
  invisible(outputs)
}

if (sys.nframe() == 0L) run_gam_diagnostics()
