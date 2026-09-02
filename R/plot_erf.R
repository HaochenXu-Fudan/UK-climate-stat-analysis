# Plot the SSP5-8.5 effective radiative forcing series used by the smoothers.

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
library(ggplot2)
data(SSP585data)

plot_erf <- function() {
  out_dir <- ensure_dir(file.path(output_root, "ERF"))
  fig_path <- file.path(out_dir, "SSP585_NetERF.png")

  p <- ggplot(ERF585, aes(x = Year, y = NetERF)) +
    geom_line(color = "dodgerblue", linewidth = 1.2) +
    scale_x_continuous(breaks = seq(1950, 2090, by = 20)) +
    labs(x = "Year", y = expression(W ~ m^-2)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold")
    )

  ggsave(fig_path, plot = p, width = 10, height = 4, dpi = 300, bg = "white")
  cat("Saved:", fig_path, "\n")
  invisible(fig_path)
}

if (sys.nframe() == 0L) plot_erf()
