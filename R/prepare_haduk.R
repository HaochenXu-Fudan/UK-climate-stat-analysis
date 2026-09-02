# Create the 1950-2021 HadUK subset used by the analysis scripts.

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

prepare_haduk <- function(min_year = 1950, max_year = 2021) {
  input <- file.path(data_dir, "HadUK.csv")
  output <- file.path(data_dir, "HadUK_1950_2021.csv")
  if (!file.exists(input)) stop("Cannot find ", input)

  dat <- read.csv(input, na.strings = c("  NA", "NA", "", "NaN"))
  required <- c("Year", "Season", "Region", "fwd", "wsmax")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("HadUK.csv is missing columns: ", paste(missing, collapse = ", "))
  }

  dat$Year <- as_numeric_clean(dat$Year)
  dat <- dat[dat$Year >= min_year & dat$Year <= max_year, , drop = FALSE]
  dat <- dat[order(dat$Region, dat$Season, dat$Year), , drop = FALSE]
  write.csv(dat, output, row.names = FALSE, na = "NA")
  cat("Saved", nrow(dat), "rows to", output, "\n")
  invisible(output)
}

if (sys.nframe() == 0L) prepare_haduk()
