# Shared configuration and data helpers for the climate analysis scripts.
# Runnable scripts set project_dir before sourcing this file.

if (!exists("project_dir", inherits = FALSE)) {
  stop("Set project_dir before sourcing R/climate_common.R.")
}

data_dir <- file.path(project_dir, "PrecipData")
output_root <- file.path(project_dir, "outputs")

climate_indices <- c("fwd", "wsmax")
climate_seasons <- c("Annual", "DJF", "MAM", "JJA", "SON")
season_labels <- c(
  Annual = "Annual",
  DJF = "Winter",
  MAM = "Spring",
  JJA = "Summer",
  SON = "Autumn"
)

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

validate_index <- function(index) {
  if (length(index) != 1L || !(index %in% climate_indices)) {
    stop("index must be either 'fwd' or 'wsmax'.")
  }
  invisible(index)
}

index_units <- function(index) {
  validate_index(index)
  if (index == "fwd") "proportion" else "days"
}

to_model_scale <- function(x, index) {
  validate_index(index)
  x <- as_numeric_clean(x)
  if (index == "wsmax") pmax(x, 0)^(1 / 3) else x
}

to_report_scale <- function(x, index) {
  validate_index(index)
  if (index == "wsmax") pmax(x, 0)^3 else x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  normalizePath(path, mustWork = TRUE)
}

read_haduk <- function(prefer_filtered = TRUE) {
  candidates <- if (prefer_filtered) {
    c(
      file.path(data_dir, "HadUK_1950_2021.csv"),
      file.path(data_dir, "HadUK.csv")
    )
  } else {
    c(
      file.path(data_dir, "HadUK.csv"),
      file.path(data_dir, "HadUK_1950_2021.csv")
    )
  }

  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) stop("Cannot find HadUK.csv or HadUK_1950_2021.csv.")

  dat <- read.csv(path, na.strings = c("  NA", "NA", "", "NaN"))
  required <- c("Year", "Season", "Region", climate_indices)
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("HadUK file is missing columns: ", paste(missing, collapse = ", "))
  }
  dat$Year <- as_numeric_clean(dat$Year)
  attr(dat, "source_file") <- path
  dat
}

find_ensemble_files <- function(source) {
  source <- match.arg(source, c("UKCP18", "CORDEX"))
  pattern <- if (source == "UKCP18") {
    "^UKCP18_[0-9]+\\.csv$"
  } else {
    "^CORDEX_.*\\.csv$"
  }

  files <- sort(list.files(data_dir, pattern = pattern, full.names = TRUE))
  if (source == "CORDEX" && length(files) == 0L) {
    subdir <- file.path(data_dir, "CORDEX-data")
    if (dir.exists(subdir)) {
      files <- sort(list.files(subdir, pattern = pattern, full.names = TRUE))
    }
  }
  if (length(files) == 0L) stop("No ", source, " CSV files found.")
  files
}

parse_cordex_member <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  matched <- regexec(
    "^CORDEX_(.+)_([^_]+)_(r[0-9]+i[0-9]+p[0-9]+)$",
    stem
  )
  parts <- regmatches(stem, matched)[[1]]
  if (length(parts) != 4L) {
    stop("Unexpected CORDEX filename: ", basename(path))
  }
  list(gcm = parts[2], rcm = parts[3], run = parts[4], name = stem)
}

ensemble_member_name <- function(path, source) {
  source <- match.arg(source, c("UKCP18", "CORDEX"))
  if (source == "UKCP18") {
    sub("\\.csv$", "", basename(path))
  } else {
    parse_cordex_member(path)$name
  }
}

read_index_series <- function(dat, region, season, index, model_scale = FALSE) {
  validate_index(index)
  required <- c("Year", "Season", "Region", index)
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("Input data is missing columns: ", paste(missing, collapse = ", "))
  }

  ans <- dat[
    dat$Region == region & dat$Season == season,
    c("Year", index)
  ]
  names(ans) <- c("Year", "value")
  ans$Year <- as_numeric_clean(ans$Year)
  ans$value <- if (model_scale) {
    to_model_scale(ans$value, index)
  } else {
    as_numeric_clean(ans$value)
  }
  ans <- ans[order(ans$Year), , drop = FALSE]
  rownames(ans) <- NULL
  ans
}

build_ensemble_matrix <- function(
    haduk,
    files,
    source,
    region,
    season,
    index,
    model_scale = FALSE,
    year_min = 1950,
    year_max = 2080
) {
  source <- match.arg(source, c("UKCP18", "CORDEX"))
  obs <- read_index_series(haduk, region, season, index, model_scale)
  names(obs)[2] <- "HadUK"
  merged <- obs

  for (path in files) {
    member <- read.csv(path, na.strings = c("  NA", "NA", "", "NaN"))
    sub <- read_index_series(member, region, season, index, model_scale)
    names(sub)[2] <- ensemble_member_name(path, source)
    merged <- merge(merged, sub, by = "Year", all = TRUE)
  }

  merged <- merged[
    merged$Year >= year_min & merged$Year <= year_max,
    ,
    drop = FALSE
  ]
  merged <- merged[order(merged$Year), , drop = FALSE]
  rownames(merged) <- NULL
  merged
}

align_erf <- function(years, erf) {
  if (is.vector(erf) && !is.data.frame(erf)) {
    if (length(erf) != length(years)) {
      stop("ERF vector length does not match the requested years.")
    }
    return(as_numeric_clean(erf))
  }

  erf <- as.data.frame(erf)
  year_col <- grep("^year$|year", names(erf), ignore.case = TRUE, value = TRUE)[1]
  value_col <- grep("^NetERF$|net.*erf|erf", names(erf), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(year_col) || is.na(value_col)) {
    stop("Cannot identify Year and NetERF columns in ERF data.")
  }

  ans <- as_numeric_clean(erf[[value_col]])[
    match(years, as_numeric_clean(erf[[year_col]]))
  ]
  if (anyNA(ans)) {
    stop("Missing ERF values for years: ", paste(years[is.na(ans)], collapse = ", "))
  }
  ans
}

safe_file_label <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "", gsub(" ", "", x))
}
