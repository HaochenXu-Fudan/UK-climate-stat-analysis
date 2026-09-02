# Run the two-way EuroCORDEX smoother for both precipitation indices.
# fwd stays on proportion scale; wsmax is fitted on cube-root scale and
# converted back to days in plots and summaries.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else if (basename(getwd()) == "R") {
  normalizePath(getwd())
} else {
  normalizePath(file.path(getwd(), "R"))
}
project_dir <- normalizePath(file.path(script_dir, ".."))

indices_to_run <- c("fwd", "wsmax")

run_cordex_smoothing <- function() {
  core_file <- file.path(script_dir, "smooth_cordex_index.R")
  if (!file.exists(core_file)) stop("Cannot find ", core_file)

  all_results <- setNames(vector("list", length(indices_to_run)), indices_to_run)
  for (index in indices_to_run) {
    run_env <- new.env(parent = globalenv())
    run_env$project_dir <- project_dir
    run_env$var_name <- index
    sys.source(core_file, envir = run_env)
    all_results[[index]] <- run_env$result
  }

  out_root <- file.path(project_dir, "outputs", "CORDEX_smooth")
  if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)
  saveRDS(
    all_results,
    file.path(out_root, "CORDEX_London_fwd_wsmax_all_results.rds")
  )

  cat("\nAll requested EuroCORDEX postprocessing runs finished.\n")
  cat("Main output folder:", out_root, "\n")
  invisible(all_results)
}

if (sys.nframe() == 0L) run_cordex_smoothing()
