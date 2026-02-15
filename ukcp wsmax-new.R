library(TimSPEC)

## ============================
## 0) Read observations (HadUK)
## ============================
## 1950-2021
haduk <- read.csv("/Users/persevere/Downloads/压缩/PrecipData/HadUK_1950_2021.csv", na.strings = c("  NA", "NA", ""))

## 1892-2021
# haduk <- read.csv("/Users/persevere/Downloads/压缩/PrecipData/HadUK.csv", na.strings = c("  NA", "NA", ""))

## 这一步：UKCP 的 wsmax
var_name <- "wsmax"

## ============================
## 1) Find UKCP files
## UKCP18_##.csv
## ============================
csv_files <- list.files(pattern = "\\.csv$")
csv_ukcp  <- csv_files[grep("^UKCP18_\\d+\\.csv$", csv_files)]

cat("Number of UKCP files found:", length(csv_ukcp), "\n")
if (length(csv_ukcp) == 0) stop("No UKCP18_*.csv files found in the working directory.")

## ============================
## 2) Output settings (paper-friendly)
## ============================
out_dir <- file.path(path.expand("~"), "Downloads", "HadUK_byRegion_2figs")
# out_dir <- "figures_UKCP_wsmax_byRegion"
if (!dir.exists(out_dir)) dir.create(out_dir)

fig_width  <- 10   # 5 panels -> wider
fig_height <- 7
fig_res    <- 300

seasons <- c("Annual", "DJF", "MAM", "JJA", "SON")
sea_label <- c(
  Annual = "Annual",
  DJF    = "winter",
  MAM    = "spring",
  JJA    = "summer",
  SON    = "autumn"
)
regions <- sort(unique(haduk$Region))

## ============================
## 3) Helper: build merged Data for one (region, season)
## ============================
build_ukcp_merged <- function(reg, sea, var_name, haduk, csv_ukcp) {
  
  ## --- Obs ---
  sub_obs <- subset(haduk, Region == reg & Season == sea)
  if (nrow(sub_obs) == 0) return(NULL)
  
  Data_obs <- data.frame(
    year = sub_obs$Year,
    obs  = as.numeric(sub_obs[[var_name]])
  )
  
  ## --- UKCP ensemble members ---
  ensemble_list <- list()
  
  for (file in csv_ukcp) {
    
    dat <- read.csv(file, na.strings = c("  NA", "NA", ""))
    
    dat_sub <- subset(dat, Region == reg & Season == sea)
    if (nrow(dat_sub) == 0) next
    
    ## 解析 UKCP 文件名：UKCP18_01.csv -> id = "01"
    id <- strsplit(file, "[_.]")[[1]][2]
    colname <- paste0("ens", id)
    
    df <- data.frame(
      year  = dat_sub$Year,
      value = as.numeric(dat_sub[[var_name]])
    )
    names(df)[2] <- colname
    
    ensemble_list[[colname]] <- df
  }
  
  ## 如果该 region-season 没有 UKCP members，就只返回 obs
  if (length(ensemble_list) == 0) {
    return(Data_obs)
  }
  
  ## Merge obs + all ensembles
  Merged <- Reduce(
    function(x, y) merge(x, y, by = "year", all = TRUE),
    c(list(Data_obs), ensemble_list)
  )
  
  return(Merged)
}

## ============================
## 4) Batch: One Region -> One File (5 panels)
## Using: for loop + paste0 + x11 + dev.copy
## ============================
for (reg in regions) {
  
  reg_clean <- gsub(" ", "", reg)
  fig_name <- paste0("UKCP18_", reg_clean, "_", var_name, "_Obs_Ensemble_5seasons.png")
  fig_path <- file.path(out_dir, fig_name)
  
  x11(width = fig_width, height = fig_height)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  for (sea in seasons) {
    
    Merged <- build_ukcp_merged(reg, sea, var_name, haduk, csv_ukcp)
    show_sea <- sea_label[[sea]]
    
    if (is.null(Merged)) {
      plot.new()
      title(main = paste0(reg, " ", show_sea, " (no obs)"))
      next
    }
    
    ## 只有 obs
    if (ncol(Merged) <= 2) {
      plot(Merged$year, Merged$obs, type = "l",
           xlab = "Year", ylab = var_name,
           main = paste0(reg, " ", show_sea, " (Obs only)"))
      next
    }
    
    ## Colours：obs=black, ensMean=red, members=coral3
    n_ens <- ncol(Merged) - 2
    Lwd <- c(
      0.8,                  # obs (black)
      0.8,                  # ensMean (red)
      rep(0.4, n_ens)        # members (coral3)
    )
    ens_cols <- rep("coral3", n_ens)
    Colours <- c("black", "red", ens_cols)
    
    PlotEnsTS(
      Merged,
      Colours = Colours,
      Lwd = Lwd,
      EnsTransp = 0.15,
      Units = "days",
      main = paste0(reg, " ", show_sea, " – ", var_name)
    )
  }
  
  plot.new()  # 空 panel
  
  mtext(paste0("UKCP18 ensemble vs HadUK observations – ", reg, " (", var_name, ")"),
        outer = TRUE, cex = 1.1)
  
  dev.copy(png, filename = fig_path,
           width  = fig_width  * fig_res,
           height = fig_height * fig_res,
           res    = fig_res,
           bg     = "white")
  
  dev.off()  # close png
  dev.off()  # close x11
  
  cat("Saved:", fig_path, "\n")
}


