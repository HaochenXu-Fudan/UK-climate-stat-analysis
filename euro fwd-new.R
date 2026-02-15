library(TimSPEC)

## ============================
## 0) Read observations (HadUK)
## ============================
## 1950-2021
haduk <- read.csv("/Users/persevere/Downloads/压缩/PrecipData/HadUK_1950_2021.csv", na.strings = c("  NA", "NA", ""))

## 1892-2021
# haduk <- read.csv("/Users/persevere/Downloads/压缩/PrecipData/HadUK.csv", na.strings = c("  NA", "NA", ""))

## 你这一步是 EuroCORDEX 的 fwd
var_name <- "fwd"

## ============================
## 1) Find EuroCORDEX files
## CORDEX_$$$_%%%_@@@.csv
## ============================
setwd("/Users/persevere/Downloads/压缩/PrecipData")
csv_files <- list.files(pattern = "\\.csv$")
csv_cordex <- csv_files[grep("^CORDEX_.*\\.csv$", csv_files)]

cat("Number of CORDEX files found:", length(csv_cordex), "\n")
if (length(csv_cordex) == 0) stop("No CORDEX_*.csv files found in the working directory.")

## ============================
## 2) Output settings (paper-friendly)
## ============================
out_dir <- file.path(path.expand("~"), "Downloads", "HadUK_byRegion_2figs")
# out_dir <- "figures_CORDEX_fwd_byRegion"
if (!dir.exists(out_dir)) dir.create(out_dir)

fig_width  <- 10   # 一张图放 5 个panel，建议更宽一点
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

## 如果你想先测试一两个地区，临时用：
# regions <- c("London", "South East")

## ============================
## 3) Helper: build merged Data for one (region, season)
## ============================
build_cordex_merged <- function(reg, sea, var_name, haduk, csv_cordex) {
  
  ## --- Obs ---
  sub_obs <- subset(haduk, Region == reg & Season == sea)
  
  if (nrow(sub_obs) == 0) return(NULL)
  
  Data_obs <- data.frame(
    year = sub_obs$Year,
    obs  = as.numeric(sub_obs[[var_name]])
  )
  
  ## --- CORDEX ensemble members ---
  ensemble_list <- list()
  
  for (file in csv_cordex) {
    
    dat <- read.csv(file, na.strings = c("  NA", "NA", ""))
    
    dat_sub <- subset(dat, Region == reg & Season == sea)
    if (nrow(dat_sub) == 0) next
    
    ## 解析文件名：CORDEX_GCM_RCM_RUN.csv
    parts <- strsplit(gsub("\\.csv$", "", file), "_")[[1]]
    if (length(parts) < 4) next
    
    gcm <- parts[2]
    rcm <- parts[3]
    run <- parts[4]
    
    colname <- paste0("ens_", gcm, "_", rcm, "_", run)
    
    df <- data.frame(
      year  = dat_sub$Year,
      value = as.numeric(dat_sub[[var_name]])
    )
    names(df)[2] <- colname
    
    ensemble_list[[colname]] <- df
  }
  
  ## 若该 region-season 在 CORDEX 没有任何成员，就只返回 obs
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
  
  ## 论文友好文件名：去空格
  reg_clean <- gsub(" ", "", reg)
  fig_name <- paste0("CORDEX_", reg_clean, "_", var_name, "_Obs_Ensemble_5seasons.png")
  fig_path <- file.path(out_dir, fig_name)
  
  ## 打开固定尺寸设备（导师建议）
  x11(width = fig_width, height = fig_height)
  
  ## 5 panels：2x3 布局，留一个空位
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  for (sea in seasons) {
    
    Merged <- build_cordex_merged(reg, sea, var_name, haduk, csv_cordex)
    show_sea <- sea_label[[sea]] 
    
    if (is.null(Merged)) {
      ## 如果观测都没有，就画空图占位
      plot.new()
      title(main = paste0(reg, " ", show_sea, " (no obs)"))
      next
    }
    
    ## 如果只有两列(year, obs)，说明没 CORDEX members
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
      Lwd = Lwd,            # 或者试试 lwd = Lwd
      EnsTransp = 0.15,
      Units = "Proportion",
      main = paste0(reg, " ", show_sea, " – ", var_name)
    )
  }
  
  ## 最后一个 panel 空出来（因为 2x3=6格，我们只用5格）
  plot.new()
  
  ## 外标题（整张图的title）
  mtext(paste0("EuroCORDEX ensemble vs HadUK observations – ", reg, " (", var_name, ")"),
        outer = TRUE, cex = 1.1)
  
  ## 保存：dev.copy -> png（导师建议）
  dev.copy(png, filename = fig_path,
           width  = fig_width  * fig_res,
           height = fig_height * fig_res,
           res    = fig_res,
           bg     = "white")
  
  dev.off()  # close png
  dev.off()  # close x11
  
  cat("Saved:", fig_path, "\n")
}



