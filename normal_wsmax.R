library(MASS)
# read table
haduk <- read.csv("/Users/persevere/Downloads/normal/PrecipData/HadUK_1950_2021.csv", na.strings = c("  NA", "NA", ""))

# file out
out_dir <- file.path(path.expand("~"), "Downloads", "HadUK_hist_wsmax")
if (!dir.exists(out_dir)) dir.create(out_dir)

var_name <- "wamax"

fig_width  <- 10   # 5 panels -> wider
fig_height <- 7
fig_res    <- 300

# class
seasons <- c("Annual", "DJF", "MAM", "JJA", "SON")
sea_label <- c(
  Annual = "Annual",
  DJF    = "winter",
  MAM    = "spring",
  JJA    = "summer",
  SON    = "autumn"
)
regions <- sort(unique(haduk$Region))

for (reg in regions){
  reg_clean <- gsub(" ", "", reg)
  fig_name <- paste0("hist_", reg_clean, "_", var_name, ".png")
  fig_path <- file.path(out_dir, fig_name)
  
  x11(width = fig_width, height = fig_height)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  for (sea in seasons){
    sub_obs <- subset(haduk, Region == reg & Season == sea)
    show_sea <- sea_label[[sea]]
    if (nrow(sub_obs) == 0){
      plot.new()
      title(main = paste0(reg, " ", show_sea, " (no obs)"))
      next
    }
    x <- na.omit(sub_obs$fwd)
    
    if (length(x) == 0) {
      plot.new()
      title(main = paste0(reg, " ", show_sea, " (no valid data)"))
    } else {
      hist(x,
           probability = TRUE,
           breaks = 35,
           main = paste0(reg, " ", show_sea, " - fwd density"),
           xlab = "fwd",
           col = "lightblue",
           border = "white")
      
      mu <- mean(x)
      sigma <- sd(x)
      
      curve(dnorm(x, mean = mu, sd = sigma),
            add = TRUE,
            col = "red",
            lwd = 2)
    }
  }
  plot.new()
  
  mtext(paste0("Histogram of HadUK observations – ", reg, " (", var_name, ")"),
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