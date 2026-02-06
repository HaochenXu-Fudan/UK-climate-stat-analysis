# version 1
# 读取
df <- read.csv("/Users/persevere/Downloads/压缩/PrecipData/HadUK.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 第一列当作年份列（不管它列名叫啥）
year_col <- names(df)[1]

# 确保是数值型（如果是字符也能转）
df[[year_col]] <- as.integer(df[[year_col]])

# 过滤：只保留 year >= 1950 的行（等价于删掉 <=1949）
df2 <- df[df[[year_col]] >= 1950, ]

# 保存为新文件
write.csv(df2, "/Users/persevere/Downloads/压缩/PrecipData/HadUK_1950_2021.csv", row.names = FALSE)

###############################
# version 2

# library(readr)
# library(dplyr)

# df2 <- read_csv("HadUK.csv", show_col_types = FALSE) %>%
#   filter(as.integer(.data[[names(.)[1]]]) >= 1950)

# write_csv(df2, "/Users/persevere/Downloads/压缩/PrecipData/HadUK_1950_2021_new.csv")
