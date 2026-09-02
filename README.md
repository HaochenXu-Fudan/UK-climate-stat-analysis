# UK 区域降水指数：数据与集合后处理

本项目分析英国 16 个区域的两个降水指数，比较 HadUK 观测与 UKCP18、EuroCORDEX
集合模拟，并完成数据诊断、趋势平滑、后验预测和集合权重诊断。

维护中的 R 代码统一放在 `R/`，输入数据放在 `PrecipData/`，新生成的图、表和
RDS 对象统一写入 `outputs/`。

> **数据安全说明**：常规分析脚本只读取 `PrecipData/`。只有
> `R/prepare_haduk.R` 会写入该目录；它会根据 `HadUK.csv` 重新生成并覆盖
> `HadUK_1950_2021.csv`，不会修改 `HadUK.csv`、UKCP18 或 CORDEX 原始文件。

## 1. 两个降水指数

| 变量 | 含义 | 湿日定义 | 建模与展示方式 |
|---|---|---|---|
| `fwd` | 区域内各网格点湿日比例的平均值 | 日降水量至少 1 mm | 在原始比例尺度上建模和展示 |
| `wsmax` | 区域内各网格点“期内开始的最长湿润连续过程”的平均持续时间 | 连续湿日，且前后均为干日 | 先作立方根变换后建模，再还原为天数作图和汇总 |

`fwd` 和 `wsmax` 共用读数、季节循环和绘图框架，但不是同一个统计量。代码为它们
保留了不同的尺度变换、先验和初始方差，不能通过简单改列名互换。

## 2. 时间、季节与区域

### 年份定义

数据中的一年按 **上一年 12 月至当年 11 月** 计算。例如，`Year = 1970`
表示 1969 年 12 月至 1970 年 11 月。

### 季节代码

- `Annual`：全年
- `DJF`：冬季
- `MAM`：春季
- `JJA`：夏季
- `SON`：秋季

### 16 个区域

Channel Islands、East Midlands、East Scotland、East of England、Isle of Man、
London、North East England、North Scotland、North West England、Northern Ireland、
South East England、South West England、Wales、West Midlands、West Scotland、
Yorkshire and Humber。

HadUK 中 Channel Islands 的 `fwd` 和 `wsmax` 全部缺失；UKCP18 和 CORDEX 文件
中包含该区域的模拟值。

## 3. 目录结构

```text
climate/
├── README.md                  # 本说明
├── CODE_OVERVIEW.md           # 早期代码整理摘要
├── R/                         # 当前维护的 R 代码
├── PrecipData/                # 输入数据与原始数据说明
├── outputs/                   # 新运行产生的图、CSV 和 RDS
├── Interpreting time series from ensembles of climate—EnsemblePostprocessing.pdf
├── Interpreting time series from ensembles of climate—EnsemblePostprocessing_Supplement.pdf
└── Thesis_Template_for_University_College_London__1_.pdf
```

`PrecipData/figures/`、`PrecipData/figures_UKCP18/`、压缩包和测试图片是已有的历史
图形或辅助文件，不是当前维护脚本的输入。当前脚本的新结果一律写入 `outputs/`。

## 4. 数据文件

| 文件 | 数量 | 时间范围 | 说明 |
|---|---:|---|---|
| `PrecipData/HadUK.csv` | 1 | 1892–2021 | 由 HadUK 格点数据计算得到的区域观测值 |
| `PrecipData/HadUK_1950_2021.csv` | 1 | 1950–2021 | 当前分析常用的 HadUK 子集，可由准备脚本重建 |
| `PrecipData/UKCP18_##.csv` | 12 | 1981–2080 | UKCP18 集合成员；`##` 是成员编号 |
| `PrecipData/CORDEX_<GCM>_<RCM>_<run>.csv` | 64 | 1981–2078/2079/2080 | EuroCORDEX 成员；文件名记录 GCM、RCM 和运行编号 |
| `PrecipData/regionmask-region_osgb.nc` | 1 | 不适用 | HadUK 网格上的 16 区域掩膜，维度为 82 × 112 × 16 |

UKCP18 和 CORDEX 投影在 2015 年之前使用历史温室气体排放，之后使用 RCP8.5
情景。当前数据中的 12 个 UKCP18 成员均覆盖至 2080；64 个 CORDEX 成员中，
47 个覆盖至 2080、9 个覆盖至 2079、8 个覆盖至 2078。

所有 CSV 使用同一组字段：

| 字段 | 内容 |
|---|---|
| `Year` | 按上一年 12 月至当年 11 月定义的年份 |
| `Season` | `Annual`、`DJF`、`MAM`、`JJA` 或 `SON` |
| `Region` | 16 个区域之一 |
| `fwd` | 湿日比例 |
| `wsmax` | 最长湿润连续过程的平均持续时间，单位为天 |

### 读取区域掩膜

区域掩膜不参与当前 CSV 后处理脚本；需要检查原始网格区域时，可单独使用
`RNetCDF`：

```r
library(RNetCDF)

nc_file <- open.nc("PrecipData/regionmask-region_osgb.nc")
region_names <- var.get.nc(nc_file, "geo_region")
region_masks <- var.get.nc(nc_file, "region_mask")
close.nc(nc_file)

map_data <- apply(
  region_masks,
  MARGIN = 1:2,
  FUN = function(x) sum(x * seq_along(region_names))
)
map_data[map_data == 0] <- NA
image(map_data, col = hcl.colors(length(region_names), palette = "Set 2"))
```

## 5. 软件依赖

需要 R，以及下列 R 包：

- `TimSPEC`：集合 EBM 平滑、后验预测、集合时序图和 SSP5-8.5 NetERF 数据；
- `mgcv`：GAM 拟合与残差诊断；
- `ggplot2`：NetERF 曲线；
- `RNetCDF`：仅在读取区域掩膜时需要。

项目目前没有自动安装依赖的脚本。运行前请先确认以上包在所用 R 环境中可用。

## 6. 代码入口

请从项目根目录运行脚本，例如：

```sh
cd /Users/persevere/Downloads/climate
Rscript R/plot_ensemble.R
```

### 可直接运行的脚本

| 脚本 | 功能 | 主要输出 |
|---|---|---|
| `R/prepare_haduk.R` | 从完整 HadUK 数据生成 1950–2021 子集 | `PrecipData/HadUK_1950_2021.csv` |
| `R/plot_ensemble.R` | 绘制 HadUK、集合均值和所有成员的五季节原始时序图 | `outputs/ensemble_raw/` |
| `R/diagnose_gam.R` | 对 HadUK、每个 UKCP18 成员和每个 CORDEX 成员作 GAM 残差 QQ 图；另作 HadUK 分布图 | `outputs/diagnostics/` |
| `R/smooth_observations.R` | 对 HadUK 两个指数分别作 observation-only EBM 趋势平滑 | `outputs/observation_smooth/` |
| `R/smooth_ukcp.R` | 对 UKCP18 两个指数作集合平滑、两类后验预测和权重诊断 | `outputs/UKCP_smooth/` |
| `R/smooth_cordex.R` | 对 CORDEX 两个指数作两层集合平滑、两类后验预测和权重诊断 | `outputs/CORDEX_smooth/` |
| `R/plot_erf.R` | 绘制平滑模型使用的 SSP5-8.5 NetERF 协变量 | `outputs/ERF/` |

### 内部模块

以下文件由入口脚本载入，不应单独运行：

- `R/climate_common.R`：路径、索引变换、数据读取、成员合并和 ERF 对齐；
- `R/smooth_helpers.R`：后验权重、分位数和 `wsmax` 天数尺度绘图；
- `R/smooth_cordex_index.R`：由 `R/smooth_cordex.R` 分别以 `fwd`、`wsmax`
  调用的单指数 CORDEX 实现。

## 7. 推荐运行顺序

如果 `HadUK_1950_2021.csv` 已存在且不需要更新，跳过第 1 步。

```sh
# 1. 可选：重建 HadUK 子集；此命令会覆盖同名子集文件
Rscript R/prepare_haduk.R

# 2. 查看输入数据与模型协变量
Rscript R/plot_ensemble.R
Rscript R/plot_erf.R

# 3. 检查分布与 GAM 残差
Rscript R/diagnose_gam.R

# 4. 观测数据平滑
Rscript R/smooth_observations.R

# 5. 集合后处理
Rscript R/smooth_ukcp.R
Rscript R/smooth_cordex.R
```

`diagnose_gam.R` 会为每个集合成员分别生成诊断图，因此文件较多。
`smooth_ukcp.R` 和 `smooth_cordex.R` 默认对两个指数、五个季节分别拟合并进行
后验抽样，是整个流程中耗时较长的部分。

## 8. 默认分析设置

- 区域：London；
- 季节：Annual、DJF、MAM、JJA、SON；
- 观测平滑年份：1950–2021；
- 集合后处理年份：1950–2080；
- 初始观测时段：1950–1979；
- 每个季节的后验预测样本数：1,000；
- 后验抽样随机种子：2,000。

脚本当前不解析命令行参数。需要更换区域、指数、年份或抽样数时，修改对应脚本
顶部的 `*_to_run`、`region_to_use`、`year_min`、`year_max` 或 `N_pps`；CORDEX
的具体设置位于内部文件 `R/smooth_cordex_index.R`。正式结果建议保留
`N_pps = 1000`；仅作代码调试时可暂时降低。

> **更换区域时请注意**：UKCP18 和 CORDEX 平滑脚本的部分输出文件名目前直接包含
> `London`。除修改 `region_to_use` 外，还应同步修改文件名模板，避免结果内容与文件名
> 不一致。

## 9. 集合处理约定

### UKCP18

- 12 个成员按可交换集合处理，`Groups = NULL`；
- 不人为拆分 UKCP18_A、UKCP18_B 等不存在的分组；
- 如果同一年只有部分成员缺失，脚本会停止并报告年份，不会自动用行均值填补。

### EuroCORDEX

- 从文件名解析 GCM、RCM 和运行编号；
- 当前 64 个成员包含 6 个 GCM、10 个 RCM 和 4 种运行编号；
- 平滑模型保留 RCM 与 GCM 两个分组维度，且不加入二者交互项；
- 某年仅有部分成员缺失时，用该年其余可用成员的平均值补齐；全体成员都缺失时
  不作这种填补。这个规则主要处理部分 CORDEX 成员在 2078 或 2079 后结束的情况。

## 10. 平滑分析保留的结果

UKCP18 和 CORDEX 对每个指数均保留下列不同用途的结果：

1. MAP 平滑趋势图；
2. 基于 Laplace 近似的后验预测图；
3. 基于重要性抽样的后验预测图；
4. 五个季节各自的重要性权重诊断图；
5. 一张包含五个季节的累计权重诊断图；
6. 初始值、MAP 参数和 2080 后验汇总 CSV；
7. 每个指数的完整 RDS，以及同时包含两个指数的合并 RDS。

`wsmax` 的模型内部结果保留在立方根尺度；面向解释的图和 2080 汇总会转换回天数。
凡文件名带 `_days` 的 `wsmax` 图，纵轴均为还原后的天数。MAP、Laplace 近似、
重要性抽样和权重诊断回答的问题不同，不应互相替代或删除。

## 11. 输出目录

```text
outputs/
├── ERF/
├── ensemble_raw/
│   ├── UKCP18/{fwd,wsmax}/
│   └── CORDEX/{fwd,wsmax}/
├── diagnostics/
│   ├── HadUK/{fwd,wsmax}/
│   ├── UKCP18/{fwd,wsmax}/
│   └── CORDEX/{fwd,wsmax}/
├── observation_smooth/{fwd,wsmax}/
├── UKCP_smooth/
│   ├── UKCP18_fwd/
│   └── UKCP18_wsmax/
└── CORDEX_smooth/
    ├── CORDEX_fwd/
    └── CORDEX_wsmax/
```

## 12. 数据来源说明

数据文件的原始说明由 Richard Chandler 于 2025 年 11 月 19 日提供，现保存在
`PrecipData/README.docx` 和 `PrecipData/README.txt`。本 README 在保留原始变量、
文件命名和数据来源含义的基础上，补充了当前代码结构、实际文件覆盖范围、运行方法、
输出位置以及 `fwd` 与 `wsmax` 的不同处理方式。
