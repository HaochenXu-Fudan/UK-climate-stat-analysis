# Climate analysis code

All maintained R code is in `R/`. Input data remains in `PrecipData/`, and all
new figures, tables and RDS files are written to `outputs/`.

## Variable treatment

- `fwd`: mean wet-day fraction. It is fitted and reported on the original
  proportion scale.
- `wsmax`: mean maximum wet-spell duration. It is fitted on cube-root scale to
  improve the residual distribution, then transformed back to days for plots
  and reported summaries.

The two indices share data-handling and plotting code, but keep separate prior
settings, initial variances and reporting transformations where required.

## Runnable scripts

Run commands from the project directory, for example:

```sh
Rscript R/plot_ensemble.R
```

- `R/prepare_haduk.R`: recreate `PrecipData/HadUK_1950_2021.csv` from
  `HadUK.csv`.
- `R/plot_ensemble.R`: raw five-season HadUK-versus-ensemble plots for UKCP18
  and CORDEX, preserving observations, ensemble mean and individual members.
- `R/diagnose_gam.R`: HadUK, UKCP18-member and CORDEX-member GAM residual QQ
  plots, plus HadUK distribution plots.
- `R/smooth_observations.R`: observation-only EBM smoothing for both indices.
- `R/smooth_ukcp.R`: final UKCP18 analysis for both indices. UKCP18 members are
  treated as an exchangeable ensemble (`Groups = NULL`).
- `R/smooth_cordex.R`: final CORDEX analysis for both indices. GCM and RCM are
  retained as the two ensemble grouping dimensions.
- `R/plot_erf.R`: plot the SSP5-8.5 NetERF covariate.

`R/climate_common.R` contains shared data and transformation helpers, while
`R/smooth_helpers.R` contains the ERF alignment, posterior summaries and
day-scale smooth plotting shared by UKCP18 and CORDEX. `R/smooth_cordex_index.R`
is an internal implementation sourced by
`R/smooth_cordex.R`; do not run it directly.

## Ensemble smoother outputs retained

Both final ensemble smoothers retain these distinct outputs for each index:

1. MAP smoothed trend plot.
2. Posterior predictive plot using the Laplace approximation.
3. Posterior predictive plot using importance sampling.
4. One importance-weight diagnostic per season.
5. A five-season cumulative-weight diagnostic.
6. Initial-value, MAP-parameter and 2080 posterior-summary CSV tables.
7. Per-index and combined RDS result objects.

Exact compatibility copies of the same `wsmax` figure were removed. The
remaining `wsmax` filenames use `_days` where the plotted values have been
converted back from cube-root scale.

## Default scope

The maintained scripts default to London and the five periods Annual, DJF,
MAM, JJA and SON. The settings near the top of each runnable script can be
changed when other regions or a smaller development run are needed. The final
ensemble smoothers use 1,000 posterior predictive samples by default.
