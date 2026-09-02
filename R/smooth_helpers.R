# Shared posterior and plotting helpers for UKCP18 and CORDEX smoothers.

to_num <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

get_erf_vector <- function(years) {
  if (!exists("ERF585")) {
    stop("ERF585 was not loaded by data(SSP585data).")
  }

  erf <- ERF585
  if (is.vector(erf) && !is.data.frame(erf)) {
    if (length(erf) != length(years)) {
      stop("ERF585 is a vector but its length does not match the data years.")
    }
    return(as.numeric(erf))
  }

  erf <- as.data.frame(erf)
  year_col <- grep("^year$|year", names(erf), ignore.case = TRUE, value = TRUE)[1]
  net_col <- grep("^NetERF$|net.*erf|erf", names(erf), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(year_col) || is.na(net_col)) {
    stop("Cannot identify Year and NetERF columns in ERF585.")
  }

  Xt <- to_num(erf[[net_col]])[match(years, to_num(erf[[year_col]]))]
  if (length(Xt) != length(years) || anyNA(Xt)) {
    stop(
      "ERF alignment failed. Missing ERF for years: ",
      paste(years[is.na(Xt)], collapse = ", ")
    )
  }
  Xt
}

weighted_quantile <- function(
    x,
    w = NULL,
    probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
) {
  if (is.null(w)) w <- rep(1, length(x))
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0L) return(rep(NA_real_, length(probs)))

  w <- w / sum(w)
  ord <- order(x)
  x <- x[ord]
  cw <- cumsum(w[ord])
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

get_weights <- function(pps) {
  n <- if (is.matrix(pps$Obs)) ncol(pps$Obs) else dim(pps$Obs)[3]
  w <- if (!is.null(pps$Weights) && !is.null(pps$Weights$w)) {
    pps$Weights$w
  } else {
    rep(1 / n, n)
  }
  w / sum(w)
}

extract_obs_matrix <- function(pps) {
  obs <- pps$Obs
  if (length(dim(obs)) == 3L) obs <- obs[, 1, , drop = TRUE]
  obs
}

wsmax_to_days <- function(x) {
  pmax(x, 0)^3
}

plot_wsmax_day_panel <- function(
    Data,
    Smooth = NULL,
    Samples = NULL,
    Groups = NULL,
    DatColours = c("black", "coral3"),
    DatTypes = c(1, 1),
    PredColours = c("darkblue", "darkgoldenrod"),
    alpha = c(0.6, 0.1),
    PlotMu0 = TRUE,
    PlotConsensus = TRUE,
    SmoothPIs = FALSE,
    Units = "days",
    plot.title = "",
    Legend = TRUE,
    LegPos = "topleft",
    EnsTransp = 0.2
) {
  if (!xor(is.null(Smooth), is.null(Samples))) {
    stop("Supply exactly one of Smooth or Samples.")
  }

  years <- Data[, 1]
  PlotData <- Data
  PlotData[, -1] <- as.data.frame(
    wsmax_to_days(as.matrix(Data[, -1, drop = FALSE]))
  )

  pred_cols <- ci_cols <- t(col2rgb(PredColours)) / 255
  pred_cols <- rgb(pred_cols, alpha = alpha[1])
  ci_cols <- rgb(ci_cols, alpha = alpha[2])

  if (is.null(Samples)) {
    predictions <- dlm.ObsPred(Smooth)
    mu0_model <- Smooth$Smooth$s[-1, 1]
    se0_model <- predictions$SE[, 1]
    mu0_hat <- wsmax_to_days(mu0_model)
    ci_limits <- rbind(
      wsmax_to_days(mu0_model - 1.96 * se0_model),
      wsmax_to_days(mu0_model + 1.96 * se0_model)
    )

    if (PlotConsensus) {
      if (!("Consensus" %in% names(Smooth$Model))) {
        stop("Cannot plot ensemble consensus without Smooth$Model$Consensus.")
      }
      consensus_matrix <- matrix(Smooth$Model$Consensus, nrow = 1)
      consensus_model <- as.vector(Smooth$Smooth$s[-1, ] %*% t(consensus_matrix))
      consensus <- wsmax_to_days(consensus_model)
    }
  } else {
    weights <- NULL
    if (!is.null(Samples$Weights) && !is.null(Samples$Weights$w)) {
      weights <- Samples$Weights$w / sum(Samples$Weights$w)
    }

    mu0_samples <- wsmax_to_days(Samples$States[, , 1])
    obs_samples <- wsmax_to_days(extract_obs_matrix(Samples))

    if (is.null(weights)) {
      mu0_hat <- colMeans(mu0_samples, na.rm = TRUE)
      ci_limits <- apply(
        obs_samples,
        1,
        quantile,
        probs = c(0.025, 0.975),
        na.rm = TRUE
      )
    } else {
      mu0_hat <- apply(
        mu0_samples,
        2,
        weighted.mean,
        w = weights,
        na.rm = TRUE
      )
      ci_limits <- apply(
        obs_samples,
        1,
        function(x) weighted_quantile(x, weights, c(0.025, 0.975))
      )
    }

    if (PlotConsensus) {
      if (!("Consensus" %in% names(Samples$Model))) {
        stop("Cannot plot ensemble consensus without Samples$Model$Consensus.")
      }
      consensus_model <- apply(
        Samples$States,
        1:2,
        function(x) as.numeric(x %*% Samples$Model$Consensus)
      )
      consensus_days <- wsmax_to_days(consensus_model)
      consensus <- if (is.null(weights)) {
        colMeans(consensus_days, na.rm = TRUE)
      } else {
        apply(consensus_days, 2, weighted.mean, w = weights, na.rm = TRUE)
      }
    }

    if (SmoothPIs) {
      ci_limits <- t(apply(
        ci_limits,
        1,
        function(x) {
          position <- seq_along(x)
          pmax(loess(x ~ position)$fitted, 0)
        }
      ))
    }
  }

  y_values <- c(
    as.numeric(as.matrix(PlotData[, -1, drop = FALSE])),
    mu0_hat,
    ci_limits,
    if (PlotConsensus) consensus else NULL
  )
  ylim <- range(y_values[is.finite(y_values)], na.rm = TRUE)
  padding <- 0.05 * diff(ylim)
  if (!is.finite(padding) || padding == 0) padding <- 1
  ylim <- ylim + c(-padding, padding)

  PlotEnsTS(
    PlotData,
    Colours = DatColours,
    Types = DatTypes,
    Groups = Groups,
    EnsTransp = EnsTransp,
    Units = Units,
    ylim = ylim
  )

  year_limits <- par("usr")[1:2]
  abline(h = pretty(par("usr")[3:4], n = 10), col = grey(0.8), lty = 2)
  abline(v = pretty(year_limits, n = diff(year_limits) / 10), col = grey(0.8), lty = 2)
  if (PlotConsensus) lines(years, consensus, lwd = 5, col = pred_cols[2])
  polygon(
    c(years, rev(years)),
    c(ci_limits[1, ], rev(ci_limits[2, ])),
    border = NA,
    col = ci_cols[1]
  )
  if (PlotMu0) lines(years, mu0_hat, lwd = 5, col = pred_cols[1])
  title(main = plot.title)

  if (Legend) {
    legend_colours <- c(DatColours[c(1, 2, 2)], PredColours[c(1, 2)], NA)
    legend_widths <- c(5, 5, 1, 5, 5, NA)
    legend_fill <- c(rep(NA, 5), ci_cols[1])
    legend_text <- c(
      "Observations",
      "Ensemble mean",
      "Ensemble members",
      expression(hat(mu)[0](t)),
      "Ensemble consensus",
      expression("95% interval for " * Y[0](t))
    )
    wanted <- rep(TRUE, length(legend_colours))
    if (!PlotConsensus) wanted[5] <- FALSE
    if (!PlotMu0) wanted[4] <- FALSE
    if (ncol(PlotData) < 3) wanted[c(2, 3, 5)] <- FALSE
    legend(
      LegPos,
      col = legend_colours[wanted],
      lwd = legend_widths[wanted],
      fill = legend_fill[wanted],
      border = NA,
      ncol = 2,
      bg = "white",
      legend = legend_text[wanted]
    )
  }

  result <- data.frame(
    Time = years,
    Mu0 = mu0_hat,
    LL95 = ci_limits[1, ],
    UL95 = ci_limits[2, ]
  )
  if (PlotConsensus) result$Consensus <- consensus
  invisible(result)
}
