# =============================================================================
# Interval diagnostics for the fixed S&P 500 stock universe
#
# The rolling-result directory is used only to recover:
#   1. the exact stock universe from selected_stock_manifest.csv;
#   2. the original sample range and return multiplier from
#      rolling_analysis_design.csv.
#
# No white-noise test is rerun. No method, rolling midpoint, ranking, or
# top-variable selection is used. For a user-specified calendar interval the
# script produces:
#   1. a time plot of all stock returns;
#   2. one lag-q autocorrelation for each stock;
#   3. all p^2 directed lag-q correlations, with a cross-stock-only plot;
#   4. descriptive summaries of the two correlation groups.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# ----------------------------------------------------------------------------
# Utilities and reconstruction of the original stock universe
# ----------------------------------------------------------------------------

save_plot_pair <- function(plot, stem, width, height, dpi = 180L) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    paste0(stem, ".png"), plot,
    width = width, height = height, dpi = dpi, bg = "white"
  )
  ggsave(
    paste0(stem, ".pdf"), plot,
    width = width, height = height, device = "pdf"
  )
  invisible(plot)
}

read_interval_context <- function(
    results_dir,
    stock_file = "./realdata/kp9ncdlcmmbmoaby.csv") {
  required <- c(
    "selected_stock_manifest.csv",
    "rolling_analysis_design.csv"
  )
  missing <- required[!file.exists(file.path(results_dir, required))]
  if (length(missing)) {
    stop("The result directory is missing files: ", paste(missing, collapse = ", "))
  }
  if (!file.exists(stock_file)) stop("Original stock document not found:", stock_file)
  
  manifest <- fread(file.path(results_dir, "selected_stock_manifest.csv"))
  design <- fread(file.path(results_dir, "rolling_analysis_design.csv"))
  if (!nrow(manifest) || !nrow(design)) {
    stop("The stock list or rolling design file is empty.")
  }
  if (!"PERMNO" %in% names(manifest)) {
    stop("selected_stock_manifest.csv is missing PERMNO.")
  }
  
  universe_start <- as.IDate(design$universe_start[1L])
  universe_end <- as.IDate(design$universe_end[1L])
  return_multiplier <- as.numeric(design$return_multiplier[1L])
  if (anyNA(c(universe_start, universe_end, return_multiplier))) {
    stop("The sample range or return multiple of rolling_analysis_design.csv is invalid.")
  }
  
  selected_permnos <- sort(unique(as.integer(manifest$PERMNO)))
  raw <- fread(
    stock_file,
    select = c("DlyCalDt", "PERMNO", "DlyRet"),
    showProgress = TRUE
  )
  raw[, `:=`(
    DlyCalDt = as.IDate(DlyCalDt),
    PERMNO = as.integer(PERMNO)
  )]
  raw <- raw[
    PERMNO %in% selected_permnos &
      DlyCalDt >= universe_start & DlyCalDt <= universe_end
  ]
  if (anyDuplicated(raw[, .(DlyCalDt, PERMNO)])) {
    stop("The raw data contains duplicate DlyCalDt x PERMNO.")
  }
  
  wide <- dcast(raw, DlyCalDt ~ PERMNO, value.var = "DlyRet")
  setorder(wide, DlyCalDt)
  permno_columns <- as.character(selected_permnos)
  missing_permnos <- setdiff(permno_columns, names(wide))
  if (length(missing_permnos)) {
    stop("The raw data is missing the selected stock:", paste(missing_permnos, collapse = ", "))
  }
  setcolorder(wide, c("DlyCalDt", permno_columns))
  
  X <- as.matrix(wide[, ..permno_columns])
  storage.mode(X) <- "double"
  X <- X * return_multiplier
  colnames(X) <- paste0("PERMNO_", permno_columns)
  if (anyNA(X) || any(!is.finite(X))) {
    stop("The reconstructed common stock return matrix contains NA/Inf.")
  }
  
  if ("start_ticker" %in% names(manifest)) {
    ticker <- as.character(manifest$start_ticker)
  } else {
    ticker <- rep(NA_character_, nrow(manifest))
  }
  ticker_map <- data.table(
    PERMNO = as.integer(manifest$PERMNO),
    ticker = ticker
  )
  ticker_map[
    is.na(ticker) | !nzchar(ticker),
    ticker := paste0("PERMNO ", PERMNO)
  ]
  ticker_map[, variable := paste0("PERMNO_", PERMNO)]
  ticker_map[, display_label := paste0(ticker, " [", PERMNO, "]")]
  
  expected_p <- suppressWarnings(as.integer(design$p_common[1L]))
  if (!is.na(expected_p) && ncol(X) != expected_p) {
    stop("Refactor the number of stocks to ", ncol(X), ", but the design document records it as ", expected_p, ".")
  }
  
  list(
    results_dir = normalizePath(results_dir, mustWork = TRUE),
    manifest = manifest,
    design = design,
    universe_start = universe_start,
    universe_end = universe_end,
    dates = as.IDate(wide$DlyCalDt),
    X = X,
    ticker_map = ticker_map,
    return_multiplier = return_multiplier
  )
}


# ----------------------------------------------------------------------------
# Plots and lagged correlations
# ----------------------------------------------------------------------------

plot_all_stock_returns <- function(X, dates, title = NULL) {
  p <- ncol(X)
  colors <- adjustcolor(hcl.colors(p, "Dark 3"), alpha.f = 0.55)
  long <- data.table(
    date = rep(as.Date(dates), times = p),
    variable = rep(colnames(X), each = nrow(X)),
    value = as.vector(X)
  )
  
  ggplot(long, aes(date, value, group = variable, color = variable)) +
    geom_line(linewidth = 0.25, alpha = 0.58) +
    scale_color_manual(values = setNames(colors, colnames(X)), guide = "none") +
    scale_x_date(expand = expansion(mult = c(0.01, 0.045))) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey35") +
    labs(title = title, x = "Date", y = "Daily return (%)") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 7, r = 22, b = 7, l = 7)
    )
}

compute_lag_correlation_tables <- function(X, ticker_map, lag = 1L) {
  lag <- as.integer(lag)
  if (length(lag) != 1L || is.na(lag) || lag < 1L || lag >= nrow(X)) {
    stop("lag must be a positive integer less than the interval sample size.")
  }
  
  fit <- stats::acf(
    X,
    lag.max = lag,
    type = "correlation",
    plot = FALSE,
    demean = TRUE
  )
  lag_matrix <- fit$acf[lag + 1L, , , drop = TRUE]
  p <- ncol(X)
  if (!identical(dim(lag_matrix), c(p, p))) {
    lag_matrix <- matrix(lag_matrix, nrow = p, ncol = p)
  }
  
  map <- ticker_map[match(colnames(X), variable)]
  if (anyNA(map$variable)) stop("The stock name cannot be fully matched with ticker_map.")
  
  stock_acf <- data.table(
    stock_index = seq_len(p),
    variable = colnames(X),
    PERMNO = map$PERMNO,
    ticker = map$ticker,
    display_label = map$display_label,
    lag = lag,
    acf = diag(lag_matrix)
  )
  
  all_pairs <- as.data.table(
    which(matrix(TRUE, nrow = p, ncol = p), arr.ind = TRUE)
  )
  setnames(all_pairs, c("row", "col"), c("current_index", "lagged_index"))
  all_pairs[, `:=`(
    pair_index = seq_len(.N),
    lag = lag,
    current_variable = colnames(X)[current_index],
    lagged_variable = colnames(X)[lagged_index],
    current_ticker = map$ticker[current_index],
    lagged_ticker = map$ticker[lagged_index],
    correlation = lag_matrix[cbind(current_index, lagged_index)],
    is_autocorrelation = current_index == lagged_index
  )]
  all_pairs[, abs_correlation := abs(correlation)]
  setcolorder(
    all_pairs,
    c(
      "pair_index", "lag", "current_index", "lagged_index",
      "current_variable", "lagged_variable",
      "current_ticker", "lagged_ticker",
      "correlation", "abs_correlation", "is_autocorrelation"
    )
  )
  
  list(stock_acf = stock_acf, all_pairs = all_pairs)
}

plot_stock_acf_scatter <- function(stock_acf, n_obs) {
  pointwise_bound <- 1.96 / sqrt(n_obs)
  
  ggplot(stock_acf, aes(stock_index, acf)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey30") +
    geom_hline(
      yintercept = c(-pointwise_bound, pointwise_bound),
      linewidth = 0.45,
      linetype = "dashed",
      color = "#B2182B"
    ) +
    geom_point(size = 1.15, alpha = 0.72, color = "#2166AC") +
    scale_x_continuous(breaks = NULL) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      x = "Stocks",
      y = paste0("Lag-", unique(stock_acf$lag), " sample ACF")
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

plot_cross_stock_lag_scatter <- function(all_pairs, n_obs) {
  cross_pairs <- copy(all_pairs[is_autocorrelation == FALSE])
  cross_pairs[, cross_pair_index := seq_len(.N)]
  pointwise_bound <- 1.96 / sqrt(n_obs)
  
  ggplot(cross_pairs, aes(cross_pair_index, correlation)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey30") +
    geom_hline(
      yintercept = c(-pointwise_bound, pointwise_bound),
      linewidth = 0.4,
      linetype = "dashed",
      color = "#B2182B"
    ) +
    geom_point(size = 0.22, alpha = 0.42, color = "#2166AC") +
    scale_x_continuous(breaks = NULL) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      x = "Stock pairs",
      y = paste0(
        "Lag-", unique(cross_pairs$lag), " cross-correlation"
      )
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

summarize_lag_correlations <- function(all_pairs, n_obs) {
  pointwise_bound <- 1.96 / sqrt(n_obs)
  
  all_pairs[, .(
    n_values = .N,
    mean_correlation = mean(correlation),
    median_correlation = median(correlation),
    mean_abs_correlation = mean(abs_correlation),
    median_abs_correlation = median(abs_correlation),
    q90_abs_correlation = as.numeric(quantile(abs_correlation, 0.90)),
    q95_abs_correlation = as.numeric(quantile(abs_correlation, 0.95)),
    q99_abs_correlation = as.numeric(quantile(abs_correlation, 0.99)),
    max_abs_correlation = max(abs_correlation),
    pointwise_bound = pointwise_bound,
    n_above_pointwise_bound = sum(abs_correlation > pointwise_bound),
    fraction_above_pointwise_bound = mean(abs_correlation > pointwise_bound)
  ), by = .(
    correlation_type = fifelse(
      is_autocorrelation,
      "same_stock_autocorrelation",
      "cross_stock_lag_correlation"
    )
  )]
}


# ----------------------------------------------------------------------------
# Main interval analysis
# ----------------------------------------------------------------------------

analyze_return_interval <- function(
    start_date,
    end_date,
    lag = 1L,
    results_dir = "./realdata_sp500_gfc_2005_2010_rolling_pvalues",
    stock_file = "./realdata/kp9ncdlcmmbmoaby.csv",
    output_dir = NULL,
    dpi = 180L,
    save_all_pair_pdf = FALSE,
    returns_title = NULL) {
  context <- read_interval_context(results_dir, stock_file)
  start_date <- as.IDate(start_date)
  end_date <- as.IDate(end_date)
  if (length(start_date) != 1L || length(end_date) != 1L ||
      is.na(start_date) || is.na(end_date) || start_date > end_date) {
    stop("start_date and end_date must be valid dates, and start_date <= end_date.")
  }
  if (start_date < context$universe_start || end_date > context$universe_end) {
    stop(
      "The specified interval must lie within the original experimental range ", context$universe_start,
      " to ", context$universe_end, "."
    )
  }
  
  keep <- context$dates >= start_date & context$dates <= end_date
  if (!any(keep)) stop("No trading day data in the specified interval.")
  X_interval <- context$X[keep, , drop = FALSE]
  dates_interval <- context$dates[keep]
  if (nrow(X_interval) <= as.integer(lag)) {
    stop("The specified interval is too short to calculate this lag.")
  }
  
  actual_start <- min(dates_interval)
  actual_end <- max(dates_interval)
  if (is.null(output_dir)) {
    output_dir <- file.path(
      results_dir,
      "interval_diagnostics",
      paste0(actual_start, "_to_", actual_end, "_lag", as.integer(lag))
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  metadata <- data.table(
    requested_start = start_date,
    requested_end = end_date,
    actual_start = actual_start,
    actual_end = actual_end,
    n_trading_days = nrow(X_interval),
    n_stocks = ncol(X_interval),
    lag = as.integer(lag),
    return_multiplier = context$return_multiplier,
    stock_universe_start = context$universe_start,
    stock_universe_end = context$universe_end,
    results_dir = context$results_dir
  )
  fwrite(metadata, file.path(output_dir, "interval_metadata.csv"))
  
  returns_plot <- plot_all_stock_returns(
    X_interval, dates_interval, title = returns_title
  )
  save_plot_pair(
    returns_plot,
    file.path(output_dir, "all_stock_returns"),
    width = 11.5, height = 6.8, dpi = dpi
  )
  
  lag_results <- compute_lag_correlation_tables(
    X_interval, context$ticker_map, lag = lag
  )
  fwrite(
    lag_results$stock_acf,
    file.path(output_dir, paste0("lag", lag, "_stock_acf.csv"))
  )
  fwrite(
    lag_results$all_pairs,
    file.path(output_dir, paste0("lag", lag, "_all_pair_correlations.csv"))
  )
  
  correlation_summary <- summarize_lag_correlations(
    lag_results$all_pairs, nrow(X_interval)
  )
  fwrite(
    correlation_summary,
    file.path(output_dir, paste0("lag", lag, "_correlation_summary.csv"))
  )
  
  acf_scatter <- plot_stock_acf_scatter(
    lag_results$stock_acf, nrow(X_interval)
  )
  save_plot_pair(
    acf_scatter,
    file.path(output_dir, paste0("lag", lag, "_stock_acf_scatter")),
    width = 11, height = 6.2, dpi = dpi
  )
  
  cross_scatter <- plot_cross_stock_lag_scatter(
    lag_results$all_pairs, nrow(X_interval)
  )
  cross_stem <- file.path(
    output_dir, paste0("lag", lag, "_all_pair_correlation_scatter")
  )
  ggsave(
    paste0(cross_stem, ".png"), cross_scatter,
    width = 11.5, height = 6.5, dpi = dpi, bg = "white"
  )
  if (isTRUE(save_all_pair_pdf)) {
    ggsave(
      paste0(cross_stem, ".pdf"), cross_scatter,
      width = 11.5, height = 6.5, device = "pdf"
    )
  }
  
  message(
    "Interval diagnosis completed:", actual_start, " to ", actual_end,
    ";T=", nrow(X_interval), "; p=", ncol(X_interval),
    ";Self ACF=", ncol(X_interval),
    ";Cross-stock pair=", ncol(X_interval) * (ncol(X_interval) - 1L),
    ";Results director:", normalizePath(output_dir, mustWork = TRUE)
  )
  
  invisible(list(
    metadata = metadata,
    stock_acf = lag_results$stock_acf,
    all_pair_correlations = lag_results$all_pairs,
    correlation_summary = correlation_summary,
    plots = list(
      returns = returns_plot,
      stock_acf = acf_scatter,
      cross_stock = cross_scatter
    )
  ))
}


# Edit the two dates below for each interval of interest.
diagnostics <- analyze_return_interval(
  start_date = "2005-02-23",
  end_date = "2006-06-09",
  lag = 1,
  results_dir = "./result/realdata_sp500_gfc_2005_2010_rolling_pvalues",
  stock_file = "./realdata/kp9ncdlcmmbmoaby.csv"
)
