# =============================================================================
# Rolling-window p-values for the S&P 500 GFC versus COVID-19 application
#
# This script:
#   * runs selected methods one method at a time;
#   * parallelizes windows within the current method on Unix-like servers;
#   * checkpoints after every batch of windows;
#   * writes one CSV and one PNG/PDF figure per method;
#   * rebuilds a combined CSV, summary, and all-method figure from every
#     completed method file already present in the output directory.
#
# Recommended first run (Proposed only, defaults n_w=252, s=5, alpha=.10):
#   Rscript realdata_sp500_rolling_pvalues.R
#
# Run another method into the same output directory:
#   ROLLING_METHODS=Tsay2020 Rscript realdata_sp500_rolling_pvalues.R
#   ROLLING_METHODS=Li2019 Rscript realdata_sp500_rolling_pvalues.R
#   ROLLING_METHODS=Feng2022_T_FC Rscript realdata_sp500_rolling_pvalues.R
#   ROLLING_METHODS=Chen2025_rho Rscript realdata_sp500_rolling_pvalues.R
#
# Override the design explicitly:
#   ROLLING_WINDOW_SIZE=252 ROLLING_STEP_SIZE=5 ROLLING_ALPHA=0.10 \
#   ROLLING_K_VALUES=1,2 ROLLING_WORKERS=4 ROLLING_METHODS=our_resampling \
#   Rscript realdata_sp500_rolling_pvalues.R
#
# Multiple selected methods are allowed, but are deliberately run sequentially;
# only the windows inside the current method are parallelized.
# =============================================================================

suppressPackageStartupMessages({
  library(Rcpp)
  library(HDTSA)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

if (!file.exists("./src/tool.cpp")) {
  stop("Please run from the project root directory; ./src/tool.cpp cannot be found in the current directory.")
}

sourceCpp("./src/tool.cpp")
source("./src/test_fun_statistic.R")
source("./src/wnSVAR_statistic.R")
source("./src/HDWHtest.R")


# ----------------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------------

ALL_ROLLING_METHODS <- c(
  "Li2019",
  "Li2019_correction",
  "Chang2017",
  "Wang2022",
  "Tsay2020",
  "our",
  "our_resampling",
  "Feng2022_T_SUM",
  "Feng2022_T_MAX",
  "Feng2022_T_FC",
  "Chen2025_rho",
  "Chen2025_tau",
  "Chen2025_D",
  "Chen2025_R",
  "Chen2025_tau_star"
)

METHOD_LABELS <- c(
  Li2019 = "LLYY",
  Li2019_correction = "LLYY correction",
  Chang2017 = "CYZ",
  Wang2022 = "WKX",
  Tsay2020 = "Tsay",
  our = "Proposed (asymptotic)",
  our_resampling = "Proposed",
  Feng2022_T_SUM = "FLM SUM",
  Feng2022_T_MAX = "FLM MAX",
  Feng2022_T_FC = "FLM FC",
  Chen2025_rho = "CSF L-rho",
  Chen2025_tau = "CSF L-tau",
  Chen2025_D = "CSF L-D",
  Chen2025_R = "CSF L-R",
  Chen2025_tau_star = "CSF L-tau*"
)

# Plotmath versions used only for figure labels.  These render the method
# names in the same mathematical style as the LaTeX tables, without requiring
# an additional LaTeX/ggtext dependency.
METHOD_MATH_LABELS <- c(
  Li2019 = "plain(LLYY)",
  Li2019_correction = "plain(LLYY)~plain('(correction)')",
  Chang2017 = "plain(CYZ)",
  Wang2022 = "plain(WKX)",
  Tsay2020 = "plain(Tsay)",
  our = "plain(Proposed)~plain('(asymptotic)')",
  our_resampling = "plain(Proposed)",
  Feng2022_T_SUM = "plain(FLM)[plain(SUM)]",
  Feng2022_T_MAX = "plain(FLM)[plain(MAX)]",
  Feng2022_T_FC = "plain(FLM)[plain(FC)]",
  Chen2025_rho = "plain(CSF)[L[rho]]",
  Chen2025_tau = "plain(CSF)[L[tau]]",
  Chen2025_D = "plain(CSF)[L[D]]",
  Chen2025_R = "plain(CSF)[L[R]]",
  Chen2025_tau_star = "plain(CSF)[L[tau^'*']]"
)

DEFAULT_COMBINED_PLOT_METHODS <- c(
  "our_resampling",
  "Li2019",
  "Chang2017",
  "Wang2022",
  "Tsay2020",
  "Feng2022_T_FC",
  # "Chen2025_rho",
  # "Chen2025_tau",
  "Chen2025_D",
  "Chen2025_R",
  "Chen2025_tau_star"
)

PERIODS <- data.table(
  period_id = c("gfc_2007_2010"),
  period_label = c("The global financial crisis period"),
  period_start = as.IDate(c("2007-01-01")),
  period_end = as.IDate(c("2010-12-31")),
  event_date = as.IDate(c("2008-09-15"))
)

env_integer <- function(name, default, lower = 1L) {
  text <- trimws(Sys.getenv(name, ""))
  if (!nzchar(text)) return(as.integer(default))
  value <- suppressWarnings(as.integer(text))
  if (length(value) != 1L || is.na(value) || value < lower) {
    stop(name, " must be an integer greater than or equal to ", lower, ".")
  }
  value
}

env_numeric <- function(name, default) {
  text <- trimws(Sys.getenv(name, ""))
  if (!nzchar(text)) return(as.numeric(default))
  value <- suppressWarnings(as.numeric(text))
  if (length(value) != 1L || !is.finite(value)) {
    stop(name, " must be a finite numeric value.")
  }
  value
}

env_logical <- function(name, default) {
  text <- tolower(trimws(Sys.getenv(name, "")))
  if (!nzchar(text)) return(isTRUE(default))
  if (text %in% c("1", "true", "yes", "y")) return(TRUE)
  if (text %in% c("0", "false", "no", "n")) return(FALSE)
  stop(name, " must be one of TRUE/FALSE, 1/0, yes/no.")
}

env_csv <- function(name, default, allowed = NULL, integer = FALSE) {
  text <- trimws(Sys.getenv(name, ""))
  values <- if (nzchar(text)) {
    unique(trimws(strsplit(text, ",", fixed = TRUE)[[1L]]))
  } else {
    default
  }
  if (integer) {
    values <- suppressWarnings(as.integer(values))
    if (!length(values) || anyNA(values) || any(values < 1L)) {
      sstop(name, " must be a comma-separated positive integer.")
    }
    return(sort(unique(values)))
  }
  if (!is.null(allowed)) {
    unknown <- setdiff(values, allowed)
    if (length(unknown)) {
      stop(name, " Contains unknown method:", paste(unknown, collapse = ", "))
    }
  }
  values
}

window_size <- env_integer("ROLLING_WINDOW_SIZE", 252L, lower = 10L)
step_size <- env_integer("ROLLING_STEP_SIZE", 5L)
alpha <- env_numeric("ROLLING_ALPHA", 0.10)
if (alpha <= 0 || alpha >= 1) stop("ROLLING_ALPHA must be in (0,1).")
K_values <- env_csv("ROLLING_K_VALUES", c(1L, 2L), integer = TRUE)

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!length(detected_cores) || is.na(detected_cores)) detected_cores <- 2L
default_workers <- max(1L, max(4L, detected_cores - 4L))
workers <- env_integer("ROLLING_WORKERS", default_workers)

alpha_tag <- gsub("\\.", "p", format(alpha, scientific = FALSE, trim = TRUE))
default_output_dir <- sprintf(
  "./realdata_sp500_rolling_results_n%d_s%d_a%s",
  window_size, step_size, alpha_tag
)

CONFIG <- list(
  stock_file = "./realdata/kp9ncdlcmmbmoaby.csv",
  output_dir = Sys.getenv("ROLLING_OUTPUT_DIR", default_output_dir),
  periods = copy(PERIODS),
  universe_start = min(PERIODS$period_start),
  universe_end = max(PERIODS$period_end),
  return_multiplier = env_numeric("ROLLING_RETURN_MULTIPLIER", 100),
  window_size = window_size,
  step_size = step_size,
  alpha = alpha,
  K_values = K_values,
  methods = env_csv(
    "ROLLING_METHODS", "our_resampling", allowed = ALL_ROLLING_METHODS
  ),
  plot_methods = env_csv(
    "ROLLING_PLOT_METHODS",
    DEFAULT_COMBINED_PLOT_METHODS,
    allowed = ALL_ROLLING_METHODS
  ),
  our_B = env_integer("ROLLING_OUR_B", 1000L),
  our_boot_type = env_integer("ROLLING_OUR_BOOT_TYPE", 1L),
  chang_B = env_integer("ROLLING_CHANG_B", 1000L),
  chang_pre = env_logical("ROLLING_CHANG_PRE", TRUE),
  chang_kernel = Sys.getenv("ROLLING_CHANG_KERNEL", "QS"),
  chang_control_PCA = list(),
  wang_B = env_integer("ROLLING_WANG_B", 1250L),
  parallel = env_logical("ROLLING_PARALLEL", TRUE),
  workers = workers,
  checkpoint_batch_size = env_integer(
    "ROLLING_BATCH_SIZE", max(8L, 2L * workers)
  ),
  seed = env_integer("ROLLING_SEED", 20260821L),
  plot_width = env_numeric("ROLLING_PLOT_WIDTH", 12),
  plot_height = env_numeric("ROLLING_PLOT_HEIGHT", 7),
  plot_dpi = env_integer("ROLLING_PLOT_DPI", 180L)
)

if (CONFIG$window_size <= max(CONFIG$K_values)) {
  stop("ROLLING_WINDOW_SIZE must be greater than the maximum K.")
}
if (!CONFIG$our_boot_type %in% 1:3) {
  stop("ROLLING_OUR_BOOT_TYPE can only be 1, 2, or 3.")
}
CONFIG$chang_kernel <- match.arg(CONFIG$chang_kernel, c("QS", "Par", "Bart"))


# ----------------------------------------------------------------------------
# 2. Balanced universe and rolling-window construction
# ----------------------------------------------------------------------------

prepare_balanced_universe <- function(config = CONFIG) {
  if (!file.exists(config$stock_file)) {
    stop("Stock file not found:", config$stock_file)
  }
  
  message(
    "Construct a common stock pool:", config$universe_start, " to ",
    config$universe_end, "……"
  )
  z_all <- fread(
    config$stock_file,
    select = c(
      "MbrStartDt", "MbrEndDt", "PERMNO", "Ticker", "PrimaryExch",
      "DlyCalDt", "DlyCap", "DlyRet"
    ),
    showProgress = TRUE
  )
  z_all <- z_all[
    DlyCalDt >= config$universe_start & DlyCalDt <= config$universe_end
  ]
  if (!nrow(z_all)) stop("No stock records within the common screening span.")
  
  trading_dates <- sort(unique(z_all$DlyCalDt))
  first_date <- min(trading_dates)
  audit <- z_all[, .(
    membership_start = min(MbrStartDt),
    membership_end = max(MbrEndDt),
    continuous_membership = any(
      MbrStartDt <= config$universe_start & MbrEndDt >= config$universe_end
    ),
    n_dates = uniqueN(DlyCalDt),
    n_finite_returns = sum(is.finite(DlyRet)),
    start_cap = {
      value <- DlyCap[DlyCalDt == first_date]
      if (length(value) && is.finite(value[[1L]])) value[[1L]] else NA_real_
    },
    start_ticker = {
      value <- Ticker[DlyCalDt == first_date]
      if (length(value)) value[[1L]] else NA_character_
    },
    start_exchange = {
      value <- PrimaryExch[DlyCalDt == first_date]
      if (length(value)) value[[1L]] else NA_character_
    }
  ), by = PERMNO]
  
  audit[, complete_panel :=
          continuous_membership &
          n_dates == length(trading_dates) &
          n_finite_returns == length(trading_dates)]
  selected_permnos <- sort(audit[complete_panel == TRUE, PERMNO])
  if (length(selected_permnos) < 2L) {
    stop("Within the common screening span, there are fewer than two stocks with complete and continuous data.")
  }
  
  z <- z_all[
    PERMNO %in% selected_permnos,
    .(DlyCalDt, PERMNO, DlyRet)
  ]
  if (anyDuplicated(z[, .(DlyCalDt, PERMNO)])) {
    stop("The common stock pool contains duplicate DlyCalDt × PERMNO records.")
  }
  wide <- dcast(z, DlyCalDt ~ PERMNO, value.var = "DlyRet")
  setorder(wide, DlyCalDt)
  X <- as.matrix(wide[, -"DlyCalDt"])
  storage.mode(X) <- "double"
  colnames(X) <- paste0("PERMNO_", names(wide)[-1L])
  X <- X * config$return_multiplier
  if (anyNA(X) || any(!is.finite(X))) {
    stop("The common payoff matrix still contains NA/Inf.")
  }
  
  list(
    X = X,
    dates = as.IDate(wide$DlyCalDt),
    audit = audit[order(PERMNO)],
    selected_permnos = selected_permnos
  )
}

extract_period <- function(universe, period_row) {
  keep <- universe$dates >= period_row$period_start &
    universe$dates <= period_row$period_end
  if (!any(keep)) stop("There are no trading days within period ", period_row$period_id, ".")
  list(
    X = universe$X[keep, , drop = FALSE],
    dates = universe$dates[keep]
  )
}

make_window_index <- function(dates, period_row, config = CONFIG) {
  n_period <- length(dates)
  if (n_period < config$window_size) {
    stop(
      "period ", period_row$period_id, " only have ", n_period,
      " trading days, less than the window length ", config$window_size, "。"
    )
  }
  starts <- seq.int(
    from = 1L,
    to = n_period - config$window_size + 1L,
    by = config$step_size
  )
  ends <- starts + config$window_size - 1L
  mids <- floor((starts + ends) / 2)
  data.table(
    period_id = period_row$period_id,
    period_label = period_row$period_label,
    window_index = seq_along(starts),
    window_id = sprintf("%s_w%04d", period_row$period_id, seq_along(starts)),
    row_start = starts,
    row_end = ends,
    window_start = dates[starts],
    window_end = dates[ends],
    window_mid = dates[mids],
    n_window = config$window_size
  )
}

center_panel <- function(X) sweep(X, 2L, colMeans(X), FUN = "-")

stable_seed <- function(base_seed, ...) {
  text <- paste(..., collapse = "|")
  ints <- utf8ToInt(text)
  increment <- if (length(ints)) {
    sum((seq_along(ints) %% 1009L + 1L) * ints)
  } else 0
  as.integer((as.double(base_seed) + increment) %% 2147483646) + 1L
}

make_chen_input <- function(X, seed) {
  tied <- which(vapply(
    seq_len(ncol(X)),
    function(j) anyDuplicated(X[, j]) > 0L,
    logical(1L)
  ))
  if (!length(tied)) return(X)
  set.seed(seed)
  for (j in tied) {
    scale_j <- max(sd(X[, j]), 1)
    X[, j] <- X[, j] + runif(nrow(X), -1, 1) * scale_j * 1e-10
  }
  X
}

# Tsay-style time plot of all stock-return series in one period.  X is the
# uncentered percentage-return matrix produced by prepare_balanced_universe().
save_period_return_plot <- function(
    X, dates, period_id, period_label, ylim, config = CONFIG) {
  figure_dir <- file.path(config$output_dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  colors <- grDevices::adjustcolor(
    grDevices::hcl.colors(ncol(X), palette = "Dark 3"),
    alpha.f = 0.62
  )
  y_label <- if (isTRUE(all.equal(config$return_multiplier, 100))) {
    "Daily return (%)"
  } else {
    "Daily return"
  }
  
  draw_plot <- function() {
    graphics::par(mar = c(4.5, 4.8, 3.2, 1.2), las = 1)
    graphics::plot(
      as.Date(dates), X[, 1L],
      type = "n",
      ylim = ylim,
      xlab = "Year",
      ylab = y_label,
      main = period_label,
      # main = paste0(
      #   period_label, ": all ", ncol(X), " stock-return series"
      # ),
      bty = "o"
    )
    graphics::matlines(
      as.Date(dates), X,
      type = "l",
      lty = 1,
      lwd = 0.38,
      col = colors
    )
    graphics::abline(h = 0, col = "#555555", lwd = 0.45)
    graphics::box(lwd = 1)
  }
  
  stem <- file.path(figure_dir, paste0("all_stock_returns_", period_id))
  grDevices::png(
    paste0(stem, ".png"),
    width = 2200,
    height = 1350,
    res = 190,
    bg = "white"
  )
  draw_plot()
  grDevices::dev.off()
  
  grDevices::pdf(
    paste0(stem, ".pdf"),
    width = 11.5,
    height = 7.1,
    useDingbats = FALSE
  )
  draw_plot()
  grDevices::dev.off()
  invisible(stem)
}


# ----------------------------------------------------------------------------
# 3. One-method test dispatcher
# ----------------------------------------------------------------------------

result_rows <- function(method, K_values, reject, p_value, statistic, cv) {
  data.table(
    method = method,
    K = as.integer(K_values),
    reject = as.logical(reject),
    p_value = pmin(1, pmax(0, as.numeric(p_value))),
    statistic = as.numeric(statistic),
    critical_value = as.numeric(cv)
  )
}

gumbel_upper_p <- function(score, multiplier = 1 / sqrt(pi)) {
  lambda <- multiplier * exp(-as.numeric(score) / 2)
  -expm1(-lambda)
}

# CYZ implementation used by the user's earlier HDTSA simulations.  This is
# the old L_inf path expanded explicitly, rather than a call to a newer
# HDTSA::WN_test wrapper.  The *_fun object names are intentional because the
# underlying C/R helpers are unexported namespace objects.
chang_wn_old_with_cv <- function(
    Y,
    lag.k = 2L,
    B = 1000L,
    kernel.type = c("QS", "Par", "Bart"),
    pre = FALSE,
    alpha = 0.05,
    control.PCA = list()) {
  n <- nrow(Y)
  p <- ncol(Y)
  kernel.type <- match.arg(kernel.type)
  ken_type <- switch(kernel.type, QS = 1L, Par = 2L, Bart = 3L)
  
  if (isTRUE(pre)) {
    con <- list(
      lag.k = 5L,
      thresh = FALSE,
      delta = 2 * sqrt(log(ncol(Y)) / nrow(Y)),
      opt = 1L
    )
    con[(namc <- names(control.PCA))] <- control.PCA
    segmentTS_fun <- get("segmentTS", asNamespace("HDTSA"))
    X_pre <- segmentTS_fun(
      Y,
      lag.k = con$lag.k,
      thresh = con$thresh,
      delta = con$delta,
      opt = con$opt,
      control = control.PCA
    )
    Y <- X_pre$Z
  }
  
  WN_teststat_fun <- get("WN_teststatC", asNamespace("HDTSA"))
  WN_ftC_fun <- get("WN_ftC", asNamespace("HDTSA"))
  bandwith_fun <- get("bandwith", asNamespace("HDTSA"))
  WN_bootc_fun <- get("WN_bootc", asNamespace("HDTSA"))
  
  Tn_list <- WN_teststat_fun(Y, n, p, lag.k)
  statistic <- as.numeric(Tn_list$Tn)
  sigma_zero <- Tn_list$sigma_zero
  X_mean <- Tn_list$X_mean
  ft <- WN_ftC_fun(n, lag.k, p, Y, X_mean)
  bn <- bandwith_fun(ft, lag.k, p, p, ken_type)
  boot_nomal <- matrix(rnorm(B * (n - lag.k)), B, n - lag.k)
  bootstrap_statistics <- WN_bootc_fun(
    n, lag.k, p, B, bn, ken_type, ft, Y, sigma_zero, boot_nomal
  )
  bootstrap_statistics <- as.numeric(bootstrap_statistics)
  
  list(
    statistic = statistic,
    p.value = mean(bootstrap_statistics > statistic),
    critical_value = unname(quantile(
      bootstrap_statistics,
      probs = 1 - alpha,
      type = 1,
      names = FALSE
    ))
  )
}

run_test_method <- function(X, method, config = CONFIG, chen_seed = NULL) {
  X <- center_panel(X)
  K_values <- config$K_values
  Kmax <- max(K_values)
  
  if (method %in% c("Li2019", "Li2019_correction")) {
    ans <- LiWn_test(X, tau = Kmax, alpha = config$alpha)
    corrected <- identical(method, "Li2019_correction")
    statistic <- if (corrected) ans$statistic_correct else ans$statistic
    reject <- if (corrected) ans$res_correct else ans$res
    return(result_rows(
      method, K_values, reject[K_values],
      pnorm(statistic[K_values], lower.tail = FALSE),
      statistic[K_values], ans$critical_value[K_values]
    ))
  }
  
  if (method == "Chang2017") {
    # Use the explicitly expanded old-HDTSA implementation above, matching the
    # version used for the simulation study.
    ans <- lapply(K_values, function(K) {
      chang_wn_old_with_cv(
        X,
        lag.k = K,
        B = config$chang_B,
        kernel.type = config$chang_kernel,
        pre = config$chang_pre,
        alpha = config$alpha,
        control.PCA = config$chang_control_PCA
      )
    })
    p_value <- vapply(ans, function(z) as.numeric(z$p.value), numeric(1L))
    statistic <- vapply(ans, function(z) as.numeric(z$statistic), numeric(1L))
    cv <- vapply(ans, function(z) as.numeric(z$critical_value), numeric(1L))
    return(result_rows(
      method, K_values, p_value < config$alpha,
      p_value, statistic, cv
    ))
  }
  
  if (method == "Wang2022") {
    ans <- wnSVAR(
      X, K = Kmax, Boot = config$wang_B, alpha = config$alpha
    )
    return(result_rows(
      method, K_values,
      ans$pValue2[K_values] < config$alpha,
      ans$pValue2[K_values], ans$statistic2[K_values],
      ans$criticalValue2[K_values]
    ))
  }
  
  if (method == "Tsay2020") {
    ans <- lapply(K_values, function(K) {
      HDWNtest(X, lag = K, alpha = config$alpha, output = FALSE)
    })
    p_value <- vapply(ans, function(z) as.numeric(z$pv), numeric(1L))
    statistic <- vapply(ans, function(z) as.numeric(z$Test), numeric(1L))
    cv <- vapply(ans, function(z) as.numeric(z$CVm), numeric(1L))
    return(result_rows(
      method, K_values, p_value < config$alpha, p_value, statistic, cv
    ))
  }
  
  if (method %in% c("our", "our_resampling")) {
    use_resampling <- identical(method, "our_resampling")
    ans <- newWn_test(
      X,
      tau = Kmax,
      B = config$our_B,
      alpha = config$alpha,
      resampling = use_resampling,
      boot_type = config$our_boot_type
    )
    reject <- if (use_resampling) ans$res_resampling else ans$res
    p_value <- if (use_resampling) ans$p_value_resampling else ans$p_value
    cv <- if (use_resampling) {
      ans$critical_value_resampling
    } else ans$critical_value
    return(result_rows(
      method, K_values, reject[K_values], p_value[K_values],
      ans$statistic[K_values], cv[K_values]
    ))
  }
  
  if (method %in% c(
    "Feng2022_T_SUM", "Feng2022_T_MAX", "Feng2022_T_FC"
  )) {
    what <- switch(
      method,
      Feng2022_T_SUM = "sum",
      Feng2022_T_MAX = "max",
      Feng2022_T_FC = "fc"
    )
    ans <- FengWn_test(
      X, tau = Kmax, alpha = config$alpha, what = what
    )
    item <- switch(
      method,
      Feng2022_T_SUM = list(
        reject = ans$res_sum,
        p = ans$p_value_sum,
        statistic = ans$statistic_sum,
        cv = ans$critical_value_sum
      ),
      Feng2022_T_MAX = list(
        reject = ans$res_max,
        p = ans$p_value_max,
        statistic = ans$statistic_max,
        cv = ans$critical_value_max
      ),
      Feng2022_T_FC = list(
        reject = ans$res_fc,
        p = pchisq(ans$statistic_fc, df = 4, lower.tail = FALSE),
        statistic = ans$statistic_fc,
        cv = ans$critical_value_fc
      )
    )
    return(result_rows(
      method, K_values, item$reject[K_values], item$p[K_values],
      item$statistic[K_values], item$cv[K_values]
    ))
  }
  
  if (startsWith(method, "Chen2025_")) {
    what <- sub("^Chen2025_", "", method)
    if (is.null(chen_seed)) chen_seed <- config$seed
    X_chen <- make_chen_input(X, chen_seed)
    ans <- ChenWn_test(
      X_chen,
      tau = Kmax,
      alpha = config$alpha,
      what = what,
      check_ties = FALSE
    )
    result_name <- paste0("res_", what)
    statistic_name <- paste0("score_", what)
    statistic <- ans[[statistic_name]]
    degenerate <- what %in% c("D", "R", "tau_star")
    if (degenerate) {
      p_value <- gumbel_upper_p(statistic, multiplier = 2.467 / sqrt(pi))
      cv <- rep(ans$critical_value_degenerate, Kmax)
    } else {
      p_value <- gumbel_upper_p(statistic)
      cv <- rep(ans$critical_value_linear, Kmax)
    }
    return(result_rows(
      method, K_values, ans[[result_name]][K_values], p_value[K_values],
      statistic[K_values], cv[K_values]
    ))
  }
  
  stop("Scrolling method not yet implemented: ", method)
}


# ----------------------------------------------------------------------------
# 4. Parallel rolling execution and checkpoints
# ----------------------------------------------------------------------------

method_boot_count <- function(method, config = CONFIG) {
  if (method == "our_resampling") return(config$our_B)
  if (method == "Chang2017") return(config$chang_B)
  if (method == "Wang2022") return(config$wang_B)
  NA_integer_
}

method_signature <- function(method, config = CONFIG) {
  paste(
    "rolling_v1",
    method,
    config$window_size,
    config$step_size,
    format(config$alpha, digits = 16),
    paste(config$K_values, collapse = "-"),
    method_boot_count(method, config),
    config$our_boot_type,
    config$chang_pre,
    config$chang_kernel,
    config$seed,
    sep = "|"
  )
}

run_window_task <- function(task_row, period_data, method, config = CONFIG) {
  started <- proc.time()[["elapsed"]]
  task_row <- as.list(task_row)
  current_method <- method
  current_B <- method_boot_count(current_method, config)
  current_signature <- method_signature(current_method, config)
  seed <- stable_seed(config$seed, method, task_row$window_id)
  set.seed(seed)
  X_period <- period_data[[task_row$period_id]]$X
  X_window <- X_period[
    task_row$row_start:task_row$row_end, , drop = FALSE
  ]
  
  answer <- tryCatch(
    list(
      value = run_test_method(
        X_window, method, config,
        chen_seed = stable_seed(seed, "chen_ties")
      ),
      error = NULL
    ),
    error = function(e) list(value = NULL, error = conditionMessage(e))
  )
  elapsed <- proc.time()[["elapsed"]] - started
  
  if (is.null(answer$error)) {
    rows <- answer$value
    rows[, `:=`(status = "ok", error = "")]
  } else {
    rows <- data.table(
      method = method,
      K = config$K_values,
      reject = NA,
      p_value = NA_real_,
      statistic = NA_real_,
      critical_value = NA_real_,
      status = "error",
      error = answer$error
    )
  }
  
  rows[, `:=`(
    period_id = task_row$period_id,
    period_label = task_row$period_label,
    window_index = task_row$window_index,
    window_id = task_row$window_id,
    window_start = as.IDate(task_row$window_start),
    window_end = as.IDate(task_row$window_end),
    window_mid = as.IDate(task_row$window_mid),
    n_window = task_row$n_window,
    p = ncol(X_window),
    alpha = config$alpha,
    step_size = config$step_size,
    B = current_B,
    seed = seed,
    run_signature = current_signature,
    elapsed_seconds = elapsed
  )]
  setcolorder(rows, c(
    "period_id", "period_label", "window_index", "window_id",
    "window_start", "window_end", "window_mid", "n_window", "p",
    "method", "K", "alpha", "reject", "p_value", "statistic",
    "critical_value", "B", "seed", "step_size", "run_signature",
    "status", "error", "elapsed_seconds"
  ))
  rows
}

completed_window_ids <- function(saved, method, config = CONFIG) {
  if (!nrow(saved)) return(character())
  required <- c("window_id", "K", "status", "run_signature")
  if (!all(required %in% names(saved))) return(character())
  signature <- method_signature(method, config)
  saved[
    status == "ok" & run_signature == signature & K %in% config$K_values,
    .(complete = setequal(K, config$K_values)),
    by = window_id
  ][complete == TRUE, window_id]
}

append_checkpoint <- function(saved, new_rows, checkpoint) {
  out <- rbindlist(list(saved, new_rows), fill = TRUE)
  setorder(out, window_id, K)
  out <- unique(out, by = c("window_id", "method", "K"), fromLast = TRUE)
  fwrite(out, checkpoint)
  out
}

run_one_method_rolling <- function(
    method, task_table, period_data, config = CONFIG) {
  method_dir <- file.path(config$output_dir, "method_results")
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- file.path(method_dir, paste0("rolling_", method, ".csv"))
  saved <- if (file.exists(checkpoint)) fread(checkpoint) else data.table()
  
  if (nrow(saved) && "run_signature" %in% names(saved)) {
    incompatible <- unique(saved$run_signature)
    expected <- method_signature(method, config)
    if (any(incompatible != expected)) {
      stop(
        "Existing ", method, " checkpoints' design parameters differ from the current settings.",
        "Please change ROLLING_OUTPUT_DIR to avoid mixing different rolling designs."
      )
    }
  }
  
  completed <- completed_window_ids(saved, method, config)
  pending <- task_table[!window_id %in% completed]
  message(
    "[", method, "] Total window =", nrow(task_table),
    ", completed=", length(completed), ", pending execution=", nrow(pending)
  )
  
  if (nrow(pending)) {
    batches <- split(
      seq_len(nrow(pending)),
      ceiling(seq_len(nrow(pending)) / config$checkpoint_batch_size)
    )
    for (batch_number in seq_along(batches)) {
      batch <- pending[batches[[batch_number]]]
      indices <- seq_len(nrow(batch))
      worker_fun <- function(i) {
        run_window_task(batch[i], period_data, method, config)
      }
      
      use_fork <- isTRUE(config$parallel) && config$workers > 1L &&
        .Platform$OS.type != "windows"
      pieces <- if (use_fork) {
        parallel::mclapply(
          indices,
          worker_fun,
          mc.cores = min(config$workers, length(indices)),
          mc.preschedule = FALSE,
          mc.set.seed = FALSE
        )
      } else {
        if (isTRUE(config$parallel) && .Platform$OS.type == "windows") {
          warning("Windows does not support fork; falling back to serial execution this time.")
        }
        lapply(indices, worker_fun)
      }
      
      new_rows <- rbindlist(pieces, fill = TRUE)
      saved <- append_checkpoint(saved, new_rows, checkpoint)
      ok_windows <- uniqueN(new_rows[status == "ok", window_id])
      error_windows <- uniqueN(new_rows[status == "error", window_id])
      message(
        "[", method, "] batch ", batch_number, "/", length(batches),
        " complete：ok=", ok_windows, "，error=", error_windows
      )
    }
  }
  saved
}


# ----------------------------------------------------------------------------
# 5. Figures and summaries
# ----------------------------------------------------------------------------

method_color <- function(method) {
  palette <- c(
    Li2019 = "#1B9E77", Li2019_correction = "#66A61E",
    Chang2017 = "#7570B3", Wang2022 = "#E6AB02",
    Tsay2020 = "#D95F02", our = "#E7298A",
    our_resampling = "#C51B7D", Feng2022_T_SUM = "#1F78B4",
    Feng2022_T_MAX = "#6A3D9A", Feng2022_T_FC = "#E31A1C",
    Chen2025_rho = "#00A6A6", Chen2025_tau = "#33A02C",
    Chen2025_D = "#B15928", Chen2025_R = "#A6761D",
    Chen2025_tau_star = "#666666"
  )
  unname(palette[[method]] %||% "#333333")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

base_rolling_plot <- function(
    data, config = CONFIG,
    y_breaks = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 1)) {
  data <- copy(data)
  method_order <- unique(c(config$plot_methods, ALL_ROLLING_METHODS))
  # Saved checkpoints may contain older display labels.  Recover the current
  # label from the stable period_id so changing PERIODS$period_label never
  # collapses both periods into an NA facet.
  period_label_lookup <- setNames(
    config$periods$period_label,
    config$periods$period_id
  )
  current_period_label <- unname(
    period_label_lookup[as.character(data$period_id)]
  )
  if (anyNA(current_period_label)) {
    stop("The drawing data contains a period_id that does not exist in config$periods.")
  }
  data[, `:=`(
    window_mid = as.Date(window_mid),
    period_label = factor(
      current_period_label,
      levels = config$periods$period_label
    ),
    K_panel = factor(paste0("K = ", K), levels = paste0("K = ", config$K_values)),
    method_label = factor(
      unname(METHOD_LABELS[method]),
      levels = unname(METHOD_LABELS[method_order])
    ),
    method_math_label = factor(
      unname(METHOD_MATH_LABELS[method]),
      levels = unname(METHOD_MATH_LABELS[method_order])
    )
  )]
  events <- copy(config$events)
  
  events[, period_label := factor(
    unname(
      period_label_lookup[as.character(period_id)]
    ),
    levels = config$periods$period_label
  )]
  
  events[, event_date := as.Date(event_date)]
  
  if (anyNA(events$period_label)) {
    stop("There are period_ids in EVENTS that cannot be matched to PERIODS.")
  }
  
  ggplot(data, aes(x = window_mid, y = p_value)) +
    geom_hline(
      yintercept = config$alpha,
      linewidth = 0.55,
      linetype = "dashed",
      color = "#777777"
    ) +
    {
      if (abs(config$alpha - 0.05) > 1e-12) {
        geom_hline(
          yintercept = 0.05,
          linewidth = 0.45,
          linetype = "solid",
          color = "#333333"
        )
      }
    } +
    geom_vline(
      data = events,
      aes(
        xintercept = event_date,
        linetype = event_linetype
      ),
      inherit.aes = FALSE,
      color = "#888888",
      linewidth = 0.5
    ) +
    scale_linetype_identity() +
    facet_grid(K_panel ~ period_label, scales = "free_x") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = y_breaks
    ) +
    labs(
      # x = "Rolling-window midpoint",
      x = NULL,
      y = "p-value",
      caption = paste0(
        "Window size = ", config$window_size,
        " trading days; step size= ", config$step_size,
        # ; dashed horizontal line: alpha = ", 0.05,
        " trading days; dashed horizontal line: ", config$alpha,
        if (abs(config$alpha - 0.05) > 1e-12) "; dotted horizontal line: 0.05" else ""
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#F2F2F2", color = "#777777")
    )
}

save_method_plot <- function(data, method, config = CONFIG) {
  target_method <- method
  good <- data[
    status == "ok" & is.finite(p_value) & get("method") == target_method
  ]
  if (!nrow(good)) return(invisible(NULL))
  label <- unname(METHOD_LABELS[[method]])
  plot <- base_rolling_plot(good, config) +
    geom_line(
      aes(group = interaction(period_id, K)),
      linewidth = 0.65,
      color = method_color(method),
      na.rm = TRUE
    ) +
    # labs(title = paste0(label, ": rolling-window white-noise p-values")) +
    guides(color = "none")
  
  figure_dir <- file.path(config$output_dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(figure_dir, paste0("rolling_", method))
  ggsave(
    paste0(stem, ".png"), plot,
    width = config$plot_width, height = config$plot_height,
    dpi = config$plot_dpi, bg = "white"
  )
  ggsave(
    paste0(stem, ".pdf"), plot,
    width = config$plot_width, height = config$plot_height,
    device = "pdf"
  )
  invisible(plot)
}

longest_true_run <- function(x) {
  x <- as.logical(x)
  x[is.na(x)] <- FALSE
  runs <- rle(x)
  if (!any(runs$values)) return(0L)
  max(runs$lengths[runs$values])
}

combine_existing_methods <- function(config = CONFIG) {
  method_dir <- file.path(config$output_dir, "method_results")
  files <- list.files(
    method_dir,
    pattern = "^rolling_.*\\.csv$",
    full.names = TRUE
  )
  if (!length(files)) return(invisible(NULL))
  combined <- rbindlist(lapply(files, fread), fill = TRUE)
  combined[, `:=`(
    window_start = as.IDate(window_start),
    window_end = as.IDate(window_end),
    window_mid = as.IDate(window_mid)
  )]
  setorder(combined, method, period_id, window_index, K)
  fwrite(combined, file.path(config$output_dir, "combined_rolling_results.csv"))
  
  good <- combined[status == "ok" & is.finite(p_value)]
  if (!nrow(good)) return(invisible(combined))
  summary <- good[order(window_index), .(
    n_windows = .N,
    rejection_rate = mean(reject),
    min_p_value = min(p_value),
    median_p_value = median(p_value),
    longest_rejection_run = longest_true_run(reject)
  ), by = .(period_id, period_label, method, K, alpha)]
  fwrite(summary, file.path(config$output_dir, "rolling_summary.csv"))
  
  plot_good <- good[method %in% config$plot_methods]
  if (!nrow(plot_good)) return(invisible(combined))
  
  figure_dir <- file.path(config$output_dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  
  # One combined figure for each K.  Within a figure, rows are methods and
  # columns are the two periods; methods are therefore no longer overlaid.
  available_K <- intersect(config$K_values, sort(unique(plot_good$K)))
  for (K_now in available_K) {
    good_K <- plot_good[K == K_now]
    n_methods <- uniqueN(good_K$method)
    combined_height <- max(config$plot_height, 1.15 * n_methods + 2.2)
    
    # In the short method panels, labeling 0.00, 0.05, and 0.10 at once
    # causes overlap.  Keep the p=0.05 dashed line, but omit only its tick
    # label in the combined figure.  Individual-method figures are unchanged.
    plot <- base_rolling_plot(
      good_K, config,
      y_breaks = c(0, 0.10, 0.25, 0.50, 0.75, 1)
    ) +
      geom_line(
        aes(
          color = method_label,
          group = interaction(method, period_id)
        ),
        linewidth = 0.55,
        alpha = 0.9,
        na.rm = TRUE
      ) +
      facet_grid(
        method_math_label ~ period_label,
        scales = "free_x",
        labeller = labeller(method_math_label = label_parsed)
      ) +
      scale_color_manual(
        values = setNames(
          vapply(unique(good_K$method), method_color, character(1L)),
          unname(METHOD_LABELS[unique(good_K$method)])
        ),
        drop = TRUE
      ) +
      # labs(
      #   title = paste0(
      #     "Rolling-window white-noise p-values (K = ", K_now, ")"
      #   ),
      #   color = "Method"
      # ) +
      guides(color = "none") +
      theme(
        legend.position = "none",
        strip.text.y = element_text(angle = 0, hjust = 0)
      )
    
    stem <- file.path(
      figure_dir, paste0("rolling_all_methods_K", K_now)
    )
    ggsave(
      paste0(stem, ".png"), plot,
      width = config$plot_width, height = combined_height,
      dpi = config$plot_dpi, bg = "white"
    )
    ggsave(
      paste0(stem, ".pdf"), plot,
      width = config$plot_width, height = combined_height,
      device = "pdf"
    )
  }
  invisible(combined)
}


# ----------------------------------------------------------------------------
# 6. Main entry point
# ----------------------------------------------------------------------------

run_rolling_application <- function(config = CONFIG) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  universe <- prepare_balanced_universe(config)
  fwrite(universe$audit, file.path(config$output_dir, "stock_selection_audit.csv"))
  fwrite(
    universe$audit[PERMNO %in% universe$selected_permnos],
    file.path(config$output_dir, "selected_stock_manifest.csv")
  )
  
  period_data <- list()
  task_parts <- list()
  for (i in seq_len(nrow(config$periods))) {
    period_row <- config$periods[i]
    panel <- extract_period(universe, period_row)
    period_data[[period_row$period_id]] <- panel
    task_parts[[period_row$period_id]] <- make_window_index(
      panel$dates, period_row, config
    )
  }
  task_table <- rbindlist(task_parts)
  
  # Draw the two Tsay-style raw-return figures before any rolling tests.  A
  # common y-axis range is used so that volatility is comparable across periods.
  common_return_ylim <- range(
    unlist(lapply(period_data, function(z) z$X), use.names = FALSE),
    finite = TRUE
  )
  for (i in seq_len(nrow(config$periods))) {
    period_row <- config$periods[i]
    panel <- period_data[[period_row$period_id]]
    save_period_return_plot(
      panel$X,
      panel$dates,
      period_row$period_id,
      period_row$period_label,
      common_return_ylim,
      config
    )
  }
  
  design <- data.table(
    universe_start = config$universe_start,
    universe_end = config$universe_end,
    p_common = ncol(universe$X),
    window_size = config$window_size,
    step_size = config$step_size,
    alpha = config$alpha,
    K_values = paste(config$K_values, collapse = ","),
    methods_this_run = paste(config$methods, collapse = ","),
    methods_in_combined_plot = paste(config$plot_methods, collapse = ","),
    parallel = config$parallel,
    workers = config$workers,
    return_multiplier = config$return_multiplier,
    seed = config$seed
  )
  fwrite(design, file.path(config$output_dir, "rolling_analysis_design.csv"))
  fwrite(task_table, file.path(config$output_dir, "rolling_window_index.csv"))
  
  message(
    "rolling design: p=", ncol(universe$X),
    ", n_w=", config$window_size,
    ", s=", config$step_size,
    ", alpha=", config$alpha,
    ", K=", paste(config$K_values, collapse = ","),
    ", workers=", if (config$parallel) config$workers else 1L
  )
  
  for (method in config$methods) {
    result <- run_one_method_rolling(method, task_table, period_data, config)
    save_method_plot(result, method, config)
    combine_existing_methods(config)
  }
  
  message("Rolling analysis completed or has continued to the latest checkpoint:", config$output_dir)
  invisible(list(
    universe = universe,
    windows = task_table,
    combined = combine_existing_methods(config)
  ))
}

if (sys.nframe() == 0L) {
  run_rolling_application()
}


# 0925
PERIODS <- data.table(
  period_id = "gfc_2005_2010",
  period_label = "The global financial crisis period",
  period_start = as.IDate("2005-01-01"),
  period_end = as.IDate("2010-12-31"),
  event_date = as.IDate("2009-10-22")
)
EVENTS <- data.table(
  period_id = c(
    "gfc_2005_2010",
    "gfc_2005_2010"
  ),
  event_date = as.IDate(c(
    "2007-03-13",
    "2009-10-22"
  )),
  event_linetype = c(
    "dotted",
    "dashed"
  )
)
cfg <- CONFIG

cfg$periods <- copy(PERIODS)
cfg$events <- copy(EVENTS)
cfg$universe_start <- min(cfg$periods$period_start)
cfg$universe_end <- max(cfg$periods$period_end)

cfg$output_dir <- "./result/realdata_sp500_gfc_2005_2010_rolling_pvalues"

cfg$methods <- c(
  "our_resampling",
  "Li2019",
  "Chang2017",
  "Wang2022",
  "Tsay2020",
  "Feng2022_T_FC",
  "Chen2025_D",
  "Chen2025_R",
  "Chen2025_tau_star"
)

cfg$plot_methods <- c(
  "our_resampling",
  "Li2019",
  "Chang2017",
  "Wang2022",
  "Tsay2020",
  "Feng2022_T_FC",
  "Chen2025_D",
  "Chen2025_R",
  "Chen2025_tau_star"
)

cfg$alpha <- 0.05
cfg$K_values <- 1L

run_rolling_application(cfg)
