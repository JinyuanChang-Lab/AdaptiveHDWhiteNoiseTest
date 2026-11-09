### our method
newWn_test <- function(X, tau=1, B=1000, alpha=0.05, resampling = FALSE, boot_type = 1){
  #-----------------------------------------------------#
  # X: a n*p sample matrix
  # n: sample size
  # p: dimension
  # tau: time lag for autocovariance matrices
  # B: bootstrap time for resampling method, uesd when resampling = TRUE
  # alpha: test level
  # resampling: using resampling methon 
  # boot_type: 1 means Rademacher, 2 means normal, 3 means mammen
  #-----------------------------------------------------#
  
  n <- nrow(X)
  p <- ncol(X)
  
  Hn <- rep(0, tau)
  const <- sqrt(2) / (n-1)
  tmp <- X %*% t(X)
  diag(tmp) <- 0
  sigma_n1 <- const * sum(tmp^2)
  for(lag in c(1:tau)){
    
    Gn_1 <- sum(tmp * tmp[c(c((n-lag+1):n),c(1:(n-lag))),c(c((n-lag+1):n),c(1:(n-lag)))])
    
    
    if(lag > 1){
      Hn[lag] <- Hn[lag-1] + Gn_1/sigma_n1
    }
    else{
      Hn[lag] <- Gn_1/sigma_n1
    }
    
  }
  
  
  #--------------------calculate cv for Hn----------------------------------
  cv_vec <- qnorm(1-alpha,0,sqrt(c(1:tau)))
  p_value <- pnorm(
    Hn,
    mean = 0,
    sd = sqrt(seq_len(tau)),
    lower.tail = FALSE
  )
  
  if(resampling==TRUE){
    Hn_B <- resampling(X,n,p,B,tau,boot_type)
    cv_vec_resampling <- Hn_B[,floor(B*(1-alpha))]
    p_value_resampling <- rowMeans(
      sweep(Hn_B, 1L, Hn, FUN = ">")
    )
    return(list(
      res = Hn > cv_vec,
      res_resampling = Hn > cv_vec_resampling,
      statistic = Hn,
      p_value = p_value,
      p_value_resampling = p_value_resampling,
      critical_value = cv_vec,
      critical_value_resampling = cv_vec_resampling
    ))
  }
  
  return(list(
    res = Hn > cv_vec,
    statistic = Hn,
    p_value = p_value,
    critical_value = cv_vec
  ))
  
}



LiWn_test <- function(X, tau=1, alpha=0.05){
  #-----------------------------------------------------#
  # X: a n*p sample matrix
  # n: sample size
  # p: dimension
  # tau: time lag for autocovariance matrices
  # alpha: test level
  #-----------------------------------------------------#
  n <- nrow(X)
  p <- ncol(X)
  
  cp <- p / n
  Gq <- rep(0, tau)
  Gn <- rep(0, tau)
  if(n >= p){
    Sigma <- (1/n) * t(X) %*% X
    s1 <- sum( diag( Sigma ) )/ p
    s2 <- sum( diag( Sigma%*%Sigma )) /p
    #(1/(n))^2* sum( diag(  (t(X) %*% X) )^2) / p
    for (q in c(1:tau)){
      if(q>1){
        # Gq[q] <- Gq[q-1] + sum(diag( (1/(n^2)) * (t( X[1:(n-q), ])  %*% X[(q+1):n, ] ) %*% (t(X[(q+1):n, ]) %*% X[1:(n-q), ]) ) )
        E1 <- rbind(X[(n-q+1):n,], X[1:(n-q),])
        Gq[q] <- Gq[q-1] + sum(diag( (1/(n^2)) * (t(E1)  %*% X ) %*% (t(X)%*%E1 ) ) )
      }
      else{
        E1 <- rbind(X[(n-q+1):n,], X[1:(n-q),])
        Gq[q] <- sum(diag( (1/(n^2)) * (t(E1)  %*% X ) %*% (t(X)%*%E1 ) ) )
        # Gq[q] <- sum(diag( (1/(n^2)) * (t(X[1:(n-q), ]) %*% X[(q+1):n, ] ) %*% (t(X[(q+1):n, ]) %*% X[1:(n-q), ]) ) )
      }
      
    }
    Gn <- Gq - c(1:tau) * n* (cp^2) * (s1^2)
  }
  else{
    Sigma <- (1/n) * X %*% t(X)
    s1 <- sum( diag( Sigma ) )/ p
    s2 <- sum( diag( Sigma%*%Sigma )) /p
    
    for (q in c(1:tau)){
      if(q>1){
        # Gq[q] <- Gq[q-1] + sum(diag( (1/(n^2)) * (X[(q+1):n, ] %*% t(X[(q+1):n, ])) %*% (X[1:(n-q), ] %*% t(X[1:(n-q), ])) ) )
        
        E1 <- rbind(X[(n-q+1):n,], X[1:(n-q),])
        Gq[q] <- Gq[q-1] + sum(diag( (1/(n^2)) * (X  %*% t(X) ) %*% (E1%*%t(E1) ) ) )
      }
      else{
        # Gq[q] <- sum(diag( (1/(n^2)) * (X[(q+1):n, ] %*% t(X[(q+1):n, ])) %*% (X[1:(n-q), ] %*% t(X[1:(n-q), ])) ) )
        E1 <- rbind(X[(n-q+1):n,], X[1:(n-q),])
        Gq[q] <- sum(diag( (1/(n^2)) * (X  %*% t(X) ) %*% (E1%*%t(E1) ) ) )
      }
      
    }
    Gn <- Gq - c(1:tau) * n* (cp^2) * (s1^2)
  }
  s2_tilde <- s2 - cp * s1^2 
  #cv
  varConj <- 2* (cp^2) * (s2_tilde^2) * c(1:tau)  
  
  #correction
  tmp <- X %*% t(X)
  diag(tmp) <- 0
  sigma_n1 <- 2/ (p^2) * (1/(n*(n-1)) * sum(tmp^2) )^2
  tmp1 <- X^2
  norm <- rowSums(tmp1)
  norm_mean <- norm - mean(norm)
  sigma_n2 <- 1/ (n^4*p^2) * sum(norm_mean^4) +
    ((n-2-3*(1/n))/(p^2*n^2*(n-1)^2)) * (sum(norm_mean^2))^2
  
  Gn_tau <- Gn * (n/p)
  correction_sigma <- c(1:tau) * (sigma_n1+sigma_n2)
  
  statistic <- Gn / sqrt(varConj)
  statistic_correct <- Gn_tau / sqrt(correction_sigma)
  critical_value <- rep(qnorm(1 - alpha), tau)
  
  return(list(
    res = statistic > critical_value,
    res_correct = statistic_correct > critical_value,
    Gn = Gn,
    varConj = varConj,
    statistic = statistic,
    statistic_correct = statistic_correct,
    critical_value = critical_value
  ))
}



### Feng, Liu & Ma (2022)
FengWn_test <- function(X, tau=1, alpha=0.05, what=c("sum","max","fc")){
  #-----------------------------------------------------#
  # X: a n*p sample matrix
  # tau: time lag for autocovariance matrices (need tau <= n-2)
  # alpha: test level
  # what: A character vector, containing at least one of "sum", "max", "fc", determines which statistics to compute.
  #-----------------------------------------------------#
  if(length(what) == 0L) stop("at least one of 'sum','max','fc' must be in `what`")
  what <- match.arg(what, c("sum","max","fc"), several.ok = TRUE)
  do_sum <- "sum" %in% what
  do_max <- "max" %in% what
  do_fc  <- "fc"  %in% what
  compute_sum <- do_sum || do_fc
  compute_max <- do_max || do_fc
  need_pmat  <- compute_max
  
  n <- nrow(X)
  p <- ncol(X)
  K <- c(1:tau)
  Np2 <- K * p^2
  qalpha <- -log(pi) - 2*log(-log(1-alpha))
  
  #-------------------- sum-type --------------------#
  tmp <- X %*% t(X)
  diag(tmp) <- 0
  S <- sum(tmp^2)
  Gsum <- rep(0, tau)
  
  #-------------------- max-type --------------------#
  Z <- NULL
  Tnk <- NULL
  if(need_pmat){
    Z <- sweep(X, 2, sqrt(colMeans(X^2)), "/")
    Tnk <- rep(0, tau)
  }
  
  for(lag in 1:tau){
    # Keep exactly the same non-circular lag construction as test_fun.R.
    m <- n - lag
    Gn_1 <- sum(
      tmp[1:m, 1:m, drop = FALSE] *
        tmp[(lag+1):n, (lag+1):n, drop = FALSE]
    )
    Gsum[lag] <- if(lag > 1) Gsum[lag-1] + Gn_1 else Gn_1
    
    if(need_pmat){
      rho_hat <- t(Z[(lag+1):n,]) %*% Z[1:(n-lag),] / n
      Tnk[lag] <- sqrt(n) * max(abs(rho_hat))
    }
  }
  
  #-------------------- T_SUM p-value --------------------#
  TSUM <- Gsum * sqrt(n*(n-1)) / (S * sqrt(2*K))
  pSUM <- pnorm(TSUM, lower.tail = FALSE)
  res_sum <- if(do_sum) (TSUM > qnorm(1-alpha)) else NULL
  
  #-------------------- T_MAX p-value --------------------#
  res_max <- NULL
  pMAX <- NULL
  if(need_pmat){
    TMAX <- cummax(Tnk)
    y_max <- TMAX^2 - 2*log(Np2) + log(log(Np2))
    pMAX <- 1 - exp(-exp(-y_max/2)/sqrt(pi))
    if(do_max) res_max <- (y_max > qalpha)
  }
  
  #-------------------- Fisher Combination test T_FC --------------------#
  res_fc <- NULL
  if(do_fc){
    TFC <- -2*log(pMAX) - 2*log(pSUM)
    res_fc <- TFC > qchisq(1-alpha, df = 4)
  }
  
  res_default <- if(!is.null(res_fc)) res_fc else if(!is.null(res_max)) res_max else res_sum
  
  return(list(
    res = res_default,
    res_sum = res_sum,
    res_max = res_max,
    res_fc = res_fc,
    statistic_sum = if(compute_sum) TSUM else NULL,
    statistic_max = if(need_pmat) y_max else NULL,
    statistic_fc = if(do_fc) TFC else NULL,
    critical_value_sum = if(compute_sum) rep(qnorm(1-alpha), tau) else NULL,
    critical_value_max = if(need_pmat) rep(qalpha, tau) else NULL,
    critical_value_fc = if(do_fc) rep(qchisq(1-alpha, df=4), tau) else NULL,
    p_value_sum = if(compute_sum) pSUM else NULL,
    p_value_max = if(need_pmat) pMAX else NULL
  ))
}

# ============================================================================
# ChenWn_test: Examples 1--5
#
# what:
#   "rho"      = Spearman's rho (Example 1)
#   "tau"      = Kendall's tau (Example 2)
#   "D"        = Hoeffding's D (Example 3)
#   "R"        = Blum--Kiefer--Rosenblatt's R (Example 4)
#   "tau_star" = Bergsma--Dassios--Yanagimoto's tau* (Example 5)
#
# Example 5 requires：install.packages("TauStar")
# The current implementation is written for the paper's continuous-margin, no exact ties case.
# ============================================================================


# Exact O(m^2) implementation of the original 6th-order R U-statistic (U-statistic) in Example 4 of the paper.
# The inputs Rx and Ry must be two tie-free 1:m rank permutations.
paper_R_fast_from_ranks <- function(Rx, Ry, check = TRUE) {
  m <- length(Rx)
  if (length(Ry) != m || m < 6L) stop("Rx and Ry must have the same length and be at least 6.")
  Rx <- as.integer(Rx)
  Ry <- as.integer(Ry)
  if (check &&
      (!identical(sort(Rx), seq_len(m)) || !identical(sort(Ry), seq_len(m)))) {
    stop("Rx and Ry must both be tie-free 1:m permutations.")
  }
  
  rank_table <- matrix(0, nrow = m, ncol = m)
  rank_table[cbind(Rx, Ry)] <- 1
  
  # prefix[a,b] = #{i: Rx_i<=a, Ry_i<=b}。
  prefix <- rank_table
  if (m >= 2L) {
    for (a in 2:m) prefix[a, ] <- prefix[a, ] + prefix[a - 1L, ]
    for (b in 2:m) prefix[, b] <- prefix[, b] + prefix[, b - 1L]
  }
  
  n11 <- prefix[Rx, Ry, drop = FALSE]
  Rx_e <- matrix(Rx, nrow = m, ncol = m)
  Ry_f <- matrix(Ry, nrow = m, ncol = m, byrow = TRUE)
  n10 <- Rx_e - n11
  n01 <- Ry_f - n11
  n00 <- m - Rx_e - Ry_f + n11
  
  e_in_11 <- outer(Ry, Ry, FUN = "<=")
  n11 <- n11 - e_in_11
  n10 <- n10 - (!e_in_11)
  
  f_in_11 <- outer(Rx, Rx, FUN = ">=")
  n11 <- n11 - f_in_11
  n01 <- n01 - (!f_in_11)
  
  A <- n00 * n11
  B <- n01 * n10
  row_sq <- n00 * n11^2 + n11 * n00^2 +
    n01 * n10^2 + n10 * n01^2
  S <- 4 * (A - B)^2 - 4 * row_sq + 4 * (A + B)
  diag(S) <- 0
  
  sum(S) / (32 * choose(m, 6L))
}


ChenWn_test <- function(
    X,
    tau = 1L,
    alpha = 0.05,
    what = c("rho", "tau", "D"),
    check_ties = TRUE,
    parallel = FALSE,
    workers = 2L) {
  
  choices <- c("rho", "tau", "D", "R", "tau_star")
  if (length(what) == 0L) {
    stop("what must contain at least one of the following names: ", paste(choices, collapse = ", "))
  }
  what <- match.arg(what, choices, several.ok = TRUE)
  
  do_rho <- "rho" %in% what
  do_tau <- "tau" %in% what
  do_D <- "D" %in% what
  do_R <- "R" %in% what
  do_tau_star <- "tau_star" %in% what
  do_rank <- do_rho || do_D || do_R
  
  if (!is.matrix(X) || !is.numeric(X)) stop("X must be an n*p numeric matrix.")
  if (anyNA(X) || any(!is.finite(X))) stop("X cannot contain NA, NaN, or Inf.")
  
  n <- nrow(X)
  p <- ncol(X)
  if (n < 2L || p < 1L) stop("X requires at least 2 rows and 1 column.")
  
  if (length(tau) != 1L || is.na(tau) || tau != as.integer(tau) || tau < 1L) {
    stop("tau must be a positive integer.")
  }
  tau <- as.integer(tau)
  if (length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0,1).")
  }
  
  min_m <- 2L
  if (do_D) min_m <- max(min_m, 5L)
  if (do_R) min_m <- max(min_m, 6L)
  if (do_tau_star) min_m <- max(min_m, 4L)
  if (tau > n - min_m) {
    stop("Currently, what is required is m=n-lag >=", min_m,
         ", therefore tau must be <= ", n - min_m, "。")
  }
  
  if (do_tau_star && !requireNamespace("TauStar", quietly = TRUE)) {
    stop("Before using tau_star, please run install.packages(\"TauStar\") first.")
  }
  
  if (check_ties) {
    tied_columns <- which(vapply(seq_len(p), function(j) {
      anyDuplicated(X[, j]) > 0L
    }, logical(1L)))
    if (length(tied_columns) > 0L) {
      stop("Exact ties detected, in the following columns:", paste(tied_columns, collapse = ", "),
           ". The current formula is implemented for the continuous case without ties.")
    }
  }
  
  # The parallel switch currently only applies to Example 4 (R).
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE.")
  }
  if (!is.numeric(workers) || length(workers) != 1L ||
      is.na(workers) || !is.finite(workers) || workers < 1 ||
      workers > .Machine$integer.max || workers != floor(workers)) {
    stop("workers must be a positive integer.")
  }
  
  R_row_max <- function(i, Rx, Ry, R_fun) {
    max(vapply(seq_len(ncol(Ry)), function(j) {
      R_fun(Rx[, i], Ry[, j], check = FALSE)
    }, numeric(1L)))
  }
  environment(R_row_max) <- baseenv()
  R_fun <- if (do_R) paper_R_fast_from_ranks else NULL
  
  K <- seq_len(tau)
  Np2 <- K * p^2
  qalpha <- -log(pi) - 2 * log(-log(1 - alpha))
  kappaD <- 2.467
  QDalpha <- log(kappaD^2 / pi) - 2 * log(-log(1 - alpha))
  
  Mrho <- if (do_rho) numeric(tau) else NULL
  Mtau <- if (do_tau) numeric(tau) else NULL
  MD <- if (do_D) numeric(tau) else NULL
  MR <- if (do_R) numeric(tau) else NULL
  Mtau_star <- if (do_tau_star) numeric(tau) else NULL
  
  max_raw_D <- if (do_D) numeric(tau) else NULL
  max_raw_R <- if (do_R) numeric(tau) else NULL
  max_raw_tau_star <- if (do_tau_star) numeric(tau) else NULL
  
  rank_columns <- function(A) {
    out <- vapply(seq_len(ncol(A)), function(j) {
      rank(A[, j], ties.method = "first")
    }, numeric(nrow(A)))
    matrix(out, nrow = nrow(A), ncol = ncol(A))
  }
  
  for (lag in K) {
    m <- n - lag
    Xearly <- X[seq_len(m), , drop = FALSE]
    Xlate <- X[(lag + 1L):n, , drop = FALSE]
    
    if (do_rank) {
      Rx <- rank_columns(Xearly)
      Ry <- rank_columns(Xlate)
    }
    
    # Example 1: Spearman's rho.
    if (do_rho) {
      rbar <- (m + 1) / 2
      rho_mat <- (12 / (m * (m^2 - 1))) *
        crossprod(Rx - rbar, Ry - rbar)
      Mrho[lag] <- m * max(rho_mat^2)
    }
    
    # Example 2: Kendall's tau.
    if (do_tau) {
      tau_mat <- cor(Xearly, Xlate, method = "kendall")
      ck <- 9 * m * (m - 1) / (2 * (2 * m + 5))
      Mtau[lag] <- ck * max(tau_mat^2)
    }
    
    # Example 3: D_natural_mat = paper D_hat / 30.
    if (do_D) {
      fD <- function(r) (r - 1) * (r - 2)
      D2_mat <- crossprod(fD(Rx), fD(Ry))
      D1_mat <- matrix(0, p, p)
      D3_mat <- matrix(0, p, p)
      
      for (tt in seq_len(m)) {
        Bx_t <- sweep(Xearly, 2L, Xearly[tt, ], FUN = "<=")
        By_t <- sweep(Xlate, 2L, Xlate[tt, ], FUN = "<=")
        Qt <- crossprod(Bx_t, By_t)
        D1_mat <- D1_mat + (Qt - 1) * (Qt - 2)
        D3_mat <- D3_mat +
          outer(Rx[tt, ] - 2, Ry[tt, ] - 2) * (Qt - 1)
      }
      
      denominator_D <- m * (m - 1) * (m - 2) * (m - 3) * (m - 4)
      D_natural_mat <- ((m - 2) * (m - 3) * D1_mat + D2_mat -
                          2 * (m - 2) * D3_mat) / denominator_D
      max_raw_D[lag] <- 30 * max(D_natural_mat)
      MD[lag] <- pi^4 * (m - 1) * max(D_natural_mat)
    }
    
    # Example 4: The paper's sixth-order R itself; each pair of variables is O(m^2).
    if (do_R) {
      # R_mat <- matrix(NA_real_, p, p)
      # for (i in seq_len(p)) {
      #   for (j in seq_len(p)) {
      #     R_mat[i, j] <- paper_R_fast_from_ranks(
      #       Rx[, i], Ry[, j], check = FALSE
      #     )
      #   }
      # }
      # max_raw_R[lag] <- max(R_mat)
      
      row_max <- if (parallel && workers > 1L && p > 1L) {
        out <- parallel::mclapply(
          seq_len(p),
          R_row_max,
          Rx = Rx,
          Ry = Ry,
          R_fun = R_fun,
          mc.cores = min(as.integer(workers), p),
          mc.preschedule = TRUE,
          mc.set.seed = FALSE
        )
        
        valid <- vapply(out, function(x) {
          is.numeric(x) && length(x) == 1L && is.finite(x)
        }, logical(1L))
        if (!all(valid)) {
          stop("Some R parallel tasks failed. Please check the subprocess errors or the server memory.")
        }
        
        unlist(out, use.names = FALSE)
      } else {
        vapply(
          seq_len(p), R_row_max, numeric(1L),
          Rx = Rx, Ry = Ry, R_fun = R_fun
        )
      }
      max_raw_R[lag] <- max(row_max)
      MR[lag] <- pi^4 * (m - 1) * max_raw_R[lag] / 90
    }
    
    # Example 5: The standard TauStar output must be multiplied by 3/2 to obtain the paper scale.
    if (do_tau_star) {
      tau_star_mat <- matrix(NA_real_, p, p)
      for (i in seq_len(p)) {
        for (j in seq_len(p)) {
          tau_star_mat[i, j] <- 1.5 * TauStar::tStar(
            Xearly[, i], Xlate[, j],
            vStatistic = FALSE,
            method = "fastest"
          )
        }
      }
      max_raw_tau_star[lag] <- max(tau_star_mat)
      Mtau_star[lag] <- pi^4 * (m - 1) * max_raw_tau_star[lag] / 54
    }
  }
  
  score_rho <- if (do_rho) {
    cummax(Mrho) - 2 * log(Np2) + log(log(Np2))
  } else NULL
  score_tau <- if (do_tau) {
    cummax(Mtau) - 2 * log(Np2) + log(log(Np2))
  } else NULL
  score_D <- if (do_D) {
    cummax(MD) - 2 * log(Np2) + log(log(Np2)) + pi^4 / 36
  } else NULL
  score_R <- if (do_R) {
    cummax(MR) - 2 * log(Np2) + log(log(Np2)) + pi^4 / 36
  } else NULL
  score_tau_star <- if (do_tau_star) {
    cummax(Mtau_star) - 2 * log(Np2) + log(log(Np2)) + pi^4 / 36
  } else NULL
  
  res_rho <- if (do_rho) score_rho > qalpha else NULL
  res_tau <- if (do_tau) score_tau > qalpha else NULL
  res_D <- if (do_D) score_D > QDalpha else NULL
  res_R <- if (do_R) score_R > QDalpha else NULL
  res_tau_star <- if (do_tau_star) score_tau_star > QDalpha else NULL
  
  res_default <- if (!is.null(res_rho)) {
    res_rho
  } else if (!is.null(res_tau)) {
    res_tau
  } else if (!is.null(res_D)) {
    res_D
  } else if (!is.null(res_R)) {
    res_R
  } else {
    res_tau_star
  }
  
  list(
    res = res_default,
    res_rho = res_rho,
    res_tau = res_tau,
    res_D = res_D,
    res_R = res_R,
    res_tau_star = res_tau_star,
    score_rho = score_rho,
    score_tau = score_tau,
    score_D = score_D,
    score_R = score_R,
    score_tau_star = score_tau_star,
    max_component_rho = Mrho,
    max_component_tau = Mtau,
    max_component_D = MD,
    max_component_R = MR,
    max_component_tau_star = Mtau_star,
    max_raw_D = max_raw_D,
    max_raw_R = max_raw_R,
    max_raw_tau_star = max_raw_tau_star,
    critical_value_linear = qalpha,
    critical_value_degenerate = QDalpha
  )
}
