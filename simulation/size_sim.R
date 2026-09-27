library(parallel)


# ==========================================
# Simulation Settings
# ==========================================
n_seq <- c(50, 100, 200,300,400)
p_types <- c("n", "n^2")
reps <- 10000   
df_seq <- c(3, 5, 10, Inf)

dir.create("./result/plot_for_size", recursive = TRUE, showWarnings = FALSE)
out_file <- "./result/plot_for_size/figure_simulation_results_new.csv"
if (!file.exists(out_file)) {
  cat("n,p,p_type,df_t,empirical_size\n", file = out_file)
}


num_cores <- max(1, detectCores() - 1)
RNGkind("L'Ecuyer-CMRG")
cat(sprintf("using cores: %d\n", num_cores))

run_single_sim <- function(time) {
  # 1. Generate the Gram matrix of high-dimensional data W = X^T X
  if(p >= n){
    if (is.infinite(df_t)) {
      # normal
      W <- rWishart(1, df = p, Sigma = diag(n))[,,1] 
    } else {
      # Multivariate t-distribution
      u <- rchisq(n, df = df_t)
      scale_factor <- sqrt((df_t) / u) 
      Y_TY <- rWishart(1, df = p, Sigma = diag(n))[,,1]
      W <- t(t(Y_TY) * scale_factor) * scale_factor 
    }
  }
  else {
    Z <- matrix(rnorm(p * n), nrow = p, ncol = n)
    Z_TZ <- crossprod(Z) 
    
    if (is.infinite(df_t)) {
      W <- Z_TZ
    } else {
      u <- rchisq(n, df = df_t)
      scale_factor <- sqrt((df_t) / u)
      W <- t(t(Z_TZ) * scale_factor) * scale_factor 
    }
  }
  
  tr_S0 <- sum(diag(W)) / n
  tr_S0_sq <- sum(W^2) / n^2
  Gq <- sum(W * W[idx1, idx1])/n^2 + sum(W * W[idx2, idx2])/n^2 - (2/n) * tr_S0^2
  
  s2_tilde <- (1 / p) * tr_S0_sq - (1 / (n * p)) * (tr_S0^2)
  cv_sqrt <- (2 *p* s2_tilde) / n

  return(Gq / cv_sqrt)
}

for (n in n_seq) {
  
  idx1 <- c(2:n, 1)        
  idx2 <- c(3:n, 1:2)      
  
  for (p_type in p_types) {
    
    if (p_type == "log(n)") {
      p <- max(2, ceiling(n^(1/8)))
    } else if (p_type == "n") {
      p <- n
    } else if (p_type == "n^2") {
      p <- n^2
    }
    
    for (df_t in df_seq) {
      
      cat(sprintf("\n[running] n = %d, p = %d (%s), df = %s ...\n", 
                  n, p, p_type, as.character(df_t)))
      
      set.seed(123) 
      
      run_time <- system.time({
        res_list <- mclapply(1:reps, run_single_sim, mc.cores = num_cores)
      })

      res <- do.call(rbind, res_list)
      size_1 <- mean(res[, 1] > 1.645) * 100
      
      cat(sprintf("  -> time: %.2f sec | Size: %.2f%%\n", run_time["elapsed"], size_1))
      
      current_result <- data.frame(
        n = n,
        p = p,
        p_type = p_type,
        df_t = df_t,
        empirical_size = size_1
      )
      
      write.table(current_result, file = out_file, sep = ",", 
                  append = TRUE, col.names = FALSE, row.names = FALSE)
      
    }
  }
}

cat("\nAll simulation matrices have finished running! The results have been incrementally saved in:", out_file, "\n")