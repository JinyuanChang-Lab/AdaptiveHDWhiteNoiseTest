library(Rcpp)
library(HDTSA)
library(parallel)

sourceCpp("./src/tool.cpp")
source("./src/test_fun.R")
source("./src/datamodel_add.R")
source("./src/wnSVAR.R")
source("./src/HDWHtest.R")


run_simulation_model78_only <- function(
    model_func,
    model_name,
    start_ind = 1L,
    end_ind = 1L,
    reptime = 5000L,
    cores = 200L,
    seed_offset = 12345L) {
  
  reported_lags <- c(2L, 4L, 6L, 8L)
  k <- max(reported_lags)
  alpha <- 0.05
  
  col_names <- c(
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
  row_names <- c(
    "n50p5", "n100p5", "n200p5",
    "n50p7", "n100p10", "n200p14",
    "n50p50", "n100p100", "n200p200",
    "n50p2500", "n100p10000", "n200p40000"
  )
  
  new_result_matrix <- function() {
    ans <- matrix(
      NA_real_,
      nrow = length(row_names),
      ncol = length(col_names),
      dimnames = list(row_names, col_names)
    )
    ans
  }
  
  result1 <- new_result_matrix()
  result2 <- new_result_matrix()
  result3 <- new_result_matrix()
  result4 <- new_result_matrix()
  
  n_vec <- c(50L, 100L, 200L, 50L, 100L, 200L,
             50L, 100L, 200L, 50L, 100L, 200L)
  p_vec <- c(5L, 5L, 5L,
             floor(sqrt(n_vec[1:3])),
             50L, 100L, 200L,
             n_vec[1:3]^2)
  
  start_ind <- as.integer(start_ind)
  end_ind <- as.integer(end_ind)
  reptime <- as.integer(reptime)
  cores <- as.integer(cores)
  seed_offset <- as.integer(seed_offset)
  
  if (start_ind < 1L || end_ind > length(n_vec) || start_ind > end_ind) {
    stop("Must satisfy 1 <= start_ind <= end_ind <= 12.")
  }
  if (reptime < 1L || cores < 1L) stop("reptime and cores must be positive integers.")
  
  output_dir <- file.path("./result/", paste0(model_name, "_sim78"))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (ind in seq.int(start_ind, end_ind)) {
    ind_start_time <- Sys.time()
    n <- n_vec[ind]
    p <- p_vec[ind]
    
    cores_now <- cores
    if (ind == 9L) cores_now <- min(cores_now, 200L)
    if (ind == 10L) cores_now <- min(cores_now, 200L)
    if (ind == 11L) cores_now <- min(cores_now, 200L)
    
    cat(
      "\n=== running: ", model_name,
      " | ind=", ind,
      " | n=", n,
      " | p=", p,
      " | cores=", cores_now,
      " | start=", format(ind_start_time, "%Y-%m-%d %H:%M:%S"),
      " ===\n",
      sep = ""
    )
    
    R_results <- mclapply(
      seq_len(reptime),
      function(time) {
        set.seed(time + seed_offset)
        
        if (time %% cores_now == 0L) {
          cat(
            "model=", model_name,
            " ind=", ind,
            " time=", time,
            "\n",
            file = file.path(output_dir, "log.txt"),
            append = TRUE,
            sep = ""
          )
        }
        
        X_model <- model_func(n, p, burnin=0)
        
        res <- matrix(
          NA_real_,
          nrow = length(reported_lags),
          ncol = length(col_names),
          dimnames = list(paste0("lag", reported_lags), col_names)
        )
        # Our newWn_test
        if("our" %in% col_names){
          tmp <- newWn_test(X_model, k, resampling = TRUE)
          res[, "our_resampling"] <- tmp$res_resampling[c(2,4,6,8)]
        }
        # LLYY test
        tmp <- LiWn_test(X_model, k)
        res[, "Li2019"] <- tmp$res[c(2,4,6,8)]
        
        res
      },
      mc.cores = cores_now
      # mc.preschedule = FALSE
    )
    
    mean_for_lag <- function(lag_row) {
      values <- do.call(rbind, lapply(R_results, function(z) z[lag_row, ]))
      colMeans(values, na.rm = FALSE) * 100
    }
    
    result1[ind, ] <- mean_for_lag(1L)
    result2[ind, ] <- mean_for_lag(2L)
    result3[ind, ] <- mean_for_lag(3L)
    result4[ind, ] <- mean_for_lag(4L)
    
    file_suffix <- paste0("_", model_name, "_n", n, "p", p, ".csv")
    write.csv(
      result1,
      file.path(output_dir, paste0("res_q1", file_suffix)),
      na = "NA"
    )
    write.csv(
      result2,
      file.path(output_dir, paste0("res_q2", file_suffix)),
      na = "NA"
    )
    write.csv(
      result3,
      file.path(output_dir, paste0("res_q3", file_suffix)),
      na = "NA"
    )
    write.csv(
      result4,
      file.path(output_dir, paste0("res_q4", file_suffix)),
      na = "NA"
    )
    ind_end_time <- Sys.time()
    ind_elapsed <- as.numeric(
      difftime(ind_end_time, ind_start_time, units = "secs")
    )
    cat(
      "=== Completed: ", model_name,
      " | ind=", ind,
      " | n=", n,
      " | p=", p,
      " | end=", format(ind_end_time, "%Y-%m-%d %H:%M:%S"),
      " | elapsed=", sprintf(
        "%02d:%02d:%06.2f",
        floor(ind_elapsed / 3600),
        floor((ind_elapsed %% 3600) / 60),
        ind_elapsed %% 60
      ),
      " ===\n",
      sep = ""
    )
    cat("result save to: ", output_dir, "\n", sep = "")
    
    
    rm(R_results)
    invisible(gc())
  }
  
  invisible(list(
    q1 = result1,
    q2 = result2,
    q3 = result3,
    q4 = result4,
    output_dir = output_dir
  ))
}



# run_simulation_model78_only(
#   data_model7, "model7", start_ind=1, end_ind=12, reptime=5000, cores=360
# )
# run_simulation_model78_only(
#   data_model8, "model8", start_ind=1, end_ind=12, reptime=5000, cores=360
# )


