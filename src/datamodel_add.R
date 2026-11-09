
data_model7_fast <- function(n, p, burnin = 100) {
  n_all <- n + burnin
  
  gamma0 <- runif(p, 0.25, 0.50)
  gamma1 <- runif(p, 0, 0.25)
  beta1  <- runif(p, 0, 0.25)
  
  result <- matrix(0, nrow = n, ncol = p)
  result_row <- 1L
  
  sigma2_prev <- gamma0 / (1 - gamma1 - beta1)
  u_prev <- sqrt(sigma2_prev) * rnorm(p)
  
  if (1L > burnin) {
    result[result_row, ] <- u_prev
    result_row <- result_row + 1L
  }
  
  if (n_all >= 2L) {
    for (tt in 2:n_all) {
      sigma2_now <-
        gamma0 +
        gamma1 * u_prev^2 +
        beta1 * sigma2_prev
      
      u_now <- sqrt(sigma2_now) * rnorm(p)
      
      if (tt > burnin) {
        result[result_row, ] <- u_now
        result_row <- result_row + 1L
      }
      
      u_prev <- u_now
      sigma2_prev <- sigma2_now
    }
  }
  
  result
}

data_model8_fast <- function(n, p, burnin = 100) {
  n_all <- n + burnin
  k0 <- max(1L, p %/% 4L)
  active <- seq_len(k0)
  
  # AR(1) coefficients of the first k0 components
  a_vec <- rep(0.4, k0)
  
  gamma0 <- runif(p, 0.25, 0.50)
  gamma1 <- runif(p, 0, 0.25)
  beta1  <- runif(p, 0, 0.25)
  
  result <- matrix(0, nrow = n, ncol = p)
  result_row <- 1L
  
  sigma2_prev <- gamma0 / (1 - gamma1 - beta1)
  e_prev <- sqrt(sigma2_prev) * rnorm(p)
  
  x_prev <- numeric(p)
  x_now <- e_prev
  # x_now[active] <-
  #   drop(A0 %*% x_prev[active]) + e_prev[active]
  x_now[active] <-
    a_vec * x_prev[active] + e_prev[active]
  
  if (1L > burnin) {
    result[result_row, ] <- x_now
    result_row <- result_row + 1L
  }
  
  x_prev <- x_now
  
  if (n_all >= 2L) {
    for (tt in 2:n_all) {
      sigma2_now <-
        gamma0 +
        gamma1 * e_prev^2 +
        beta1 * sigma2_prev
      
      e_now <- sqrt(sigma2_now) * rnorm(p)
      
      # When j > k0, the entire row corresponding to A is all zeros, so x_tj = e_tj
      x_now <- e_now
      
      # Update only the first k0 non-white-noise components
      # x_now[active] <-
      #   drop(A0 %*% x_prev[active]) + e_now[active]
      x_now[active] <-
        a_vec * x_prev[active] + e_now[active]
      
      
      if (tt > burnin) {
        result[result_row, ] <- x_now
        result_row <- result_row + 1L
      }
      
      e_prev <- e_now
      sigma2_prev <- sigma2_now
      x_prev <- x_now
    }
  }
  
  result
}
data_model7 <- data_model7_fast
data_model8 <- data_model8_fast


