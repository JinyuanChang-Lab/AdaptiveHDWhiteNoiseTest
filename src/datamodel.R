#model1
data_model1 <- function(n,p){
  zt <- matrix(rgamma(n*p,4,2), n, p)-2
  # sigma0 <- diag(1,p,p)
  return(zt)
}

#model2
data_model2 <- function(n, p) {
  nu <- 9
  Y <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  
  v <- rchisq(n, df = nu)
  
  scale_factor <- sqrt(7 / v)
  
  Z <- Y * scale_factor
  return(Z)
}
#model3
data_model3 <- function(n, p) {
  norm <- matrix(rnorm(n*p), n, p)
  zt1 <- sqrt(0.2) * norm
  zt2 <- sqrt(1.8) * norm
  et1 <- rbinom(n,1,0.5)
  et2 <- 1- et1
  zt <- et1*zt1 + et2*zt2
  return(zt)
}

#model4
data_model4 <- function(n, p) {
  zt <- matrix(rgamma(n * p, 4, 2), p, n) - 2
  yt <- matrix(0, p, n)
  yt[, 1] <- zt[, 1]
  for(t in 2:n) {
    yt[, t] <- 0.1 * yt[, t-1] + zt[, t]
  }
  return(t(yt))
}

#model5
# data_model5 <- function(n,p){
#   theta <- mat.or.vec(p,p)
#   index <- which(abs(row(theta) - col(theta)) <= 1, arr.ind = TRUE)
#   random_probs <- runif(nrow(index))
#   random_values <- runif(nrow(index), -0.95, 0.95)
#   theta[index] <- ifelse(random_probs > 0.3, 0, random_values)
#   
#   nu <- 9
#   Y <- matrix(rnorm((n+1) * p, mean = 0, sd = 1), nrow = n+1, ncol = p)
#   v <- rchisq(n+1, df = nu)
#   scale_factor <- sqrt(7 / v)
#   et <- t(Y * scale_factor)
#   
#   X <- mat.or.vec(p,n)
#   # et <- matrix(rnorm((n+1)*p),p,n+1)
#   for(i in c(1:n))
#     X[, i] <- et[, i+1] - theta %*% et[, i]
#   return(t(X))
# }

#model6
data_model6 <- function(n, p){
  norm <- matrix(rnorm(n * p), n, p)
  zt1 <- sqrt(0.2) * norm
  zt2 <- sqrt(1.8) * norm
  et1 <- rbinom(n, 1, 0.5)
  et2 <- 1 - et1
  zt <- t(et1 * zt1 + et2 * zt2)
  xt <- zt + 0.07 * cbind(0, zt[, 1:(n-1)])
  return(t(xt))
}


data_model5_lowmem <- function(n, p) {
  n <- as.integer(n)
  p <- as.integer(p)
  
  # ------------------------------------------------------------
  # 1. Generate the 3p-2 candidate elements of the tridiagonal theta.
  # ------------------------------------------------------------
  n_theta <- 3L * p - 2L
  
  random_probs <- runif(n_theta)
  random_values <- runif(n_theta, -0.95, 0.95)
  
  theta_values <- ifelse(
    random_probs > 0.3,
    0,
    random_values
  )
  
  theta_diag <- theta_values[3L * seq_len(p) - 2L]
  
  if (p > 1L) {
    id <- seq_len(p - 1L)
    theta_lower <- theta_values[3L * id - 1L]
    theta_upper <- theta_values[3L * id]
  } else {
    theta_lower <- numeric(0L)
    theta_upper <- numeric(0L)
  }
  
  rm(random_probs, random_values, theta_values)
  
  # ------------------------------------------------------------
  # 2. Generate e_1,...,e_{n+1} all at once.
  # Each row of E is a p-dimensional innovation.
  # ------------------------------------------------------------
  nu <- 9
  
  E <- matrix(
    rnorm((n + 1L) * p),
    nrow = n + 1L,
    ncol = p
  )
  
  v <- rchisq(n + 1L, df = nu)
  scale_factor <- sqrt(7 / v)
  
  E <- E * scale_factor
  
  # ------------------------------------------------------------
  # 3. Compute X_t = e_{t+1} - theta e_t
  # ------------------------------------------------------------
  X <- matrix(0, nrow = n, ncol = p)
  
  for (tt in seq_len(n)) {
    e_current <- E[tt, ]
    e_next <- E[tt + 1L, ]
    
    theta_e <- theta_diag * e_current
    
    if (p > 1L) {
      # Upper diagonal
      theta_e[seq_len(p - 1L)] <-
        theta_e[seq_len(p - 1L)] +
        theta_upper * e_current[2L:p]
      
      # Lower diagonal 
      theta_e[2L:p] <-
        theta_e[2L:p] +
        theta_lower * e_current[seq_len(p - 1L)]
    }
    
    X[tt, ] <- e_next - theta_e
  }
  
  X
}
data_model5 <- data_model5_lowmem
