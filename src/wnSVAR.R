library(ncvreg)

VARest <- function(X, ARorder=1, lambda=-1, method)
{
  
  if (!is.matrix(X)) X = matrix(X, ncol=1)
  
  np = dim(X);
  n = np[1]
  p = np[2]
  
  
  Y = matrix(X[(ARorder+1):n,], ncol=p)
  XL = c()
  for (ilag in 1:ARorder)
    XL = cbind(XL,  X[ilag:(n-ARorder+ilag-1),])
  
  X = XL
  E = Y
  AR = E
  
  A = matrix(0, p, p*ARorder)
  if (max(lambda) < 0)
    lambda = rep(-1, p)
  
  
  for (i in 1:p)
  {
    y = Y[,i]
    if (toupper(method) == "LASSO")
    {
      
      penalty = rep(1, ncol(X))
      #      penalty[1] = 0
      
      if (lambda[i] < 0)
      {
        #        cvfit <- cv.glmnet(X, y)
        cvfit <- cv.ncvreg(X, y)
        lambda[i] = cvfit$lambda.min
      }
      beta <- ncvreg(X, y, lambda = lambda[i])$beta
      E[,i] = y - (beta[1] + X%*%matrix(beta[-1], ncol=1))
      A[i,] = as.numeric(beta[-1])
    }else
    {
      fit = lm(y~X)
      A[i,] = fit$coefficients[-1]
      E[,i] = fit$residuals
    }
  }
  
  return(list(E=E, lambda=lambda, A = A))
  
}


Macf <- function(y,K)
{
  np = dim(y);
  n = np[1]
  p = np[2]
  
  
  y = y - t(matrix(rep(colMeans(y), n), p, n));
  
  ACF2 = c();
  ACF0 = c();
  ACFs1 = c();
  ACFs10 = c();
  
  S = 1/sqrt(colSums(y^2))
  if (p > 1) S = diag(S, p)
  
  if ((p > 1)& (p<n/5))    ######################
  {
    S = t(y)%*%y;
    S = svd(S)
    S = S$u%*%diag(1/(sqrt(S$d)+1/n^3))%*%t(S$u)
  }
  
  acf0 = 0;
  acf2 = 0;
  
  acf_s1 = 0;
  acf_s10 = 0;
  
  rA = 0
  #  sn = min(floor(n/log(p)^2), p^2)
  sn = p
  
  for (i in 1:K)
  {
    rk = abs(S%*% (t(y[1:(n-i), ])%*%y[(i+1):n, ]) %*% S)/(n-i)
    rA = rA + rk
    si = rowSums(rk)
    
    srk = sort(rA, decreasing = TRUE) 
    acf0 = sum(srk[1]);    
    acf2 = acf2 + sum(si);    
    
    
    rk = sort(rk, decreasing = TRUE)
    acf_s1 =  max(acf_s1, sum(srk[1:sn]));
    acf_s10 = acf_s10 + sum(rk[1:sn]);
    
    
    ACF2 = c(ACF2, acf2);
    ACF0 = c(ACF0, acf0);
    ACFs1 = c(ACFs1, acf_s1);
    ACFs10 = c(ACFs10, acf_s10);
  }
  return(list(ACF2=ACF2, ACF0=ACF0, ACFs1=ACFs1, ACFs10=ACFs10))
}


wnSVAR <- function(X, ARorder=0, K= 10, method= "LASSO", Boot=1250)
{
  
  p = dim(X)[2]
  for (i in 1:p)
    X[,i] = (X[,i]-mean(X[,i]))/(sd(X[,i])+1.0e-10)
  
  JB2 = matrix(0,Boot, K);
  JB0 = matrix(0,Boot, K);
  JBs1 = matrix(0,Boot, K);
  JBs10 = matrix(0,Boot, K);
  
  
  if (ARorder == 0)
  {
    np = dim(X);
    n <<- np[1]
    p <<- np[2]
    K <<- K
    
    LB0 = Macf(X, K);
    jb2 <<- LB0$ACF2 
    jb0 <<- LB0$ACF0 
    jbs1 <<- LB0$ACFs1
    jbs10 <<- LB0$ACFs10
    
    for (i in 1:Boot)
    {
      e = (runif(n)>0.5)*2-1;
      LBi = Macf(X*(matrix(rep(e, p), n, p)), K);
      JB2[i,] = LBi$ACF2 - jb2;
      JB0[i,] = LBi$ACF0 - jb0;    
      JBs1[i,] = LBi$ACFs1 - jbs1;
      JBs10[i,] = LBi$ACFs10 - jbs10;    
    }
    
    JB2 = t(JB2)
    JB0 = t(JB0)
    JBs1 = t(JBs1)
    JBs10 = t(JBs10)
    
    A = NULL
  }else
  {
    lambda = rep(-1, p)
    AR0 = VARest(X, ARorder=ARorder, lambda=lambda, method=method)
    A = AR0$A
    B = matrix(0, p*ARorder, p*ARorder)
    B[1:p, ] = A
    if (ARorder > 1)
      B[(p+1):(ARorder*p), 1:((ARorder-1)*p)] = diag(1, (ARorder-1)*p)
    val = max(abs(eigen(B)$values))
    if (val >= 1)
      A = AR0$A/(val+0.1)
    
    A = A*(abs(A)>1.0e-12);
    
    
    lambda = AR0$lambda
    E = AR0$E
    n = nrow(E)
    
    k0 = max(ARorder) + 1   ###################
    LB0 = Macf(E[k0:n,], K);   
    jb2 = LB0$ACF2 
    jb0 = LB0$ACF0 
    jbs1 = LB0$ACFs1
    jbs10 = LB0$ACFs10
    
    
    n0 = 200
    rep = 1
    JB = c()
    #    JB <- foreach(rep=1:Boot, .combine='cbind') %dopar%
    for (rep in 1:Boot)
    {
      
      #        source('Macf.R')
      #        source('VARest.R')
      
      
      Y = matrix(0, p,n+n0)
      for (i in (ARorder+1):n0)
      {
        e = E[sample(1:n, 1),] 
        xi = c()
        for (ilag in 1:ARorder)
          xi = c(xi, Y[,i-ilag])
        Y[,i] = A%*%matrix(xi, ncol=1) + e;
      }
      for (i in (n0+1):(n+n0))
      {
        e = E[i-n0, ]*(2*(runif(1)>0.5)-1) 
        xi = c()
        for (ilag in 1:ARorder)
          xi = c(xi, Y[,i-ilag])
        Y[,i] = A%*%matrix(xi, ncol=1) + e;
      }
      Y = t(Y[,(n0+1):(n0+n)])
      for (i in 1:p)
        Y[,i] = (Y[,i]-mean(Y[,i]))/(sd(Y[,i])+1.0e-10)
      
      Eb = VARest(Y, ARorder=ARorder, lambda=lambda, method=method)$E
      n1 = nrow(Eb)
      LBi = Macf(Eb[(n1-n+k0):n1,], K);
      
      JB = cbind(JB, c(LBi$ACF0, LBi$ACF2,  LBi$ACFs1, LBi$ACFs10))
    }
    
    JB0 = JB[1:K,] - matrix(jb0, K, Boot);
    JB2 = JB[(K+1):(2*K),] - matrix(jb2, K, Boot);
    JBs1 = JB[(2*K+1):(3*K),] - matrix(jbs1, K, Boot);
    JBs10 = JB[(3*K+1):(4*K),] - matrix(jbs10, K, Boot);
    
    
  }
  
  pValue2 = rowMeans(JB2>0);
  pValue0 = rowMeans(JB0>0);
  pValues1 = rowMeans(JBs1>0);
  pValues10 = rowMeans(JBs10>0);
  
  
  return(list(pValue0=pValue0, pValue2=pValue2, 
              pValue0sn=pValues1, pValue2sn=pValues10, A = A))
}


