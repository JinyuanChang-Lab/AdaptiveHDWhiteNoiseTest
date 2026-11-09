"HDWNtest" <- function(xt,lag=10,alpha=c(0.1,0.05,0.01),subfr=0.75,output=TRUE){
### High-dimensional white-noise test
### modified on April 17, 2019 by Ruey Tsay (based on the revision of the paper)
#### xt: n-by-p data matrix of time series with n observations and p dimensions
#### lag: number of lags used in the joint test.
#### alpha: type-I error, the default is 10%, 5% and 1%.
#### subfr: fraction of dimension selected for testing. 
####            Specifically, the subset used for testing when p >= n is floor(n*subfr). Default is 75%.
####
if(!is.matrix(xt))xt <- as.matrix(xt)
p <- ncol(xt); n <- nrow(xt)
#
Tst <- NULL; Test <- NULL; Sel <- NULL
if(p < n){
 m1 <- HDwnT(xt,lag,alpha=alpha,output=output)
   }else{
  sdim <- floor(n*subfr)
  # cat("Dimension of sub-samples: ",sdim,"\n")
  wgt <- Selwgt(xt,lag)
  mm <- sort(wgt,index=TRUE)
  Sel <- mm$ix[p:(p-sdim+1)]
  zt <- xt[,Sel]
  m1 <- HDwnT(zt,lag,alpha=alpha,output=output)
 }
 Tst <- m1$Tst; Test <- m1$Test; cndn <- m1$cndn; cmdm <- m1$cmdm
 CVi <- m1$CVi; CVm <- m1$CVm
 pv <- m1$pv; pvi <- m1$pvi
##
HDwnTest <- list(Tst=Tst,Test=Test,cndn=cndn,cmdm=cmdm,CVi=CVi,CVm=CVm,Sel=Sel,pv=pv,pvi=pvi)
}

"HDwnT" <- function(xt,lag,lags,alpha,output=TRUE){
 n <- nrow(xt); p <- ncol(xt)
 nnm1 <- n*(n-1)
 tmp <- c((n-1):(n-lag))*5 - 4
 tmp <- tmp/(5*n*n)
 eff <- 1/tmp
 rootT <- sqrt(eff)
 nrhovar <- n*(n^2-1)/12
 rankmean <- (n+1)/2
 bias <- 1/n
 prob <- 1-(alpha/2)
 xo <- - log(-log(prob))
### Transformation
 m1 <- princomp(xt)
 yt <- m1$scores
 yt <- apply(yt,2,rank)-rankmean
 efb <- p*p
### compute location and normalization factor
 nsize <- efb
 cn <- 1/sqrt(2*log(nsize))
 dn <- sqrt(2*log(nsize))-(log(4*pi)+log(log(nsize)))/(2*sqrt(2*log(nsize)))
 cndn <- c(cn,dn)
 xi <- cn*xo+dn
###
### Perform the test
 Tst <- NULL; pvalue <- NULL
 for (i in 1:lag){
 bias <- (n-i)/nnm1
 G2 <- crossprod(yt[(i+1):n,],yt[1:(n-i),])/nrhovar+diag(rep(bias,p))
  tmp1 <- max(G2); tmp2 <- -min(G2)
  Tst1 <- max(tmp1,tmp2)*rootT[i]
#### Compute the asymptotic p-value
  tt <- (Tst1-dn)/cn
  pv <- 1 - (exp(-exp(-tt)))^2
 Tst <- c(Tst, Tst1)
 pvalue <- c(pvalue,pv)
 }
 Tst <- round(Tst,3)
 pvalue <- round(pvalue,3)
 if(output){
  cat("Test results using individual lag: ","\n")
  cat("alpha-level: ",alpha,"\n")
  cat("Critical value: ",xi,"\n")
  Lag <- c(1:lag)
  tbl <- cbind(Lag,Tst,pvalue)
  print(tbl)
  cat(" ","\n")
 }
### Joint test
  m <- lag
  cnm <- 1/sqrt(2*log(nsize*m))
  dnm <- sqrt(2*log(m*nsize))-(log(4*pi)+log(log(nsize*m)))/(2*sqrt(2*log(nsize*m)))
  cmdm <- c(cnm,dnm)
  Test <- max(Tst)
  tt <- (Test-dnm)/cnm
  pv <- 1 - (exp(-exp(-tt)))^2
  xm <- cnm*xo+dnm
 if(output){
  cat("Joint test for the first ",lag," lags: ","\n")
  cat("alpha: ",alpha,"\n")
  cat("Critical value: ",xm,"\n")
  cat("Test statistic and p-value: ",round(c(Test,pv),3),"\n")
 }
#
HDwnT <- list(Tst=Tst,Test=Test,CVi=xi,CVm=xm,cndn=cndn,cmdm=cmdm,pv=pv,pvi=pvalue)
}

##### computing sub-sampling probabilities
"Selwgt" <- function(x,lag=10){
### select top sub-series
if(!is.matrix(x))x <- as.matrix(x)
n <- nrow(x); p <- ncol(x)
###
Select <- NULL
wgt <- rep(1,p)/p
for (i in 1:lag){
 rho <- cor(x[(i+1):n,],x[1:(n-i),],method="spearman")
 rho <- abs(rho)
 wgt <- wgt + apply(rho,1,max)
 }
wgt <- wgt/sum(wgt)
##cat("prob: ","\n")
###print(round(wgt,4))
Selwgt <- wgt
}

#### Use simulation to compute the empircal critical values of joint test if needed.
###
"HDcriVa" <- function(p=100,n=300,iter=1000,df=3,lag=10,subfr=0.75,alpha=c(0.1,0.05,0.01)){
###
Tst <- Test <- NULL
for (it in 1:iter){
  if(df > 25){
   xt <- matrix(rnorm(p*n),n,p)
   }else{
    xt <- matrix(rt(n*p,df=df),n,p)
  }
  m1 <- HDwnTest(xt,lag=lag,alpha=alpha,subfr=subfr,output=FALSE)
  Tst <- rbind(Tst,m1$Tst)
  Test <- c(Test,m1$Test)
 }

crit <- quantile(Test,prob=(1-alpha))
cat("Critical values for tail prob: ",alpha,"\n")
print(round(crit,2))

HDcritVa <- list(Tst=Tst,Test=Test,cv=crit)
}