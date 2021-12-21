#' @title Performance of lasso
#' @description A function to measure the performance of lasso in high-dimensioanl linear regression model
#' @param X the independent variables in the linear regression model
#' @param beta the coefficient 
#' @param sigma the standard error of epsilon
#' @param N the number of simulation
#' @return the average number of regressors missed from the true model,the average number of regressors selected outside the model,relative average empirical risk
#' @import glmnet
#' @importFrom stats rnorm coef qnorm
#' @examples 
#' \dontrun{
#' library(MASS)
#' n=100
#' p=500
#' sigma <- matrix(numeric(p^2),nrow=p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'     if(i==j){sigma[i,j]<-1}
#'     else{sigma[i,j]<- (1/2)^abs(i-j)}
#'   }
#' }
#' mean <- numeric(p)
#' set.seed(1234)
#' x <- mvrnorm(n,mean,sigma)
#' beta <- numeric(p)
#' beta[1:5]=1
#' avermiss(x,beta,1,100)
#'  }
#' @export



avermiss <- function(X,beta,sigma,N){
  n=nrow(X);p=ncol(X);
  m=length(beta[which(beta>0.0001)])
  En <- t(X)%*%X/n  
  x0=X[,1:m]
  A=t(x0)%*%x0
  nmiss=0
  nfalse=0
  rrisk=0
  alpha=0.05
  c=1.1
  s1=c/10*sigma*qnorm(1-alpha/(2*p))
  for(i in 1:N){
    epsilon=rnorm(n)
    yn=X%*%beta+sigma*epsilon
    bet <- solve(A)%*%t(x0)%*%yn
    beta0 <- numeric(p)
    beta0[1:m]=bet
    delta1 <- beta0-beta
    so=t(delta1)%*%En%*%delta1
    
    oracler <- sqrt(so)
    fit <- glmnet(X,yn,intercept = FALSE)
    cofficient <- coef(fit,s=s1)
    cofficient <- cofficient[-1]
    index <- which(cofficient>0.0001)
    betan <- cofficient[index]
    delta <- cofficient-beta
    sumtrue=0
    sumfalse=0
    if(length(index)==0){miss=m; sumfalse=0}
    else {for(i in  1:length(index)){
      if(0 < index[i] & index[i]< m+1) {sumtrue=sumtrue+1}
      else {sumfalse=sumfalse+1}
    }
      miss=m-sumtrue}
    ab <- c(miss,sumfalse)
    s=t(delta)%*%En%*%delta
    nrisk <- sqrt(s)
    nmiss=nmiss+ab[1]
    nfalse=nfalse+ab[2]
    rrisk=rrisk+nrisk/oracler
  }
  
  b <- c(nmiss/N,nfalse/N,rrisk/N)
  return(b)
}


