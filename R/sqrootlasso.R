#' @title Square-root lasso
#' @description A function to estimate high-dimensional sparse linear regression model 
#' @param x the independent variables in the linear regression model
#' @param y the response variable in the linear regression model
#' @param lambda the penalty factor
#' @param maxIter the maximum number of iterations
#' @return beta a estimate of regression coefficient
#' @examples 
#' \dontrun{
#' library(MASS)
#' n=100
#' p=500
#' sigma <- matrix(numeric(p^2),nrow=p)
#'  for(i in 1:p){
#'  for(j in 1:p){
#'  if(i==j){sigma[i,j]<-1}
#'  else{sigma[i,j]<- (1/2)^abs(i-j)}}}
#'  mean <- numeric(p)
#'  set.seed(1234)
#'   x <- mvrnorm(n,mean,sigma)
#'   beta <- numeric(p)
#'   beta[1:5]=1
#'   epsilon=rnorm(n) 
#'   y=x%*%beta+epsilon  
#'   b <- sqrootlasso(x,y,42.7,1e4)
#'   }
#' @export


sqrootlasso <- function(x,y,lambda,maxIter) {
  n= nrow(x);p = ncol(x);
  maxerror = 1e-6;
  tmp <- error <- vector();
  #initial beta 
  xx = t(x)%*%x
  xy = t(x)%*%y
  M  = diag(1,p,p)*lambda;
  beta=solve((xx + M))%*%xy; 
  xx = xx/n;    #Gram matrix             
  xy = xy/n;  
  error = y - x%*%beta;     #residuals
  qhat = sum(error^2)/n;  #average of squared residuals
  Iter = 0;
  while (Iter < maxIter){ 
    beta1 = beta;
    Iter = Iter + 1;
    for(j in 1:p){ 
      s0 = sum(xx[j,]%*%beta)-xx[j,j]*beta[j]-xy[j];
      #error = y - x%*%beta +  x[,j]%*%beta[j];
      if (abs(beta[j])>0){
        error = error+x[,j]*beta[j];
        qhat = sum(error^2)/n;
      }
      tmp = (lambda/n)*sqrt(qhat);
      qhat1=qhat-(s0^2)/xx[j,j];
      if(qhat1 <0){qhat1=0;}
      if (s0 > tmp) {
        beta[j] = ((lambda/sqrt(n^2- lambda^2/xx[j,j]))*sqrt(qhat1)-s0 )/xx[j,j];
        error = error-x[,j]*beta[j];
      }
      if (s0< -tmp){
        beta[j] = (-(lambda/sqrt(n^2- lambda^2/xx[j,j]))*sqrt(qhat1)-s0)/xx[j,j];
        error = error-x[,j]*beta[j];
      }
      
      if (abs(s0) <= tmp){ 
        beta[j] = 0;
      }
    } 
    if (norm(beta-beta1,"1") < maxerror ) {
      Iter = maxIter + 1;
    } 
  } 
  return(beta);
}

