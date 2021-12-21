#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp 
//' @description A Gibbs sampler using Rcpp
//' @param n the number of sample
//' @param n1 the number of between-sample random numbers
//' @param a a constant
//' @param b a constant
//' @return a random sample of size n
//' @importFrom stats rbinom rbeta 
//' @examples
//' \dontrun{
//' n <- 1e4 #length of chain
//' a <- 2
//' b <- 2
//' n1=10
//' mtrix <- matrix(0,n,2)
//' mtrix <- gibbsC(n,n1,a,b)
//' } 
//' @export

// [[Rcpp::export]]
NumericMatrix gibbsC(int n,int n1,double a, double b){
  NumericMatrix mtsample(n, 2);
  mtsample(0,0)=1;mtsample(0,1)=0.1;
  for(int i = 1;i < n ; i++) {
    double y = mtsample(i-1, 1);
    mtsample(i, 0) = rbinom(1, n1, y)[0];
    double x =mtsample(i, 0);
    mtsample(i, 1) = rbeta(1,x+a,n1-x+b)[0];
  }
  return(mtsample);
}
