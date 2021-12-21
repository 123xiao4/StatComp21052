## -----------------------------------------------------------------------------
knitr::kable(summary(women))

## -----------------------------------------------------------------------------
plot(women$height,women$weight,xlab = "Height(in inches)",ylab = "Weight(in lbs)")

## -----------------------------------------------------------------------------
plot(women$height,women$weight,xlab = "Height(in inches)",ylab = "Weight(in lbs)")
mymodel <- lm(weight~height,data=women)
abline(mymodel)
knitr::kable(summary(mymodel)$coef)

## -----------------------------------------------------------------------------
# dev.new()
# par(mfrow=c(2,2))
# plot(mymodel)

## -----------------------------------------------------------------------------
mymodel2 <- lm(weight~height+I(height^2),data=women)
plot(women$height,women$weight,xlab = "Height(in inches)",ylab = "Weight(in lbs)")
mymodel <- lm(weight~height,data=women)
lines(women$height,fitted(mymodel2))
knitr::kable(summary(mymodel2)$coef)

## -----------------------------------------------------------------------------
# dev.new()
# par(mfrow=c(2,2))
# plot(mymodel2)

## -----------------------------------------------------------------------------
# library(forecast)
# library(tseries)
# plot(Nile)
# ndiffs(Nile)
# dnile <- diff(Nile)
# plot(dnile)
# adf.test(dnile)


## -----------------------------------------------------------------------------
# Acf(dnile)
# pacf(dnile)


## -----------------------------------------------------------------------------
# mymodel3 <-  arima(Nile,order=c(0,1,1))
# mymodel3
# knitr::kable(accuracy(mymodel3))

## -----------------------------------------------------------------------------
# qqnorm(mymodel3$residuals)
# qqline(mymodel3$residuals)
# Box.test(mymodel3$residuals,type="Ljung-Box")

## -----------------------------------------------------------------------------
m <- 10000
v <- runif(m)
sigma <- c(0.1,0.5,1,1.5,2,5)
# dev.new()
# par(mfrow=c(3,2))
for(i in 1:length(sigma)){
y <- sqrt(2)*sigma[i]*log(1/(1-v))
za <- seq(0,range(y)[2],0.01)
hist(y,prob=TRUE, main = expression(f(y) == frac(y,sigma^2)*e^(-frac(y^2,sigma^2))))
lines(za,(za/(sigma[i]^2))*exp(-(za^2)/(sigma[i]^2)))
}

## -----------------------------------------------------------------------------
m <- 1000
U1 <- rnorm(m,mean=0,sd=1)
U2 <- rnorm(m,mean=3,sd=1)
v <- runif(m)
t <- as.integer(v>0.75)
U <- t*U1+(1-t)*U2
hist(U,prob=TRUE)

## -----------------------------------------------------------------------------
# l=1e4
# dev.new()
# par(mfrow=c(3,2))
# f <- c(0.25,0.35,0.45,0.55,0.65,0.85)
# for(i in 1:length(f)){
#   U1 <- rnorm(l,mean=0,sd=1)
# U2 <- rnorm(l,mean=3,sd=1)
# w <- runif(l)
# p <- as.integer(w>f[i])
# U <- p*U1+(1-p)*U2
# hist(U,prob=TRUE)
#   }


## -----------------------------------------------------------------------------
lambda <- 2
#set the value of lambda
to <- 2
#set the value of t
an <- numeric(1000)
for(i in 1:1000){ 
  tn <- rexp(1000,lambda) #interarrival times
sn <- cumsum(tn)#arrival times
m <- min(which(sn > to))#arrivals+1 in [0,t]
an[i] <- m-1
}
mean(an)
yn <- rgamma(mean(an),1,1)
# according to the arrivals,generate corresponding number of random samples from gamma distribution,Y_i
xn <- sum(yn) #sum the Y_i,get the value of compound process 
print(xn)

## -----------------------------------------------------------------------------
lambda <- 1
that <- 10
alph <- 1
scal <- 1
xn <- numeric(100)
wn <- numeric(100)
for(j in 1:100){ 
 for(i in 1:100){
 tn <- rexp(100,lambda)
sn <- cumsum(tn)
m <- min(which(sn > that))
xn[i] <- m-1 #generate the value of N(t)
y <- rgamma(mean(xn),alph,scal)#generate y according to the mean of N(t)
wn[j] <- sum(y)
}}
c(mean(wn),lambda*that*alph*scal)#E(X(t))=lambda*t*E(Y),E(Y)=alpha*scale
c(var(wn),lambda*that*alph*scal^2)#Var(X(t))=lambda*t*Var(Y),Var(Y)=alpha*scale^2

## ----echo=FALSE---------------------------------------------------------------
lambda <- 2
that <- 10
alph <- 1
scal <- 1
xn <- numeric(100)
wn <- numeric(100)
for(j in 1:100){ 
 for(i in 1:100){
 tn <- rexp(100,lambda)
sn <- cumsum(tn)
m <- min(which(sn > that))
xn[i] <- m-1 #generate the value of N(t)
y <- rgamma(mean(xn),alph,scal)#generate y according to the mean of N(t)
wn[j] <- sum(y)
}}
c(mean(wn),lambda*that*alph*scal)#E(X(t))=lambda*t*E(Y),E(Y)=alpha*scale
c(var(wn),lambda*that*alph*scal^2)#Var(X(t))=lambda*t*Var(Y),Var(Y)=alpha*scale^2


## ----echo=FALSE---------------------------------------------------------------
lambda <- 3
that <- 10
alph <- 1
scal <- 1
xn <- numeric(100)
wn <- numeric(100)
for(j in 1:100){ 
 for(i in 1:100){
 tn <- rexp(100,lambda)
sn <- cumsum(tn)
m <- min(which(sn > that))
xn[i] <- m-1 #generate the value of N(t)
y <- rgamma(mean(xn),alph,scal)#generate y according to the mean of N(t)
wn[j] <- sum(y)
}}
c(mean(wn),lambda*that*alph*scal)#E(X(t))=lambda*t*E(Y),E(Y)=alpha*scale
c(var(wn),lambda*that*alph*scal^2)#Var(X(t))=lambda*t*Var(Y),Var(Y)=alpha*scale^2


## -----------------------------------------------------------------------------
cumugama <- function(t){
  l <- 1e5
  x <- runif(l,min=0,max=t)
  gamahat <- 30*t*mean((x^2)*((1-x)^2))
  print(gamahat)
}
u <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
gamahat <- numeric(9)
gama <- numeric(9)
for(i in 1:length(u)){
  gamahat[i]=cumugama(u[i])
  gama[i]=pbeta(u[i],3,3)
}
error <- abs(gamahat-gama)
result <- cbind(u,gamahat,gama,error)
knitr::kable(result)

## -----------------------------------------------------------------------------

mcrandom <- function(m,s,antithetic=TRUE){
  w <- runif(m/2)
  if(!antithetic) wq <- runif(m/2)
  else  wq <- 1-w 
  xq <- numeric(m/2)
  x <- sqrt(-2*log(w))*s
  xq <- sqrt(-2*log(wq))*s
  return(list(x=x,xq=xq))
}

l=1e4
simu <- c(1,2,3,4,5)
reduction <- numeric(length(simu))
for(i in 1:length(simu)){
mca <- mcrandom(l,simu[i])
mcb <- mcrandom(l,simu[i],antithetic = FALSE)
v1 <- var(mca$x+mca$xq)
v2 <- var(mcb$x+mcb$xq)
reduction[i] <- (v2-v1)/v2
}
print(reduction)


## ----echo=FALSE,fig.width=10--------------------------------------------------
    x <- seq(1, 100, 0.1)

    g <- x^2 *(1/sqrt(2*pi))*exp(-x^2/2)
    f1 <- 1/sqrt(2*pi)*exp(-x^2/2)
    f2 <- 1/2*(x^2)*exp(-x)
    gs <- c(expression(g(x)==1/sqrt(2*pi)*x^2*e^(-x^2/2)),
            expression(f[1](x)==1/sqrt(2*pi)*e^(-x^2/2)),
            expression(f[2](x)==1/2*x^2*e^(-x)))
    plot(x, g, type = "l", ylab = "",
         ylim = c(0,0.5), col=1)
    lines(x, f1, lty = 2, col=2)
    lines(x, f2, lty = 3, col=3)
    legend("topright", legend = gs,lty=1:3,col=1:3)

## -----------------------------------------------------------------------------
qu <- 1e5
set.seed(121)
w1 <- rnorm(qu)
g1 <- function(x){
  x^2*(x>1)
}
est1 <- mean(g1(w1))
sd1 <- sd(g1(w1))
set.seed(21)
w2 <- rgamma(qu,3,1)
g2 <- function(x){
  2/(sqrt(2*pi))*exp(x-1/2*(x^2))*(x>1)
}
est2 <- mean(g2(w2))
sd2 <- sd(g2(w2))
est1
est2
sd1
sd2

## -----------------------------------------------------------------------------
set.seed(11)
l=1e4
n=20
a=0.05
u=2
m <- numeric(l)
for(j in 1:l){
  x <- rchisq(n,df=2)
  c1 <- mean(x)-sd(x)*qt(1-a/2,df=n-1)/sqrt(n)
  c2 <- mean(x)-sd(x)*qt(a/2,df=n-1)/sqrt(n)
  if(u<=c2 & u>=c1) m[j]=1
  else m[j]=0
}
mean(m)



## -----------------------------------------------------------------------------
##(1) chis-square distribution
cmu <- 1
n=10
l=1e2
hat.cmu <-se.cmu <- numeric()
pva1 <- numeric()
set.seed(123)
for(j in 1:l){
  x <- rchisq(n,1)
 hat.cmu[j] <- mean(x)
 se.cmu[j] <- sd(x)
 test.t <- t.test(x,mu=cmu)
 pva1[j] <- test.t$p.value
}
pva2 <- 2*(1-pt(abs(sqrt(n)*(hat.cmu-cmu)/se.cmu),df=n-1))
mean(pva1<=0.05)
mean(pva2<=0.05)
mean(pva1<=0.01)
mean(pva2<=0.01)

## -----------------------------------------------------------------------------
##(2) uniform distribution on [0,2]
umu <- 1
n=10
l=1e2
hat.umu <-se.umu <- numeric()
pvb1 <- numeric()
set.seed(111)
for(j in 1:l){
  x <- runif(n,min=0,max=2)
 hat.umu[j] <- mean(x)
 se.umu[j] <- sd(x)
 test.t <- t.test(x,mu=umu)
 pvb1[j] <- test.t$p.value
}
pvb2 <- 2*(1-pt(abs(sqrt(n)*(hat.umu-umu)/se.umu),df=n-1))
mean(pvb1<=0.05)
mean(pvb2<=0.05)
mean(pvb1<=0.01)
mean(pvb2<=0.01)

## -----------------------------------------------------------------------------
##(3) Exponential(rate=1)
emu <- 1
n=10
l=1e2
hat.emu <-se.emu <- numeric()
pvc1 <- numeric()
set.seed(123)
for(j in 1:l){
  x <- runif(n,min=0,max=2)
 hat.emu[j] <- mean(x)
 se.emu[j] <- sd(x)
 test.t <- t.test(x,mu=emu)
 pvc1[j] <- test.t$p.value
}
pvc2 <- 2*(1-pt(abs(sqrt(n)*(hat.emu-emu)/se.emu),df=n-1))
mean(pvc1<=0.05)
mean(pvc2<=0.05)
mean(pvc1<=0.01)
mean(pvc2<=0.01)



## -----------------------------------------------------------------------------
library(MASS)
#write a function to calculate the value of the multivariate skewness statistic and judge if it reject the null hypothesis
bd <- function(x){
  ns <- nrow(x)
  d <- ncol(x)
  y <- x
  for(i in 1:d){
    y[,i]=x[,i]-mean(x[,i])
  }
  sigmamax <- t(y)%*%y/ns
  nsigmahat <- solve(sigmamax)
  smatrix <- y%*%nsigmahat%*%t(y)
  bvalue <- sum(rowSums(smatrix^3))/(ns^2)
  freedom <- d*(d+1)*(d+2)/6
  bound <- qchisq(0.95,freedom)
  newb=ns*bvalue/6
  ifreject <- as.integer(newb>bound)
  return(ifreject)
}
# set the parameter for sample of multiple-normal distribution
#I set d=3
setmean <- c(0,0,0)
setcov <- matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,nrow=3)
nt=1e3#the number of repetition
ns <- c(30,50,100)#the vector of the number of sample 
nns <- length(ns)
probreject <- numeric(nns)
set.seed(80512)
for(j in 1:nns){
  prreject <- numeric(nt)
for(k in 1:nt){
samplex <- mvrnorm(ns[j],setmean,setcov)
prreject[k] <- bd(samplex)
}
probreject[j]=mean(prreject)
}
ns <- as.integer(ns)
probreject 
table <- rbind(as.integer(ns),round(probreject,3))
table

## -----------------------------------------------------------------------------

library(MASS)
set.seed(202110)
setmean <- c(0,0)
proportion <- c(seq(0,.15,.01),seq(.16,1,0.05))
np <- length(proportion)
proreject <- numeric(np)
nt=1e2
ns=100
setcov1 <- matrix(c(1,0,0,1),ncol=2,nrow=2)
setcov2 <- matrix(c(100,0,0,100),ncol=2,nrow=2)
for(k in 1:np){
  bdt <- numeric(nt)
  pro <- proportion[k]
      for(i in 1:nt){
        U1 <- mvrnorm(ns,setmean,setcov1)
        U2 <- mvrnorm(ns,setmean,setcov2)
        v <- runif(ns)
        t <- as.integer(v<1-pro)
        samplex <- t*U1+(1-t)*U2
        bdt[i] <- bd(samplex)
  }
proreject[k] <- mean(bdt)
}
proreject
#plot power vs epsilon
plot(proportion,proreject,type = 'b',xlab=bquote(epsilon),ylab="power",ylim=c(0,1))
abline(h=0.1,lty=3)
hatse <- sqrt(proreject*(1-proreject)/nt)
lines(proportion, proreject+hatse,lty=3)
lines(proportion,proreject-hatse,lty=3)

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor,package = "bootstrap")# load the data
n <- nrow(scor)
sigma <- (n-1)/n*cov(scor) #the MLE of covariance matrix
lambda <- eigen(sigma)$values
theta_hat <- lambda[1]/sum(lambda)
print(theta_hat)

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor,package = "bootstrap")# load the data
m <- 1000
theta.b <- numeric(m)
set.seed(1029)
for (i in 1:m){
  j <- sample(1:n,size=n,replace=TRUE)
  mech <- scor$mec[j]
  vect <- scor$vec[j]
  alge <- scor$alg[j]
  anal <- scor$ana[j]
  stat <- scor$sta[j]
  da <- data.frame(mech,vect,alge,anal,stat)
  bsigma <- (n-1)/n*cov(da)
  blambda <- eigen(bsigma)$values
  theta.b[i] <- blambda[1]/sum(blambda)
}
bbias <- mean(theta.b)-theta_hat
bbias
bsdhat <- sd(theta.b)
bsdhat


## -----------------------------------------------------------------------------
theta.jack <- numeric(n)
for(i in 1:n){
  da <- scor[-i,]
  newsigma <- (n-1)/n*cov(da)
  newlambda <- eigen(newsigma)$values
  theta.jack[i] <- newlambda[1]/sum(newlambda)
}
jbias <- (n-1)*(mean(theta.jack)-theta_hat)
print(jbias)
sdhat <- ((n-1)/sqrt(n))*sd(theta.jack)
print(sdhat)



## -----------------------------------------------------------------------------
library(bootstrap)
data(scor,package = "bootstrap")
n=nrow(scor)
m <- 1e2
theta.b <- numeric(m)
set.seed(1029)
km=1e2
percent1 <- percent2 <- bca1 <- bca2 <- numeric(km)
for(k in 1:km){
for (i in 1:m){
  j <- sample(1:n,size=n,replace=TRUE)
  mech <- scor$mec[j]
  vect <- scor$vec[j]
  alge <- scor$alg[j]
  anal <- scor$ana[j]
  stat <- scor$sta[j]
  da <- data.frame(mech,vect,alge,anal,stat)
  bsigma <- (n-1)/n*cov(da)
  blambda <- eigen(bsigma)$values
  theta.b[i] <- blambda[1]/sum(blambda)
}
##percentile c.i.
percent1[k]=quantile(theta.b,0.25)
percent2[k]=quantile(theta.b,0.975)

## bca c.i.
a1=sum((mean(theta.jack)-theta.jack)^3)
thetanew<-(mean(theta.jack)-theta.jack)^2
a2=sum(thetanew^1.5)
a=a1/a2
b=mean(theta.b<theta_hat)
hatz0=qnorm(b)
z=qnorm(0.25)
c1=z+hatz0
c2=hatz0-z
d1=1-a*c1
d2=1-a*c2
a1value=pnorm(hatz0+c1/d1)
a2value=pnorm(hatz0+c2/d2)
bca1[k]=quantile(theta.b,a1value)
bca2[k]=quantile(theta.b,a2value)
}
mean(percent1)
mean(percent2)
mean(bca1)
mean(bca2)

## -----------------------------------------------------------------------------
library(boot)
skness <- function(x,i){
  barx <- mean(x[i,])
  hatm2 =mean((x[i,]-barx)^2)
  hatm3=mean((x[i,]-barx)^3)
  return(hatm3/hatm2^1.5)
}
#normal population
sknorm=0
l=100
times=100
cinormb<-cinormbasic<-cinormperc<-matrix(NA,times,2)
set.seed(1122)
for(i in 1:times){
spnorm <- as.matrix(rnorm(l))
spnorm.b <- boot(spnorm,statistic=skness,R=100)
normb.ci <- boot.ci(spnorm.b,type=c("norm","basic","perc"))
cinormb[i,]<-normb.ci$norm[2:3]
cinormbasic[i,]<-normb.ci$basic[4:5]
cinormperc[i,]<-normb.ci$percent[4:5]
}
norm.norm =mean(cinormb[,1]<=sknorm & cinormb[,2]>=sknorm)
basic.norm =mean(cinormbasic[,1]<=sknorm & cinormbasic[,2]>=sknorm)
percentile.norm=mean(cinormperc[,1]<=sknorm & cinormperc[,2]>=sknorm)
cat('norm.norm =',norm.norm,'basic.norm =',basic.norm,'percentile.norm=',percentile.norm)

missnorml=mean(cinormb[,1]>sknorm)
missnormr=mean(cinormb[,2]<sknorm)
missbcal=mean(cinormbasic[,1]>sknorm)
missbcar=mean(cinormbasic[,2]<sknorm)
misspercl=mean(cinormperc[,1]>sknorm)
misspercr=mean(cinormperc[,2]<sknorm)
missleft.norm <- c(missnorml,missbcal,misspercl)
missright.norm <- c(missnormr,missbcar,misspercr)
miss.norm <- cbind(missleft.norm,missright.norm)
rownames(miss.norm)<-c("norm","bca","percentile")
miss.norm
##chisquare population
skchi=sqrt(8/5)
l=100
times=100
cichib<-cichibasic<-cichiperc<-matrix(NA,times,2)
set.seed(1022)
for(i in 1:times){
spchi <- as.matrix(rchisq(l,5))
spchi.b <- boot(spchi,statistic=skness,R=100)
chib.ci <- boot.ci(spchi.b,type=c("norm","basic","perc"))
cichib[i,]<-chib.ci$norm[2:3]
cichibasic[i,]<-chib.ci$basic[4:5]
cichiperc[i,]<-chib.ci$percent[4:5]
}
norm.chi =mean(cichib[,1]<=skchi & cichib[,2]>=skchi)
basic.chi =mean(cichibasic[,1]<=skchi & cichibasic[,2]>=skchi)
percentile.chi=mean(cichiperc[,1]<=skchi & cichiperc[,2]>=skchi)
cat('norm.chi =',norm.chi,
    'basic.chi =',basic.chi,
    'percentile.chi=',percentile.chi)

missnorml1=mean(cichib[,1]>skchi)
missnormr1=mean(cichib[,2]<skchi)
missbcal1=mean(cichibasic[,1]>skchi)
missbcar1=mean(cichibasic[,2]<skchi)
misspercl1=mean(cichiperc[,1]>skchi)
misspercr1=mean(cichiperc[,2]<skchi)
missleft.chi<- c(missnorml1,missbcal1,misspercl1)
missright.chi <- c(missnormr1,missbcar1,misspercr1)
miss.chi <- cbind(missleft.chi,missright.chi)
rownames(miss.chi)<-c("norm","bca","percentile")
miss.chi

difference <- c(norm.norm-norm.chi,basic.norm-basic.chi,percentile.norm-percentile.chi)
difference

## -----------------------------------------------------------------------------
mpg2=c(21.0,22.8,21.4 ,18.7, 18.1, 14.3, 24.4, 19.2, 17.3 ,15.2 ,10.4, 14.7, 32.4 ,30.4,33.9,21.5, 15.5)
wt2<-c(2.620,2.320,3.215 ,3.440, 3.460, 3.570, 3.190,3.780, 4.070, 3.730 ,5.280, 5.424,5.345, 2.200, 1.615, 1.835, 2.465)
ns=length(mpg2)
initial <- cor(mpg2,wt2,method="spearman")
nt=1e3
simulation <- numeric(nt)
set.seed(1234)
for(i in 1:nt){
  newmpg <- sample(mpg2,size=ns,replace=F)
  newwt<- sample(wt2,size=ns,replace=F)
  simulation[i]=cor(newmpg,newwt,method="spearman")
}
simu=mean(abs(c(initial,simulation))>=abs(initial))
exsit=cor.test(mpg2,wt2,method="spearman")$p.value
round(c(simu,exsit),6)


## -----------------------------------------------------------------------------
library(Ball)
library(energy)
library(RANN)
library(boot)
statn <- function(x, ix, ns,k) {
  la <- ns[1]; lb <- ns[2]; l <- la + lb
  if(is.vector(x)) x <- data.frame(x,0);
  x <- x[ix, ];
  NN <- nn2(data=x, k=k+1)
  regiona <- NN$nn.idx[1:la,-1] 
  regionb <- NN$nn.idx[(la+1):l,-1] 
  ia <- sum(regiona < la+.5);ib <- sum(regionb >lb+.5) 
  (ia+ib)/(k*l)
}

test.nn<- function(z,ns,k){
  simu.boot <- boot(data=z,statistic=statn,R=R,
  sim = "permutation", ns = ns,k=k)
  st <- c(simu.boot$t0,simu.boot$t)
  pvalue <- mean(st>=st[1])
  list(statistic=st[1],p.value=pvalue)
}
nsimu <- 1e2; k<-3; p<-2; 
set.seed(1199)
nx<-ny <- 30;n <- nx+ny; N = c(nx,ny)
nt<-100; R=100
pvalue <- matrix(NA,nsimu,3)
for(i in 1:nsimu){
  xsample <- matrix(rnorm(nx*p,0,1.5),ncol=p);
  ysample <- cbind(rnorm(ny),rnorm(ny));
  ztest <- rbind(xsample,ysample)
  pvalue[i,1] <- test.nn(ztest,N,k)$p.value
  pvalue[i,2] <- eqdist.etest(ztest,sizes=N,R=nt)$p.value
  pvalue[i,3] <- bd.test(x=xsample,y=ysample,num.permutations=100,seed=i*1122)$p.value
}
alpha <- 0.1; 
getpower1 <- colMeans(pvalue<alpha)
getpower <- as.data.frame(getpower1,nrow=3,row.names=c("NN",'Energy','Ball'))
colnames(getpower)="p-value"
knitr::kable(getpower)

## -----------------------------------------------------------------------------
nsimu2 <- 100; k<-3; p<-2; u=0.5;set.seed(1199)
nx <- n2 <- 10; nt<-100; n <- nx+ny; Nc = c(nx,ny)
pvalue <- matrix(NA,nsimu2,3)
for(i in 1:nsimu2){
  xsample <- matrix(rnorm(nx*p,0,1.2),ncol=p);
  ysample <- cbind(rnorm(ny),rnorm(ny,u));
  ztest <- rbind(xsample,ysample)
  pvalue[i,1] <- test.nn(ztest,Nc,k)$p.value
  pvalue[i,2] <- eqdist.etest(ztest,sizes=Nc,R=nt)$p.value
  pvalue[i,3] <- bd.test(x=xsample,y=ysample,num.permutations=1e2,seed=i*2120)$p.value
}
alpha <- 0.1; 
getpower2 <- colMeans(pvalue<alpha)
getpower <- as.data.frame(getpower2,nrow=3,row.names=c("NN",'Energy','Ball'))
colnames(getpower)="p-value"
knitr::kable(getpower)

## -----------------------------------------------------------------------------
nsimu3 <- 100; k<-3; p<-2; u=0.5;set.seed(1199)
n1 <- n2 <- 10; nt<-1e2; n <- n1+n2; N = c(n1,n2)
pvalue <- matrix(NA,nsimu3,3)
for(i in 1:nsimu3){
  xsample <- matrix(rt(n1*p,1),ncol=p);
 U1 <- rnorm(n2,mean=0,sd=1)
 U2 <- rnorm(n2,mean=3,sd=1)
 v1 <- runif(n2)
 t1 <- as.integer(v1<0.7)
 U <- t1*U1+(1-t1)*U2
 U3 <- rnorm(n2,mean=0,sd=1)
 U4 <- rnorm(n2,mean=3,sd=1)
 v2 <- runif(n2)
 t2 <- as.integer(v2<0.7)
 V <- t2*U3+(1-t2)*U4
 ysample <- cbind(U,V);
  ztest <- rbind(xsample,ysample)
  pvalue[i,1] <- test.nn(ztest,N,k)$p.value
  pvalue[i,2] <- eqdist.etest(ztest,sizes=N,R=nt)$p.value
  pvalue[i,3] <- bd.test(x=xsample,y=ysample,num.permutations=1e3,seed=i*1122)$p.value
}
alpha <- 0.1; 
getpower3 <- colMeans(pvalue<alpha)
getpower <- as.data.frame(getpower3,nrow=3,row.names=c("NN",'Energy','Ball'))
colnames(getpower)="p-value"
knitr::kable(getpower)

## -----------------------------------------------------------------------------
nsimu4 <- 100; k<-3; p<-2; u=1.2;set.seed(190)
nx <- 30; ny <- 50; nt<-1e3; n <- nx+ny; Nc = c(nx,ny)
pvalue <- matrix(NA,nsimu4,3)
R=999
for(i in 1:nsimu4){
  xsample <- matrix(rnorm(nx*p,1.2,1),ncol=p);
  ysample <- cbind(rnorm(ny),rnorm(ny,u));
  ztest <- rbind(xsample,ysample)
  pvalue[i,1] <- test.nn(ztest,Nc,k)$p.value
  pvalue[i,2] <- eqdist.etest(ztest,sizes=Nc,R=nt)$p.value
  pvalue[i,3] <- bd.test(x=xsample,y=ysample,num.permutations=1e3,seed=i*1522)$p.value
}
alpha <- 0.1; 
getpower4 <- colMeans(pvalue<alpha)
getpower <- as.data.frame(getpower4,nrow=3,row.names=c("NN",'Energy','Ball'))
colnames(getpower)="p-value"
knitr::kable(getpower)

## -----------------------------------------------------------------------------
#write a function to generate the chain
samplemh <- function(ns,xinit){
  mhx <- numeric(ns)
  mhx[1] <- xinit
  ut <- runif(ns)
   for (j in 2:ns) {
        xt <- mhx[j-1]
        y <- rnorm(1,xt,1)
        pv <- dcauchy(y)*dnorm(xt,y,1)/(dcauchy(xt)*dnorm(y,xt,1))
        if (ut[j] <= pv){ mhx[j] = y} 
        else {mhx[j]=xt} 
   }
 return(mhx)
}
set.seed(20211)
xit <- rnorm(1,0,1)#set initial value
ns=2e4
mhy <- samplemh(ns,xit)
plot(mhy[15000:15500], type="l", main="Time series of simulation", ylab="x")

#compare deciles
prb <- seq(0.1,0.9,0.1)
dq1 <- quantile(mhy[1001:ns],probs = prb)
dq2 <- qcauchy(prb)
dq1-dq2

## -----------------------------------------------------------------------------
#write the function to compute G-R statistic 
GR <- function(phi) {
        phi <- as.matrix(phi)
        m <- ncol(phi)
        n <- nrow(phi)
phimeans <- rowMeans(phi)     
        B <- m*var(phimeans) 
    phiw <- apply(phi, 1, "var")
        w <- mean(phiw)        
        vhat <- w*(m-1)/m + (B/m)
        Rhat <- vhat/w  
        return(Rhat)
}

## -----------------------------------------------------------------------------
    k <- 3         
    ns <- 1e4     
    bn <- 1e3   
    set.seed(1321)
    xint <- rnorm(k)
    spx <- matrix(0, nrow=k, ncol=ns)
    for (i in 1:k)
        spx[i, ] <- samplemh(ns,xint[i])

    #compute diagnostic statistics
    phi <- t(apply(spx, 1, cumsum))
    for (i in 1:nrow(phi))
        phi[i,] <- phi[i,] / (1:ncol(phi))


    for (i in 1:k)
      if(i==1){
        plot((bn+1):ns,phi[i, (bn+1):ns], type="l",ylim=c(-1.5,1),
            xlab='Index', ylab=bquote(phi))
      }else{
        lines(phi[i, (bn+1):ns], col=i)
    }
  
 #plot the sequence of R-hat statistics
  rhat <- rep(0, ns)
    for (j in (bn+1):ns)
        rhat[j] <- GR(phi[,1:j])
    plot(rhat[(bn+1):ns], type="l", xlab="", ylab="R")
    abline(h=1.2, lty=2)


## -----------------------------------------------------------------------------
#write a function to generate the chain by using Gibbs sample
bichain <- function(nt,n,a,b){
  sample <- matrix(0, nt, 2)
  x1 <- runif(1,min=0,max=n)
  y1 <- runif(1)
  sample[1, ] <- c(x1,y1) #initialize
  for (j in 2:nt) {
  y <- sample[j-1, 2]
  sample[j, 1] <- rbinom(1, n, y)
  x <- sample[j, 1]
  sample[j, 2] <- rbeta(1,x+a,n-x+b)
  }
  return(sample)
}
#initialize constants and parameters
nt <- 1e4 #length of chain
a <- 2
b <- 2
n=10
burn <- 1e3 #burn-in length
bn <- burn+1
#generate a chain with target joint density f(x,y)
spxy <- bichain(nt,n,a,b)[(bn:nt),]


## -----------------------------------------------------------------------------
set.seed(1321)
spxy1 <- bichain(nt,n,a,b)
spxy2 <- bichain(nt,n,a,b)
spxy3 <- bichain(nt,n,a,b)

sspx=rbind(spxy1[,1],spxy2[,1],spxy3[,1])
sspy=rbind(spxy1[,2],spxy2[,2],spxy3[,2]) 
phi1 <- t(apply(sspx, 1, cumsum))
phi2 <- t(apply(sspy, 1, cumsum))
for (i in 1:nrow(phi1))
{
phi1[i,] <- phi1[i,] / (1:ncol(phi1))
phi2[i,] <- phi2[i,] / (1:ncol(phi2))

}
print(c(GR(phi1),GR(phi2)))
r.bino <- numeric(nt)
r.beta <-numeric(nt)

for(l in(bn):nt){
r.bino[l]<-GR(phi1[,1:l])
r.beta[l]<-GR(phi2[,1:l])
}
par(mfrow=c(1,2)) 
plot(r.bino[(bn+1):nt],type="l",xlab="N",ylab="R")
abline(h=1.2,lty=2)
plot(r.beta[(bn+1):nt],type="l",xlab="N",ylab="R")
abline(h=1.2,lty=2)


## -----------------------------------------------------------------------------
# the function to compute the kth term
termk <- function(a,k){
  d=length(a)
  odist <- norm(a,type = "2")
  s=0
  for(i in 1:k){s=s+log(i)}
  tk <- (-1)^k *exp((k+1)*log(odist^2)+lgamma((d+1)/2)+lgamma(k+1.5)-s-k*log(2)-log(2*k+1)-log(2*k+2)-lgamma(k+d/2+1))
  return(tk)
}
a=c(1,2)
k=100
termk(a,k)


## -----------------------------------------------------------------------------
#(b) the function that it computes and returns the sum
fsum <- function(a){
  od <- norm(a,type="2")
  n=length(a)
  initn=od^2*gamma((n+1)/2)*gamma(1.5)
  initd=2*gamma(n/2+1)
  init=initn/initd
  for(k in 1:1e4){
    if(abs(termk(a,k))<1e-24){init=init}
    else{init <- init+termk(a,k)
    }
  }
    return(init)
}

#(c)
a=c(1,2)
fsum(a)


## -----------------------------------------------------------------------------


upbd <- function(a,k){sqrt(a^2*k/(k+1-a^2))}
encore <- function(x,k){(1+x^2/k)^(-(k+1)/2)}
coe <- function(k){2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))}
sofunc <- function(a,k){
  if(!is.vector(a)) a <- c(a)
  n=length(a)
  out <- numeric(n)
  for(i in 1:n){
  inte1<-integrate(encore,lower=0,upper=upbd(a[i],k-1),k=k-1)$value
inte2 <- integrate(encore,lower=0,upper=upbd(a[i],k),k=k)$value
out[i] <- coe(k)*inte2-coe(k-1)*inte1
}
out
}
result1 <- numeric(15)
for(k in 4:20){
result1[k-3] <- uniroot(sofunc,interval = c(1e-5,(k-1)^0.5), k=k)$root
}
result1

## -----------------------------------------------------------------------------
#11.4
ga <- function(a,k){
  pt(upbd(a,k),df=k)-pt(upbd(a,k-1),df=k-1)
}
result2 <- numeric(15)
for(k in 4:20){
  result2[k-3] <- uniroot(ga,interval = c(1e-6,(k-1)^0.5), k=k)$root
}
result1-result2

## -----------------------------------------------------------------------------

y<-c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
sum(y)
length(y)


## -----------------------------------------------------------------------------

trims <- c(0,0.1,0.2,0.5)
x <- rcauchy(100)
lapply(trims,function(trim) mean(x,trim=trim))
lapply(trims,mean,x=x)


## -----------------------------------------------------------------------------
#extract r_2 in exercise 4
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
formu <- lapply(seq_along(formulas), function(i) lm(formulas[[i]],data=mtcars))
rsq <- function(mod) summary(mod)$r.squared
lapply(seq_along(formu), function(i) rsq(formu[[i]]))


##extract r_2 in exercise 5
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars.new <- mtcars[rows, ]
lm(mpg~disp,data=mtcars.new)
})

lapply(1:10, function(i) rsq(bootstraps[[i]]))

## -----------------------------------------------------------------------------
#(a)

da <- matrix(c(1:20),nrow=5)
da <- as.data.frame(da)
vapply(da,sd,FUN.VALUE = 0) 
#existing data frame
vapply(mtcars, sd, FUN.VALUE = 0)

## -----------------------------------------------------------------------------
#(b)
#when the non-numeric variable is easy to drop
vapply(iris,is.numeric,logical(1))
iris.new <- iris[,-5]
vapply(iris.new,sd,FUN.VALUE = 0)

#write a function to come true
sd.mix <- function(x) vapply(x[vapply(x,is.numeric,logical(1))], sd, FUN.VALUE = 0)

sd.mix(iris)


## -----------------------------------------------------------------------------
library(parallel)
mcsapply <- function(cluster = NULL, X, FUN, ..., simplify = TRUE,USE.NAMES = TRUE){
   FUN <- match.fun(FUN)
  answer <- parallel::parLapply(cluster,X = X, fun = FUN, ...)
if (USE.NAMES && is.character(X) && is.null(names(answer))) names(answer) <- X
if (!identical(simplify, FALSE) && length(answer))
      simplify2array(answer, higher = (simplify == "array"))
  else answer
} 
cores <- detectCores()
cluster <- makePSOCKcluster(cores)
mcsapply(cluster, 1:20, get("+"),3)
da <- matrix(1:8*1024,nrow=8)
da <- as.data.frame(da)
system.time(mcsapply(cluster,da,sd))

## -----------------------------------------------------------------------------
###########Write a Rcpp function for Exercise 9.8
library(Rcpp)
#direct_cpp <- "D:/Rcpp/"
#sourceCpp(paste0(direct_cpp, "gibbs_cpp.cpp")) 
#initialize constants and parameters
cppFunction('NumericMatrix gibbsC(int n,int n1,double a, double b){
  NumericMatrix mtsample(n, 2);
    mtsample(0,0)=1;mtsample(0,1)=0.1;
    for(int i = 1;i < n ; i++) {
 double y = mtsample(i-1, 1);
  mtsample(i, 0) = rbinom(1, n1, y)[0];
  double x =mtsample(i, 0);
  mtsample(i, 1) =rbeta(1,x+a,n1-x+b)[0];
  }
    
    return(mtsample);
}
'
)
n <- 1e4 #length of chain
a <- 2
b <- 2
n1=10
mtrix <- matrix(0,n,2)
mtrix <- gibbsC(n,n1,a,b)

###################
bichain <- function(nt,n,a,b){
  sampleb <- matrix(0, nt, 2)
  x1 <- runif(1,min=0,max=n)
  y1 <- runif(1)
  sampleb[1, ] <- c(x1,y1) #initialize
  for (j in 2:nt) {
  y <- sampleb[j-1, 2]
  sampleb[j, 1] <- rbinom(1, n, y)
  x <- sampleb[j, 1]
  sampleb[j, 2] <- rbeta(1,x+a,n-x+b)
  }
  return(sampleb)
}
n <- 1e4 #length of chain
a <- 2
b <- 2
n1=10
mtrix <- matrix(0,n,2)
mtrix <- gibbsC(n,n1,a,b)
mtrix1 <- matrix(0,n,2)
mtrix1 <- bichain(n,n1,a,b)
qqplot(mtrix[,1],mtrix1[,1],xlab = "quantiles of gibbC ",ylab = "quantiles of gibbR",main="qq-plot of x")
qqplot(mtrix[,2],mtrix1[,2],xlab = "quantiles of gibbC ",ylab = "quantiles of gibbR",main="qq-plot of y",xlim=c(0,1),ylim=c(0,1))

## -----------------------------------------------------------------------------
library(microbenchmark)
nt <- c(100,500,1e3,5e3,1e4)
mt =length(nt)
ts <- vector("list", length=mt)
sumts <- vector("list", length=mt)
for(i in 1:mt){
ts[[i]] <- microbenchmark(gibbC=gibbsC(nt[i],10,2,2), 
                         gibbR=bichain(nt[i],10,2,2))
    sumts[[i]]<-summary(ts[[i]])[,c(1,3,5,6)]
}
sumts

