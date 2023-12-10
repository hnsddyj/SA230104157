## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
cat(x,'\n')
cat(y,'\n')

n<-length(x)
m<-length(y)
cvmtest <- function(x,y){
#样本x1到xn的ecdf
fn<-function(a){
  flag1<-0
  for(i in 1:length(x)){
    if(x[i]<=a)
      flag1 <- flag1+1
  }
  return(flag1/length(x))
}
#样本y1到ym的ecdf
gm<-function(b){
  flag2<-0
  for(j in 1:length(y)){
    if(y[j]<=b)
      flag2 <- flag2+1
  }
  return(flag2/length(y))
}
#统计量W2
sum1<-sum2<-0
for(i in 1:n){
  sum1<-sum1+(fn(x[i])-gm(x[i]))^2
}
for(j in 1:m){
  sum2<-sum2+(fn(y[j])-gm(y[j]))^2
}
W2<-(sum1+sum2)*m*n/(m+n)^2  
return(W2)
}

R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:26
W2 <- numeric(R) #storage for replicates
W20 <- cvmtest(x,y)
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = 14, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
W2[i] <- cvmtest(x1,y1)
}
p <- mean(c(W20, W2) >= W20)
p

## -----------------------------------------------------------------------------
set.seed(1234)
# Count 5 test
COUNT5test = function(x, y) {
X = x - mean(x)
Y = y - mean(y)
OUTx = sum(X > max(Y)) + sum(X < min(Y))#统计极值点
OUTy = sum(Y > max(X)) + sum(Y < min(X))
return(as.integer(max(c(OUTx, OUTy)) > 5))# 返回1 则拒绝等方差的原假设，返回0则不能拒绝等方差的假设
}

## -----------------------------------------------------------------------------
# Count 5 test permutation
COUNT5test_permutation = function(a) {
n = length(a)
x = a[1:(n/2)]
y = a[-(1:(n/2))]
X = x - mean(x)
Y = y - mean(y)
OUTx = sum(X > max(Y)) + sum(X < min(Y)) 
OUTy = sum(Y > max(X)) + sum(Y < min(X))# 返回1 则拒绝等方差的原假设，返回0则不能拒绝等方差的假设
return(as.integer(max(c(OUTx, OUTy)) > 5))
}

permutation = function(a,R) {
  n = length(a)
  result = numeric(R)
  for (r in 1: R){
      i = sample(1:n ,n ,replace = FALSE)
      result[r] = COUNT5test_permutation(z[i])
  }
  sum(result)/R
}              

## -----------------------------------------------------------------------------
n1 <- 20
n2 <- 50
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m = 10000

Alphahat1 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean
y = y - mean(y)
COUNT5test(x, y)
}))

Alphahat2 = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean 
y = y - mean(y)
z = c(x,y)
permutation(z,1000) 
})<0.05)

round(c(COUNT5test=Alphahat1,COUNT5test_permutation=Alphahat2),5)

