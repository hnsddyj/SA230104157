## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
U <- c(11,8,27,13,16,0,23,10,24,2)
V <- c(12,9,28,14,17,1,24,11,25,3)
l <- function(Lambda){
  l <- 0
  for (i in 1:10){
    n <- max(1.5e-5,pexp(V[i],rate = Lambda)-pexp(U[i],rate = Lambda))
    l <- l + log(n)
  }
  return(l)
}

optimize(l,lower=1e-3,upper=5,maximum=T)

EX <- function(u,v,lambda){
  b = u*exp(-lambda*u)-v*exp(-lambda*v)+(exp(-lambda*u)-exp(-lambda*v))/lambda
  a = pexp(v,lambda)-pexp(u,lambda)
  return(b/a)
}

EM <- function(u,v,max.it=1e4){
  Lambda1 = 10
  Lambda2 = 1
  y = numeric(10)
  while (abs(Lambda1-Lambda2) >= 1.5e-5) {
    Lambda1 = Lambda2
    for (i in 1:10) {
      y[i] = EX(u[i],v[i],Lambda1)
    }
    Lambda2 = 10/sum(y)
  }
  return(Lambda2)
}
EM(U,V)

## -----------------------------------------------------------------------------
solve.game <- function(A) {
  min_A <- min(A)
  A <- A - min_A #so that v >= 0
  max_A <- max(A)
  A <- A / max(A)
  m <- nrow(A)
  n <- ncol(A)
  it <- n^3
  a <- c(rep(0, m), 1) #objective function
  A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
  b1 <- rep(0, n)
  A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
  b3 <- 1
  sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
                maxi=TRUE, n.iter=it)
  #the ’solution’ is [x1,x2,...,xm | value of game]
  #
  #minimize v subject to ...
  #let y strategies 1:n, with v as extra variable
  a <- c(rep(0, n), 1) #objective function
  A1 <- cbind(A, rep(-1, m)) #constraints <=
  b1 <- rep(0, m)
  A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
  b3 <- 1
  sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
                maxi=FALSE, n.iter=it)
  soln <- list("A" = A * max_A + min_A,
               "x" = sx$soln[1:m],
               "y" = sy$soln[1:n],
               "v" = sx$soln[m+1] * max_A + min_A)
soln
}
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
               2,0,0,0,-3,-3,4,0,0,
               2,0,0,3,0,0,0,-4,-4,
               -3,0,-3,0,4,0,0,5,0,
               0,3,0,-4,0,-4,0,5,0,
               0,3,0,0,4,0,-5,0,-5,
               -4,-4,0,0,0,5,0,0,6,
               0,0,4,-5,-5,0,0,0,6,
               0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) #needed for simplex function
S <- solve.game(A+2)
round(cbind(S$x, S$y), 7)
S$v
#通过输入收益矩阵A运行结果符合11.15这条结论，收益为0；B<-A+2运行结果符合11.15这条结论，收益为2；



