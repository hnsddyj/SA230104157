## -----------------------------------------------------------------------------
SIS <- function(x,y,d,xtrue){
  SIS_flag <- 0
  n=nrow(x)
  p=ncol(x)
  R=c(rep(0,p))
  for(i in 1:p){
    R[i]=(sum(((y-mean(y))*x[,i]))/sd(x[,i]))^2
  }
  r=order(R,decreasing=T)
  if(all( xtrue %in% r[1:d]))
    SIS_flag <- 1
  return(list(select=r[1:d],SIS_flag=SIS_flag))
}

SIS.SCAD <- function(x,y,d,xtrue){
  p <- ncol(x)
  beta_SIS.SCAD <- c(rep(0,p))
  m <- SIS(x,y,d,c(1,2,3))$select
  fit_scad <- cv.ncvreg(x[,m],y,penalty="SCAD")
  beta_SIS.SCAD[m]<-as.vector(coef(fit_scad))[2:(d+1)]
  A <- which(abs(beta_SIS.SCAD)>0)
  size <- length(A)
  return(list(A=A,size=size))
}

ISIS <- function(x,y,d,xtrue){
  size <- 0
  ISIS_flag <- 0
  AA <- c()
  yy=y
  xx=x
  while(size<d)
  {
    B <- SIS.SCAD(x,y,d,xtrue)
    m <- B$size
    if(m>0){
      ISIS_data <- data.frame(y,x[,B$A])
      model <-lm(y~x[,B$A],data = ISIS_data)
      AA <- c(B$A,AA)
      size <- size + m
      x <- x[,-B$A]
      residual <- residuals(model)
      y <- residual
    }
    else break
  }
  x=xx
  y=yy
  if(all( xtrue %in% AA))
    ISIS_flag <- 1
  return(list(ISIS_flag=ISIS_flag,AA=AA))
}

## -----------------------------------------------------------------------------
set.seed(1234)
library(MASS)
library(glmnet)
library(ncvreg)
SIS_flag1  <- 0
ISIS_flag1  <- 0
SIS_result1 <- 0
ISIS_result1 <- 0
p <- 100
n <- 20
rho <- 0.1
k <- 50
xtrue<-c(1,2,3)
d <- 2.5*floor(n/log(n))
for(i in 1:k){
  beta <- c(5,5,5,rep(0,(p-3)))
  mean <- c(rep(0,p))
  sigma <- matrix(c(rep(0.5,p*p)), nrow=p, ncol=p)+diag(p)*(1-rho)
  x <- mvrnorm(n, mean, sigma)
  y <- x%*%beta+rnorm(n)
  #SIS
  SIS_flag1 <- SIS_flag1+SIS(x,y,d,c(1,2,3))$SIS_flag
  #ISIS
  ISIS_flag1 <- ISIS_flag1+ISIS(x,y,d,c(1,2,3))$ISIS_flag
}
SIS_result1 <- SIS_flag1/k
ISIS_result1 <- ISIS_flag1/k

SIS_result1
ISIS_result1

