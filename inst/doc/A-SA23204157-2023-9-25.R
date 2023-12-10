## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 2e5
x1 <- runif(m)
theta_hat_1 <- exp(x1)

## -----------------------------------------------------------------------------
set.seed(123)
m <- 2e5
x2 <- runif(m);
theta_hat_2 <- (exp(x2) + exp(1-x2)) / 2

## -----------------------------------------------------------------------------
a <- 1-var(theta_hat_2) / var(theta_hat_1)
a

## -----------------------------------------------------------------------------
#编写一个函数来求不同ρ下的方差，这里不妨设d=1
my.simulation <- function(a){
set.seed(12345)
l <- a
d<- 1
m <- 1e6
j <- 1
k<- 100
pihat<-numeric(k)
while(j<=k)
{X <- runif(m,0,d/2)
Y <- runif(m,0,pi/2)
pihat[j] <- 2*l/d/mean(l/2*sin(Y)>X)
j <-j+1
}
return(var(pihat))
}
##求ρ=0.3、0.6和1的情况
my.simulation(0.3)
my.simulation(0.6)
my.simulation(1)

