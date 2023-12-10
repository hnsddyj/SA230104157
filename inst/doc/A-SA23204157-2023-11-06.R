## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(1234)
f<-function(N,b1,b2,b3,f0){
  x1<-rpois(N, 1)
  x2<-rexp(N,1)
  x3<-sample(0:1,N,replace = TRUE)
  g<-function(alpha){
    tmp<-exp(-alpha-b1*x1-b2*x2-b3*x3)
  p<-1/(1+tmp)
  mean(p)-f0}
  result1<-uniroot(g,c(-15,0))
  return(unlist(result1))
}
a1<-f(10^6,0,1,-1,0.1)[1]
a2<-f(10^6,0,1,-1,0.01)[1]
a3<-f(10^6,0,1,-1,0.001)[1]
a4<-f(10^6,0,1,-1,0.0001)[1]
a<-cbind(a1,a2,a3,a4)
f0<-cbind(0.1,0.01,0.001,0.0001)
plot(log(f0),a)
rownames(a) = "output a"
colnames(a) = paste("f0",f0)
knitr::kable(a)

## -----------------------------------------------------------------------------
set.seed(1234)
f_Lap<-function(x){
  y<-exp(-abs(x))/2
}

Metropolis<-function(Sigma, x0, n)
{
flag<-0
x<-numeric(n)
x[1]<-x0
u<-runif(n)
for(j in 2:n){
  y = rnorm(1, x[j-1], Sigma)
  if (u[j] <= (f_Lap(y) / f_Lap(x[j-1]))) x[j] = y 
  else {
  x[j] = x[j-1]
  flag<-flag +1
  }
}
  return(list(x=x,flag=flag))
}

Sigma = c(.05, .5, 2, 16)
RW1 = Metropolis(Sigma[1],25,1000)
RW2 = Metropolis(Sigma[2],25,1000)
RW3 = Metropolis(Sigma[3],25,1000)
RW4 = Metropolis(Sigma[4],25,1000)
print(c(RW1$flag, RW2$flag, RW3$flag, RW4$flag))
REJ = cbind(RW1$flag, RW2$flag, RW3$flag, RW4$flag)
ACC = round((1000-REJ)/1000,4)
ACC
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))  #画图
    RW = cbind(RW1$x, RW2$x, RW3$x,  RW4$x)
    for (i in 1:4) {
        plot(RW[,i], type="l",
             xlab=bquote(sigma == .(round(Sigma[i],4))),
             ylab="X", ylim=range(RW[,i]))
    }

## -----------------------------------------------------------------------------
# initialize constants and parameters
n <- 5000 # length of chain
Burn <- 1000 # burn-in length
X <- matrix(0, n, 2) # the chain, a bivariate sample
rhou <- 0.9 # correlation
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rhou^2)*sigma1
s2 <- sqrt(1-rhou^2)*sigma2

X[1, ] <- c(mu1, mu2) #initialize
for (i in 2:n) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rhou * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rhou * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
  }
b <- Burn + 1
x <- X[b:n, ]

# compare sample statistics to parameters
colMeans(x)
cov(x)
cor(x)
plot(x, main="", cex=.5, xlab=bquote(X[1]),ylab=bquote(X[2]), ylim=range(x[,2]))


## -----------------------------------------------------------------------------
# Fit a simple linear regression model
datalm<-as.data.frame(X)
plot(X[,1],X[,2],pch=20,xlab="x4",ylab="y4")
model <- lm(X[,1] ~ X[,2], data = datalm) #回归拟合
summary(model) 
anova(model)
abline(model, col = 2, lty = 2)

## -----------------------------------------------------------------------------
Gelman_Rubin <- function(psi) {
psi <- as.matrix(psi)
k <- nrow(psi)
n <- ncol(psi)
psi_means <- rowMeans(psi) #row means
B <- n * var(psi_means) #between variance est.
psi_w <- apply(psi, 1, "var") #within variances
W <- mean(psi_w) #within est.
v_hat <- W*(n-1)/n + (B/n) #upper variance est.
R_hat <- v_hat / W #G-R statistic
return(R_hat)
}

Normal.Chain <- function(sigma, N, X1) {
x <- rep(0, N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
  xt <- x[i-1]
  y <- rnorm(1, xt, sigma) 
  r1 <- dnorm(y, 0, 1) * dnorm(xt, y, sigma)
  r2 <- dnorm(xt, 0, 1) * dnorm(y, xt, sigma)
  r <- r1 / r2
  if (u[i] <= r) 
    x[i] <- y 
  else
    x[i] <- xt}
return(x)
}

sigma <- .2 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 15000 #length of chain
sb <- 1000 #burn-in length
#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- Normal.Chain(sigma, n, x0[i])
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman_Rubin(psi))
par(mfrow=c(2,2))
for (i in 1:k)
  par(mar=c(1,1,1,1))
  plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1))

#plot the sequence of R-hat statistics
r.hat <- rep(0, n)
for (j in (b+1):n)
  r.hat[j] <- Gelman_Rubin(psi[,1:j])
plot(r.hat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)

## -----------------------------------------------------------------------------
library(coda)
Normal.Chain <- function(sigma, N, X1) {
x <- rep(0, N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
  xt <- x[i-1]
  y <- rnorm(1, xt, sigma) #candidate point
  r1 <- dnorm(y, 0, 1) * dnorm(xt, y, sigma)
  r2 <- dnorm(xt, 0, 1) * dnorm(y, xt, sigma)
  r <- r1 / r2
  if (u[i] <= r) 
    x[i] <- y 
  else
    x[i] <- xt}
return(x)
}
x0 <- c(-5, 5, 10)
for (i in 1:3)
  X[i, ] <- Normal.Chain(0.2, 15000, x0[i])
List = mcmc.list(mcmc(X[1, ]), mcmc(X[2, ]),mcmc(X[3, ]))
gelman.diag(List) 
gelman.plot(List)

