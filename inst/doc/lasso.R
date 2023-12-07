## -----------------------------------------------------------------------------
Lasso <- function(x,y,xtrue){
  lasso_flag <-0
  lasso.mod <- glmnet(x, y, alpha = 1)
  cv.out <- cv.glmnet(x, y,group=FALSE, alpha = 1)
  fit_lasso <- glmnet(x,y,lambda = cv.out$lambda.min,family = "gaussian",standardize=TRUE)
  beta_lasso  <- as.vector(fit_lasso$beta)
  A <- which(abs(beta_lasso)>0)
  if(all( xtrue %in% A))
    lasso_flag <- 1
  return(list(lasso_flag=lasso_flag,beta_lasso=beta_lasso))
}

## -----------------------------------------------------------------------------
set.seed(1234)
library(glmnet)
library(ncvreg)
library(MASS)
lasso_flag1 <- 0
lasso_result1 <- 0
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
  #lasso回归
  lasso_flag1 <- lasso_flag1+Lasso(x,y,c(1,2,3))$lasso_flag
}
lasso_result1 <- lasso_flag1/k

lasso_result1

