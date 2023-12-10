## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(boot)
set.seed(1234)
aircondit <- as.matrix(aircondit)
BOOT.mean <- function(x,j) mean(x[j])
BOOT.obj <- boot(aircondit, statistic=BOOT.mean, R=2000)
print(boot.ci(BOOT.obj, type = c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(1234)
#计算theta的估计值
datacov<-cov(scor)
Lambdahat<-eigen(datacov)$values
Thetahat<-Lambdahat[1]/sum(Lambdahat)
#利用刀切法
n<-nrow(scor)
Thetajack <- numeric(n)
for(i in 1:n){
  datanew<-cov(scor[-i,])
  Lambdanew<-eigen(datanew)$values
  Thetajack[i]<-Lambdanew[1]/sum(Lambdanew)
}

bias.jack <- (n-1)*(mean(Thetajack)-Thetahat)
se.jack <- sqrt((n-1)*mean((Thetajack-Thetahat)^2))
round(c(original=Thetahat,bias.jack=bias.jack,se.jack=se.jack),4)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
k <- length(magnetic)
res1 <- res2 <- res3 <- res4 <- numeric(k*(k-1)/2)
# for n-fold cross validation
# fit models on leave-two-out samples
a <- 0
for (i in 1:(k-1))
  for (j in (i+1):k) {
    a <- a+1
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    Q1 <- lm(y~x)
    y11 <- Q1$coef[1]+chemical[i]*Q1$coef[2]  
    y12 <- Q1$coef[1]+chemical[j]*Q1$coef[2] 
    res1[a] <- (magnetic[i]-y11)^2+(magnetic[j]-y12)^2
    
    Q2 <- lm(y~x+I(x^2))
    y21 <- Q2$coef[1] + Q2$coef[2] * chemical[i] + Q2$coef[3] * chemical[i]^2
    y22 <- Q2$coef[1] + Q2$coef[2] * chemical[j] + Q2$coef[3] * chemical[j]^2
    res2[a] <- (magnetic[i]-y21)^2+(magnetic[j]-y22)^2
    
    Q3 <- lm(log(y)~x)
    logy31 <- exp(Q3$coef[1] + Q3$coef[2] * chemical[i])
    logy32 <- exp(Q3$coef[1] + Q3$coef[2] * chemical[j])
    res3[a] <- (magnetic[i]-logy31)^2+(magnetic[j]-logy32)^2
    
    Q4 <- lm(log(y)~log(x))
    logy41 <- exp(Q4$coef[1] + Q4$coef[2] * log(chemical[i]))
    logy42 <- exp(Q4$coef[1] + Q4$coef[2] * log(chemical[j]))
    res4[a] <- (magnetic[i]-logy41)^2+(magnetic[j]-logy42)^2
  }

res = c(mean(res1)/2,mean(res2)/2,mean(res3)/2,mean(res4)/2)
res
#table<-cbind(prediction error,Quadratic,Exponential,Log-Log)
#knitr::kable(table,align="c")
matrix(res, nrow=1, dimnames=list("prediction error", c("Linear","Quadratic"," Exponential","Log-Log")))


