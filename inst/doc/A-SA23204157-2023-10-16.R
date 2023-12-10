## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
M <- 1000
p1<-runif(950)
p2<-rbeta(50,0.1,1)
p3<-c(p1,p2)
p<-sort(p3)
p.adjust1<-p.adjust(p,method="bonferroni",n=length(p))
p.adjust2<-p.adjust(p,method="fdr",n=length(p))
#计算bonferroni检验的FWER,FDR(V/R),TPR(S/M-M0)
R1=0
m0=0
v1=0
for(i in 1:length(p)){
  if(p.adjust1[i]<0.1)
    R1 <- R1+1
  if(p[i]<0.1)
    m0 <- m0+1
  if(p.adjust1[i]<0.1&&p[i]<0.1)
    v1 <- v1+1
}
S1=R1-v1
FWER1=v1/M
FDR1=v1/R1
TPR1=S1/(M-m0)
#计算BH检验的FWER,FDR(V/R),TPR(S/M-M0)
R2=0
m02=0
v2=0
for(i in 1:length(p)){
  if(p.adjust2[i]<0.1)
    R2 <- R2+1
  if(p[i]<0.1)
    m02 <- m02+1
  if(p.adjust2[i]<0.1&&p[i]<0.1)
    v2 <- v2+1
}
S2=R2-v2
FWER2=v2/M
FDR2=v2/R2
TPR2=S2/(M-m02)

## -----------------------------------------------------------------------------
set.seed(1234)
n <- c(5,10,20)
B <- 1000
m <- 1000
j<-1
#模拟一千次
while(j<=m){
for(i in 1:length(n)){
  x <- rexp(n[i],2)
  lambda<-mean(x)
  lambdastar <- numeric(B)
  for(b in 1:B)
  {
   xstar<-sample(x,replace=TRUE)
   lambdastar[b]<-mean(xstar)
  }
y<-c(bias=mean(lambdastar)-lambda,se.boot=sd(lambdastar),se.samp=sd(x)/sqrt(length(x)))
j<-j+1
}
}

## -----------------------------------------------------------------------------
set.seed(1234)
n <- c(5,10,20)
bias<-numeric(3)
standarderror<-numeric(3)
lambda<-2
for (i in 1:length(n)){
  bias[i]<- 2/(n[i]-1)
  standarderror[i]<- 2*n[i]/((n[i]-1)*sqrt(n[i]-2))
}
bias
standarderror

## -----------------------------------------------------------------------------
BOOT.t.ci <- function(x, B = 500, R = 100, level = .95, statistic){
#计算 bootstrap的t置信区间
x <- as.matrix(x)
n <- nrow(x)
se <- numeric(B)
STAT <- numeric(B)
#用局部函数去计算bootstrap
BOOT.se <- function(x, R, f) {
m <- nrow(x)
x <- as.matrix(x)
TH <- replicate(R, expr = {
  i <- sample(1:m, size = m, replace = TRUE)
  f(x[i, ])
})
return(sd(TH))
}
for (b in 1:B) {
j <- sample(1:n, size = n, replace = TRUE)
y <- x[j, ]
STAT[b] <- statistic(y)
stat1 <- statistic(x)
se[b] <- BOOT.se(y, R = R, f = statistic)
}
T.stats <- (STAT - stat1) / se
Alpha <- 1 - level
SE0 <- sd(STAT)
QT <- quantile(T.stats, c(Alpha/2, 1-Alpha/2), type = 1)
names(QT) <- rev(names(QT))
CI <- rev(stat1 - QT * SE0)
}

## -----------------------------------------------------------------------------
LSAT<-c(576,635,558,578,666,580,555,661,651,605,653,575,545,572,594)
GPA<-c(339,330,281,303,344,307,300,343,336,313,312,274,276,288,296)
data <- cbind(LSAT, GPA)
stat <- function(data) {
cor(data[,1],data[,2])}
CI <- BOOT.t.ci(data, statistic = stat, B=2000, R=200)
print(CI)

