## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
#绘制四个图像
set.seed(1234)
x<-seq(0,10,0.01)
g<-x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
f1<-x*exp(-x/2)/4*(x>0)
f2<-4*(x^2)*exp(-2*x)*(x>0)
f3<-x*exp((-x)/2)/sqrt(2*pi)
f4<-x*exp((-x^2+1)/2)/sqrt(2*pi)
gs<-c(expression(g(x)==x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)),
            expression(f[1](x)==x*exp(-x/2)/4*(x>0)),
            expression(f[2](x)==4(x^2)*exp(-2*x)*(x>0)),
            expression(f[3](x)==x*exp((-x)/2)/sqrt(2*pi)),
            expression(f[4](x)==x*exp((-x^2+1)/2)/sqrt(2*pi)))

plot(x, g,type = "l",main="g(x)和fi(x)的图像",ylab='',ylim=c(0,0.6),lwd = 1.5,col=1)
lines(x,f1,lty = 2, lwd = 2,col=2)
lines(x,f2,lty = 2, lwd = 2,col=3)
lines(x,f3,lty = 2, lwd = 2,col=4)
lines(x,f4,lty = 2, lwd = 2,col=5)
legend("topright", legend = gs, lty = 1:5, lwd = 2, inset = 0.02,col=1:5)

## -----------------------------------------------------------------------------
#绘制g(x)/fi(x)的图像
plot(x, g/f1, type = "l",main="g(x)/fi(x)的图像", ylab = "",
        ylim = c(0,15), lwd = 2, lty = 2,col=2)
    lines(x, g/f2, lty = 2, lwd = 2,col=3)
    lines(x, g/f3, lty = 2, lwd = 2,col=4)
    lines(x, g/f4, lty = 2, lwd = 2,col=5)
    legend("topright", legend = gs,
           lty = 2, lwd = 2, inset = 0.02,col=1:5)

## -----------------------------------------------------------------------------
set.seed(1234)
n<-10000
h1<-function(x){
  g<-x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
  f1<-x*exp(-x/2)/4*(x>0)
  h<-g/f1
  return(h)
}
h2<-function(x){
  g<-x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
  f2<-4*(x^2)*exp(-2*x)*(x>0)
  h<-g/f2
  return(h)
}
#生成服从f1的随机数
v<-rgamma(n,2,1/2)
y1<-h1(v)
mean1<-mean(y1)
sd1<-sd(y1)
mean1
sd1
#生成服从f2的随机数
u<-rgamma(n,3,2)
y2<-h2(u)
mean2<-mean(y2)
sd2<-sd(y2)
mean2
sd2


## -----------------------------------------------------------------------------
M <- 10000
k <- 10#k为分割的大小
r <- M/k 
N <- 50 
T2 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(x)exp(-1/x^2)/(x^4)*sqrt(2*pi)*(x>0)*(x<1)
for (i in 1:N) {
  est[i, 1] <- mean(g(runif(M)))
  for(j in 1:k)T2[j]<-mean(g(runif(M/k,(j-1)/k,j/k)))
  est[i, 2] <- mean(T2)}
  round(apply(est,2,mean),4)
  round(apply(est,2,sd),5)

## -----------------------------------------------------------------------------
M <- 10000
N <- 50 
k <- 5#k为分割的大小 
r <- M/k 
T1 <- numeric(k)
est <- matrix(0, N, 2)
g<-function(x,a,b) exp(-x)/(1+x^2)*(x>a)*(x<b)
gf<-function(x,a,b) g(x,a,b)/(exp(-x)/(exp(-a)-exp(-b)))
l<-function(u,a,b) -log(exp(-a)-u*(exp(-a)-exp(-b)))
for (i in 1:N) {
  u.imp<-runif(M)
  u.str<-runif(M/k)
  est[i, 1] <- mean(gf(l(u.imp,0,1),0,1))
  for(j in 1:k) 
    T1[j]<- mean(gf(l(u.str,(j-1)/k,j/k),(j-1)/k,j/k))
  est[i, 2] <- sum(T1)
}
round(apply(est,2,mean),5)
round(apply(est,2,sd),5)

## -----------------------------------------------------------------------------
set.seed(1234)
n<-20
m<-1000
x<-numeric(m)
mu<-numeric(m)
sigma<-numeric(m)
flag=0
for(i in 1:m){
  x<-rchisq(n,df=2)
  mu[i]<-mean(x)
  sigma[i]<-sd(x)
  if(mu[i]-2.093*sigma[i]/sqrt(20)<=2&mu[i]+2.093*sigma[i]/sqrt(20)>=2){
  flag<-flag+1
  }
  cp<-flag/m
}
  cp

## -----------------------------------------------------------------------------
set.seed(1234)
n<-20
m<-1000
x<-numeric(m)
mu<-numeric(m)
sigma<-numeric(m)
flag=0
for(i in 1:m){
  x<-rlnorm(n,0,1)
  mu[i]<-mean(x)
  sigma[i]<-sd(x)
  if(mu[i]-2.093*sigma[i]/sqrt(20)<=exp(1/2)&mu[i]+2.093*sigma[i]/sqrt(20)>=exp(1/2)){
  flag<-flag+1
  }
  cp<-flag/m
}
  cp

## -----------------------------------------------------------------------------
set.seed(1234)
m<-1000
n<-20
ucl<-numeric(m)
y<-numeric(m)
for(i in 1:m){
  x<-rchisq(n,df=2)
  ucl[i]<-(n-1)*var(x)/qchisq(0.05,df=n-1)
  y[i]<-(ucl[i]>=4)
}
cp<-mean(y)
cp

## -----------------------------------------------------------------------------
n<-c(10,100,1000)#样本容量
m<-1000#试验次数
set.seed(1234)
T0<-numeric(m)
Tie1<-Tie2<-Tie3<-numeric(3)
p.val1<-p.val2<-p.val3<-numeric(m)
for(i in 1:length(n)){
  flag1<-flag2<-flag3<-0
  #1.自由度为1的卡方分布
  for(j in 1:m){
    x<-rchisq(n[i],1)
    T0[j]<-(mean(x)-1)*sqrt(n[i])/sd(x)#计算$T^o$
    p.val1[j]<-2*(1-pt(abs(T0[j]),n[i]-1))
   if(p.val1[j]<=0.05)
     flag1<-flag1+1
  }
  Tie1[i]<-flag1/m
  #2.均匀分布Uniform(0,2)
  for(j in 1:m){
    x<-runif(n[i],0,2)
    T0[j]<-(mean(x)-1)*sqrt(n[i])/sd(x)#计算$T^o$
    p.val2[j]<-2*(1-pt(abs(T0[j]),n[i]-1))
    if(p.val2[j]<=0.05)
     flag2<-flag2+1
  }
  Tie2[i]<-flag2/m
  #3.指数分布Exponential(rate=1)
  for(j in 1:m){
    x<-rexp(n[i],1)
    T0[j]<-(mean(x)-1)*sqrt(n[i])/sd(x)#计算$T^o$
    p.val3[j]<-2*(1-pt(abs(T0[j]),n[i]-1))
    if(p.val3[j]<=0.05)
     flag3<-flag3+1
  }
  Tie3[i]<-flag3/m
}

Tie1
Tie2
Tie3

