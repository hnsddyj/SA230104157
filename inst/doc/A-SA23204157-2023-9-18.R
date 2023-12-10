## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
#使用逆变换方法从该分布生成大小为1000的随机样本
n <- 1000
u <-runif(n)
x <- abs(log(12/u))
x
#样本与目标分布进行比较
hist(x, prob = TRUE, main = expression(f(x)==12*exp(-abs(x))), xlab = "value")
y<-seq(0,10,0.01)
lines(y,12*exp(-abs(y)))

## -----------------------------------------------------------------------------
#编写函数
Myfunction <- function(m,n){
  N <- 1000
  j <- 0
  i <- 0
  y <- numeric(N)
  while (i < N)
  {u <- runif(1)
  j <- j + 1
  x <- runif(1)
  if (x^(m-1) * (1-x)^(n-1)> u)
  #选择上界c为1/B(a,b)，g(x)为U(0,1)均匀分布的密度函数
  {
    i <- i + 1
    y[i] <- x
  }
  j
  }
  return(y)
}
Myfunction(3,2)
#用理论Beta(3,2)密度叠加的样本直方图。
hist(Myfunction(3,2), freq = FALSE, breaks = seq(0,1,0.01), main = "histogram of the sample with the theoretical Beta(3,2) density", xlab = "value")
ff <- function(x) {
  x^2*(1-x)^{2-1}/beta(3,2)}
curve(ff, 0, 1,col = 2, add = TRUE)   
legend('topleft',"real density", col = 2, lwd = 1, cex=0.7)

## -----------------------------------------------------------------------------
UU1 <- runif(10000, min = -1, max = 1)    
UU2 <- runif(10000, min = -1, max = 1)
UU3 <- runif(10000, min = -1, max = 1)
UU <- ifelse((abs(UU3)>abs(UU2) & abs(UU3)>abs(UU1)), UU2, UU3)
hist(UU, freq = FALSE, breaks = seq(-1,1,0.01), main = "histogram density estimate of a large simulated random sample", xlab = "value")
ff <- function(X) {3/4*(1-X^2)}
curve(ff, -1, 1, col = 2, add = TRUE)   
legend('topleft',"real density", col = 2, lwd = 1, cex=0.6)   


## -----------------------------------------------------------------------------
#首先考虑默认情况，相当于随机的有放回抽样。要抽取的数组a所对应的编号相当于离散型的随机变量x，生成的均匀分布的随机数为y。因此可以用逆变换来生成对应的随机数编号。
my.sample1 <- function(a,b){
c <-length(a)
   x <- seq(1, c, 1)#数组a产生的编号
   p <- rep(1:1,c)/c#每个编号所对应的概率
   cp <- c(0,cumsum(p)) #累加得到分布函数
   y1 = numeric(b)
  for(j in 1:b)#一共进行b次，抽取b个数据
  { u=runif(1)#产生随机数
    for(i in 1:(length(cp)-1))
     { if(u>cp[i] & u<=cp[i+1])#利用逆变换的到x（即对应的编号）
        {y1[j]=x[i]
          break}
    }
   }
  return(a[y1])
}
#举例
a<-c(1,4,5,2,7,9,6,46,0,8)
my.sample1(a,3)
my.sample1(a,2)
#考虑，replace等于TRUE的情况，相当于随机的无放回抽样。
my.sample2 <- function(a,b){
  c <-length(a)
  x <- seq(1, c, 1)#数组a产生的编号
  p <- rep(1:1,c)/c#每个编号所对应的概率
  cp <- c(0,cumsum(p)) #累加得到分布函数
  y1 = numeric(b)
  j=1
  while(j<=b)#一共进行b次，抽取b个数据
  { u=runif(1)#产生随机数
    for(i in 1:(length(cp)-1))
    if(u>cp[i] & u<=cp[i+1])#利用逆变换的到x（即对应的编号）
    {xx=x[i]
      break
    }
    m=1
    flag=1#用于标记是否符合不重复的条件
    while(m<=j){
     if(y1[m]==xx)#判断是否与之前的有重复
      { flag=0
      break
      }
    m=m+1
    }
    if(flag==1)#不与之前的重复，则放入y1中
    { y1[j]=xx
       j=j+1
    }
  }
 return(a[y1])
}

#举例
b<-c(32,34,1,5,65,74,89)
my.sample2(b,4)
my.sample2(b,2)

