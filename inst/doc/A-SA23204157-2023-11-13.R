## -----------------------------------------------------------------------------
# 数值型数据的dataframe
numeric_df <- data.frame(
  A = c(1, 2, 3, 4),
  B = c(5, 6, 7, 8),
  C = c(9, 10, 11, 12)
)
# 用于计算单个向量的标准差的函数
compute_sd <- function(x) {
  sd(x)
}
# (a) Compute the standard deviation of every column in a numeric data frame
result_a <- vapply(numeric_df, compute_sd, numeric(1))
print(result_a)

# 混合数据的dataframe
mixed_df <- data.frame(
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 22),
  Score = c(95, 88, 75)
)

# (b) Compute the standard deviation of every numeric column in a mixed data frame
numeric_columns <- sapply(mixed_df, is.numeric)
result_b <- vapply(mixed_df[numeric_columns], compute_sd, numeric(1))
print(result_b)

## -----------------------------------------------------------------------------
set.seed(1234)
library(Rcpp)
library(microbenchmark)
#编写一个r函数
a <- 2
b <- 2
n <- 10
gxfunc<-function(m,N){
  GbN<-N
  gx<-matrix(0,ncol=GbN,nrow=2)
  gx[,1]<-m
  for(i in 2:GbN){
  x<-gx[,i-1][1]
  gx[,i][2]<-rbeta(1,shape1=x+a,shape2=n-x+b)
  y<-gx[,i][2]
  gx[,i][1]<-rbinom(1,size=n,prob=y)
}
  return(gx)
}
gxm <- c(1,0.1)
gx_r <- gxfunc(gxm,1000)

#编写一个rcpp函数
cppFunction(
'NumericMatrix mygibbs(double x,double y,int N) {
  NumericMatrix gx(2,N);
  int n=10,a=2,b=2;
  gx(0,0)=x;
  gx(1,0)=y;  
  for(int i=1;i<N;i++){
    double x0=gx(0,i-1);
    gx(1,i)=rbeta(1,x0+a,n-x0+b)[0];
    double y0=gx(1,i);
    gx(0,i)=rbinom(1,n,y0)[0];
  }
  return gx;
}')
gx_Rcpp<-mygibbs(1,0.1,1000)

#比较两个函数的计算时间
Time <- microbenchmark(Gibbs1=mygibbs(1,0.1,1000),Gibbs2=gxfunc(gxm,1000))
summary(Time)[,c(1,3,5,6)]


