#' @title About sure independent screening methods
#' @description SIS and ISIS using R

#' @export
#' @importFrom glmnet glmnet
#' @importFrom ncvreg ncvreg
#' @importFrom glmnet cv.glmnet
#' @importFrom ncvreg cv.ncvreg
#' @importFrom MASS mvrnorm
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom stats residuals
#' @importFrom stats coef
#' @importFrom Rcpp evalCpp
#' @param x The design matrix, of dimensions n * p, without an intercept.
#' @param y The response vector of dimension n * 1.
#' @param xtrue The size of true model
#' @return Returns an object with \item{ISIS_flag}{A sign of whether ISIS contains a true value.}
#' \item{AA}{Variables for ISIS screening.}
#' @examples
#' \dontrun{
#' a=Lasso(x,y,xtrue)
#' }

#' set.seed(1234)
#' p <- 100
#' n <- 20
#' rho <- 0.1
#' k <- 50
#' d <- 2.5*floor(n/log(n))
#' for(i in 1:k){
#' beta <- c(5,5,5,rep(0,(p-3)))
#' mean <- c(rep(0,p))
#' sigma <- matrix(c(rep(0.5,p*p)), nrow=p, ncol=p)+diag(p)*(1-rho)
#' @importFrom MASS mvrnorm
#' x <-mvrnorm(n,mean,sigma)
#' y <- x%*%beta+rnorm(n)
#' #SIS
#' SIS_flag1 <- SIS_flag1+SIS(x,y,d,c(1,2,3))$SIS_flag
#' #lasso回归
#' lasso_flag1 <- lasso_flag1+Lasso(x,y,c(1,2,3))$lasso_flag
#' #ISIS
#' ISIS_flag1 <- ISIS_flag1+ISIS(x,y,d,c(1,2,3))$ISIS_flag
#' }
#' SIS_result1 <- SIS_flag1/k
#' lasso_result1 <- lasso_flag1/k
# ISIS_result1 <- ISIS_flag1/k

Lasso <- function(x,y,xtrue){
  lasso_flag <-0
  lasso.mod <- glmnet(x, y, alpha = 1)
  cv.out <- cv.glmnet(x, y,grouped=FALSE, alpha = 1)
  fit_lasso <- glmnet(x,y,lambda = cv.out$lambda.min,family = "gaussian",standardize=TRUE)
  beta_lasso  <- as.vector(fit_lasso$beta)
  A <- which(abs(beta_lasso)>0)
  if(all( xtrue %in% A))
    lasso_flag <- 1
  return(list(lasso_flag=lasso_flag,beta_lasso=beta_lasso))
}

#' Title
#'
#' @param x x
#' @param y y
#' @param d d
#' @param xtrue xtrue
#'
#' @return listdata
#' @export
#'
#' @examples
#' \dontrun{
#' b=SIS(x,y,d,xtrue)
#' }


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

#' Title
#'
#' @param x x
#' @param y y
#' @param d d
#' @param xtrue  xtrue
#'
#' @return listdata
#' @export
#'
#' @examples
#' \dontrun{
#' a=SIS.SCAD(x,y,d,xtrue)
#' }


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
#' Title
#'
#' @param x x
#' @param y y
#' @param d d
#' @param xtrue  xtrue
#'
#' @return listdata
#' @export
#'
#' @examples
#' \dontrun{
#' a=ISIS(x,y,d,xtrue)
#' }


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

