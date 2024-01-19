# Lasso regression

LASSO <- function(XX, Xy, p=ncol(XX), b=rep(0,p), lambda=1, b0=rep(0,p), lambda0=0, nIter=50, returnPath=FALSE) {
  B=array(dim=c(p,nIter,length(lambda)))
  for (h in 1:length(lambda)) {
    B[,1,h]=b
    for (i in 2:ncol(B)) {
      for (j in 1:p) {
        Q=(Xy[j]-XX[j,-j] %*% B[,i-1,h][-j])/XX[j,j]
        if (Q+lambda/XX[j,j] > lambda0*b0[j]) {
          B[,i,h][j]=Q+lambda/XX[j,j]
        } else if (Q-lambda/XX[j,j] < lambda0*b0[j]) {
          B[,i,h][j]=Q-lambda/XX[j,j]
        } else {
          B[,i,h][j]=lambda0*b0[j]
        }
      }
    }
  }
  iterations=paste0('iter_',1:nIter)
  dimnames(B)=list(colnames(XX),iterations,paste0('lambda_',lambda))
  if (returnPath) {
    return(B)
  } else {
    return(B[,nIter,])
  }
}

OLS<- function(XX, Xy, p=ncol(XX), b=rep(0,p),lambda=0,nIter=100,returnPath=TRUE) {
  B=array(dim=c(p,nIter,length(lambda)))
  bIni=b
  
  for (h in 1:length(lambda)) {
    B[,1,h]=bIni
    for (i in 2:nIter) {
       for (j in 1:p) {
        offset=sum(XX[,j]*b)-XX[j,j]*b[j]
        bOLS=(Xy[j]-offset)/XX[j,j]
        b[j]=bOLS
       }
       B[,i,h]=b
      }
    }
    return(B[,,,drop=TRUE])
}

## Testing OLS (it Works!)
if(FALSE){
  n=1000
  p=10
  X=matrix(nrow=n,ncol=p,rnorm(n*p))
  b=rexp(p)
  signal=X%*%b
  error=rnorm(sd=sd(signal),n=n)
  y=signal+error

  fm0=lm(y~X-1)

  fm=OLS(XX=crossprod(X),Xy=crossprod(X,y),nIter=1000)
  plot(coef(fm0),fm[,1000])
}

## Now LASSO (using Coordinate Descent)

LASSO.CD<- function(XX, Xy, p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda=NULL,nIter=100,returnPath=TRUE) {

  ## Notes
  #  To match glment, XX and Xy should be crossprod(X)/n and crossprod(X,y)/n
  #  With both X and y centered
  ##

  ## Default lambda for LASSO
  if(is.null(lambda)){
    bOLS=Xy/diag(XX)
    lambda.max=max(abs(bOLS-b0))+1e-4
    lambda.min=lambda.max/100
    lambda=exp(seq(from=log(lambda.max),to=log(lambda.min),length=10))# glmnet uses length=100
  }

  
  B=array(dim=c(p,nIter,length(lambda)))
  bIni=b

  D=diag(XX)
  
  for (h in 1:length(lambda)) {
 
    B[,1,h]=bIni
    for (i in 2:nIter) {
       for (j in 1:p) {
        offset=sum(XX[,j]*b)-XX[j,j]*b[j]
        bOLS=(Xy[j]-offset)/XX[j,j]
        if(abs(bOLS)>lambda[h]/D[j]){
          b[j]=bOLS-sign(bOLS)*lambda[h]/D[j]
        }else{
          b[j]=0
        }
       }
       B[,i,h]=b
      }
    }
    return(B[,,,drop=TRUE])
}


## Testing LASSO.CD (It seems to work) with scaled predictors
if(FALSE){
  n=1000
  p=10
  QTL=c(2,4,6,8)

  ## With this scaling the diagonals of X'X are all equalt to n
  X=scale(matrix(nrow=n,ncol=p,rnorm(n*p)),center=TRUE,scale=TRUE)*sqrt((n)/(n-1))
 
  b=rep(0,p)
  b[QTL]=1
  signal=X%*%b
  error=rnorm(sd=sd(signal),n=n)
  y=signal+error

  fm0=lm(y~X-1)
  library(glmnet)
  fmL=glmnet(y=y,x=X,standarize=FALSE)

  ## Note: to match we need to multiply lambda by n...
  fm=LASSO.CD(XX=crossprod(X),Xy=crossprod(X,y),nIter=5000,lambda=fmL$lambda*n)
  par(mfrow=c(3,3))
  for(i in seq(from=2,to=18,by=2)){
   plot(fmL$beta[,i],fm[,1000,i],col=4,cex=1.5);abline(a=0,b=1)
  }
  
}


## Now LASSO with initial beta_0 (using Gauss-Seidel type algorithm)
LASSO.GS2<- function(XX, Xy, p=ncol(XX), b=rep(0,p), lambda=0, b0=rep(0,p), lambda0=0, nIter=100,returnPath=TRUE) {
  B=array(dim=c(p,nIter,length(lambda)))
  bIni=b
  
  for (h in 1:length(lambda)) {
    B[,1,h]=bIni
    for (i in 2:nIter) {
      for (j in 1:p) {
        offset=sum(XX[,j]*b)-XX[j,j]*b[j]
        bOLS=(Xy[j]-offset)/XX[j,j]
        if(abs(bOLS-b0[j])>lambda[h]){
          b[j]=bOLS-sign(bOLS)*lambda[h]
        }else{
          b[j]=b0[j]
        }
      }
      B[,i,h]=b
    }
  }
  iterations=paste0('iter_',1:nIter)
  dimnames(B)=list(colnames(XX),iterations,paste0('lambda_',lambda))
  if (returnPath) {
    return(B[,,,drop=TRUE])
  } else {
    return(B[,nIter,,drop=TRUE])
  }
}


if(FALSE){
  n=1000
  p=10
  QTL=c(2,4,6,8)
  X=matrix(nrow=n,ncol=p,rnorm(n*p))
  b=rep(0,p)
  b[QTL]=1
  signal=X%*%b
  error=rnorm(sd=sd(signal),n=n)
  y=signal+error

  fm0=lm(y~X-1)
  library(glmnet)
  fmL=glmnet(y=y,x=X)
  fm=LASSO.GS(XX=crossprod(X),Xy=crossprod(X,y),nIter=1000,lambda=fmL$lambda)
  fm2=LASSO.GS2(XX=crossprod(X),Xy=crossprod(X,y),nIter=1000,lambda=fmL$lambda,b0=b)

  par(mfrow=c(3,3))
  for(i in seq(from=2,to=18,by=2)){
   plot(fmL$beta[,i],fm[,1000,i],col=4,cex=1.5);abline(a=0,b=1)
  }
}

