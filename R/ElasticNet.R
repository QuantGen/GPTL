# Elastic Net

# C code implemented
ElasticNet<- function(XX, Xy, p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda=NULL,alpha=0.5,nIter=100,returnPath=TRUE) {
  
  B=array(dim=c(p,nIter,length(lambda)))
  bIni=b
  
  for (h in 1:length(lambda)) {
 
    B[,1,h]=bIni
    for (i in 2:nIter) {
      B[,i,h]=.Call("ElasticNet",XX, Xy, B[,i-1,h],p, 1, lambda[h], b0, alpha)
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

  XX=crossprod(X)
  Xy=crossprod(X,y)

  library(glmnet)
  fmL=glmnet(y=y,x=X, alpha=0.5)
  fm_EN=ElasticNet(XX, Xy, p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda=fmL$lambda,nIter=50,returnPath=T, alpha=0.5)


  par(mfrow=c(3,3))
  for(i in seq(from=2,to=18,by=2)){
    plot(fmL$beta[,i],fm_EN[,50,i],col=4,cex=1.5);abline(a=0,b=1)
  }

}
