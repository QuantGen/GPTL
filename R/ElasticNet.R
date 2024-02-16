# Elastic Net

# C code implemented
ElasticNet<- function(XX, Xy, p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda=NULL,nLambda=30,alpha=0.5,nIter=100,returnPath=TRUE) {
  #alpha=0 -> Ridge; alpha=1 -> Lasso
  
  if(is.null(lambda)){
    lambda.max=max(abs(Xy))+1e-4
    lambda.min=lambda.max/30
    lambda=seq(from=lambda.max,to=lambda.min,length=nLambda)
  }
  
  B=array(dim=c(p,nIter,length(lambda)))
  dimnames(B)=list(colnames(XX),paste0('iter_',1:nIter),paste0('lambda_',round(lambda,4)))
  
  lambda1=lambda*alpha
  lambda2=lambda*(1-alpha)*0.5
  
  for (h in 1:length(lambda)) {
    B[,1,h]=b
    for (i in 2:nIter) {
      B[,i,h]=.Call("ElasticNet",XX, Xy, B[,i-1,h], p, 1, lambda1[h], lambda2[h], b0)
    }
  }

  if (returnPath) {
    return(list(B=B[,,,drop=TRUE], lambda=lambda, alpha=alpha))
  } else {
    return(list(B=B[,nIter,,drop=TRUE], lambda=lambda, alpha=alpha))
  }
}

