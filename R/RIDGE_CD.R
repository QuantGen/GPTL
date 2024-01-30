# Ridge regression

# C code implemented
RIDGE.CD<- function(XX, Xy, p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda=NULL,nIter=100,returnPath=TRUE) {
  
  B=array(dim=c(p,nIter,length(lambda)))
  bIni=b
  
  for (h in 1:length(lambda)) {
 
    B[,1,h]=bIni
    for (i in 2:nIter) {
      B[,i,h]=.Call("RIDGE_CD",XX, Xy, B[,i-1,h],p, 1, lambda[h], b0)
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
