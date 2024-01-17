# Lasso regression

LASSO <- function(XX, Xy, p=ncol(XX), b=rep(0,p), lambda=1, b0=rep(0,p), lambda0=0, nIter=50, returnPath=FALSE) {
  B=array(dim=c(p,nIter,length(lambda)))
  for (h in 1:length(lambda)) {
    B[,1,h]=b
    for (i in 2:ncol(B)) {
      for (j in 1:p) {
        Q=(Xy[j]-XX[j,-j] %*% B[,i-1,h][-j])/XX[j,j]
        if (Q+lambda/XX[j,j] < lambda0*b0[j]) {
          B[,i,h][j]=Q+lambda/XX[j,j]
        } else if (Q-lambda/XX[j,j] > lambda0*b0[j]) {
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
