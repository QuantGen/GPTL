# Lasso regression
# ! Still Developing

LASSO <- function(XX, Xy, p=ncol(XX), b=rep(0,p), lambda=1, b0=rep(0,p), lambda0=1, nIter=50, returnPath=FALSE) {
  B=array(dim=c(p,nIter,length(lambda)))
  for (h in 1:length(lambda)) {
    B[,1,h]=b
    for (i in 2:ncol(B)) {
      for (j in 1:p) {
        B2j=(Xy[j]-XX[j,-j] %*% B[,i-1,h][-j])/XX[j,j]
        if (abs(B2j-lambda0*b0[j]) < lambda) {
          B[,i,h][j]=lambda0*b0[j]
        } else {
          B[,i,h][j]=B2j-sign(B2j-lambda0*b0[j]-(lambda/XX[j,j]))*(lambda/XX[j,j])
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
