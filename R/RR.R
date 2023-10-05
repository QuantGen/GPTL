# A wrappper for ridge regression, possibly shrinking towards b0
# For now we just have ridge regression.. we will extend to LASSOand Elastic Net
 RR <- function(XX, Xy,lambda,p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda0=0, active=1:p, RSS=10000, maxIter=1000, tol=1e-5) {
     # Loss function is (y-Xb)'(y-Xb)+ lambda|b-lambda0*b0 |^2
    # adding the shrinkage parameter to the left-hand-side
    diag(XX)=diag(XX)+lambda

     # adding prior mean to the right-hand-side
    if(lambda0>0){
      Xy=Xy+lambda*lambda0*b0
    }
    
    ans<-fitLSYS(XX, Xy, b, active, RSS, maxIter, tol)[[1]] 
   
    return(ans)
 }
