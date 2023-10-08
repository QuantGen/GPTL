
# A wrapper to the C-function that performs Gradient Descent in a system of linear equations

GD<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/10,lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE){
    if(lambda>0){   diag(XX)=diag(XX)+lambda }

    if(lambda>0 & lambda0>0 ){    Xy=Xy+lambda*lambda0*b0 }

    learning_rate=learning_rate/mean(diag(XX))
    
    if(returnPath){
     B=matrix(nrow=ncol(XX),ncol=nIter+1,NA)
     B[,1]=b
     for(i in 2:ncol(B)){
        B[,i]=.Call("GRAD_DESC",XX, Xy, B[,i-1],p, 1, learning_rate)
     }
     rownames(B)=rownames(XX)
     colnames(B)=paste0('iter_',0:(nIter-1))
    }else{
         B <- .Call("GRAD_DESC",XX, Xy, b, p, nIter, learning_rate)
    }
    return(B)
}
