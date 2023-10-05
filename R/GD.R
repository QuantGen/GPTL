
# A wrapper to the C-function that performs Gradient Descent in a system of linear equations

GD<- function(XX,Xy,b=rep(0,ncol(XX)),active=1:ncol(XX), RSS=1,nIter=10,learning_rate=1/5,lambda=0,b0=rep(0,ncol(XX)),lambda0=1,returnPath=FALSE){
    
    diag(XX)=diag(XX)+lambda
    learning_rate=learning_rate/mean(diag(XX))
    Xy=Xy+lambda*lambda0*b0
    active <- active - 1L # for the 0-based index
    if(returnPath){
     B=matrix(nrow=ncol(XX),ncol=nIter+1,NA)
     B[,1]=b
     for(i in 2:ncol(B)){
        B[,i]=.Call("GRAD_DESC",XX, Xy, B[,i-1], active, 1, learning_rate)[[1]]
     }
     rownames(B)=rownames(XX)
     colnames(B)=paste0('iter_',0:(nIter-1))
    }else{
         B <- .Call("GRAD_DESC",XX, Xy, b, active, nIter, learning_rate)[[1]]
    }
    return(B)
}
