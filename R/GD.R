
# A wrapper to the C-function that performs Gradient Descent in a system of linear equations

GD0<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/10,lambda=0,b0=rep(0,p),
              lambda0=1,returnPath=FALSE){

    
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
     colnames(B)=paste0('iter_',0:(nIter))
    }else{ 
        B=.Call("GRAD_DESC",XX, Xy, b+0,p, nIter, learning_rate)
    }
    #if(!returnPath){
    #  B=B[,ncol(B),drop=TRUE]
    #}
    return(B)
}

GD<- function( XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/10,
               lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE,sortGradient=TRUE){
    if(sortGradient){  
      tmp=order(abs(Xy),decreasing=TRUE)
      XX=XX[tmp,tmp]
      Xy=Xy[tmp]
    }
    previous_lambda=0
    B=array(dim=c(p,ifelse(returnPath,nIter,1),length(lambda)))

    for(h in 1:length(lambda)){
        diag(XX)=diag(XX)+(lambda[h]- previous_lambda)
        LR=learning_rate/mean(diag(XX))
  
        if( lambda0>0 ){    
            Xy=Xy+(lambda[h]-previous_lambda)*lambda0*b0 
        }
        previous_lambda=lambda[h]

        if(returnPath){
            B[,1,h]=b
            for(i in 2:ncol(B)){
              B[,i,h]=.Call("GRAD_DESC",XX, Xy, B[,i-1,h],p, 1, LR)
            }
         }else{ 
            B[,1,h]=.Call("GRAD_DESC",XX, Xy, b+0,p, nIter, LR)
        }
    }
    if(returnPath){
      iterations=paste0('iter_',1:nIter)
    }else{
      iterations=paste0('iter_',nIter)
    }
  
    dimnames(B)=list(colnames(XX),iterations,paste0('lambda_',lambda))
  
    B=B[,,,drop=TRUE]
    return(B)
}

# Old R implementation
GD.R<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/10,lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE){
    learning_rate=learning_rate/mean(diag(XX))
    B=matrix(nrow=ncol(XX),ncol=nIter+1,NA)
     B[,1]=b
     for(i in 2:ncol(B)){
        b=B[,i-1]
        for(j in 1:p){
            gradient=sum(XX[,j]*b)+lambda*b[j]-Xy[j]-lambda*lambda0*b0[j]
            b[j]=B[j,i-1]-learning_rate*gradient
        }
         B[,i]=b
     }
     rownames(B)=rownames(XX)
     colnames(B)=paste0('iter_',0:(nIter))
   if(!returnPath){
      B=as.vector(B[,ncol(B)])
   }
    return(B)
}
