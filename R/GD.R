GD.SS<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/50,
               lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE){
  
  	if(!(is(XX,"matrix") | is(XX,"dgCMatrix"))) stop("XX must be a matrix or dgCMatrix\n")
  	
    previous_lambda=0
    B=array(dim=c(p,ifelse(returnPath,nIter,1),length(lambda)))

    for(h in 1:length(lambda))
    {
    	if(is(XX,"dgCMatrix"))
    	{	
    	    #Sparse matrix
        	Matrix::diag(XX)<-Matrix::diag(XX) + (lambda[h]- previous_lambda)
        	LR=learning_rate/mean(Matrix::diag(XX))
        }else{
        	#Dense matrix
        	diag(XX)=diag(XX)+(lambda[h]- previous_lambda)
        	LR=learning_rate/mean(diag(XX))
        }
  
        if( lambda0>0 ){    
            Xy=Xy+(lambda[h]-previous_lambda)*lambda0*b0 
        }
        previous_lambda=lambda[h]

        if(returnPath)
        {
            B[,1,h]=b
            
            if(is(XX,"dgCMatrix"))
            {
            	#Sparse matrix
            	for(i in 2:ncol(B))
            	{
            		B[,i,h]=.Call("GRAD_DESC_sparse",XX@x,XX@p,XX@i,Xy, B[,i-1,h],p, 1, LR)
            	}
            }else{
            	#Dense matrix
            	for(i in 2:ncol(B))
            	{
            		B[,i,h]=.Call("GRAD_DESC",XX, Xy, B[,i-1,h],p, 1, LR)
            	}
            }            
         }else{
         	if(is(XX,"dgCMatrix"))
         	{
         		#Sparse matrix
             	B[,1,h]=.Call("GRAD_DESC_sparse",XX@x,XX@p,XX@i,Xy, b+0,p, nIter, LR)
            }else{
            	#Dense matrix
            	B[,1,h]=.Call("GRAD_DESC",XX, Xy, b+0,p, nIter, LR)
            }
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
GD.R<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/50,lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE){
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


# A wrapper to the C-function that performs Gradient Descent in a system of linear equations
GD0<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/50,lambda=0,b0=rep(0,p),
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

# Stable GD function before implementing sparse sufficient statistics
GD1<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/50,
               lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE){
    # disabled for now because orderBack may not be working properly
    if(FALSE){ 
      gradient0=Xy
      if(any(abs(b)>.Machine$double.eps)){ gradient0=gradient0-XX%*%b }
      tmp=order(abs(gradient0),decreasing=TRUE)
      XX=XX[tmp,tmp]
      Xy=Xy[tmp]
      orderBack=order(tmp)
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

    # sortGradient disabled for now because orderBack may not be working properly
    if(FALSE){ 
      B=B[orderBack,,,drop=FALSE] 
    }
  
    B=B[,,,drop=TRUE]

  
    return(B)
}


# First sparse version
GD_sparse<- function(XX,Xy,p=ncol(XX),b=rep(0,p), nIter=10,learning_rate=1/50,
               lambda=0,b0=rep(0,p),lambda0=1,returnPath=FALSE){
    # disabled for now because orderBack may not be working properly
    if(FALSE){ 
      gradient0=Xy
      if(any(abs(b)>.Machine$double.eps)){ gradient0=gradient0-XX%*%b }
      tmp=order(abs(gradient0),decreasing=TRUE)
      XX=XX[tmp,tmp]
      Xy=Xy[tmp]
      orderBack=order(tmp)
    }
  
    previous_lambda=0
    B=array(dim=c(p,ifelse(returnPath,nIter,1),length(lambda)))

    for(h in 1:length(lambda))
    {
        diag(XX)=diag(XX)+(lambda[h]- previous_lambda)
        LR=learning_rate/mean(diag(XX))
  
        if( lambda0>0 ){    
            Xy=Xy+(lambda[h]-previous_lambda)*lambda0*b0 
        }
        previous_lambda=lambda[h]

        if(returnPath)
        {
        	#Inefficient, but it is just for testing purposes
            C<-as(XX,"dgCMatrix")
            B[,1,h]=b
            for(i in 2:ncol(B))
            {
              B[,i,h]=.Call("GRAD_DESC_sparse",C@x,C@p,C@i,Xy, B[,i-1,h],p, 1, LR)
              #B[,i,h]=.Call("GRAD_DESC",XX, Xy, B[,i-1,h],p, 1, LR)
            }
         }else{
         	 #Inefficient, but it is just for testing purposes
             C<-as(XX,"dgCMatrix")
             B[,1,h]=.Call("GRAD_DESC_sparse",C@x,C@p,C@i,Xy, b+0,p, nIter, LR)
            #B[,1,h]=.Call("GRAD_DESC",XX, Xy, b+0,p, nIter, LR)
        }
    }
    if(returnPath){
      iterations=paste0('iter_',1:nIter)
    }else{
      iterations=paste0('iter_',nIter)
    }
  
    dimnames(B)=list(colnames(XX),iterations,paste0('lambda_',lambda))

    # sortGradient disabled for now because orderBack may not be working properly
    if(FALSE){ 
      B=B[orderBack,,,drop=FALSE] 
    }
  
    B=B[,,,drop=TRUE]

  
    return(B)
}
