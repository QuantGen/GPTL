GD<- function(XX, Xy, b, nIter=10, learning_rate=1/50, lambda=0, lambda0=1, returnPath=FALSE, verbose=TRUE){

    if(!(is(XX,"matrix") | is(XX,"dgCMatrix"))) stop("XX must be a matrix or dgCMatrix\n")

    if(is.null(rownames(XX)) | is.null(colnames(XX))){
        stop('XX must have variant IDs as row and column names\n')
    }
    
    if (!all(rownames(XX) == colnames(XX))) stop("Row and column names in XX not match\n")

    if(is.null(names(Xy))){
        stop('Xy must have variant IDs as names\n')
    }

    if (is.vector(b)) {
        if (is.null(names(b))) {
            stop("The prior estimates vector (b) must have variant IDs as names\n")
        }
        b=as.data.frame(b)
    } else if (is.matrix(b) | is.data.frame(b)) {
        if (is.null(rownames(b))) {
            stop("The prior estimates matrix (b) must have variant IDs as row names\n")
        }
        b=as.data.frame(b)
    } else {
        stop("b must be in one of these formats: vector, matrix, data.frame\n")
    }
    
    snp_list=Reduce(intersect, list(rownames(XX),names(Xy),rownames(b)))

    if (length(snp_list) == 0){ 
        stop("No matched variants in XX, Xy, and prior\n")
    }else{
        if(verbose){
            message(' There were ',length(snp_list), ' variants in common between XX, Xy, and the prior.\n')
        }
    }
    
    XX=XX[snp_list,snp_list,drop = FALSE]
    Xy=Xy[snp_list,,drop = FALSE]
    b=b[snp_list,,drop = FALSE]

    p=length(Xy)
    b0=rep(0,p)
  	
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


GD.ld<- function(ld, gwas, b, nIter=10, learning_rate=1/50, lambda=0, lambda0=1, returnPath=FALSE){
    
    if (!all(rownames(ld) == colnames(ld))) stop("Rowname and colname in LD not match\n")

    if (!all(c('id', 'beta', 'se', 'n', 'allele_freq') %in% colnames(gwas))) stop("Must provide GWAS results that consist of columns: id (variant IDs), beta (variant effects), se (variant standard errors), n (sample sizes for GWAS), allele_freq (variant allele frequency)\n")
    
    snp_list=Reduce(intersect, list(rownames(ld),gwas$id,names(b)))
    if (length(snp_list) == 0) stop("No matched SNPs in LD, GWAS, and prior\n")
    ld=ld[snp_list,snp_list,drop = FALSE]
    gwas=gwas[gwas$id %in% snp_list,,drop = FALSE]
    b=b[snp_list]

    p=nrow(gwas)
    b0=rep(0,p)
  
    allele_freq=gwas$allele_freq
    beta=gwas$beta
    n_gwas=gwas$n
    sd=sqrt(2 * allele_freq * (1-allele_freq))
    sd_diag=Diagonal(x=sd)
    XX=(n_gwas-1) * (sd_diag %*% ld %*% sd_diag)
    Xy=beta * Matrix::diag(XX)
  	
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
