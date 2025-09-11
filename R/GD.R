GD<- function(XX, Xy, b, nIter=10, learning_rate=1/50, lambda=0, lambda0=1, returnPath=FALSE, verbose=TRUE){

    if(!(is(XX,"matrix") | is(XX,"dgCMatrix"))) stop("XX must be a matrix or dgCMatrix\n")

    if(is.null(rownames(XX)) | is.null(colnames(XX))){
        stop('XX must have variant IDs as row and column names\n')
    }
    
    if (!all(rownames(XX) == colnames(XX))) stop("Row and column names in XX not match\n")

    if (is.vector(Xy)) {
        if (is.null(names(Xy))) {
            stop("Xy must have variant IDs as names\n")
        }
    } else if (is.matrix(Xy) | is.data.frame(Xy)) {
        if (is.null(rownames(Xy))) {
            stop("Xy must have variant IDs as row names\n")
        }
        XyName=rownames(Xy)
        Xy=as.vector(as.matrix(Xy))
        names(Xy)=XyName
    } else {
        stop("Xy must be in one of these formats: vector, matrix, data.frame\n")
    }

    if (is.vector(b)) {
        if (is.null(names(b))) {
            stop("The prior estimates vector (b) must have variant IDs as names\n")
        }
    } else if (is.matrix(b) | is.data.frame(b)) {
        if (is.null(rownames(b))) {
            stop("The prior estimates matrix (b) must have variant IDs as row names\n")
        }
        bName=rownames(b)
        b=as.vector(as.matrix(b))
        names(b)=bName
    } else {
        stop("b must be in one of these formats: vector, matrix, data.frame\n")
    }
    
    snp_list=Reduce(intersect, list(rownames(XX),names(Xy),names(b)))

    if (length(snp_list) == 0){ 
        stop("No matched variants in XX, Xy, and prior\n")
    }else{
        if(verbose){
            message(' There were ',length(snp_list), ' variants in common between XX, Xy, and the prior.\n')
        }
    }
    
    XX=XX[snp_list,snp_list,drop = FALSE]
    Xy=Xy[snp_list]
    b=b[snp_list]

    p=nrow(XX)
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

