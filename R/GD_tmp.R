GD_tmp<- function(XX, Xy, b=NULL, nIter=10, learningRate=1/50, lambda=0, returnPath=FALSE, verbose=TRUE){

    if(!(is(XX,"matrix") | is(XX,"dgCMatrix") | (nrow(XX)==ncol(XX)))) stop("XX must be a square matrix or dgCMatrix\n")
    if(!(is(Xy,"vector") | is(Xy,"matrix") | is(Xy,"data.frame"))) stop("Xy must be in one of these formats: vector, matrix or data.frame with single column\n")
    
    if(is.null(b)){
        b=rep(0, nrow(XX))
        names(b)=rownames(XX)
    }

    if(!(is(b,"vector") | is(b,"matrix") | is(b,"data.frame"))) stop("The prior estimates (b) must be in one of these formats: vector, matrix or data.frame with single column\n")

    Xy=as.matrix(Xy)
    b=as.matrix(b)

    nameWarningFlag=0
    if(is.null(rownames(XX)) | is.null(colnames(XX)) | is.null(rownames(Xy)) | is.null(rownames(b))){
        warning('Variant IDs are missing in one or more of inputs: XX, Xy, or b\n')
        nameWarningFlag=1
    }

    snp_list=Reduce(intersect, list(rownames(XX),rownames(Xy),rownames(b)))
    if (length(snp_list) == 0){ 
        warning('Variant IDs are not matching between inputs: XX, Xy, or b\n')
        nameWarningFlag=1
    }

    if ((nrow(XX)!=nrow(Xy) | nrow(XX)!=nrow(b)) & (nameWarningFlag==1)){
        stop('Distinct number of variants detected in inputs: XX, Xy, and b, while variant IDs are missing in one or more of inputs: XX, Xy, or b\n')
    }

    if (nameWarningFlag==0){
        XX=XX[snp_list,snp_list,drop = FALSE]
        Xy=Xy[snp_list,,drop = FALSE]
        b=b[snp_list,,drop = FALSE]
        if(verbose){
            message(length(snp_list), ' variants in common between XX, Xy, and b are retained\n')
        }
    } else {
        if(verbose){
            message(nrow(XX), ' variants are retained\n')
        }
    }

    p=nrow(XX)
    b0=rep(0,p)
  	
    previous_lambda=0
    B=array(dim=c(p,ifelse(returnPath,nIter+1,1),length(lambda)))
    RSS=numeric(nIter+1)
    RSS[1]=-2*t(b)%*%Xy+t(b)%*%XX%*%b

    for(h in 1:length(lambda))
    {
    	if(is(XX,"dgCMatrix"))
    	{	
    	    #Sparse matrix
        	Matrix::diag(XX)<-Matrix::diag(XX) + (lambda[h]- previous_lambda)
        	LR=learningRate/mean(Matrix::diag(XX))
        }else{
        	#Dense matrix
        	diag(XX)=diag(XX)+(lambda[h]- previous_lambda)
        	LR=learningRate/mean(diag(XX))
        }
  
        #if( lambda0>0 ){    
        #    Xy=Xy+(lambda[h]-previous_lambda)*lambda0*b0 
        #}
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
                    RSS[i]=-2*t(B[,i,h])%*%Xy+t(B[,i,h])%*%XX%*%B[,i,h]
            	}
            }
         }else{
         	if(is(XX,"dgCMatrix"))
         	{
         		#Sparse matrix
             	B[,1,h]=.Call("GRAD_DESC_sparse",XX@x,XX@p,XX@i,Xy, b+0,p, nIter+1, LR)
            }else{
            	#Dense matrix
            	B[,1,h]=.Call("GRAD_DESC",XX, Xy, b+0,p, nIter+1, LR)
            }
        }
    }
    
    
    if(returnPath){
      iterations=paste0('iter_',0:nIter)
    }else{
      iterations=paste0('iter_',nIter)
    }
    
    dimnames(B)=list(colnames(XX),iterations,paste0('lambda_',lambda))

    if(returnPath){
      B=B[,-1,,drop=TRUE]
      RSS=RSS[-1]
    }else{
      B=B[,,,drop=TRUE]
    }
    
    return(list(B=B,RSS=RSS))
    return(B)
}


fm=GD(XX, Xy, b=NULL, nIter=100, learningRate=1/50, lambda=0, returnPath=TRUE, verbose=TRUE)














