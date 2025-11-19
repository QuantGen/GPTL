PR <- function(XX, Xy, b=NULL, lambda=NULL, nLambda=30, alpha=0, convThreshold=1e-4, maxIter=500, returnPath=FALSE, verbose=TRUE) {
  #alpha=0 -> Ridge; alpha=1 -> Lasso

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
  b0=b
  
  diagXX=as.vector(Matrix::diag(XX))
  
  if(is.null(lambda)){
    if (alpha == 0) {
      grid.max=10000
      grid.min=1
      grid=exp(seq(from=log(grid.max),to=log(grid.min),length=nLambda))
      K=mean(diagXX)
      lambda=K*grid
    } else {
      # Adjust for initial beta later...
      lambda.max=(max(abs(Xy))+1e-5)/alpha
      lambda.min=lambda.max/100
      lambda=exp(seq(from=log(lambda.max),to=log(lambda.min),length=nLambda))
    }
  }
  
  B=array(dim=c(p,maxIter,length(lambda)))
  dimnames(B)=list(colnames(XX),paste0('iter_',1:maxIter),paste0('lambda_',round(lambda,4)))
  
  lambda1=lambda*alpha
  lambda2=lambda*(1-alpha)*0.5

  conv_iter=rep(maxIter, length(lambda))
  
  for (h in 1:length(lambda)) {
    B[,1,h]=b
    
    if(is(XX,"dgCMatrix"))
    {
    	#Sparse matrix
    	for (i in 2:maxIter) {
      		B[,i,h]=.Call("ElasticNet_sparse",XX@x,XX@p,XX@i,diagXX,
      		              Xy, B[,i-1,h], p, 1, lambda1[h], lambda2[h], b0)
      		if (max(abs(B[,i,h]-B[,i-1,h])) < convThreshold) 
      		{
        		conv_iter[h]=i
        		break
      		}
    	}	
    }else{
    	#Dense matrix
    	for (i in 2:maxIter) {
      		B[,i,h]=.Call("ElasticNet",XX, Xy, B[,i-1,h], p, 1, lambda1[h], lambda2[h], b0)
      		if (max(abs(B[,i,h]-B[,i-1,h])) < convThreshold) 
      		{
        		conv_iter[h]=i
        		break
      		}
    	}	
    }
  }

  num_not_conv=sum(conv_iter == maxIter)

  if (num_not_conv > 0) {
    message(' There were ', num_not_conv, ' models reached maximum iterations (', maxIter, ') and not converged, check output for details.\n')
  }
  
  if (returnPath) {
    return(list(B=B[,,,drop=TRUE], lambda=lambda, alpha=alpha, conv_iter=conv_iter))
  } else {
    B_last=B[,conv_iter[1],1,drop=TRUE]
    if (length(lambda)>1) {
      for (h in 2:length(lambda)) {
        B_last=cbind(B_last, B[,conv_iter[h],h,drop=TRUE])
      }
      colnames(B_last)=dimnames(B)[[3]]
    }
    return(list(B=B_last, lambda=lambda, alpha=alpha, conv_iter=conv_iter))
  }
}

