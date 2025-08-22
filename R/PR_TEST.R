PR_TEST <- function(XX, Xy, b, lambda=NULL, nLambda=30, alpha=0, conv_threshold=1e-4, maxIter=500, returnPath=TRUE, verbose=TRUE) {
  #alpha=0 -> Ridge; alpha=1 -> Lasso

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
  b0=b
  
  diagXX=as.vector(Matrix::diag(XX))
  
  if(is.null(lambda)){
    if (alpha == 0) {
      grid.min=1
      grid.max=10000
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
      		if (max(abs(B[,i,h]-B[,i-1,h])) < conv_threshold) 
      		{
        		conv_iter[h]=i
        		break
      		}
    	}	
    }else{
    	#Dense matrix
    	for (i in 2:maxIter) {
      		B[,i,h]=.Call("ElasticNet",XX, Xy, B[,i-1,h], p, 1, lambda1[h], lambda2[h], b0)
      		if (max(abs(B[,i,h]-B[,i-1,h])) < conv_threshold) 
      		{
        		conv_iter[h]=i
        		break
      		}
    	}	
    }
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
