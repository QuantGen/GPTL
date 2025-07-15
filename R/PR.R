PR.SS<- function(XX, Xy, b, lambda=NULL, nLambda=30, alpha=0, conv_threshold=1e-4, maxIter=500, returnPath=TRUE) {
  #alpha=0 -> Ridge; alpha=1 -> Lasso

  if(!(is(XX,"matrix") | is(XX,"dgCMatrix"))) stop("XX must be a matrix or dgCMatrix\n")
  if (rownames(XX) != colnames(XX)) stop("Rowname and colname in XX not match\n")
  snp_list=Reduce(intersect, list(rownames(XX),rownames(Xy),rownames(b)))
  if (length(snp_list) == 0) stop("No matched SNPs in XX, Xy, and prior\n")
  XX=XX[snp_list,snp_list]
  Xy=Xy[snp_list,]
  b=b[snp_list,]

  p=ncol(XX)
  b0=rep(0,p)
  
  diagXX=as.vector(Matrix::diag(XX))
  
  if(is.null(lambda)){
    if (alpha == 0) {
      h2.min=0.01
      h2.max=0.9
      h2.grid=exp(seq(from=log(h2.min),to=log(h2.max),length=nLambda))
      K=mean(diagXX)
      lambda=K*(1-h2.grid)/h2.grid
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


# Stable PR function before implementing sparse sufficient statistics
PR1<- function(XX, Xy, p=ncol(XX), b=rep(0,p),b0=rep(0,p),lambda=NULL,nLambda=30,alpha=0,conv_threshold=1e-4,maxIter=500,returnPath=TRUE) {
  #alpha=0 -> Ridge; alpha=1 -> Lasso
  
  if(is.null(lambda)){
    if (alpha == 0) {
      h2.min=0.01
      h2.max=0.9
      h2.grid=exp(seq(from=log(h2.min),to=log(h2.max),length=nLambda))
      K=mean(diag(XX))
      lambda=K*(1-h2.grid)/h2.grid
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
    for (i in 2:maxIter) {
      B[,i,h]=.Call("ElasticNet",XX, Xy, B[,i-1,h], p, 1, lambda1[h], lambda2[h], b0)
      if (max(abs(B[,i,h]-B[,i-1,h])) < conv_threshold) {
        conv_iter[h]=i
        break
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
