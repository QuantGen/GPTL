BMM_Block=function(XX, Xy, B, my, vy, n, nIter=1200, burnIn=200, thin=5, R2=0.25,
	        nComp=matrix(ncol(B)), K=1/nComp, df0.E=5, S0.E=vy*(1-R2)*df0.E, df0.b=rep(10,nComp), 
	        priorProb=rep(1/nComp,nComp), priorCounts=rep(2*nComp,nComp), 
			fixVarE=FALSE, fixVarB=rep(FALSE,ncol(B)), verbose=TRUE, ...) {
	
	if(!(is(XX,"matrix") | is(XX,"dgCMatrix") | (nrow(XX)==ncol(XX)))) stop("XX must be a square matrix or dgCMatrix\n")
    if(!(is(Xy,"vector") | is(Xy,"matrix") | is(Xy,"data.frame"))) stop("Xy must be in one of these formats: vector, matrix or data.frame with single column\n")
    
    if(is.null(B)){
        B=rep(0, nrow(XX))
        names(B)=rownames(XX)
    }

    if(!(is(B,"vector") | is(B,"matrix") | is(B,"data.frame"))) stop("The prior estimates (B) must be in one of these formats: vector, matrix or data.frame with single column\n")

    Xy=as.matrix(Xy)
    B=as.matrix(B)

    nameWarningFlag=0
    if(is.null(rownames(XX)) | is.null(colnames(XX)) | is.null(rownames(Xy)) | is.null(rownames(B))){
        warning('Variant IDs are missing in one or more of inputs: XX, Xy, or B\n')
        nameWarningFlag=1
    }

    snp_list=Reduce(intersect, list(rownames(XX),rownames(Xy),rownames(B)))
    if (length(snp_list) == 0){ 
        warning('Variant IDs are not matching between inputs: XX, Xy, or B\n')
        nameWarningFlag=1
    }

    if ((nrow(XX)!=nrow(Xy) | nrow(XX)!=nrow(b)) & (nameWarningFlag==1)){
        stop('Distinct number of variants detected in inputs: XX, Xy, and B, while variant IDs are missing in one or more of inputs: XX, Xy, or B\n')
    }

    if (nameWarningFlag==0){
        XX=XX[snp_list,snp_list,drop = FALSE]
        Xy=Xy[snp_list,,drop = FALSE]
        B=B[snp_list,,drop = FALSE]
        if(verbose){
            message(length(snp_list), ' variants in common between XX, Xy, and B are retained\n')
        }
    } else {
        if(verbose){
            message(nrow(XX), ' variants are retained\n')
        }
    }

	diagXX=as.vector(Matrix::diag(XX))
	
	ld_index=get_block_ids(XX)
	nBlocks=max(ld_index)

	b=numeric(nrow(XX))
	names(b)=rownames(XX)

	POST.PROB=matrix(NA, nrow=nrow(XX), ncol=ncol(B))

	MeanVarB=matrix(NA, nrow=nBlocks, ncol=ncol(B))
	MeanVarE=numeric(nBlocks)

	# loop over blocks
	for (blk in 1:nBlocks) {
		timeIn=proc.time()[3]

		snps=rownames(XX)[ld_index == blk]

		XXblk=XX[snps,snps]
		Bblk=B[snps,]
		Xyblk=Xy[snps,]

		diagXXblk=as.vector(Matrix::diag(XXblk))
		R2blk=R2*sum(diagXXblk)/sum(diagXX)

		if (length(snps)==1) {
			XXblk=as.matrix(XXblk)
			colnames(XXblk)=snps
            rownames(XXblk)=snps
			Bblk=matrix(Bblk, nrow=1)
			rownames(Bblk)=snps
		}

		fm=BMM(XX=XXblk, Xy=Xyblk, B=Bblk, n=n, my=my, vy=vy, nIter=nIter, burnIn=burnIn, thin=5, R2=R2blk, fixVarE=fixVarE, fixVarB=fixVarB, verbose=FALSE, ...)
		
		b[snps]=fm$b
		POST.PROB[ld_index == blk,]=fm$POST.PROB
		MeanVarB[blk,]=fm$postMeanVarB
		MeanVarE[blk]=mean(fm$samplesVarE)

		timeOut=proc.time()[3]

		if (verbose) {
			message(' Block ', blk, '/', nBlocks, ' varE=', round(MeanVarE[blk],4), '; varB=[', paste(round(MeanVarB[blk,],8),collapse=' , '),']; Time=',round(timeOut-timeIn, 4),' sec.')
		}
	}
	postProb=colMeans(POST.PROB)

	return(list(b=b, POST.PROB=POST.PROB, postProb=postProb, MeanVarB=MeanVarB, MeanVarE=MeanVarE))
}



get_block_ids = function(M) {
  p=M@p
  i=M@i
  n=ncol(M)

  MIN=integer(n)
  MAX=integer(n)

  # For each column j, get min/max row index of its nonzeros
  for (j in 1:n) {
    idx=(p[j] + 1):p[j + 1]
    if (length(idx) == 0) {
      MIN[j]=j
      MAX[j]=j
    } else {
      rows=i[idx] + 1

      MIN[j]=rows[1]
      MAX[j]=rows[length(rows)]
    }
  }

  blk=integer(n)
  b=1
  blk[1]=b
  for (j in 2:n) {
    if (MIN[j] > MAX[j - 1]) {
      b=b + 1
    }
    blk[j]=b
  }
  blk
}
