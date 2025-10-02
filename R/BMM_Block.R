BMM_Block=function(XX, Xy, B, my, vy, n, nIter=1200, burnIn=200, thin=5, R2=0.25,
	        nComp=matrix(ncol(B)), K=1/nComp, df0.E=5, S0.E=vy*(1-R2)*df0.E, df0.b=rep(10,nComp), 
	        priorProb=rep(1/nComp,nComp), priorCounts=rep(2*nComp,nComp), verbose=TRUE) {
	
	if(!(is(XX,"dgCMatrix"))) stop("XX must be a dgCMatrix\n")

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

 	if (is.vector(B)) {
 		if (is.null(names(B))) {
 			stop("The prior estimates vector (B) must have variant IDs as names\n")
 		}
 		B=as.data.frame(B)
 	} else if (is.matrix(B) | is.data.frame(B)) {
 		if (is.null(rownames(B))) {
 			stop("The prior estimates matrix (B) must have variant IDs as row names\n")
 		}
 		B=as.data.frame(B)
 	} else {
 		stop("B must be in one of these formats: vector, matrix, data.frame\n")
 	}
	
 	snp_list=Reduce(intersect, list(rownames(XX),names(Xy),rownames(B)))

 	if (length(snp_list) == 0){ 
 		stop("No matched variants in XX, Xy, and prior\n")
 	}else{
 		if(verbose){
 			message(' There were ',length(snp_list), ' variants in common between XX, Xy, and the prior.\n')
 		}
 	}
	
 	XX=XX[snp_list,snp_list,drop = FALSE]
 	Xy=Xy[snp_list]
 	B=B[snp_list,,drop = FALSE]

 	B=as.matrix(B)
	
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
		Xyblk=Xy[snps]

		if (length(snps==1)) {
			XXblk=as.matrix(XXblk)
			colnames(XXblk)=snps
            rownames(XXblk)=snps
			Bblk=matrix(Bblk, nrow=1)
			rownames(Bblk)=snps
		}

		fm=BMM(XX=XXblk, Xy=Xyblk, B=Bblk, n=n, my=my, vy=vy, nIter=nIter, burnIn=burnIn, thin=5, R2=R2/nrow(XX)*length(snps), fixVarE=FALSE, fixVarB=rep(FALSE,ncol(B)), verbose=FALSE)
		
		b[snps]=fm$b
		POST.PROB[ld_index == blk,]=fm$POST.PROB
		MeanVarB[blk,]=fm$postMeanVarB
		MeanVarE[blk]=mean(fm$samplesVarE)

		TimeOut=proc.time()[3]

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
