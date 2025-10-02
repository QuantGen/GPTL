# A Gibbs Sampler for a Bayesian Mixture Model
# XX=X'X, Xy=X'y, my=mean(y), vy=var(y),n=length(y)
# Note: we don't include an intercept, so, before computing X'y, X'X, you should center X and y arount their means.
# B, a matrix with prior estiamtes of marker effects, as many columns as prior estimates. We also suggest including one column
#     full of zeroes, to allow for just plain shrinkage towards zero, if needed (this is the default value)
# nIter, burnIn: the number of iterations and burnIn used for the Gibbs Sampler
# R2: the prior proportion of variance of y explained by the model (this is used to derive hyper-parameters for the prior variances)
# priorProb a vector of prior probabilities for the mixture compoentns (as many entries as columsn in B), we internally scale it to add up to one
# priorCounts, the number of counts associated to priorProb, using priorCounts very large (e.g., 1e8) fixes the prior probabilities. The default value is 2*ncol(B)
#  To do: sample error variance ; add prior probabilities for the mixtures (right now equivalent to 1/nClasses); priors for the variances
##

BMM=function(XX, Xy, B, my, vy, n, nIter=1200, burnIn=200, thin=5, R2=0.25,
	        nComp=matrix(ncol(B)), K=1/nComp, df0.E=5, S0.E=vy*(1-R2)*df0.E, df0.b=rep(10,nComp), 
	        priorProb=rep(1/nComp,nComp), priorCounts=rep(2*nComp,nComp), verbose=TRUE, fixVarE=FALSE,fixVarB=rep(FALSE,nComp)){

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
 # nIter=150;burnIn=50;R2=.5;nComp=matrix(ncol(B));df0.E=5;S0.E=vy*R2*df0.E;df0.b=rep(5,nComp);alpha=.1;my=mean(y); vy=var(y); B=cbind(rep(0,p),-1,1)
 p=nrow(XX)
 b=rowMeans(B)
 d=rep(1,p) # indicator variable for the group
 POST.PROB=matrix(nrow=p,ncol=nComp,0)
 
 # dividing R2/10 assumes that most of the vairance is between components.
 if(is(XX,"dgCMatrix"))
 {
 	S0.b=c(df0.b)*c(vy)*c(R2*K)/c(sum(as.vector(Matrix::diag(XX)))/n)
 }else{
 	S0.b=c(df0.b)*c(vy)*c(R2*K)/c(sum(diag(XX))/n) 
 }
 varB=(S0.b/df0.b)

 priorProb=priorProb/sum(priorProb)
 compProb=priorProb

 postMeanCompProb=rep(0,nComp)
	
 postMeanB=rep(0,p)
 postMeanVarB=rep(0,nComp)
 postProb=rep(0,nComp)

 #Need to convert to numeric because if XX is sparse, the result is an object 
 #of class dgeMatrix
if(fixVarE){
	RSS=vy*(n-1)*(1-R2)
}else{
	RSS=vy*(n-1)+crossprod(b,XX)%*%b-2*crossprod(b,Xy)
}
 RSS=as.numeric(RSS)

 varE=RSS/n
	
 counts=priorCounts/as.vector(nComp)
	
 PROBS=matrix(nrow=nComp,ncol=p)
 	
 timeEffects=0; timeProb=0; timeApply=0

 samplesVarB=matrix(nrow=nIter,ncol=nComp,NA)
 samplesB=matrix(nrow=nIter,ncol=p,NA)
 samplesVarE=rep(NA,nIter)

 weightPostMeans=1/round((nIter-burnIn)/thin)
	
 for(i in 1:nIter){        
	  ## Sampling effects
	  timeIn=proc.time()[3]
	  time0=timeIn
	  
	  if(is(XX,"dgCMatrix"))
	  {
	  	#Sparse matrix
	  	tmp=sample_effects_new_sparse(C=XX,rhs=Xy,b=b,d=d,B0=B,varE=varE,varB=varB[d],RSS=RSS)
	  }else{
	  	#Dense matrix
	  	tmp=sample_effects_new(C=XX,rhs=Xy,b=b,d=d,B0=B,varE=varE,varB=varB[d],RSS=RSS)
	  }
	  
      b=tmp[[1]]
      RSS=tmp[[2]]
	 timeEffects=timeEffects+(proc.time()[3]-timeIn)
	## End of C-code
	 samplesB[i,]=b
	 
	 ## Sampling mixture components 
	timeIn=proc.time()[3]
	 for(k in 1:nComp){
	 PROBS[k,]=dnorm(b,mean=B[,k],sd=sqrt(varB[k]))#*compProb[k]	
	 }
 	timeProb=timeProb+(proc.time()[3]-timeIn)
	timeIn=proc.time()[3] 
	  #d=apply(FUN=sample,x=1:nComp,X=PROBS,size=1,MARGIN=1,replace=TRUE)
	  # d=sampleComp(PROBS)
	  d=rMultinom(PROBS)
        timeApply=timeApply+(proc.time()[3]-timeIn) 
	
	 ## Sampling the variance and the prior probabilities of the mixture components
	 for(k in 1:nComp){
		tmp=(d==k)
		DF=sum(tmp)
		 
		if(!fixVarB[k]){
		 	SS=S0.b[k]
		 	if(DF>0){
			 bStar=b[tmp]-B[tmp,k]
			 SS=SS+sum(bStar^2)
		 	}	 
		 	varB[k]=SS/rchisq(df=DF+df0.b[k],n=1) 
		}
		counts[k]=DF
	 }
	 samplesVarB[i,]=varB
	
  # Sampling the probability of each component 
	compProb=rDirichlet(counts+priorCounts)
	
	# Sampling the error variance
  	#RSS2=vy*(n-1)+crossprod(b,XX)%*%b-2*crossprod(b,Xy)
	SS=RSS
	#print(c(RSS,RSS2)/n)
        DF=n
	if (!fixVarE) {
		 varE=SS/rchisq(df=DF,n=1)
	 }
        samplesVarE[i]=varE
  
	## computing posterior means 
	if(i>burnIn&(i%%thin==0)){
	 postMeanVarB= postMeanVarB+varB*weightPostMeans
	 postProb=postProb+compProb*weightPostMeans
	 postMeanB=postMeanB+b*weightPostMeans
	 for(k in 1:nComp){
		 tmp=(d==k)
		 POST.PROB[tmp,k]=POST.PROB[tmp,k]+weightPostMeans
	  }
        } 
	if(verbose){ 
		
		    message(' Cycle: ',i,' varE=',round(varE,4),'; varB=[', paste(round(varB,8),collapse=' , '),']; Time=',round(proc.time()[3]-time0,4),' sec.') 
	   
		   }
  } 


 return(list(b=postMeanB,POST.PROB=POST.PROB,postMeanVarB=postMeanVarB,postProb=postProb,
	     	samplesVarB=samplesVarB,samplesB=samplesB,samplesVarE=samplesVarE))
}




## A function to sample from a Dirichlet

rDirichlet=function(counts){
    x=rgamma(n=length(counts),shape=counts)
	return(x/sum(x))
}

## Two functions to sample integers with pre-specified probabilities
which.first=function(x){
	which(x)[1]
}

sampleComp=function(PROB){
 n=nrow(PROB)
 rSum=rowSums(PROB)
 PROB=apply(FUN=cumsum,X=PROB,MARGIN=1)
 PROB=sweep(PROB,FUN='/',MARGIN=2,STATS=rSum)
 u=runif(n)
 PROB=sweep(PROB,FUN='>',MARGIN=2,STATS=u)
 apply(FUN=which.first,X=PROB,MARGIN=2)
}


