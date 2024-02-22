# A Gibbs Sampler for a Bayesian Mixture Model
# C=X'X, rhs=X'y, my=mean(y), vy=var(y),n=length(y)
# Note: we don't include an intercept, so, before computing X'y, X'X, you should center X and y arount their means.
# B0, a matrix with prior estiamtes of marker effects, as many columns as prior estimates. We also suggest including one column
#     full of zeroes, to allow for just plain shrinkage towards zero, if needed (this is the default value)
# nIter, burnIn: the number of iterations and burnIn used for the Gibbs Sampler
# R2: the prior proportion of variance of y explained by the model (this is used to derive hyper-parameters for the prior variances)
# priorProb a vector of prior probabilities for the mixture compoentns (as many entries as columsn in B0), we internally scale it to add up to one
# priorCounts, the number of counts associated to priorProb, using priorCounts very large (e.g., 1e8) fixes the prior probabilities. The default value is 2*ncol(B0)
#  To do: sample error variance ; add prior probabilities for the mixtures (right now equivalent to 1/nClasses); priors for the variances
##

BMM_new=function(C,rhs,my,vy,n,B0=matrix(nrow=ncol(C),ncol=1,0),nIter=150,burnIn=50,thin=5,R2=.1,
	        nComp=matrix(ncol(B0)), df0.E=5,S0.E=vy*(1-R2)*df0.E,df0.b=rep(10,nComp), 
	        priorProb=rep(1/nComp,nComp),priorCounts=rep(2*nComp,nComp),verbose=TRUE){
 B0=as.matrix(B0)
 # nIter=150;burnIn=50;R2=.5;nComp=matrix(ncol(B0));df0.E=5;S0.E=vy*R2*df0.E;df0.b=rep(5,nComp);alpha=.1;my=mean(y); vy=var(y); B0=cbind(rep(0,p),-1,1)
 p=ncol(C) 
 b=rep(0,p)
 d=rep(1,p) # indicator variable for the group
 POST.PROB=matrix(nrow=p,ncol=nComp,0)
	
 S0.b=c(df0.b)*c(vy)*c(R2/2)/c(sum(diag(C))/n) # dividing R2/10 assumes that most of the vairance is between components.
 varB=S0.b/df0.b

 priorProb=priorProb/sum(priorProb)
 compProb=priorProb

 postMeanCompProb=rep(0,nComp)
	
 postMeanB=rep(0,p)
 postMeanVarB=rep(0,nComp)
 postProb=rep(0,nComp)

 RSS=vy*(n-1)+crossprod(b,C)%*%b-2*crossprod(b,rhs)
  
 varE=vy*(1-R2)
 counts=priorCounts/as.vector(nComp)
	
 PROBS=matrix(nrow=nComp,ncol=p)
 	
 timeEffects=0
 timeProb=0
 timeApply=0

 samplesVarB=matrix(nrow=nIter,ncol=nComp,NA)
 samplesB=matrix(nrow=nIter,ncol=p,NA)
 samplesVarE=rep(NA,nIter)

weightPostMeans=1/round((nIter-burnIn)/thin)
for(i in 1:nIter){
         
	 ## Sampling effects
	 timeIn=proc.time()[3]
	  tmp=sample_effects_new(C=C,rhs=rhs,b=b,d=d,B0=B0,varE=varE,varB=varB,RSS=RSS)
    b=tmp[1]
    RSS=tmp[2]
	 timeEffects=timeEffects+(proc.time()[3]-timeIn)
	## End of C-code
	 samplesB[i,]=b
	 
	 ## Sampling mixture components 
	timeIn=proc.time()[3]
	 for(k in 1:nComp){
	 PROBS[k,]=dnorm(b,mean=B0[,k],sd=sqrt(varB[k]))#*compProb[k]	
	 }
 	timeProb=timeProb+(proc.time()[3]-timeIn)
	tiemIn=proc.time()[3] 
	  #d=apply(FUN=sample,x=1:nComp,X=PROBS,size=1,MARGIN=1,replace=TRUE)
	  # d=sampleComp(PROBS)
	  d=rMultinom(PROBS)
        timeApply=timeApply+(proc.time()[3]-timeIn) 
	
	 ## Sampling the variance and the prior probabilities of the mixture components
	 for(k in 1:nComp){
		 tmp=(d==k)
		 DF=sum(tmp)
		 SS=S0.b[k]
		 if(DF>0){
			 bStar=b[tmp]-B0[tmp,k]
			 SS=SS+sum(bStar^2)
		 }	 
		 varB[k]=SS/rchisq(df=DF+df0.b[k],n=1) 
		 counts[k]=DF
	 }
	 samplesVarB[i,]=varB
	
  # Sampling the probability of each component 
	compProb=rDirichlet(counts+priorCounts)

	# Sampling the error variance
   SS=RSS+S0.E
   DF=n+df0.E
	 varE=SS/rchisq(df=DF,n=1)
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
	if(verbose){ print(i) }
  } 

 if(verbose){
   message('Time Effects= ', timeEffects)
   message('Time Prob= ', timeProb)
   message('Time Apply= ', timeApply)
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
 PROB=sweep(PROB,FUN='/',MARGIN=2,STAT=rSum)
 u=runif(n)
 PROB=sweep(PROB,FUN='>',MARGIN=2,STAT=u)
 apply(FUN=which.first,X=PROB,MARGIN=2)
}

