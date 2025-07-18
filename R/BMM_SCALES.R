## A test to add prior scales to the prior means

BMM_SCALES=function(XX, Xy, B,SCALES_B=NULL, my, vy, n, nIter=150, burnIn=50, thin=5, R2=0.25,
           nComp=matrix(ncol(B)), K=1/nComp, df0.E=5, S0.E=vy*(1-R2)*df0.E, df0.b=rep(10,nComp), 
           priorProb=rep(1/nComp,nComp), priorCounts=rep(2*nComp,nComp), verbose=TRUE){




##


linearIndex=function(rowIndex,colIndex,nRow){
   ans = (colIndex - 1) * nRow+rowIndex
   return(ans)
}

getVarB=function(varB,SCALES_B,d){
   VARS=rep(NA,nrow(SCALES_B))
   groups=unique(d)
   for(g in 1:length(groups)){
      tmp=(d==groups[g])
      if(any(tmp)){
         tmp2=linearIndex(rowIndex=which(tmp),colIndex=d[tmp],nRow=nrow(SCALES_B))
         VARS[tmp]=SCALES_B[tmp2]*varB[d[tmp]]
      }
   }
   return(VARS)
 }


if(FALSE){

library(GPTL)
data(mice)
X=scale(mice.X[,1:55],center=TRUE,scale=FALSE)
p=ncol(X)
n=nrow(X)
b=rep(0,p)
b[seq(from=5,to=50,by=5)]=rep(c(-1,1),each=5)
signal=X%*%b
error=rnorm(n,sd=sd(signal)*2)
y=signal+error


XX=crossprod(X)
Xy=crossprod(X,y)
my=mean(y)
vy=var(y)



 B=cbind(0,b); SCALES_B=NULL; nIter=150; burnIn=50; thin=5; R2=0.25;
           nComp=matrix(ncol(B)); K=1/nComp; df0.E=5; S0.E=vy*(1-R2)*df0.E; df0.b=rep(10,nComp);
           priorProb=rep(1/nComp,nComp); priorCounts=rep(2*nComp,nComp); verbose=TRUE
 #

 rownames(B)=colnames(X)
}
##

 B=as.matrix(B)


 if(is.null(SCALES_B)){
   SCALES_B=matrix(nrow=nrow(B),ncol=ncol(B),1)
 }else{
    SCALES_B=as.matrix(SCALES_B)
   for(i in 1:ncol(SCALES_B)){
      SCALES_B[,i]=SCALES_B[,i]/mean(SCALES_B[,i])
   }

 }

 p=length(Xy) 
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
 RSS=vy*(n-1)+crossprod(b,XX)%*%b-2*crossprod(b,Xy)
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
     
     VARS=getVarB(varB,SCALES_B,d)
     if(is(XX,"dgCMatrix"))
     {
      #Sparse matrix
      tmp=sample_effects_new_sparse(C=XX,rhs=Xy,b=b,d=d,B0=B,varE=varE,varB=VARS,RSS=RSS)
     }else{
      #Dense matrix
      tmp=sample_effects_new(C=XX,rhs=Xy,b=b,d=d,B0=B,varE=varE,varB=VARS,RSS=RSS)
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
          bStar=b[tmp]-B[tmp,k]
          SS=SS+sum(bStar^2)
       }  
       varB[k]=SS/rchisq(df=DF+df0.b[k],n=1) 
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

