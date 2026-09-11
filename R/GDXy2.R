# A function to perform Grad. Desc. that takes X and y, instead of XX and Xy

GDXy2<-function(X,y,b=NULL,learningRate=1/50,nIter=100,threshold=3/100,earlyStop=FALSE,returnPath=FALSE,verbose=FALSE){
   
   p=ncol(X)
   X=scale(X,center=TRUE,scale=FALSE)
   y=scale(y,center=TRUE)
    
   if(is.null(b)){
        b=rep(0, p)
        names(b)=colnames(X)
   }
   x2=colSums(X^2)

   learningRate=learningRate/mean(x2)
   
   e=y-X%*%b
   
   B=matrix(nrow=p,ncol=nIter+1,NA)
   B[,1]=b

   if(earlyStop){ 
      RSS=rep(NA_real_,nIter+1)
      RSS[1]=sum(e^2)
   }
   
   for(i in 2:(nIter+1)){
        .Call("GRAD_DESC_Xy", X, x2, b, e, nrow(X), ncol(X), 1, learningRate)
        B[,i]=b
        if(earlyStop){
            RSS[i]=sum((y-X%*%b)^2)
            propChange=(1-RSS[i]/RSS[i-1])
            if(verbose){ print(RSS[i]) }
            if( propChange<threshold){ 
                B=B[,1:i]
                RSS=RSS[1:i]
                break() 
            }
         }
   }
   if(returnPath){
      B=B[,-1,drop=FALSE]
   }else{
      B=B[,ncol(B)]
   }
   return(B)
}
