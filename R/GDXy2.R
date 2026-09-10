# A function to perform Grad. Desc. that takes X and y, instead of XX and Xy

GDXy2<-function(X,y,b=NULL,learningRate=1/100,nIter=100,minChange=1,earlyStop=FALSE){

   p=ncol(X)
   X=scale(X,center=TRUE,scale=FALSE)
   y=scale(y,center=TRUE)
    
   if(is.null(b)){
        b=rep(0, p)
        names(b)=colnames(X)
   }
   x2=colSums(X^2)

   learningRate=learningRate/mean(x2)

   if(earlyStop){
      minChange=learningRate*minChange
   }

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
            if( (1-(RSS[i]/RSS[i-1]))<minChange){ 
                B=B[,1:i]
                RSS=RSS[1:i]
                print(RSS)
                break() 
            }
         }
    }
   return(B)
}
