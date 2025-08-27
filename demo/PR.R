if(FALSE){
 library(remotes)
 install_github('https://github.com/quantGen/GPTL')
 }

 library(GPTL)

 rm(list=ls())

## Ridge Regression problem

 library(BGLR)
 data(wheat)

 X=scale(wheat.X,center=TRUE,scale=FALSE)
 y=wheat.Y[,1]

 XX=crossprod(X)
 Xy=crossprod(X,y)
 n=nrow(X)

 prior=rep(0,nrow(XX))
 names(prior)=rownames(XX)

 lambda=sum(diag(XX))/n


 bHat=RR(XX,Xy,lambda=lambda)
 bHat2=PR(XX=XX, Xy=Xy, b=prior, alpha=0, lambda=lambda, conv_threshold=1e-4,
          maxIter=1000, returnPath=FALSE)$B
 plot(bHat,bHat2);abline(a=0,b=1)
