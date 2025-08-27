if(FALSE){
 library(remotes)
 install_github('https://github.com/quantGen/GPTL')
 }

 library(GPTL)

 rm(list=ls())

 ## Simple OLS problem

 X=matrix(nrow=100,ncol=5,rnorm(500))
 b=rnorm(5)*2
 signal=X%*%b
 y=signal+rnorm(nrow(X),sd=sd(signal))

 bHat1=coef(lm(y~X-1))

 XX=crossprod(X)
 Xy=crossprod(X,y)

 colnames(XX)=paste0('snp_', 1:5)
 rownames(XX)=paste0('snp_', 1:5)
 rownames(Xy)=paste0('snp_', 1:5)

 prior=rep(0,5)
 names(prior)=paste0('snp_', 1:5)

 bHat2=GD(XX=XX, Xy=Xy, b=prior, learning_rate=1/50, nIter=100, returnPath=F)
 
 plot(bHat1,bHat2,cex=1.5,col=2);abline(a=0,b=1)

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
