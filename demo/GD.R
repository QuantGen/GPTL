#library(remotes)
#install_github('https://github.com/quantGen/GPTL')


library(GPTL)


 rm(list=ls())

 ## Simple OLS problem

 X=matrix(nrow=100,ncol=5,rnorm(500))
 b=rnorm(5)*2
 signal=X%*%b
 y=signal+rnorm(nrow(X),sd=sd(signal))

 bHat1=coef(lm(y~X))

 XX=crossprod(cbind(1,X),cbind(1,X))
 Xy=crossprod(cbind(1,X),y)

 bHat2=GD(XX,Xy,lambda=0,nIter=100,learning_rate=.1);plot(bHat1,bHat2,cex=1.5,col=2);abline(a=0,b=1)


## Ridge Regression problem



 library(BGLR)
 data(wheat)

 X=scale(wheat.X,center=TRUE,scale=FALSE)
 y=wheat.Y[,1]

 XX=crossprod(X)
 Xy=crossprod(X,y)
 n=nrow(X)

 lambda=sum(diag(XX))/n


 bHat=RR(XX,Xy,lambda=lambda)
 bHat2=GD(XX,Xy,lambda=lambda,nIter=100,learning_rate=.1);plot(bHat,bHat2);abline(a=0,b=1)


