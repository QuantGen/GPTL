
 
  library(remotes)
install_github('https://github.com/QuantGen/GPTL')

 library(GPTL)
 library(BGLR)
 data(wheat)
 X=wheat.X[,1:100]

 p=ncol(X)
 n=nrow(X)

 QTL=seq(from=10,to=100,by=10)
 b0=rep(0,p)
 b0[QTL]=rep(c(-1,1),each=5)

 
 signal=X%*%b0
 error=rnorm(n=n,sd=sd(signal))
 y=signal+error

 varE=var(error)

 C=crossprod(X)
 rhs=crossprod(X,y)

 B0=cbind(rep(0,p),-1,1)
 

 system.time( 
     tmp<-BMM(C=C,rhs=rhs,my=mean(y),vy=var(y),n=n,verbose=FALSE,
              B0=cbind(rep(0,p),-1,1),nIter=1000,burnIn=100) 
            )


 colQTL=ifelse(b0==0,8,ifelse(b0==1,2,4))
 pointQTL=ifelse(colQTL==8,1,19)
  par(mfrow=c(3,1))
  plot(tmp$POST.PROB[,1],col=colQTL,pch=pointQTL,ylim=c(0,1))
  plot(tmp$POST.PROB[,2],col=colQTL,pch=pointQTL,ylim=c(0,1))
  plot(tmp$POST.PROB[,3],col=colQTL,pch=pointQTL,ylim=c(0,1))

