
GD.CV=function(X,y,nTst=round(nrow(X)*0.2),nRep=10,nIter=100,seed=12345,...){

	COR=matrix(nrow=nRep,ncol=nIter)
	set.seed(seed)
	n=nrow(X)
	p=ncol(X)

	if(!exists('b')){
		b=rep(0,p)
		names(b)=colnames(X)
	}

	for(i in 1:nRep){
		tst=sample(1:n,size=nTst)
		trn=(1:n)[-tst]

		B=GDXy(X=X[trn,],y=y[trn],nIter=nIter+1,b=b,returnPath=TRUE,...)[,-1]
		for(j in 1:nIter){
			COR[i,j]=cor(y[tst],X[tst,]%*%B[,j])
		}
	}
	DF=data.frame(Interation=1:nIter,Cor=colMeans(COR),SD=apply(FUN=sd,X=COR,MARGIN=2))
	return(DF)
}

