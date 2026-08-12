# A function to perform Grad. Desc. that takes X and y, instead of XX and Xy

GDXy2<-function(X,y,b=NULL,learningRate=1/100,nIter=100)
{
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


	for(i in 2:(nIter+1))
	{
		.Call("GRAD_DESC_Xy", X, x2, b, e, nrow(X), ncol(X), 1, learningRate)
		B[,i]=b
	}
	
    #for(i in 2:(nIter+1))
	#{
    #  
    #  	for(j in 1:p)
    #  	{
    #        xj=X[,j]
    #        e=e+xj*b[j]
    #        rhs=crossprod(xj,e)
    #        dL=(x2[j]*b[j]-rhs)
    #        b[j]=b[j]-learningRate*dL
    #        e=e-xj*b[j]
    #  	}
    #  	B[,i]=b
   	#}

    return(B)

}
