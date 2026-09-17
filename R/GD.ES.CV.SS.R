# GD with early stopping using one trn/tst set
# It uses sufficient statistics
# Only accept scalar lambda

GD.ES.CV.SS<- function( X,y, trn,  centerX=TRUE,scaleX=FALSE, b=NULL, maxIter=300, 
                    learningRate=1/30,lambda=0, verbose=TRUE, ...){

    # Checking trn vector and  mapping it to integer 
    if(is.logical(trn)){ trn<-which(trn)}
    if(is.character(trn)){ 
    	if(is.null(rownames(X))){ stop('trn is charachter but X does not have rownames')}
    	trn=trn[trn%in%rownames(X)]
    	trn<-as.integer(factor(trn,levels=rownames(X)))
    }

    if(!is.integer(trn) | max(trn)>nrow(X) | min(trn)<1){stop('trn must be an integer vector with values between 1 and nrow(X)')}
    X=scale(X,center=centerX,scale=scaleX)
    y=y-mean(y)
  
    # Training/Testing data
    X_trn=X[trn,,drop=FALSE]
    X_tst=X[-trn,,drop=FALSE]
    y_trn=y[trn]
    y_tst=y[-trn]


    XX=crossprod(X_trn)
    Xy=crossprod(X_trn,y_trn)
    
    meanX2=mean(apply(FUN=function(x){sum(x^2)},X=X,MARGIN=2))
  
    if(is.null(b)){
        b=rep(rnorm(nrow(XX))/100000)
        names(b)=rownames(XX)
    }

    if(!(is(b,"vector") | is(b,"matrix") | is(b,"data.frame"))) stop("The prior estimates (b) must be in one of these formats: vector, matrix or data.frame with single column\n")

    Xy=as.matrix(Xy)
    b=as.matrix(b)

    nameWarningFlag=0
    if(is.null(rownames(XX)) | is.null(colnames(XX)) | is.null(rownames(Xy)) | is.null(rownames(b))){
        warning('Variant IDs are missing in one or more of inputs: XX, Xy, or b\n')
        nameWarningFlag=1
    }

    snp_list=Reduce(intersect, list(rownames(XX),rownames(Xy),rownames(b)))
    if (length(snp_list) == 0){ 
        warning('Variant IDs are not matching between inputs: XX, Xy, or b\n')
        nameWarningFlag=1
    }

    if ((nrow(XX)!=nrow(Xy) | nrow(XX)!=nrow(b)) & (nameWarningFlag==1)){
        stop('Distinct number of variants detected in inputs: XX, Xy, and b, while variant IDs are missing in one or more of inputs: XX, Xy, or b\n')
    }

    if (nameWarningFlag==0){
        XX=XX[snp_list,snp_list,drop = FALSE]
        Xy=Xy[snp_list,,drop = FALSE]
        b=b[snp_list,,drop = FALSE]
        if(verbose){
            message(length(snp_list), ' variants in common between XX, Xy, and b are retained\n')
        }
    } else {
        if(verbose){
            message(nrow(XX), ' variants are retained\n')
        }
    }

    p=nrow(XX)
    b0=rnorm(p)/10000
  	
    B=matrix(nrow=p,ncol=maxIter+1)
    Cor=numeric(maxIter+1)
    Cor[1]=corBootstrap(X_tst%*%b, y_tst,100)#cor(X_tst%*%b, y_tst)
    
    diag(XX)=diag(XX)+(lambda)
    LR=learningRate/meanX2

    B[,1]=b
    
    for(i in 2:ncol(B)){
        B[,i]=.Call("GRAD_DESC",XX, Xy, B[,i-1],p, 1, LR)
        Cor[i]=corBootstrap(X_tst%*%B[,i], y_tst,100)#cor(X_tst%*%B[,i], y_tst)
        if (Cor[i]<max(Cor[i-1])*0.99) {break}
    }
        
    stopIter=i-1

    XX=XX+crossprod(X_tst)
    Xy=Xy+crossprod(X_tst,y_tst)
    
    if(verbose){
    	message(' Stopped at cycle ',stopIter)
    	message('...computing final estimates.')
    }

    b=GD(XX=XX,Xy=Xy,lambda=lambda,learningRate=learningRate,nIter=stopIter)
    return(list(b=b,stopIter=stopIter))
}
