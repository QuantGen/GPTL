#GD.CVParallel=function(x,y,nRep=1,nFolds=5,acc=cor,suppress_warnings=TRUE,nIter=10,mc.cores=detectCores(),...){
#    TMP=mclapply(FUN=GD.CV,X=X,y=y,X=1:nFolds,nFolds=nFolds,acc=acc,suppress_warnings=TRUE,nIter=10,mc.cores=mc.cores,...)
#    return(TMP)
#}
#tmp<-mclapply(FUN=GD.CV,x=X,y=y,notUse=1,lambda=c(.1,1),nFolds=5,acc=cor,suppress_warnings=TRUE,nIter=10,returnPath=T,mc.cores=8,X=1:12)
## To parallelize GD.CV we can
#   - Replace the loop over reps using foreach
#   - Use mclapply with a dummy argument (notUse=1, below), the current implementation with mclappply (line 5 above) does not offer any benefit
GD.CV=function(x,y,nFolds=5,nRep=1,acc=cor,suppress_warnings=TRUE,nIter=10,lambda=001,notUse=1,...){
    print(notUse)
    DIM=c(nIter,length(lambda),nRep)
    ACC=array(dim=DIM)
    N=nrow(x)
    for(i in 1:nRep){
        folds=sample(1:nFolds,size=N,replace=TRUE)
        AVG_ACC=array(dim=c(nIter,length(lambda)),0)

        for(j in 1:nFolds){
            tst=which(folds==j)
            XTRN=x[-tst,,drop=FALSE]
            yTRN=y[-tst]
            XX=crossprod(XTRN)
            Xy=crossprod(XTRN,yTRN)
            B=GD(XX=XX,Xy=Xy,nIter=nIter,lambda=lambda,...)
            XTST=x[tst,,drop=FALSE]
            yTST=y[tst]

            
            for(Row in 1:nrow(AVG_ACC)){
                for(Col in 1:ncol(AVG_ACC)){
                    if(suppress_warnings){
                        suppressWarnings( AVG_ACC[Row,Col]<-AVG_ACC[Row,Col]+acc(yTST,XTST%*%B[,Row,Col]) )
                    }        
                }
            }
            
        }
        AVG_ACC=AVG_ACC/nFolds
        ACC[,,i]=AVG_ACC
    }
    TMP=list(lambda=lambda,nIter=nIter,ACC=apply(FUN=mean,MARGIN=c(1,2),X=ACC),SD=apply(FUN=sd,MARGIN=c(1,2),X=ACC))

    rownames(TMP$ACC)=paste0('Iter_',1:nIter)
    colnames(TMP$ACC)=paste0('Lambda_',1:length(lambda))
    rownames(TMP$SD)=paste0('Iter_',1:nIter)
    colnames(TMP$SD)=paste0('Lambda_',1:length(lambda))
    
    return(TMP)
}
