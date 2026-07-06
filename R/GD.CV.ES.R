# Comparing to GD.CV(), GD.CV.ES() will take as inputs X_trn, y_trn, X_tst, y_tst; evaluate the 
# prediction correlation each iteration; and stop the algorithm once the prediction correlation
# starts to drop. GD.CV.ES() is more time efficient as it does not require running the full path.

GD.CV.ES<- function(X_trn, y_trn, X_tst, y_tst, centerX=TRUE,scaleX=FALSE, b=NULL, maxIter=1000, learningRate=1/50, lambda=0, verbose=TRUE, ...){
GD.CV.ES<- function(XX, Xy, b=NULL, maxIter=10, learningRate=1/50, lambda=0, verbose=TRUE){

    X_trn=scale(X_trn,center=centerX,scale=scaleX)
    y_trn=scale(y_trn,center=TRUE)

    XX=crossprod(X_trn)
    Xy=crossprod(X_trn,y_trn)
    
    if(is.null(b)){
        b=rep(0, nrow(XX))
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
    b0=rep(0,p)
  	
    previous_lambda=0
    B=array(dim=c(p,maxIter+1,length(lambda)))
    Cor=numeric(maxIter+1)
    Cor[1]=cor(X_tst%*%b, y_tst)

    for(h in 1:length(lambda))
    {
    	if(is(XX,"dgCMatrix"))
    	{	
    	    #Sparse matrix
        	Matrix::diag(XX)<-Matrix::diag(XX) + (lambda[h]- previous_lambda)
        	LR=learningRate/mean(Matrix::diag(XX))
        }else{
        	#Dense matrix
        	diag(XX)=diag(XX)+(lambda[h]- previous_lambda)
        	LR=learningRate/mean(diag(XX))
        }
  
        #if( lambda0>0 ){    
        #    Xy=Xy+(lambda[h]-previous_lambda)*lambda0*b0 
        #}
        previous_lambda=lambda[h]

        B[,1,h]=b
            
        if(is(XX,"dgCMatrix"))
        {
            #Sparse matrix
            for(i in 2:ncol(B))
            {
            	B[,i,h]=.Call("GRAD_DESC_sparse",XX@x,XX@p,XX@i,Xy, B[,i-1,h],p, 1, LR)
                Cor[i]=cor(X_tst%*%B[,i,h], y_tst)
                if (Cor[i]<Cor[i-1]) {break}
            }
        }else{
            #Dense matrix
            for(i in 2:ncol(B))
            {
            	B[,i,h]=.Call("GRAD_DESC",XX, Xy, B[,i-1,h],p, 1, LR)
                Cor[i]=cor(X_tst%*%B[,i,h], y_tst)
                if (Cor[i]<Cor[i-1]) {break}
            }
        }
    }
    
    

    iterations=paste0('iter_',0:maxIter)
    
    dimnames(B)=list(colnames(XX),iterations,paste0('lambda_',lambda))

    #B=B[,i,,drop=TRUE]
    
    #return(B)
    return(list(B=B[,i,,drop=TRUE], B_path=B[,2:i,,drop=TRUE], Cor=Cor[2:i], stopIter=i-1))
}
