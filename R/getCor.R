getCor=function(XX, Xy, yy, B, verbose=FALSE){
 
 if(!(is(XX,"matrix") | is(XX,"dgCMatrix"))) stop("XX must be a matrix or dgCMatrix\n")

 if(is.null(rownames(XX)) | is.null(colnames(XX))){
 	stop('XX must have variant IDs as row and column names\n')
 }
    
 if (!all(rownames(XX) == colnames(XX))) stop("Row and column names in XX not match\n")

 if (is.vector(Xy)) {
 	if (is.null(names(Xy))) {
 		stop("Xy must have variant IDs as names\n")
 	}
 } else if (is.matrix(Xy) | is.data.frame(Xy)) {
 	if (is.null(rownames(Xy))) {
 		stop("Xy must have variant IDs as row names\n")
 	}
 	XyName=rownames(Xy)
    Xy=as.vector(as.matrix(Xy))
    names(Xy)=XyName
 } else {
 	stop("Xy must be in one of these formats: vector, matrix, data.frame\n")
 }

 if (is.vector(B)) {
 	if (is.null(names(B))) {
 		stop("The prior estimates vector (B) must have variant IDs as names\n")
 	}
 	B=as.data.frame(B)
 } else if (is.matrix(B) | is.data.frame(B)) {
 	if (is.null(rownames(B))) {
 		stop("The prior estimates matrix (B) must have variant IDs as row names\n")
 	}
 	B=as.data.frame(B)
 } else {
 	stop("B must be in one of these formats: vector, matrix, data.frame\n")
 }
	
 snp_list=Reduce(intersect, list(rownames(XX),names(Xy),rownames(B)))

 if (length(snp_list) == 0){ 
 	stop("No matched variants in XX, Xy, and prior\n")
 }else{
 	if(verbose){
 		message(' There were ',length(snp_list), ' variants in common between XX, Xy, and the prior.\n')
 	}
 }
	
 XX=XX[snp_list,snp_list,drop = FALSE]
 Xy=Xy[snp_list]
 B=B[snp_list,,drop = FALSE]

 Cor=numeric(ncol(B))
 for (i in 1:ncol(B)){
	 b=B[,i]
	 Cor[i]=crossprod(b,Xy)/sqrt(yy*t(b)%*%XX%*%b)
 }
 return(Cor)
}
