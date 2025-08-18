# This function takes:
    #  ld reference panel (ld a square symmetric correlation matrix),
    #  gwas results 
    #  And prior estimats (B)
 # it matches the arguments (based on rownames, i.e., snp ids) and computes
 # sufficient statistics (X'X and X'y) and prior estimats (B) (if provided).

getSS=function(ld,gwas,B=NULL,verbose=TRUE){
    if(is.null(rownames(ld))){
        stop('The LD reference panel (ld) must have SNP IDs as rownames')
    }

    if(is.null(rownames(gwas))){
        stop('The GWAS table (gwas) must have SNP IDs as rownames')
    }

    if (!all(c('beta', 'se', 'n', 'allele_freq') %in% colnames(gwas))) {
        stop("Must provide GWAS results that consist of columns: beta (variant effects), se (variant standard errors), n (sample sizes for GWAS), allele_freq (variant allele frequency)\n")
    }
           
    if (!is.null(B)) {
        if (is.vector(B)) {
            if (is.null(names(B))) {
                stop("The prior estimates vector (B) must have SNP IDs as names\n")
            }
            B=as.data.frame(B)
        } else if (is.matrix(B) | is.data.frame(B)) {
            if (is.null(rownames(B))) {
                stop("The prior estimates matrix (B) must have SNP IDs as rownames\n")
            }
            B=as.data.frame(B)
        } else {
            stop("B must be in one of these formats: vector, matrix, data.frame\n")
        }
    }
    
    if (is.null(B)) {
        snp_list=intersect(rownames(ld),rownames(gwas))
    } else {
        snp_list=Reduce(intersect, list(rownames(ld),rownames(gwas),rownames(B)))
    }
    
    if (length(snp_list) == 0){ 
        stop("No matched SNPs in LD, GWAS, and prior\n")
    }else{
        if(verbose){
            if (is.null(B)) {
                message(' There were ',length(snp_list), ' in common between the LD reference panel, and the GWAS.')
            } else {
                message(' There were ',length(snp_list), ' in common between the LD reference panel, the GWAS and, the prior.')
            }
        }
    }

    ld=ld[snp_list,snp_list,drop = FALSE]
    gwas=gwas[snp_list,,drop = FALSE]
    if (!is.null(B)) {
        B=B[snp_list,,drop = FALSE]
    }
    p=nrow(gwas)

    stopifnot(all(rownames(gwas)==rownames(ld)))
  
    if (!is.null(B)) {
        stopifnot(all(all(rownames(B)==rownames(gwas)),all(rownames(gwas)==rownames(ld))))
    }

    allele_freq=gwas$allele_freq
    beta=gwas$beta
    n_gwas=gwas$n
    V=(2 * allele_freq * (1-allele_freq))
    D=Diagonal(x=sqrt(V)*sqrt(gwas$n))

    XX=(D%*%ld)
    XX=XX%*%D
    rownames(XX)=rownames(ld)
    colnames(XX)=colnames(ld)

    Xy=Matrix::diag(XX)*beta
    names(Xy)=rownames(gwas)
  
    if (is.null(B)) {
        out=list(XX=XX,Xy=Xy,n=mean(gwas$n))
    } else {
        out=list(XX=XX,Xy=Xy,n=mean(gwas$n),B=B)
    }
    return(out)
}
