nextLocus=function(Z1,Dprime,shape1=shape1,shape2=shape2){
    # Simulates a 2nd allele on a gamete given the first allele (Z1) and Dprime
 
    pA=mean(Z1) ; n=length(Z1)
    pB=seq(from=0.01,to=0.99,length=100)
    x=pA*(1-pB)
    y=(1-pA)*pB
    Dmax=ifelse(x<y,x,y)
    D=Dprime*Dmax
    ready=FALSE
    # D=pAB-pApB=> pAB=D+pA*Pb
    pAB=D+pA*pB
    # Using pA=pAB+pAb, we get pAb=pA-pAB
    pAb=pA-pAB
    #Using pB=pAB+paB, we get paB=pB-pAB
    paB=pB-pAB
    pab=1-pAB-pAb-paB
    # Sampling pB constrained by pA and D
    TMP=cbind(pAB,pAb,paB,pab)
    cond=rowSums(TMP>0.005)==4
    if(sum(cond)>3){
        pB=paB+pAB
        tmp=sample(which(cond),size=1,prob=dbeta(pB[cond],shape1=shape1,shape2=shape2))
        pB=pB[tmp];  pAB=pAB[tmp]; pAb=pAb[tmp];
        paB=paB[tmp];  pab=pab[tmp]; D=D[tmp]; Dmax=Dmax[tmp]
              
        # Sampling genotype for locus 2
        x=rbinom(n=n,size=1,prob=pAB/(pA))
        # sampling gametes given Z1==0
        y=rbinom(n=n,size=1,prob=paB/(1-pA))
        Z2=ifelse(Z1==1,x,y)
    }else{
        Z2=rbinom(size=2,n=n,prob=mean(Z1)/2)
    }
    return(Z2)
}

sampleHaplotype<-function(n,nLoci,Dprime,p,shape1=shape1,shape2=shape){
    # n= number of genotypes, nLoci=number of loci in the haplotype
    # Dprime, p=allele frequency at the first locus.
    Z=as.matrix(rbinom(size=1,n=n,prob=p))
    if(nLoci>1){
        for(i in 2:nLoci){
            Z=cbind(Z,nextLocus(Z[,i-1],Dprime=Dprime[i],shape1=shape1,shape2=shape2))   
        }
    }
    return(Z)
}
 
sampleGenomes=function(n,nLoci,Dprime,shape1=2,shape2=8){
    p=rbeta(shape1=shape1,shape2=shape2,n=1)
    Z1=sampleHaplotype(n=n,nLoci=nLoci,Dprime=Dprime,p=p,shape1=shape1,shape2=shape2)
    Z2=sampleHaplotype(n=n,nLoci=nLoci,Dprime=Dprime,p=p,shape1=shape1,shape2=shape2)
    X=Z1+Z2;
    return(X)
}

nextBlock=function(x,pos,maxGaps){
    p=length(x)
    # Finding ther right-break
    if(pos+maxGaps>p){
        rightBreak=p
    }else{
        nGaps=0
        for(i in (pos+1):p){
            if(x[i]){
                nGaps=0
            }else{
            nGaps=nGaps+1
            }
            if(nGaps>maxGaps){
                break()
            }
        }
        rightBreak=i-nGaps
    }
    # Finding the left break
    if(pos-maxGaps<1){
        leftBreak=1
    }else{
        nGaps=0
        for(i in (pos-1):1){
            if(x[i]){
                nGaps=0
            }else{
                nGaps=nGaps+1
            }
            if(nGaps>maxGaps){
                break()
            }  
        }
        leftBreak=i+nGaps
    }
    return(c(left=leftBreak,right=rightBreak))
}
 
findBlocks=function(R2,threshold,maxGap,h=1/2){
    A=I(R2>threshold)
    diag(A)=FALSE
    p=nrow(R2)
    K=exp(-h*as.matrix(dist(1:p)))
    blocks=rep(NA,p)
    ready=FALSE
    nextBlockNum=1
    while(!ready){
        nextSNP=which.max(colSums(A*K))[1]
        tmp=nextBlock(A[,nextSNP],nextSNP,maxGap)
        index=tmp[1]:tmp[2]
        A[index,]=FALSE
        A[,index]=FALSE
        blocks[index]=nextBlockNum
        nextBlockNum=nextBlockNum+1
        ready=!any(A)
    }
 
    n=sum(is.na(blocks))
    blocks[is.na(blocks)]=nextBlockNum:(nextBlockNum+n-1)
    blocks=as.integer(factor(blocks,levels=unique(blocks)))
    return(blocks)
}
 
mergeBlocks=function(R2,blocks,minSize){
    counts=table(blocks)
    ready<-all( counts>=minSize)
    nBlocks=length(unique(blocks))
    counter=0
    while(!ready){
        counter=counter+1
        block=which.min(counts)
       
        if(block==1){
            blocks[blocks==block]=2
        }
        if(block==nBlocks){
            blocks[blocks==block]=nBlocks-1
        }
        if((block!=1)&(block!=nBlocks)){
            left_merge=sum(R2[blocks==block,blocks==(block-1)])
            right_merge=sum(R2[blocks==block,blocks==(block+1)])
            if(left_merge>right_merge){
                blocks[blocks==block]=block-1
            }else{
                blocks[blocks==block]=block+1
            }
        }
        blocks=as.integer(as.factor(blocks))
        counts=table(blocks)
        ready<-all( counts>=minSize)
        nBlocks=length(unique(blocks))
    }
    return(blocks)
}
