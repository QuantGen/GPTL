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
