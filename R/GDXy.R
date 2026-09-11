# A function to perform Grad. Desc. that takes X and y, instead of XX and Xy

GDXy=function(X,y,centerX=TRUE,scaleX=FALSE,earlyStop=FALSE,pctChangeRSS=0.02,...){

    X=scale(X,center=centerX,scale=scaleX)
    y=scale(y,center=TRUE)

    XX=crossprod(X)
    Xy=crossprod(X,y)
    if(earlyStop){
        B=GD.Auto.RSS(XX=XX,Xy=Xy,pctChangeRSS=pctChangeRSS,...)
    }else{
        B=GD(XX=XX,Xy=Xy,...)
    }
    return(B)
}

#GD.Auto.RSS<- function(XX, Xy, b=NULL, maxIter=10, learningRate=1/50, lambda=0, verbose=TRUE,pctChangeRSS=learningRate*2)
