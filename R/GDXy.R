# A function to perform Grad. Desc. that takes X and y, instead of XX and Xy

GDXy=function(X,y,centerX=TRUE,scaleX=FALSE,sortGradient=FALSE,...){

    X=scale(X,center=centerX,scale=scaleX)
    y=scale(y,center=TRUE)

    XX=crossprod(X)
    Xy=crossprod(X,y)

    B=GD(XX=XX,Xy=Xy,sortGradient=sortGradient,...)

    return(B)
}
