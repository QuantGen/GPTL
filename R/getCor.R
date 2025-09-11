getCor=function(XX,Xy,yy,b){
	crossprod(b,Xy)/sqrt(yy*t(b)%*%XX%*%b)
}
