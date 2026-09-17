corBootstrap=function(x,y,times=100){
  fn=function(x,y,n){
    tmp=sample(1:n,size=n, replace=TRUE)
    cor(x[tmp],y[tmp])
  }
  mean(replicate(100,fn(x,y,n)))
}

