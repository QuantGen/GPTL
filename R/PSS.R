## This function simply pools summary statistics using weights, the output can be use as inputs in both BLRCross(), RR(), or GD()

if(FALSE){
PSS=function(SSList,weights){
  # SSList should be a list of sufficient statistics with the ith element having [[i]]$XX, [[i]]$Xy, [[i]]$my, [[i]]$vy, [[i]]$n
  # weights: a vecctor of weights of the same length as SSList
  # Value: a list with the pooled sufficient staitistics XX, Xy, my, vy, and the estimated effective sample size (see https://www.jepusto.com/effective-sample-size-aggregation/)

  weights=weights/sum(weights)

  if(length(weights) != length(SSList)){ 
    stop('SSList and weights must have the same length') 
  }
  
  if(any(weights<=0)){ 
    stop('All weights must be greater than 0') 
  }
  
  N_sample=numeric(length(SSList))
  wXX=list()
  wXy=list()
  wmy=numeric(length(SSList))
  wvy=numeric(length(SSList))

  for (i in 1:length(SSList)) {
    N_sample[i]=SSList[[i]]$n
    wXX[[i]]=SSList[[i]]$XX*weights[i]
    wXy[[i]]=SSList[[i]]$Xy*weights[i]
    wmy[i]=SSList[[i]]$my*weights[i]
    wvy[i]=SSList[[i]]$vy*weights[i]
  }

  ESS=sum(N_sample*weights)^2/sum(N_sample*weights^2)
  XX=ESS*Reduce("+", wXX)
  Xy=ESS*Reduce("+", wXy)
  my=sum(wmy)
  vy=sum(wvy)

  return(list(XX=XX,Xy=Xy,my=my,vy=vy,ESS=ESS)) # ESS is the effective sample size
}
}
