## This function simply pools summary statistics using weights, the output can be use as inputs in both BLRCross(), RR(), or GD()


PSS=function(SSList,weights){
  # SSList should be a list of sufficient statistics with the ith element having [[i]]$XX, [[i]]$Xy, [[i]]$my, [[i]]$vy, [[i]]$n
  # weights: a vecctor of weights of the same length as SSList
  # Value: a list with the pooled sufficient staitistics XX, Xy, my, vy, and the estimated effective sample size (see https://www.jepusto.com/effective-sample-size-aggregation/)
  weights=weights/sum(weightS)
  if(any(weights<=0)){ stop('All weights must be greater than 0') }

  p1=ncol(SSList[[1]]$XX)
  p2=length(SSList[[1]]$Xy

  #....  

  #return(list(XX=XX,Xy=Xy,my=my,vy=vy,ESS=ESS)) # ESS is the effective sample size
}
