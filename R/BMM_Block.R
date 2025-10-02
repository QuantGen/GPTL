BMM_Block=function(XX, Xy, B, my, vy, n, nIter=1200, burnIn=200, thin=5, R2=0.25,
	        nComp=matrix(ncol(B)), K=1/nComp, df0.E=5, S0.E=vy*(1-R2)*df0.E, df0.b=rep(10,nComp), 
	        priorProb=rep(1/nComp,nComp), priorCounts=rep(2*nComp,nComp), verbose=TRUE) {
  
}

get_block_ids = function(M) {
  p=M@p
  i=M@i
  n=ncol(M)

  MIN=integer(n)
  MAX=integer(n)

  # For each column j, get min/max row index of its nonzeros
  for (j in 1:n) {
    idx=(p[j] + 1):p[j + 1]
    if (length(idx) == 0) {
      MIN[j]=j
      MAX[j]=j
    } else {
      rows=i[idx] + 1

      MIN[j]=rows[1]
      MAX[j]=rows[length(rows)]
    }
  }

  blk=integer(n)
  b=1
  blk[1]=b
  for (j in 2:n) {
    if (MIN[j] > MAX[j - 1]) {
      b=b + 1
    }
    blk[j]=b
  }
  blk
}
