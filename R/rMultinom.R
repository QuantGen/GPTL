## A wrapper for a C function that sample multinomials
## PROB is a matrix with the probabilities for different categories (rows) for each sample (columns)
## For example, if you want to sample multinomials with otucomes 1:4 for n individuals PROBS is 4xn
## The values of PROB must be non-negative, but the column values do not need to add up to 1, this normalization
## is done internally.

rMultinom=function(PROB){
   n=ncol(PROB)
   p=nrow(PROB) 
   samples=.Call("rMultinomial",PROB,n,p)
   return(samples)
}

## Example
if(FALSE){
 PROB=cbind(c(2,1.5),c(.1,.8,.1)) # 3 categories, 2 samples
 rMultinom(PROBS)

}
