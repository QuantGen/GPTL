sample_effects=function(C,rhs,b,d,B0,varE,varB){  
     d<- d - 1L # for the 0-based index
     ans <- .Call("sample_effects", C,rhs,b,d,B0,varE,varB)
     return(ans)
}

sample_effects_new=function(C,rhs,b,d,B0,varE,varB,RSS){  
     d<- d - 1L # for the 0-based index
     ans <- .Call("sample_effects_new", C,rhs,b,d,B0,varE,varB,RSS)
     return(ans)
}

sample_effects_new_sparse=function(C,rhs,b,d,B0,varE,varB,RSS)
{
	d<- d - 1L # for the 0-based index
	diagC<- Matrix::diag(C)
	ans <- .Call("sample_effects_new_sparse",C@x,C@p,C@i,diagC,ncol(C),rhs,b,d,B0,varE,varB,RSS)
	return(ans)
}
