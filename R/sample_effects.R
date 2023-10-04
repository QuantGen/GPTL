sample_effects=function(C,rhs,b,d,B0,varE,varB){  
     d<- d - 1L # for the 0-based index
     ans <- .Call("sample_effects", C,rhs,b,d,B0,varE,varB)
     return(ans)
}
