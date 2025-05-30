#include "GRAD_DESC.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>

SEXP GRAD_DESC(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP learning_rate) {

    int inc=1, j, p, niter;
    double  Cjj, Cjb, gradient, LR;
    double *pC, *prhs, *pb;
    
    p=INTEGER_VALUE(nCol);
    niter=INTEGER_VALUE(nIter);
    
    LR =NUMERIC_VALUE(learning_rate);
    
    PROTECT(C=AS_NUMERIC(C));
    pC=NUMERIC_POINTER(C);

    PROTECT(rhs=AS_NUMERIC(rhs));
    prhs=NUMERIC_POINTER(rhs);

    PROTECT(b=AS_NUMERIC(b));
    pb=NUMERIC_POINTER(b);

    inc=1;
    
    for (int iter = 0; iter < niter; iter++) { 
        for (j = 0; j < p; j++) { // loop over predictors
            Cjj = pC[j * (p + 1)]; // I was using this for LR but it gave me numerical problems
            Cjb=F77_NAME(ddot)(&p, pC+j*p, &inc, pb, &inc); //C[,j]'b
            gradient= - (prhs[j]-Cjb);   
            // dividing by Cjj (the 2nd derivative) makes step smaller when the function has more curvature
            pb[j]-=LR*gradient ; 
        }
    }
    
    UNPROTECT(3); 
    return b;
   
}

/*
C: non zero values in C=X'X, 
   @x in object of class "dgCMatrix" in package  "Matrix" (Compressed Sparse Column)
cOffset: column offset, cumulative number of non zero elements in X, 
         @p in object of class "dgCMatrix" in package "Matrix" (Compressed Sparse Column)
rIndex: row index of each element in C, 
         @i in object of class "dgCMatrix" 
         in package "Matrix" (Compressed Sparse Column)
*/

SEXP GRAD_DESC_sparse(SEXP C, SEXP cOffset, SEXP rIndex, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP learning_rate) 
{

    int inc=1, j, p, niter;
    double  Cjj, Cjb, gradient, LR;
    double *pC, *prhs, *pb;
    
    int *pcOffset, *prIndex;
    int low, up;
    int m;
    
    p=INTEGER_VALUE(nCol);
    niter=INTEGER_VALUE(nIter);
    
    LR=NUMERIC_VALUE(learning_rate);
    
    PROTECT(C=AS_NUMERIC(C));
    pC=NUMERIC_POINTER(C);
    
    PROTECT(cOffset=AS_INTEGER(cOffset));
    pcOffset=INTEGER_POINTER(cOffset);
    
    PROTECT(rIndex=AS_INTEGER(rIndex));
    prIndex=INTEGER_POINTER(rIndex);

    PROTECT(rhs=AS_NUMERIC(rhs));
    prhs=NUMERIC_POINTER(rhs);

    PROTECT(b=AS_NUMERIC(b));
    pb=NUMERIC_POINTER(b);
    
    inc=1;
    
    for (int iter = 0; iter < niter; iter++) 
    { 
        // loop over predictors
        for (j = 0; j < p; j++) 
        { 
            //Cjb=F77_NAME(ddot)(&p, pC+j*p, &inc, pb, &inc); //C[,j]'b
              
            /*index the elements of C using the column offset*/
	  		low=pcOffset[j];
          	up=pcOffset[j+1]-1;
          	
          	Cjb=0;
          	
          	/*equivalent to ddot*/
          	for(m=low;m<=up;m++)
          	{
              Cjb=Cjb+pC[m]*pb[prIndex[m]];
	  		}
            
            gradient= - (prhs[j]-Cjb);   
            pb[j]-=LR*gradient ; 
        }
    }
    
    UNPROTECT(5); 
    return b;
   
}

