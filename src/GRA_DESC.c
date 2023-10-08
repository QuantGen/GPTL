#include "GRAD_DESC.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>

SEXP GRAD_DESC(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP learning_rate) {

    int inc=1, j, p, niter;
    double Cjj, Cjb, gradient, LR;
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
