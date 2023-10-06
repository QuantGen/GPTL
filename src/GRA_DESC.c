#include "GRAD_DESC.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>

SEXP GRAD_DESC(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP learning_rate) {

    int inc=1;
    int j;
    int p; 
    int niter;
    double offset;
    double Cjj;
    double gradient;
    double sol;
    double LR;
    double *pC;
    double *prhs;
    double *pb;
    
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
    
    for (iter = 0; iter < niter; iter++) {       
        for (j = 0; j < p; j++) { // loop over predictors
            Cjj = pC[j * (p + 1)]; 
            offset=0;//F77_NAME(ddot)(&p, pC+j*p, &inc, pb, &inc);
            gradient=2*Ckk*pb[j] - (prhs[j]-offset);
            sol=pb[j]-LR*gradient ;
            pb[k] = sol;
        }
    }
    
    UNPROTECT(3); 
    return b;
   
}
