#include "LASSO_CD.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>

SEXP LASSO_CD(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP lambda, SEXP b0) {

    int inc=1, j, p, niter;
    double  Cjj, Cjb, offset, bOLS, Lambda, sign_bOLS;
    double *pC, *prhs, *pb, *pb0;
    
    p=INTEGER_VALUE(nCol);
    niter=INTEGER_VALUE(nIter);
    
    Lambda =NUMERIC_VALUE(lambda);
    
    PROTECT(C=AS_NUMERIC(C));
    pC=NUMERIC_POINTER(C);

    PROTECT(rhs=AS_NUMERIC(rhs));
    prhs=NUMERIC_POINTER(rhs);

    PROTECT(b=AS_NUMERIC(b));
    pb=NUMERIC_POINTER(b);

    PROTECT(b0=AS_NUMERIC(b0));
    pb0=NUMERIC_POINTER(b0);

    inc=1;
    
    for (int iter = 0; iter < niter; iter++) { 
        for (j = 0; j < p; j++) { // loop over predictors
            Cjj = pC[j * (p + 1)]; // C[j,j]
            Cjb=F77_NAME(ddot)(&p, pC+j*p, &inc, pb, &inc); //C[,j]'b

            offset=Cjb-Cjj*pb[j];
            bOLS=(prhs[j]-offset)/Cjj;
            // sign function
            if (bOLS-pb0[j] < 0) {
                sign_bOLS=-1;
            } else if (bOLS-pb0[j] == 0) {
                sign_bOLS=0;
            } else {
                sign_bOLS=1;
            }
            
            if (abs(bOLS - pb0[j]) > Lambda) {
                pb[j]=bOLS-sign_bOLS*Lambda;
            } else {
                pb[j]=pb0[j];
            }
        }
    }
    
    UNPROTECT(4); 
    return b;
   
}
