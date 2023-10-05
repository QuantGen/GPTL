#include "GRAD_DESC.h"

SEXP GRAD_DESC(SEXP C, SEXP rhs, SEXP b, SEXP nIter, SEXP learning_rate) {

    int p = Rf_ncols(C);
    int n_iter = Rf_asInteger(nIter);
    double LR =Rf_asReal(learning_rate);
    double *pC = REAL(C);
    double *prhs = REAL(rhs);
    b = PROTECT(Rf_duplicate(b));
    double *pb = REAL(b);
    int inc=1;
    double offset ; 
    double Cjj;
    double gradient;
    double sol;
    
    for (int iter = 0; iter < n_iter; iter++) {       
        for (int j = 0; j < p; j++) { // loop over predictors
            Cjj = pC[j * (p + 1)]; 
            offset=F77_NAME(ddot)(&p, pXX+j_global*p, &inc, pbeta, &inc);
            gradient=2*Ckk*pb[j] - (prhs[j]-offset);
            sol=pb[j]-LR*gradient ;
            pb[k] = sol;
        }
    }

    // Creating a list to return results
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 1));
    SET_VECTOR_ELT(list, 0, b);
    UNPROTECT(2); // b, list
    return list;
   
}
