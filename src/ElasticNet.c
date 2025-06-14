#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>

SEXP ElasticNet(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP lambda1, SEXP lambda2, SEXP b0) {

    int inc=1, j, p, niter;
    double  Cjj, Cjb, offset, r_tilde, bOLS, Lambda1, Lambda2, sign_bOLS;
    double *pC, *prhs, *pb, *pb0;
    
    p=INTEGER_VALUE(nCol);
    niter=INTEGER_VALUE(nIter);
    
    Lambda1 =NUMERIC_VALUE(lambda1);
    Lambda2 =NUMERIC_VALUE(lambda2);
    
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
            r_tilde=prhs[j]-offset;
            bOLS=r_tilde/Cjj;
            
            if (bOLS-pb0[j] < -Lambda1/Cjj) {
                pb[j]=(r_tilde+Lambda2*pb0[j]+Lambda1)/(Cjj+Lambda2);
            } else if (bOLS-pb0[j] > Lambda1/Cjj) {
                pb[j]=(r_tilde+Lambda2*pb0[j]-Lambda1)/(Cjj+Lambda2);
            } else {
                pb[j]=pb0[j];
            }
        }
    }
    
    UNPROTECT(4); 
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

SEXP ElasticNet_sparse(SEXP C, SEXP cOffset, SEXP rIndex, SEXP diagC, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP lambda1, SEXP lambda2, SEXP b0) {

    int inc=1, j, p, niter;
    double  Cjj, Cjb, offset, r_tilde, bOLS, Lambda1, Lambda2, sign_bOLS;
    double *pC, *prhs, *pb, *pb0;
    
    double *pdiagC;
    
    int *pcOffset, *prIndex;
    int low, up;
    int m;
    
    
    p=INTEGER_VALUE(nCol);
    niter=INTEGER_VALUE(nIter);
    
    Lambda1 =NUMERIC_VALUE(lambda1);
    Lambda2 =NUMERIC_VALUE(lambda2);
    
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

    PROTECT(b0=AS_NUMERIC(b0));
    pb0=NUMERIC_POINTER(b0);
    
    PROTECT(diagC=AS_NUMERIC(diagC));
    pdiagC=NUMERIC_POINTER(diagC);

    inc=1;
    
    for (int iter = 0; iter < niter; iter++) { 
        // loop over predictors
        for (j = 0; j < p; j++) 
        { 
            //Cjj = pC[j * (p + 1)]; // C[j,j]
            //Cjb=F77_NAME(ddot)(&p, pC+j*p, &inc, pb, &inc); //C[,j]'b
            
            Cjj=pdiagC[j];
            
            /*index the elements of C using the column offset*/
	  		low=pcOffset[j];
          	up=pcOffset[j+1]-1;
          	
          	Cjb=0;
          	
          	/*equivalent to ddot*/
          	for(m=low;m<=up;m++)
          	{
              Cjb=Cjb+pC[m]*pb[prIndex[m]];
	  		}

            offset=Cjb-Cjj*pb[j];
            r_tilde=prhs[j]-offset;
            bOLS=r_tilde/Cjj;
            
            if (bOLS-pb0[j] < -Lambda1/Cjj) {
                pb[j]=(r_tilde+Lambda2*pb0[j]+Lambda1)/(Cjj+Lambda2);
            } else if (bOLS-pb0[j] > Lambda1/Cjj) {
                pb[j]=(r_tilde+Lambda2*pb0[j]-Lambda1)/(Cjj+Lambda2);
            } else {
                pb[j]=pb0[j];
            }
        }
    }
    
    UNPROTECT(7); 
    return b;
   
}
