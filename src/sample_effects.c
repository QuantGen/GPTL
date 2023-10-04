#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>
#include "sample_effects.h"

 /* Future C code
	 for(j in 1:p){
	
		offset=sum(C[j,-j]*b[-j])
		rhs_j=(rhs[j]-offset)/varE
		rhs_j=rhs_j+B0[j,d[j]]/varB[d[j]]
 
		lhs_j=C[j,j]/varE+1/varB[d[j]]	
		sol=rhs_j/lhs_j
		b[j]=rnorm(n=1,mean=sol,sd=sqrt(1/lhs_j))
	 }
*/


SEXP sample_effects(SEXP C, SEXP rhs, SEXP b, SEXP d, SEXP B0,SEXP varE, SEXP varB) {
    
    GetRNGstate();
  
    int p = Rf_ncols(C);
 
    PROTECT(rhs=AS_NUMERIC(rhs));
    double *prhs=NUMERIC_POINTER(rhs); 

    PROTECT(C=AS_NUMERIC(C));
    double *pC=NUMERIC_POINTER(C); 
  
    PROTECT(b=AS_NUMERIC(b));
    double *pb=NUMERIC_POINTER(b); 

    PROTECT(varB=AS_NUMERIC(varB));
    double *pvarB=NUMERIC_POINTER(varB);

    PROTECT(d=AS_INTEGER(d));
    int *pd=INTEGER_POINTER(d);

    PROTECT(B0=AS_NUMERIC(B0));
    double *pB0=NUMERIC_POINTER(B0); 
  
    double VARE;
    VARE=NUMERIC_VALUE(varE);
  
    double Cjj;
    double offset;
    int inc=1;
    double rhs_offset;
    double lhs;
    double sol;
    //double z;
     
    for (int j = 0; j < p; j++) { // loop effects
       Cjj = pC[j * (p + 1)];
       offset = F77_NAME(ddot)(&p, pC+j*p, &inc, pb, &inc);
       offset -= Cjj * pb[j];
       rhs_offset = prhs[j] - offset;
       rhs_offset=(rhs_offset/VARE)+pB0[(p*pd[j] + j)]/pvarB[pd[j]];
       lhs=Cjj/VARE+1/pvarB[pd[j]];
       sol=rhs_offset/lhs;
       z=norm_rand()*sqrt(1/lhs);
       Rprintf("Z=%d\n",z);
       pb[j]=sol//+z;
    }
    // Creating a list to return results

    PutRNGstate();
      
    UNPROTECT(6); // b, d
    return b;
}
