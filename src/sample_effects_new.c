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


SEXP sample_effects_new(SEXP C, SEXP rhs, SEXP b, SEXP d, SEXP B0,SEXP varE, SEXP varB, SEXP RSS) {
    
    GetRNGstate();
  
    int p = Rf_ncols(C);
    int inc=1;

    double Cjj, offset, rhs_offset, lhs, sol, z, old_beta;
    double *pRSS;

    SEXP list;
	
    PROTECT(rhs=AS_NUMERIC(rhs));
    double *prhs=NUMERIC_POINTER(rhs); 

    PROTECT(C=AS_NUMERIC(C));
    double *pC=NUMERIC_POINTER(C); 

    PROTECT(b=AS_NUMERIC(Rf_duplicate(b)));
    double *pb=NUMERIC_POINTER(b); 

    PROTECT(varB=AS_NUMERIC(varB));
    double *pvarB=NUMERIC_POINTER(varB);

    PROTECT(d=AS_INTEGER(d));
    int *pd=INTEGER_POINTER(d);

    PROTECT(B0=AS_NUMERIC(B0));
    double *pB0=NUMERIC_POINTER(B0); 
  
    double VARE;
    VARE=NUMERIC_VALUE(varE);

    
    PROTECT(RSS=AS_NUMERIC(RSS));
    pRSS=NUMERIC_POINTER(RSS);

    for (int j = 0; j < p; j++) { // loop effects
       Cjj = pC[j * (p + 1)];
       offset = F77_NAME(ddot)(&p, pC+j*p, &inc, pb, &inc);
       offset -= Cjj * pb[j];
       rhs_offset = prhs[j] - offset;
       rhs_offset=(rhs_offset/VARE)+pB0[(p*pd[j] + j)]/pvarB[pd[j]];
       lhs=(Cjj/VARE)+(1/pvarB[pd[j]]);
       sol=rhs_offset/lhs;
       z=norm_rand()*sqrt(1/lhs);
       old_beta=pb[j];
       pb[j]=sol+z;
       //Rprintf("Z=%f\n",old_beta);
       //Rprintf("Z=%f\n",pb[j]);
	    
       pRSS[0]+=(pow(pb[j],2) - pow(old_beta,2))*Cjj  -2*(pb[j]-old_beta)*(prhs[j]-offset);

      //pRSS[0]+=(pow(pbeta[j_global],2) - pow(old_beta,2))*Cjj  -2*(pbeta[j_global]-old_beta)*(rhs-offset);
  
	    
    }

   // Creating a list with 1 vector elements:
      PROTECT(list = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(list, 0, b);
      SET_VECTOR_ELT(list, 1, RSS);
       
      UNPROTECT(8);
      
      PutRNGstate();
	
      return(list);

}
