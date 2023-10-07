#include "sample_multinomial.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>


SEXP sample_multinomial(SEXP PROB, SEXP n, SEXP p, SEXP prob) {
    SEXP samples;
    
    int i, j, nSamples, nCat;
    double sumProb, cumProb, u;
    double *p_PROB, *p_prob, *p_cumProb, *pAns;
    int *p_samples;
    
    nSamples=INTEGER_VALUE(n);
    nCat=INTEGER_VALUE(p);

    PROTECT(PROB=AS_NUMERIC(PROB));
    p_PROB=NUMERIC_POINTER(PROB);
        
    PROTECT(prob=AS_NUMERIC(prob));
    p_prob=NUMERIC_POINTER(prob);

    PROTECT(samples = NEW_INTEGER(nCat)); 
    pSamples = INTEGER_POINTER(samples); 

    GetRNGstate();
    
    for (i = 0; i < nSamples; i++){ 
      sumProb=0;

      // adding column values 
      for(j=1; j<=nCat; j++){
         sumProb+=p_PROB[(i*nCat)+j];
      }

      // normalizing
      p_prob[0]=p_prob[0]/sumProb;
      p_cumProb[0]=p_prob[0];

      // CDF
      for(j=1; j<nCat; j++){
         p_prob[j]=p_prob[j]/sumProb;
         p_cumProb[j]=p_cumProb[j-1]+p_prob[j];
      }

     // sampling
     u=rand();
     for (j = 0; j < nCat; j++) { 
         if(u<p_cumProb[j]) break;
     }
     ans[i]=j;
    }
    PutRNGstate();

    UNPROTECT(2); 
    return ans; 
}
