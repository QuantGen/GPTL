
#include "sample_multinomial.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>


SEXP rMultinomial(SEXP PROB, SEXP n, SEXP p) {
    SEXP samples;
    SEXP prob;
    SEXP cumProb;
    
    int i, j, nSamples, nCat;
    double sumProb,  u;
    double *p_PROB, *p_prob, *p_cumProb;
    int *p_samples;

    nCat = Rf_nrows(PROB);
    nSamples=INTEGER_VALUE(n);

    PROTECT(PROB=AS_NUMERIC(PROB));
    p_PROB=NUMERIC_POINTER(PROB);
        
    PROTECT(prob=NEW_NUMERIC(nCat));
    p_prob=NUMERIC_POINTER(prob);

    PROTECT(cumProb=NEW_NUMERIC(nCat));
    p_cumProb=NUMERIC_POINTER(cumProb);

    PROTECT(samples = NEW_INTEGER(nSamples)); 
    p_samples = INTEGER_POINTER(samples); 

    GetRNGstate();

    //Rprintf("nSamples=%d\n",nSamples);
    
    for (i = 0; i < nSamples; i++){ 

    //Rprintf("Iter=%d\n",i+1);
      sumProb=0;

      // assinging probs and additn to then normalize 
      for(j=0; j<nCat; j++){
         p_prob[j]=p_PROB[(i*nCat)+j];
         sumProb+=p_prob[j];
      }

      // normalizing
      p_prob[0]=p_prob[0]/sumProb;
      p_cumProb[0]=p_prob[0];

      // CDF
      for(j=1; j<nCat; j++){
         p_prob[j]=p_prob[j]/sumProb;
         p_cumProb[j]=p_cumProb[j-1]+p_prob[j];
             
        //Rprintf("Prob=%f\n",p_prob[j]);
        //Rprintf("CumProb=%f\n",p_cumProb[j]);
        
      }

     // sampling
     u=(double)rand() / RAND_MAX; 
     //Rprintf("U=%f\n",u);
     //for(j =0; j<nCat; j++){
       //Rprintf("Prob=%f\n",p_cumProb[j]);
     //}
        
     for (j = 0; j < nCat; j++) { 
         if(u<p_cumProb[j]) break;
     }
     p_samples[i]=j+1;
    }
    PutRNGstate();

    UNPROTECT(4); 
    return samples; 
}
