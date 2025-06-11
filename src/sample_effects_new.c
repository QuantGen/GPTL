#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>

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
       rhs_offset=(rhs_offset/VARE)+pB0[(p*pd[j] + j)]/pvarB[j];
       lhs=(Cjj/VARE)+(1/pvarB[j]);
       sol=rhs_offset/lhs;
       z=norm_rand()*sqrt(1/lhs);
       old_beta=pb[j];
       pb[j]=sol+z;
       //Rprintf("Z=%f\n",pb[j]);
	    
       pRSS[0]+=(pow(pb[j],2) - pow(old_beta,2))*Cjj  -2*(pb[j]-old_beta)*(prhs[j]-offset);
	    
    }

   // Creating a list with 1 vector elements:
      PROTECT(list = allocVector(VECSXP, 2));
      SET_VECTOR_ELT(list, 0, b);
      SET_VECTOR_ELT(list, 1, RSS);
       
      UNPROTECT(8);
      
      PutRNGstate();
	
      return(list);

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


SEXP sample_effects_new_sparse(SEXP C, SEXP cOffset, SEXP rIndex, SEXP diagC, SEXP nCol, SEXP rhs, SEXP b, SEXP d, SEXP B0, SEXP varE, SEXP varB, SEXP RSS) 
{
    
    GetRNGstate();
    
    int p;
    
    double Cjj, offset, rhs_offset, lhs, sol, z, old_beta;
    double *pRSS;
    
    int *pcOffset, *prIndex;
    int low, up;
    int m;
    
    double *pdiagC;
    

    SEXP list;
    
        
    p=INTEGER_VALUE(nCol);
    
    PROTECT(RSS=AS_NUMERIC(RSS));
    pRSS=NUMERIC_POINTER(RSS);
	
	
    PROTECT(rhs=AS_NUMERIC(rhs));
    double *prhs=NUMERIC_POINTER(rhs); 

    
    PROTECT(C=AS_NUMERIC(C));
    double *pC=NUMERIC_POINTER(C); 
    
    
    PROTECT(cOffset=AS_INTEGER(cOffset));
    pcOffset=INTEGER_POINTER(cOffset);
    
    
    PROTECT(rIndex=AS_INTEGER(rIndex));
    prIndex=INTEGER_POINTER(rIndex);
    
    
    PROTECT(diagC=AS_NUMERIC(diagC));
    pdiagC=NUMERIC_POINTER(diagC);
    
    
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
    
    //Rprintf("%d\n",p);
    
    for (int j = 0; j < p; j++) { // loop effects
       
       
       Cjj=pdiagC[j];
       
       low=pcOffset[j];
       up=pcOffset[j+1]-1;
          	
       offset=0;
       
       for(m=low;m<=up;m++)
       {
        	offset=offset+pC[m]*pb[prIndex[m]];
	   }
       
       offset -= Cjj * pb[j];
       rhs_offset = prhs[j] - offset;
       rhs_offset=(rhs_offset/VARE)+pB0[(p*pd[j] + j)]/pvarB[j];
       lhs=(Cjj/VARE)+(1/pvarB[j]);
       sol=rhs_offset/lhs;
       z=norm_rand()*sqrt(1/lhs);
       old_beta=pb[j];
       pb[j]=sol+z;
       
       //Rprintf("%f\n",pb[j]); 
       
       pRSS[0]+=(pow(pb[j],2) - pow(old_beta,2))*Cjj  -2*(pb[j]-old_beta)*(prhs[j]-offset);
	    
    }
    
    //Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, b);
    SET_VECTOR_ELT(list, 1, RSS);
       
    UNPROTECT(11);
    
    
    PutRNGstate();
	
    return(list);

}

