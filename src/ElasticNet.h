#define R_NO_REMAP

#include <Rinternals.h>

SEXP ElasticNet(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP lambda1, SEXP lambda2, SEXP b0);
SEXP ElasticNet_sparse(SEXP C, SEXP cOffset, SEXP rIndex, SEXP diagC, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP lambda1, SEXP lambda2, SEXP b0);
