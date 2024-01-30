#define R_NO_REMAP

#include <Rinternals.h>

SEXP ElasticNet(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP lambda, SEXP b0, SEXP alpha);
