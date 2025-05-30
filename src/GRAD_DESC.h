#define R_NO_REMAP

#include <Rinternals.h>

SEXP GRAD_DESC(SEXP C, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP learning_rate);
SEXP GRAD_DESC_sparse(SEXP C, SEXP cOffset, SEXP rIndex, SEXP rhs, SEXP b, SEXP nCol, SEXP nIter, SEXP learning_rate);
