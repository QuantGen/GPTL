#define R_NO_REMAP

#include <Rinternals.h>

SEXP sample_effects_new(SEXP C, SEXP rhs, SEXP b, SEXP d, SEXP B0,SEXP varE, SEXP varB, SEXP RSS);
SEXP sample_effects_new_sparse(SEXP C, SEXP cOffset, SEXP rIndex, SEXP diagC, SEXP nCol, SEXP rhs, SEXP b, SEXP d, SEXP B0, SEXP varE, SEXP varB, SEXP RSS);
