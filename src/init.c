#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ElasticNet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ElasticNet_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GRAD_DESC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GRAD_DESC_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LASSO_CD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RIDGE_CD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fitLSYS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rMultinomial(SEXP, SEXP, SEXP);
extern SEXP sample_effects(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sample_effects_new(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sample_effects_new_sparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ElasticNet",                (DL_FUNC) &ElasticNet,                 8},
    {"ElasticNet_sparse",         (DL_FUNC) &ElasticNet_sparse,         11},
    {"GRAD_DESC",                 (DL_FUNC) &GRAD_DESC,                  6},
    {"GRAD_DESC_sparse",          (DL_FUNC) &GRAD_DESC_sparse,           8},
    {"LASSO_CD",                  (DL_FUNC) &LASSO_CD,                   7},
    {"RIDGE_CD",                  (DL_FUNC) &RIDGE_CD,                   7},
    {"fitLSYS",                   (DL_FUNC) &fitLSYS,                    7},
    {"rMultinomial",              (DL_FUNC) &rMultinomial,               3},
    {"sample_effects",            (DL_FUNC) &sample_effects,             7},
    {"sample_effects_new",        (DL_FUNC) &sample_effects_new,         8},
    {"sample_effects_new_sparse", (DL_FUNC) &sample_effects_new_sparse, 12},
    {NULL, NULL, 0}
};

void R_init_GPTL(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

