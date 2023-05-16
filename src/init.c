#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _SIRTGP_fun_dev(SEXP, SEXP);
extern SEXP _SIRTGP_fun_mul(SEXP, SEXP);
extern SEXP _SIRTGP_G_thres_cpp(SEXP, SEXP);
extern SEXP _SIRTGP_gibbs_sample_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_loglike_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_quantile_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_rcpparma_bothproducts(SEXP);
extern SEXP _SIRTGP_rcpparma_hello_world();
extern SEXP _SIRTGP_rcpparma_innerproduct(SEXP);
extern SEXP _SIRTGP_rcpparma_outerproduct(SEXP);
extern SEXP _SIRTGP_rtrunc_Rpkg(SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_e_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_E_hat_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_eta_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_eta_hat_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_ic_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_thres1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_thres2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SIRTGP_sample_Z(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_SIRTGP_fun_dev",               (DL_FUNC) &_SIRTGP_fun_dev,                2},
    {"_SIRTGP_fun_mul",               (DL_FUNC) &_SIRTGP_fun_mul,                2},
    {"_SIRTGP_G_thres_cpp",           (DL_FUNC) &_SIRTGP_G_thres_cpp,            2},
    {"_SIRTGP_gibbs_sample_cpp",      (DL_FUNC) &_SIRTGP_gibbs_sample_cpp,      27},
    {"_SIRTGP_loglike_cpp",           (DL_FUNC) &_SIRTGP_loglike_cpp,            5},
    {"_SIRTGP_quantile_cpp",          (DL_FUNC) &_SIRTGP_quantile_cpp,           4},
    {"_SIRTGP_rcpparma_bothproducts", (DL_FUNC) &_SIRTGP_rcpparma_bothproducts,  1},
    {"_SIRTGP_rcpparma_hello_world",  (DL_FUNC) &_SIRTGP_rcpparma_hello_world,   0},
    {"_SIRTGP_rcpparma_innerproduct", (DL_FUNC) &_SIRTGP_rcpparma_innerproduct,  1},
    {"_SIRTGP_rcpparma_outerproduct", (DL_FUNC) &_SIRTGP_rcpparma_outerproduct,  1},
    {"_SIRTGP_rtrunc_Rpkg",           (DL_FUNC) &_SIRTGP_rtrunc_Rpkg,            4},
    {"_SIRTGP_sample_cpp",            (DL_FUNC) &_SIRTGP_sample_cpp,             4},
    {"_SIRTGP_sample_e_cpp",          (DL_FUNC) &_SIRTGP_sample_e_cpp,          13},
    {"_SIRTGP_sample_E_hat_cpp",      (DL_FUNC) &_SIRTGP_sample_E_hat_cpp,      14},
    {"_SIRTGP_sample_eta_cpp",        (DL_FUNC) &_SIRTGP_sample_eta_cpp,        15},
    {"_SIRTGP_sample_eta_hat_cpp",    (DL_FUNC) &_SIRTGP_sample_eta_hat_cpp,    14},
    {"_SIRTGP_sample_ic_cpp",         (DL_FUNC) &_SIRTGP_sample_ic_cpp,          5},
    {"_SIRTGP_sample_thres1",         (DL_FUNC) &_SIRTGP_sample_thres1,         13},
    {"_SIRTGP_sample_thres2",         (DL_FUNC) &_SIRTGP_sample_thres2,         10},
    {"_SIRTGP_sample_Z",              (DL_FUNC) &_SIRTGP_sample_Z,               4},
    {NULL, NULL, 0}
};

void R_init_SIRTGP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}