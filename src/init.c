//
//  init.c
//  
//
//  Created by XuZekun on 3/3/17.
//
//

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP ziphsmm_dzip(SEXP,SEXP,SEXP,SEXP);
extern SEXP ziphsmm_rzip(SEXP,SEXP,SEXP);
extern SEXP ziphsmm_zipnegloglik_nocov_cont(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP ziphsmm_grad_zipnegloglik_nocov_cont(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP ziphsmm_zipnegloglik_cov_cont(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP ziphsmm_grad_zipnegloglik_cov_cont(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
extern SEXP ziphsmm_retrieve_nocov_cont(SEXP,SEXP);
extern SEXP ziphsmm_retrieve_cov_cont(SEXP,SEXP,SEXP);
extern SEXP ziphsmm_convolution(SEXP, SEXP);

//extern SEXP RcppArmadillo_armadillo_version(SEXP);
//extern SEXP RcppArmadillo_fastLm(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ziphsmm_dzip", (DL_FUNC) &ziphsmm_dzip, 4}, //number of parms
    {"ziphsmm_rzip", (DL_FUNC) &ziphsmm_rzip, 3},
    {"ziphsmm_zipnegloglik_nocov_cont", (DL_FUNC) &ziphsmm_zipnegloglik_nocov_cont, 6},
    {"ziphsmm_grad_zipnegloglik_nocov_cont", (DL_FUNC) &ziphsmm_grad_zipnegloglik_nocov_cont, 6},
    {"ziphsmm_zipnegloglik_cov_cont", (DL_FUNC) &ziphsmm_zipnegloglik_cov_cont, 7},
    {"ziphsmm_grad_zipnegloglik_cov_cont", (DL_FUNC) &ziphsmm_grad_zipnegloglik_cov_cont, 7},
    {"ziphsmm_retrieve_nocov_cont", (DL_FUNC) &ziphsmm_retrieve_nocov_cont, 2},
    {"ziphsmm_retrieve_cov_cont", (DL_FUNC) &ziphsmm_retrieve_cov_cont, 3},
    {"ziphsmm_convolution", (DL_FUNC) &ziphsmm_convolution, 2},
    {NULL, NULL, 0}
};

void R_init_ziphsmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
