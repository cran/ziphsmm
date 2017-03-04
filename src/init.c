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



//extern SEXP RcppArmadillo_armadillo_version(SEXP);
//extern SEXP RcppArmadillo_fastLm(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
        //{"RcppArmadillo_armadillo_version",         (DL_FUNC) &RcppArmadillo_armadillo_version,         1},
    //{"RcppArmadillo_fastLm",                    (DL_FUNC) &RcppArmadillo_fastLm,                    2},
    {"ziphsmm_dzip", (DL_FUNC) &ziphsmm_dzip, 4}, //number of parms
    {"ziphsmm_rzip", (DL_FUNC) &ziphsmm_rzip, 3},
    {NULL, NULL, 0}
};

void R_init_ziphsmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
