// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// adeba_find_constants
LogicalVector adeba_find_constants(NumericMatrix x);
RcppExport SEXP adeba_adeba_find_constants(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(adeba_find_constants(x));
    return rcpp_result_gen;
END_RCPP
}
// adeba_is_constant
LogicalVector adeba_is_constant(NumericVector x);
RcppExport SEXP adeba_adeba_is_constant(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(adeba_is_constant(x));
    return rcpp_result_gen;
END_RCPP
}
// get_bandwidths
NumericMatrix get_bandwidths(NumericVector pilot, NumericVector alpha, NumericVector beta);
RcppExport SEXP adeba_get_bandwidths(SEXP pilotSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pilot(pilotSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_bandwidths(pilot, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"adeba_adeba_find_constants", (DL_FUNC) &adeba_adeba_find_constants, 1},
    {"adeba_adeba_is_constant", (DL_FUNC) &adeba_adeba_is_constant, 1},
    {"adeba_get_bandwidths", (DL_FUNC) &adeba_get_bandwidths, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_adeba(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}