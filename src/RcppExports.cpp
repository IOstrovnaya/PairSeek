// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getpairsOddXOut
NumericMatrix getpairsOddXOut(NumericMatrix MatOut, NumericMatrix MatIn, IntegerVector Yout, int pch2);
RcppExport SEXP _PairSeek_getpairsOddXOut(SEXP MatOutSEXP, SEXP MatInSEXP, SEXP YoutSEXP, SEXP pch2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type MatOut(MatOutSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type MatIn(MatInSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Yout(YoutSEXP);
    Rcpp::traits::input_parameter< int >::type pch2(pch2SEXP);
    rcpp_result_gen = Rcpp::wrap(getpairsOddXOut(MatOut, MatIn, Yout, pch2));
    return rcpp_result_gen;
END_RCPP
}
// getpairsOddXOutRandomOrder
NumericMatrix getpairsOddXOutRandomOrder(NumericMatrix MatOut, NumericMatrix MatIn, IntegerVector Yout, int pch2);
RcppExport SEXP _PairSeek_getpairsOddXOutRandomOrder(SEXP MatOutSEXP, SEXP MatInSEXP, SEXP YoutSEXP, SEXP pch2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type MatOut(MatOutSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type MatIn(MatInSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Yout(YoutSEXP);
    Rcpp::traits::input_parameter< int >::type pch2(pch2SEXP);
    rcpp_result_gen = Rcpp::wrap(getpairsOddXOutRandomOrder(MatOut, MatIn, Yout, pch2));
    return rcpp_result_gen;
END_RCPP
}
// RatioX
NumericMatrix RatioX(NumericMatrix MatIn, int pch2);
RcppExport SEXP _PairSeek_RatioX(SEXP MatInSEXP, SEXP pch2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type MatIn(MatInSEXP);
    Rcpp::traits::input_parameter< int >::type pch2(pch2SEXP);
    rcpp_result_gen = Rcpp::wrap(RatioX(MatIn, pch2));
    return rcpp_result_gen;
END_RCPP
}
// RankX
NumericMatrix RankX(NumericMatrix MatIn, int pch2);
RcppExport SEXP _PairSeek_RankX(SEXP MatInSEXP, SEXP pch2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type MatIn(MatInSEXP);
    Rcpp::traits::input_parameter< int >::type pch2(pch2SEXP);
    rcpp_result_gen = Rcpp::wrap(RankX(MatIn, pch2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PairSeek_getpairsOddXOut", (DL_FUNC) &_PairSeek_getpairsOddXOut, 4},
    {"_PairSeek_getpairsOddXOutRandomOrder", (DL_FUNC) &_PairSeek_getpairsOddXOutRandomOrder, 4},
    {"_PairSeek_RatioX", (DL_FUNC) &_PairSeek_RatioX, 2},
    {"_PairSeek_RankX", (DL_FUNC) &_PairSeek_RankX, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_PairSeek(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}