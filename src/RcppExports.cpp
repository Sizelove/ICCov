// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_internal_em_normalresponse
Rcpp::List c_internal_em_normalresponse(arma::vec Y, arma::mat X_HAZ, arma::mat X_MAIN, arma::mat X_IC, arma::mat function_values_matrix, std::vector<double> times, std::vector<double> Left, std::vector<double> Right, double stepsize, double tolerance, double pl_tolerance, arma::vec tau, unsigned int max_iter);
RcppExport SEXP _ICCov_c_internal_em_normalresponse(SEXP YSEXP, SEXP X_HAZSEXP, SEXP X_MAINSEXP, SEXP X_ICSEXP, SEXP function_values_matrixSEXP, SEXP timesSEXP, SEXP LeftSEXP, SEXP RightSEXP, SEXP stepsizeSEXP, SEXP toleranceSEXP, SEXP pl_toleranceSEXP, SEXP tauSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_HAZ(X_HAZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_MAIN(X_MAINSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_IC(X_ICSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type function_values_matrix(function_values_matrixSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type times(timesSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Left(LeftSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Right(RightSEXP);
    Rcpp::traits::input_parameter< double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< double >::type pl_tolerance(pl_toleranceSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(c_internal_em_normalresponse(Y, X_HAZ, X_MAIN, X_IC, function_values_matrix, times, Left, Right, stepsize, tolerance, pl_tolerance, tau, max_iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ICCov_c_internal_em_normalresponse", (DL_FUNC) &_ICCov_c_internal_em_normalresponse, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_ICCov(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
