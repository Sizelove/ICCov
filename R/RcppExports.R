# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

c_internal_em_normalresponse <- function(Y, X_HAZ, X_MAIN, X_IC, function_values_matrix, times, Left, Right, stepsize, tolerance, pl_tolerance, tau, max_iter) {
    .Call('_ICCov_c_internal_em_normalresponse', PACKAGE = 'ICCov', Y, X_HAZ, X_MAIN, X_IC, function_values_matrix, times, Left, Right, stepsize, tolerance, pl_tolerance, tau, max_iter)
}

