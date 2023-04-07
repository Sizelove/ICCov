#############
## File: perform_uniform_multiple_imputation
## Purpose: Implements a naive uniform multiple imputation approach.
## Author: Richard Sizelove (sizeloverich@gmail.com)
## Last edited: 03/24/2023
#############

## Usage:
  # Y : an n-length numeric vector, the response of interest
  # Y_measurement_time: a n-length numeric vector, the time at which Y is measured
  # X_MAIN: an n by p_1 matrix of main effects
  # X_IC: an n by p_2 matrix of effects corresponding to the intermediate event
  # L_i: an n-length numeric vector representing the left interval bracket
  # R_i: an n-length numeric vector representing the right interval bracket

  # Returns a list of estimates for the coefficients and their standard errors
perform_uniform_multiple_imputation = function(Y, Y_measurement_time, X_MAIN, X_IC, L_i, R_i, nimpute = 100){
  for(i in 1:nimpute){
    # treat R_i = Inf as non-incident and simulate times 
      imputed_t_bool = R_i == Inf
      t = rep(NA, length(R_i))
      t[!imputed_t_bool] = runif(n = sum(!imputed_t_bool), min = L_i[!imputed_t_bool], max = R_i[!imputed_t_bool])
      t[imputed_t_bool] = Inf
      RELU_component = pmax(Y_measurement_time - t, 0)
      
    X_IC_RELU = X_IC*RELU_component
    lm_results = lm(Y ~ -1 + X_MAIN + X_IC_RELU)
    p = length(lm_results$coefficients)
    sig_est = summary(lm_results)$sigma^2
    if(i == 1){
      estimates_matrix = c(lm_results$coefficients, sig_est)
      var_matrix = c(diag(vcov(lm_results)), (2 * sig_est^2 / (nrow(X_MAIN) - p) ) )
    }else{
      estimates_matrix = rbind(estimates_matrix,  c(lm_results$coefficients, sig_est))
      var_matrix = rbind(var_matrix,  c(diag(vcov(lm_results)), (2 * sig_est^2 / (nrow(X_MAIN) - p) ) ))
    }
  
  }
  
  estimate_means = colMeans(estimates_matrix)
  within_var = colMeans(var_matrix)
  between_var = apply(estimates_matrix, 2, var)
  total_var = within_var + (1 + (1/nimpute))*between_var
  return(list(estimates = estimate_means, var = total_var))
}