#############
## File: tabulate_functions.R
## Purpose: Provides function for summarizing results from R's lm() function.
## Author: Richard Sizelove (sizelove@live.unc.edu)
## Last edited: 03/24/2023
#############

## Usage:
  # model_results : the results from a single call of lm()
  # truth : the true model parameters
  
  # Returns a list of estimates for the coefficients, their standard errors, and coverage
tabulate_results_lm = function(model_results, truth){
  
  n = length(model_results$residuals)
  coefficients = model_results$coefficients[,1]
  see = model_results$coefficients[,2]
  cp = truth[1:6] > coefficients - qnorm(.975)*see & truth[1:6] < coefficients + qnorm(.975)*see 
  
  # sigma_sq based on chi-square confidence interval
  sigma_sq = truth[7]
  est_sigma_sq = model_results$sigma^2
  sigma_sq_se = sqrt(2*est_sigma_sq^2/(n-6))
  sigma_cp = (sample_size - 6)*est_sigma_sq/qchisq(.975, sample_size - 6) < sigma_sq & sigma_sq < (sample_size - 6)*est_sigma_sq/qchisq(.025, sample_size - 6)
  
  return(c(coefficients, est_sigma_sq, see, sigma_sq_se, cp, sigma_cp))
  
}

## Usage:
# summary_table : a table of results from each run of the simulation for the specific method
# truth : the true model parameters

# Returns a list of metrics: bias, standard errors, standard error estimates, and coverage probability
create_metric_table = function(summary_table, truth){
  
  ncoeff = ncol(summary_table) / 3
  bias = colMeans(summary_table[,1:ncoeff]) - truth
  se = apply(summary_table, 2, sd)[1:ncoeff] # true se 
  see = colMeans(summary_table[,(ncoeff+1):(2*ncoeff)]) # estimate sd
  cp = colMeans(summary_table[,(2*ncoeff+1):(3*ncoeff)])

  metric_table = cbind(bias, se, see, cp)
  rownames(metric_table) = names(truth)
  colnames(metric_table) = c("Bias", "SE", "SEE", "CP")
  return(metric_table)
}