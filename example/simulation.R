#############
## File: simulation_study_1.R
## Purpose: Performs the first simulatin study
## Author: Richard Sizelove (sizeloverich@gmail.com)
## Last edited: 03/24/2023
## Dependencies: Rcpp, bbmle, numDeriv
#############

#  set the working directory to where this file is located
setwd("")
set.seed(111)

### load required dependencies
# needed for GEl approach
library(numDeriv)
###

### source the resources
# for GEL
source("resources/gel_resources/modified_GEL.R")

# for the proposed
# for first usage, install ICCOv using the command below. May require installing devtools
#devtools::install_github("sizelove/ICCov")
library(ICCov)

# for uniform MI
source("resources/simulation_resources/uniform_MI.R")
# for the data-generating mechanism (DGM)
source("resources/simulation_resources/DGM_functions.R")
# for tabulating results
source("resources/simulation_resources/tabulate_results.R")

### simulation parameters
# a note to users: the large-scale simulations mentioned in the paper were run on
# a high-performance computing cluster. Setting sample_size or n_replicates to be high
# will require substanial computing times
n_replicates = 10 # number of replicates to run
sample_size = 200

### data parameters 
alpha = c(-1.0, 0.5, 0.3) # main effect coefficients
beta = c(0.3, 0, -0.2) # IC coefficients
gamma = c(0.5, -0.5) # covariates impact on the IC event time
sigma_sq = 3# residual variance of Y
tau = 3 # Controls spacing of intermediate event, higher values space the intermediate event

# used for calculating cp
true_parameters = c(alpha, beta, sigma_sq, gamma)
names(true_parameters) = c("alpha_0", "alpha_1", "alpha_2", "beta_0", "beta_1", "beta_2", "sigma^2", "gamma_1", "gamma_2")
true_parameters_lm = c(alpha, beta, sigma_sq)
names(true_parameters_lm) = names(true_parameters)[1:7]

for(niter in 1:n_replicates){
  
  # generate simulated data
  
  # reverse the commenting below to allow for differential missing visits
  Simulation_Data = ICCov_DGM_no_missing_visits(sample_size, alpha, beta, gamma, sigma_sq, tau)
  #Simulation_Data = ICCov_DGM_differential_missing_visits(sample_size, alpha, beta, gamma, sigma_sq, tau)
  
  
  
  # run the proposed method
  Proposed_Results = ICCov_lm(timedResponse(Y, V) ~ Cov1 + Cov2 + (Cov1 + Cov2)*intCensCov("Event", Left, Right), covariate_function = "relu", data = Simulation_Data, stepsize = 1/(5*sqrt(sample_size)), max_iter = 600, tolerance = 1e-4)
  proposed_coeff= c(Proposed_Results$Coefficients$`Main Coefficients`, Proposed_Results$Coefficients$`Interval Censored Coefficients`, Proposed_Results$Coefficients$`Residual Error`, Proposed_Results$Coefficients$`Hazard Coefficients`)
  proposed_se = c(Proposed_Results$StdErr$`Main Coefficients`, Proposed_Results$StdErr$`Interval Censored Coefficients`, Proposed_Results$StdErr$`Residual Error`, Proposed_Results$StdErr$`Hazard Coefficients`)
  proposed_cp = true_parameters > proposed_coeff - qnorm(.975)*proposed_se & true_parameters < proposed_coeff +qnorm(.975)*proposed_se
  # store the Proposed results
  proposed_vec = c(proposed_coeff, proposed_se, proposed_cp)
  names(proposed_vec) = c("alpha_0", "alpha_1", "alpha_2", "beta_0", "beta_1", "beta_2", "sigma^2", "gamma_1", "gamma_2", "se_alpha_0", "se_alpha_1", "se_alpha_2", "se_beta_0", "se_beta_1", "se_beta_2", "se_sigma^2", "se_gamma_1", "se_gamma_2", "cp_alpha_0", "cp_alpha_1", "cp_alpha_2", "cp_beta_0", "cp_beta_1", "cp_beta_2", "cp_sigma^2", "cp_gamma_1", "cp_gamma_2")
  if(niter == 1){
    Summarized_Proposed_Results = proposed_vec
  }else{
    Summarized_Proposed_Results = rbind(Summarized_Proposed_Results, proposed_vec)
  }

  
  # Run the modified GEL method on non-right-censored subjects
  Dropped_Right_Censored_Data = Simulation_Data[Simulation_Data$Right != Inf, ]
  NPMLE_support = unique(c(Dropped_Right_Censored_Data$Left, Dropped_Right_Censored_Data$Right))
  Dropped_Right_Main_Matrix = cbind(Dropped_Right_Censored_Data$Cov1, Dropped_Right_Censored_Data$Cov2)
  colnames(Dropped_Right_Main_Matrix) = c("Cov1", "Cov2")
  Dropped_Right_IC_Matrix = cbind(rep(1, nrow(Dropped_Right_Censored_Data)), Dropped_Right_Censored_Data$Cov1, Dropped_Right_Censored_Data$Cov2)
  colnames(Dropped_Right_IC_Matrix) = c("RELU", "RELU*Cov1", "RELU*Cov2")
  GEL_results = Sizelove_modified_GELalgo(Dropped_Right_Censored_Data$Y, Dropped_Right_Censored_Data$V, Dropped_Right_Main_Matrix, Dropped_Right_IC_Matrix, Dropped_Right_Censored_Data$Left, Dropped_Right_Censored_Data$Right, NPMLE_support)
  GEL_coeff = GEL_results$Theta
  GEL_se = GEL_results$SE
  GEL_cp = true_parameters_lm > GEL_coeff - qnorm(.975)*GEL_se & true_parameters_lm < GEL_coeff + qnorm(.975)*GEL_se 
  # Store the GEL results
  GEL_vec = c(GEL_coeff, GEL_se, GEL_cp)
  names(GEL_vec) = c("alpha_0", "alpha_1", "alpha_2", "beta_0", "beta_1", "beta_2", "sigma^2","se_alpha_0", "se_alpha_1", "se_alpha_2", "se_beta_0", "se_beta_1", "se_beta_2", "se_sigma^2", "cp_alpha_0", "cp_alpha_1", "cp_alpha_2", "cp_beta_0", "cp_beta_1", "cp_beta_2", "cp_sigma^2")
  if(niter == 1){
    Summarized_GEL_Results = GEL_vec
  }else{
    Summarized_GEL_Results = rbind(Summarized_GEL_Results, GEL_vec)
  }
  
  # uniform MI
  Main_Matrix = cbind(rep(1, nrow(Simulation_Data)), Simulation_Data$Cov1, Simulation_Data$Cov2)
  colnames(Main_Matrix) = c("Intercept", "Cov1", "Cov2")
  IC_Matrix = cbind(rep(1, nrow(Simulation_Data)), Simulation_Data$Cov1, Simulation_Data$Cov2)
  colnames(IC_Matrix) = c("RELU", "RELU*COV1", "RELU*Cov2")
  MI_results = perform_uniform_multiple_imputation(Simulation_Data$Y, Simulation_Data$V, Main_Matrix, IC_Matrix, Simulation_Data$Left, Simulation_Data$Right, nimpute = 100)
  MI_cp = true_parameters_lm > MI_results$estimates - qnorm(.975)*sqrt(MI_results$var) & true_parameters_lm < MI_results$estimates + qnorm(.975)*sqrt(MI_results$var)
  MI_vec = c(MI_results$estimates, sqrt(MI_results$var), MI_cp)
  if(niter == 1){
    Summarized_UMI_Results = MI_vec
  }else{
    Summarized_UMI_Results = rbind(Summarized_UMI_Results, MI_vec)
  }
  
  
  # midpoint
  mid_results = lm(Y ~ Cov1 + Cov2 + ReLU_Midpoint + Cov1*ReLU_Midpoint + Cov2*ReLU_Midpoint, data = Simulation_Data)
  mid_results = summary(mid_results)
  mid_vec = tabulate_results_lm(mid_results, true_parameters_lm)
  if(niter == 1){
    Summarized_Mid_Results = mid_vec
  }else{
    Summarized_Mid_Results = rbind(Summarized_Mid_Results, mid_vec)
  }
  
  # rightpoint
  right_results = lm(Y ~ Cov1 + Cov2 + ReLU_Rightpoint + Cov1*ReLU_Rightpoint + Cov2*ReLU_Rightpoint, data = Simulation_Data)
  right_results = summary(right_results)
  right_vec = tabulate_results_lm(right_results, true_parameters_lm)
  if(niter == 1){
    Summarized_Right_Results = right_vec
  }else{
    Summarized_Right_Results = rbind(Summarized_Right_Results, right_vec)
  }
  
  # Exact time
  true_results = lm(Y ~ Cov1 + Cov2 + ReLU_Truth + Cov1*ReLU_Truth + Cov2*ReLU_Truth, data = Simulation_Data)
  true_results = summary(true_results)
  true_vec = tabulate_results_lm(true_results, true_parameters_lm)
  if(niter == 1){
    Summarized_True_Results = true_vec
  }else{
    Summarized_True_Results = rbind(Summarized_True_Results, true_vec)
  }
  
}

create_metric_table(Summarized_Proposed_Results, true_parameters)
create_metric_table(Summarized_GEL_Results, true_parameters_lm)
create_metric_table(Summarized_Mid_Results, true_parameters_lm)
create_metric_table(Summarized_Right_Results, true_parameters_lm)
create_metric_table(Summarized_True_Results, true_parameters_lm)