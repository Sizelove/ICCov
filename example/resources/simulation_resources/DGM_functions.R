#############
## File: DGM_functions.R
## Purpose: Provides function for the data-generating mechanisms used in the two simulation studies
## Author: Richard Sizelove (sizeloverich@gmail.com)
## Last edited: 03/24/2023
#############

ICCov_DGM_no_missing_visits = function(sample_size, alpha = c(-1.0, 0.5, 0.2), beta = c(0.3, 0, -0.2), gamma = c(0.3, -0.3), sigma_sq = 3, tau = 3){
  
  X_1 = runif(sample_size, 0, 1)
  X_2 = rbinom(sample_size, 1, 0.5)
  X = as.matrix(cbind(X_1, X_2))
  colnames(X) = paste("Cov", 1:ncol(X), sep = "")
  exp_xgamma = exp(X %*% gamma)
  
  ## event time
  v = runif(sample_size, 0, 1)
  kappa = -1 * log(1 - v) / exp_xgamma
  t = (5*exp(kappa) - 5)
  t = pmax(t, 0)
  
  # generate the times
  U_1 = runif(sample_size, 0, 2)
  U_2 = U_1 + 0.1 + rexp(sample_size, 1)*tau
  U_3 = U_2 + 0.1 + rexp(sample_size, 1)*tau
  U_4 = U_3 + 0.1 + rexp(sample_size, 1)*tau
  U_5 = U_4  + 0.1 + rexp(sample_size, 1)*tau
  
  ## monitoring time for Y
  V = rep(20, sample_size) + rexp(sample_size, rate = 1)
  
  # Creates the interval-censored observation. Finds the U_j and U_{j+1} which bracket t
  visits = matrix(data = c(U_1, U_2, U_3, U_4, U_5), nrow = 5, byrow = TRUE)
  visit_data = matrix(data = rep(0, 2*sample_size), nrow = sample_size, byrow = TRUE)
  for(ii in 1:sample_size){
    a = na.omit(visits[visits[,ii] >= t[ii],ii])
    if(length(a) == 0){
      right = Inf
    }else{
      right = a[1]
    }
    visit_data[ii, 2] = right
    a = na.omit(visits[visits[,ii] <= t[ii],ii])
    if(length(a) == 0){
      left = 0
    }else{
      left = a[length(a)]
    }
    visit_data[ii, 1] = left
  }
  sort = order(visit_data[,2])
  index_sort = sort[8]

  # Simulate the response
  Y_mean = alpha[1] +  X%*%alpha[-1] + pmax((V - t), 0)*beta[1] + pmax((V - t), 0)*(X%*%beta[-1])
  Y = rnorm(n = length(t), Y_mean, sqrt(sigma_sq))
  
  # Provide the midpoint, rightpoint, and exact data as columns in the dataset
  midpoint = pmax(V - (visit_data[,2] + visit_data[,1])/2, 0.0)
  rightpoint = pmax(V - visit_data[,2], 0.0)
  truth = pmax(V - t, 0.0)
  Generated_Data = cbind(Y, V, X, visit_data, midpoint, rightpoint, truth)
  
  
  Generated_Data = as.data.frame(Generated_Data)
  colnames(Generated_Data) = c("Y", "V", colnames(X), "Left", "Right", "ReLU_Midpoint", "ReLU_Rightpoint", "ReLU_Truth")
  Generated_Data = Generated_Data[order(Generated_Data$Right), ]

  return(Generated_Data)
}

ICCov_DGM_no_missing_visits = function(sample_size, alpha = c(-1.0, 0.5, 0.2), beta = c(0.3, 0, -0.2), gamma = c(0.3, -0.3), sigma_sq = 3, tau = 3){
  
  X_1 = runif(sample_size, 0, 1)
  X_2 = rbinom(sample_size, 1, 0.5)
  X = as.matrix(cbind(X_1, X_2))
  colnames(X) = paste("Cov", 1:ncol(X), sep = "")
  exp_xgamma = exp(X %*% gamma)
  
  ## event time
  v = runif(sample_size, 0, 1)
  kappa = -1 * log(1 - v) / exp_xgamma
  t = (5*exp(kappa) - 5)
  t = pmax(t, 0)
  
  # generate the times
  U_1 = runif(sample_size, 0, 2)
  U_2 = U_1 + 0.1 + rexp(sample_size, 1)*tau
  U_3 = U_2 + 0.1 + rexp(sample_size, 1)*tau
  U_4 = U_3 + 0.1 + rexp(sample_size, 1)*tau
  U_5 = U_4  + 0.1 + rexp(sample_size, 1)*tau
  
  ## now add missingness
  missing_prob = ifelse(X_2 == 1, 0.65, 0.15)
  U_2 = ifelse(rbinom(sample_size, 1, missing_prob), NA, U_2)
  U_3 = ifelse(rbinom(sample_size, 1, missing_prob), NA, U_3)
  U_4 = ifelse(rbinom(sample_size, 1, missing_prob), NA, U_4)
  U_5 = ifelse(rbinom(sample_size, 1, missing_prob), NA, U_5)
  
  ## monitoring time for Y
  V = rep(20, sample_size) + rexp(sample_size, rate = 1)
  
  # Creates the interval-censored observation. Finds the U_j and U_{j+1} which bracket t
  visits = matrix(data = c(U_1, U_2, U_3, U_4, U_5), nrow = 5, byrow = TRUE)
  visit_data = matrix(data = rep(0, 2*sample_size), nrow = sample_size, byrow = TRUE)
  for(ii in 1:sample_size){
    a = na.omit(visits[visits[,ii] >= t[ii],ii])
    if(length(a) == 0){
      right = Inf
    }else{
      right = a[1]
    }
    visit_data[ii, 2] = right
    a = na.omit(visits[visits[,ii] <= t[ii],ii])
    if(length(a) == 0){
      left = 0
    }else{
      left = a[length(a)]
    }
    visit_data[ii, 1] = left
  }
  sort = order(visit_data[,2])
  index_sort = sort[8]
  
  # Simulate the response
  Y_mean = alpha[1] +  X%*%alpha[-1] + pmax((V - t), 0)*beta[1] + pmax((V - t), 0)*(X%*%beta[-1])
  Y = rnorm(n = length(t), Y_mean, sqrt(sigma_sq))
  
  # Provide the midpoint, rightpoint, and exact data as columns in the dataset
  midpoint = pmax(V - (visit_data[,2] + visit_data[,1])/2, 0.0)
  rightpoint = pmax(V - visit_data[,2], 0.0)
  truth = pmax(V - t, 0.0)
  Generated_Data = cbind(Y, V, X, visit_data, midpoint, rightpoint, truth)
  
  
  Generated_Data = as.data.frame(Generated_Data)
  colnames(Generated_Data) = c("Y", "V", colnames(X), "Left", "Right", "ReLU_Midpoint", "ReLU_Rightpoint", "ReLU_Truth")
  Generated_Data = Generated_Data[order(Generated_Data$Right), ]
  
  return(Generated_Data)
}