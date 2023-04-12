#'
#' 
#' 
#' 
#' Creates a "timedResponse" class variable for input into ICCov regression models
#' @name timedResponse
#' @param response A numeric response vector of interest
#' @param measurement_time A numeric vector containing the measurement time for each element in response
#' @export
timedResponse <- function(response, measurement_time){
  if(length(response) != length(measurement_time)){
    stop("response and measurement_time must be the same length")
  }
  if(any(is.na(response_time)) | any(is.na(measurement_time))){
    stop("Neither response nor measurement_time can contain NA values")
  }
  return_list = list("response" = response, "time" = measurement_time)
  return(return_list)
}

#' Fit semiparametric regression models with interval-censored covariates
#' @useDynLib ICCov
#' @param formula A object of class `"formula"`: a symbolic description of the model to be fit. See `Details` for more information.
#' @param covariate_function A `"function"` object which specifies how the intermediate event time and response are linked. See `Details` for examples.
#' @param data A `"data.frame"` object containing the variables included in the model.
#' @param data An optional object of class `"formula"` specifying how the intermediate event is linked to the covariates. If left absent, then the intermediate event will be modeled as dependent on all covariates included to the `formula` argument.
#' @param tolerance A scalar value for the convergence tolerance of the EM algorithm.
#' @param pl_tolerance A scalar value for the convergence tolerance of the profile likelihood method used for variance estimation.
#' @param stepsize A scalar value for the stepsize used in the profile likelihood method. By default, 1/(5*n^1/2) is used.
#' @export 
ICCov_lm <- function(formula, covariate_function, data, hazard_formula = NULL, tolerance = 1e-4, pl_tolerance = 1e-4, max_iter = 500, stepsize = NULL){
  # if formula isn't class formula, exit
  #print("Here inside ICCov_lm")
  if(is.null(stepsize)){
    stepsize = nrow(data)^(-1/2)
  }
  arguments <- as.list(match.call())
  callsign = match.call()
  tau = eval(arguments$response_measurement_time, data)
 # print(tau)
  
  if(class(formula) != "formula"){
    stop("Variable `formula` must have be formula class and be of the form Y ~ X_1 + X_2...")
  }
  response_name = deparse(formula[[2]])
  
  # get the data names of everything
  hazard_set_flag = FALSE
  if(!is.null(hazard_formula)){
    haz_terms = terms(as.formula(hazard_formula))
    hazard_set_flag = TRUE
  }
  formula_list = formula
  term_object = terms(formula)
  
  ## Gets just the individual variables, not interaction, and separates them
  #("Before formula list")
  formula_list = formula
  term_object = terms(formula)
  iccov_indicator = which(grepl("intCensCov\\(.*[^\\)]\\)", attr(term_object, "term.labels")))
  if(length(iccov_indicator) > 0){
    #print("Inside here")
    iccov_main = term_object[iccov_indicator[1]]
  }else{
    print("No interval censored covariate found")
  }
  ### separate out the IntCensCov variable and find the left, right, and label
  iccov_main_var = attr(iccov_main[1], "term.labels")
  iccov_label_orig = regmatches(iccov_main_var, regexpr('\\".*\\"', iccov_main_var))
  iccov_label = gsub('\\"', replacement = "", iccov_label_orig)
  iccov_split = strsplit(iccov_main_var, split = ", ")
  iccov_left_name = iccov_split[[1]][2]
  
  iccov_right_name = gsub(")", "", iccov_split[[1]][3])
  form_cov = update(old = formula_list, new = drop.terms(terms(formula_list), iccov_indicator))
  form_cov = update(form_cov, NULL ~ .)
  form_cov_terms = terms(form_cov)
  response_var = all.vars(formula)[1]
  time_var = all.vars(formula)[2]
  
  ### checking for character variables
  attr_form_cov_terms = attr(form_cov_terms, "term.labels")
  if(hazard_set_flag == FALSE){
    haz_terms = form_cov_terms
  }
  attr_haz_cov_terms = attr(haz_terms, "term.labels")
  make_character_cov = attr_form_cov_terms[which(grepl("^as.factor\\(", attr(form_cov_terms, "term.labels")) | grepl("^as.character\\(", attr(form_cov_terms, "term.labels")))]
  make_character_cov = c(make_character_cov, attr_haz_cov_terms[which(grepl("^as.factor\\(", attr(haz_terms, "term.labels")) | grepl("^as.character\\(", attr(haz_terms, "term.labels")))])
  make_character_list = gsub("^as.factor\\(", "", make_character_cov)
  make_character_list = gsub("^as.character\\(", "", make_character_list)
  make_character_list = gsub(")$", "", make_character_list)

  terms_list_complete = c(attr(form_cov_terms, "term.labels"), attr(haz_terms, "term.labels"), iccov_left_name, iccov_right_name, response_var, time_var)
  terms_list_complete = setdiff(terms_list_complete, make_character_cov)
  terms_list_complete = c(terms_list_complete, make_character_list)
  terms_list_complete %in% colnames(data)
  if(!all(terms_list_complete %in% colnames(data))){
    stop(paste("Some formula terms not found in the provided dataframe:", paste(terms_list_complete[!(terms_list_complete %in% colnames(data))], collapse = " ")))
  }
  reduced_data = data[, unique(terms_list_complete)]

  if(length(make_character_list) > 0){
    for(i in 1:length(make_character_list)){
      reduced_data[, make_character_list[i]] = as.character(reduced_data[, make_character_list[i]])
    }
  }
  data = reduced_data
  data = reduced_data[complete.cases(reduced_data), ]
  if(!isTRUE(all.equal(data, reduced_data))){
    warning("Some observations dropped for incomplete records")
  }

  iccov_left = as.numeric(unlist(data[,iccov_left_name]))
  iccov_right = as.numeric(unlist(data[,iccov_right_name]))

  ## this is just the non-IntCens effects
  form_cov = update(old = formula_list, new = drop.terms(terms(formula_list), iccov_indicator))
  form_cov = update(form_cov, NULL ~ .)
  design_matrix = model.matrix(form_cov, data=data )
  main_effect_matrix = as.matrix(design_matrix)

  if(is.null(hazard_formula)){
    main_effect_1 = main_effect_matrix == 1
    index = colSums(main_effect_1) == nrow(main_effect_matrix)
    names_hazard = colnames(main_effect_matrix)[!index]
    hazard_matrix = main_effect_matrix[, !index]
  }else{
    haz_design_matrix = model.matrix(as.formula(hazard_formula), data=data )
    haz_effect_matrix = as.matrix(haz_design_matrix)
   
    haz_effect_1 = haz_effect_matrix == 1
    index = colSums(haz_effect_1) == nrow(haz_effect_matrix)
    names_hazard = colnames(haz_effect_matrix)[!index]
    hazard_matrix = haz_effect_matrix[, !index]
    
  }

  ic_terms = terms(formula_list)
  iccov_indicator = which(!grepl("intCensCov\\(.*[^\\)]\\)", attr(ic_terms, "term.labels")))
  form_cov = update(old = formula_list, new = drop.terms(terms(formula_list), iccov_indicator))
  form_cov = update(form_cov, NULL ~ .)
  term_list = gsub(pattern = "intCensCov\\(.*[^\\)]\\):", replacement = "", attr(terms(form_cov), "term.labels"))
  term_list =  gsub(pattern = ":intCensCov\\(.*[^\\)]\\)", replacement = "", term_list)
  term_list = gsub(pattern = "intCensCov\\(.*[^\\)]\\)", replacement = "", term_list)
  ic_int_form = term_list[term_list != ""]
  ic_int_form = paste("~ ", paste(term_list, collapse = " + "), sep = "")

  if(ic_int_form == "~ "){
    ic_int_form = "~ 1"
  }

  ic_matrix = as.matrix(model.matrix(formula(ic_int_form), data = data))
  
  ## from the new object, get jut the variables other than intCensCov
  
  ## read in the interaction levels so we can create the analysis dataset

  
  ## adds in the interaction term and their names into the dataset,
  ## new variable names in cov_names
  main_cov_names = colnames(main_effect_matrix)
  int_cov_names = colnames(ic_matrix)
  haz_names = colnames(hazard_matrix)
  #complete_case_index_main = complete.cases(main_effect_matrix)
  #complete_case_index_int = complete.cases(ic_matrix)
  #complete_case_index_haz = complete.cases(hazard_matrix)
 # complete_case_index = complete_case_index_main & complete_case_index_int & complete_case_index_haz
  #data = data[complete_case_index, ]
  
  response_var = all.vars(formula)[1]
  time_var = all.vars(formula)[2]
  if(response_var %in% colnames(data)){
    response_vec = data[, response_var]
  }else{
    stop(paste("Variable", response_var, "not found in dataframe"))
  }
  response_vec = data[, response_var]
  if(time_var %in% colnames(data)){
    Tau_vec = data[, time_var]
  }else{
    stop(paste("Variable", time_var, "not found in dataframe"))
  }
  
  
  if(!is.numeric(response_vec)){
    print(response_vec)
    stop("Response must be numeric for gaussian regression")
  }
  if(length(unique(response_vec)) <= 1){
    stop("Response must have at least 2 unique values")
  }
  Left = iccov_left
  Right = iccov_right
  if(length(Left) != length(Right)){
    stop("Interval vectors left and right have different lengths")
  }
  if(!all(Left <= Right)){
    print(cbind(Left[Left > Right], Right[Left > Right]))
    stop("Some interval left endpoints are greater than their corresponding right endpoints")
  }
  
  ### reordering the variables
  order_vec = order(Right)
  response_vec = response_vec[order_vec]
  main_effect_matrix= as.matrix(main_effect_matrix[order_vec, ])
  ic_matrix = as.matrix(ic_matrix[order_vec, ])
  hazard_matrix = as.matrix(hazard_matrix)
  hazard_matrix = as.matrix(hazard_matrix[order_vec, ])
  Tau_vec = Tau_vec[order_vec]
  Left = Left[order_vec]
  Right = Right[order_vec]
  times_vec = unique(c(Left, Right))
  times_vec = sort(times_vec)
  times_vec = setdiff(times_vec, c(0, Inf))

  ### Making the user-supplied matrix of values
  modified_times = c(0, times_vec, Inf)
  user_supplied_function = covariate_function
  if(tolower(user_supplied_function) == "indicator"){
    user_supplied_function = function(Response_Timepoint, modified_times){
      ifelse(Response_Timepoint >= modified_times, 1, 0)
    }
  }else if(tolower(user_supplied_function) == "relu"){
    user_supplied_function = function(Response_Timepoint, modified_times){
      pmax(Response_Timepoint - modified_times, 0.0)
    }
  }else if(is.null(user_supplied_function)){
    stop("No function was provided!")
  }
  
    ### write later to test user supplied function
    test_vec = c(1, 2, 3)
    test_vec_2 = c(2, 3, 4, 5)
    testFunction <- function(first_argument, second_argument) {
      out <- tryCatch(
        {
          user_supplied_function(first_argument, second_argument)
        },
        error=function(cond) {
          message("Covariate function doesn't appear to be valid.")
          message("Below is the original error message returned by R:")
          message(cond)
          return(-9)
        },
        warning=function(cond) {
        },
        finally={
         
        }
      )    
      return(out)
    }
    test_value = testFunction(test_vec, test_vec_2)
    if(length(test_value) > 0){
      if(test_value == -9){
      stop("Something wrong with the user-specified covariate function.")
      }
    }
    
  value_grid = expand.grid(Tau_vec, modified_times)
  function_values_matrix = user_supplied_function(value_grid[,1], value_grid[,2])
  function_values_matrix = matrix(function_values_matrix, nrow = length(Tau_vec), ncol = length(modified_times))
  
  
  returned_results = c_internal_em_normalresponse(Y = response_vec, X_HAZ = hazard_matrix, X_MAIN = main_effect_matrix, X_IC = ic_matrix, function_values_matrix = as.matrix(function_values_matrix), times_vec, Left, Right, stepsize = stepsize, tolerance, pl_tolerance, tau = Tau_vec, max_iter)

  ## get baseline hazard
  baseline_hazard = returned_results$Lambda
  times = returned_results$times
  baseline_hazard_matrix = matrix(c(times, baseline_hazard), byrow = FALSE, ncol = 2)
  colnames(baseline_hazard_matrix) = c("t", "Lambda(t)")
  ### get the coefficients 
  alpha_list = returned_results$alpha
  names(alpha_list) = colnames(main_effect_matrix)
  beta_list = returned_results$beta
  names(beta_list) = paste(iccov_label, "*", colnames(ic_matrix))
  if(length(colnames(ic_matrix)) > 0){
    if(colnames(ic_matrix)[1] == "(Intercept)"){
      names(beta_list)[1] = iccov_label
    }
  }
  gamma_list = returned_results$gamma
  names(gamma_list) = names_hazard
  sigma_sq = returned_results$sigma_sq
  coefficients_list = list(alpha_list, beta_list, gamma_list, sigma_sq)
  names(coefficients_list) = c("Main Coefficients", "Interval Censored Coefficients", "Hazard Coefficients", "Residual Error")
  
  #### Extract the standad errors, name them, and put them into a list
  cov = returned_results$cov
  dag = diag(cov)
  gamma_se_list = sqrt(dag[1:length(gamma_list)])
  alpha_se_list = sqrt(dag[(length(gamma_list)+1):(length(gamma_list)+length(alpha_list))])
  beta_se_list = sqrt(dag[(length(gamma_list)+1+length(alpha_list)):(length(alpha_list)+length(beta_list)+length(gamma_list))])
  sigma_se = sqrt(dag[length(diag(cov))])
  
  
  names(alpha_se_list) = paste("StdErr", names(alpha_list))
  names(beta_se_list) = paste("StdErr", names(beta_list))
  names(gamma_se_list) = paste("StdErr", names(gamma_list))
  names(sigma_se) = "StdErr residual_variance"
  
  se_list = list(alpha_se_list, beta_se_list, gamma_se_list, sigma_se)
  names(se_list) = c("Main Coefficients", "Interval Censored Coefficients", "Hazard Coefficients", "Residual Error")
  
  runtime = returned_results$runtime
  post_y = returned_results$posterior_Ey
  final_list = list(coefficients_list, se_list, baseline_hazard_matrix, cov, callsign, post_y, response_vec)
  names(final_list) = c("Coefficients", "StdErr", "Baseline Hazard", "cov", "call", "fitted_values", "response")
  class(final_list) = "ICCov_LM"
  return(final_list)
}

#' @export summary.ICCov_LM
#' @export 
summary.ICCov_LM = function(object){
  if(class(object) == "myclass"){
    print(class(object))
  }
  call = object[['call']]
  est = object$Coefficients$`Main Coefficients`
  se = object$StdErr$`Main Coefficients`
  zval = est / se
  prob_z = pnorm(abs(zval), 0, 1, lower.tail = FALSE)
  main_coefficients_mat = cbind("Estimate" = est, "Std. Error" = se, 
                                                              "Z value" = zval, "Pr(>|Z|)" = prob_z)
  
  colnames(main_coefficients_mat) =  c("Estimate", "Std. Error", "z value", 
                                          "Pr(>|z|)")
  rownames(main_coefficients_mat) = names(est)
  
  est = object$Coefficients$`Interval Censored Coefficients`
  se = object$StdErr$`Interval Censored Coefficients`
  zval = est / se
  prob_z = pnorm(abs(zval), 0, 1, lower.tail = FALSE)
  expo_coefficients_mat = cbind("Estimate" = est, "Std. Error" = se, 
                                "Z value" = zval, "Pr(>|Z|)" = prob_z)
  colnames(expo_coefficients_mat) =  c("Estimate", "Std. Error", "z value", 
                                       "Pr(>|z|)")
  rownames(expo_coefficients_mat) = names(est)
  
  est = object$Coefficients$`Hazard Coefficients`
  se = object$StdErr$`Hazard Coefficients`
  zval = est / se
  prob_z = pnorm(abs(zval), 0, 1, lower.tail = FALSE)
  haz_coefficients_mat = cbind("Estimate" = est, "Std. Error" = se, 
                                "Z value" = zval, "Pr(>|Z|)" = prob_z)
  colnames(haz_coefficients_mat) =  c("Estimate", "Std. Error", "z value", 
                                       "Pr(>|z|)")
  rownames(haz_coefficients_mat) = names(est)
  
  output_object = list("call" = call, "main_coefficients_mat" = main_coefficients_mat,
                       "expo_coefficients_mat" = expo_coefficients_mat,
                       "haz_coefficients_mat" = haz_coefficients_mat)
  class(output_object) = "summary.ICCov_LM"
  print.summary_ICCov_LM(output_object)
  #return(output_object)
}

print_ICCov_lm_coeffmat = function(matrixin){
  tag =   symnum(matrixin[,4], cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 symbols = c("***", "**", "*", 
                             ".", " "))
  matrixin = round(matrixin, digits = 4)
  matrixin[matrixin[,4] < 0.001, 4] = "<0.0001"
  matrixin = cbind((matrixin), format(tag))
  print.default(matrixin, quote = FALSE, right = TRUE)
  cat("---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}
#' @export print.summary_ICCov_LM
#' @export 
print.summary_ICCov_LM = function(object){
  x = object
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", 
                         collapse = "\n"), "\n\n", sep = "")
  cat("\nCoefficients:\n")
  print_ICCov_lm_coeffmat(x$main_coefficients_mat)
  
  cat("\nInterval Censored * Coefficients: \n")
  print_ICCov_lm_coeffmat(x$expo_coefficients_mat)
  
  cat("\nHazard Coefficients: \n")
  print_ICCov_lm_coeffmat(x$haz_coefficients_mat)

}
