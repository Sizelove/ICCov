#############
## File: modified_GEL.R
## Purpose: Implements a moddified GEL (Gomez, Espinal, Langohr) method for linear regression with an interval-censored covariate.
##           changes include reformulating the covariate through the ReLU function, allowing a user-specified support, and conducting
##            variance estimation in a method similar to that described in Gomez, Espinal, and Langohr (2004)
###            and Langohr and Gomez (2014).
## Author: Richard Sizelove (sizeloverich@gmail.com)
## Last edited: 03/24/2023
## Dependencies: bbmle and numDeriv
#############

require(bbmle) ### for mle2() (required for original implementation)
require(numDeriv) ### for hessian() (required by Sizelove to compute numerical hessian)

### code edited to add "V", a n-length vector for the time at which Y is measured.
### also added "input_sj" which allows for user-control of the support of the interval-censored variable

## on April 1st, added x to be a covariate matrix rather than just a vector. tx allows for including interaction terms with the intermediate event

# change the tolerance level due to non-convergence for larger datasets

## Usage:
  # y : an n-length numeric vector, the response of interest
  # V : a n-length numeric vector, the time at which Y is measured
  # x : an n by p_1 matrix of main effects
  # tx : an n by p_2 matrix of effects corresponding to the intermediate event
  # zl:  an n-length numeric vector representing the left interval bracket
  # zr: an n-length numeric vector representing the right interval bracket
  # input_sj: a user-supplied support
  # toler: a double representing convergence tolerance 

  # Returns a list of estimates for the coefficients and their standard errors

  # A note to users: We attempted to leave this code as untouched as posisble, and only made edits where required to adequately compare the GEL method and the proposed method.
  # We have noted that at times, the mle2 function required by this function will return warnings regarding convergence failure.
  # this typically occurs within the first few iterations of the algorithm and, to our knowledge, does not seem to impact the final results.
Sizelove_modified_GELalgo <- function(y, V, x, tx, zl, zr, input_sj, toler=1e-6){
  ly <- length(y); lv <- length(V); lx <- nrow(x); ltx <- nrow(tx); ll <- length(zl); lr <- length(zr)
  
  if (min(ly, lv, lx, ltx, ll, lr) != max(ly, lv, lx, ltx, ll, lr))
    stop('Vectors must be of equal length!')
  if (any(zl>=zr))
    stop('Zl must be smaller than Zr!')
  
  # Number of observations
  n <- length(y)
  # Vector of Z's possible values
  
  #original GEL code is sj <- min(zl):max(zr)
  sj = input_sj
  
  # Length of Z's support
  m <- length(sj)
  # Auxiliary matrices of the data
  ymat <- matrix(rep(y, m), nrow=n)
  zlma <- matrix(rep(zl, m), nrow=n)
  zrma <- matrix(rep(zr, m), nrow=n)
  # Matrix of indicator variables
  smat <- matrix(rep(sj, n), byrow=T, nrow=n)
  alfas <- (zlma<=smat)*(smat<=zrma)
  
  # code added by Sizelove, Lin, Zeng.
  relu_function = function(little_v, little_sj){
    return(pmax(little_v - little_sj, 0.0))
  }
  v_by_sj = expand.grid(V, sj)
  relu_mat = matrix(relu_function(v_by_sj[,1], v_by_sj[,2]), nrow = n, ncol = length(sj))
  
  ## when the parameters are in a vector, these values give their locations
  alfIndex = 1
  betIndex = 2:(2 + ncol(x) - 1)
  gamIndex = (2 + ncol(x)):(2 + ncol(x) + ncol(tx) - 1)
  sigIndex = 2 + ncol(x) + ncol(tx)
  
  ## the edits to the code now involve replacing "smat", which gives the value of S_j
  ## and replacing it with "relu_mat", which gives the value of (V_i - S_j)_{+}
  # the indicators variable "alfas" still remains the same, as it gives us
  # the values of S_j which are plausible given Z_l and Z_r for each subject
  
  ## Initial values
  # Omega
  ome.hat <- rep(1/m, m)
  omat <- matrix(rep(ome.hat, n), byrow=T, nrow=n)
  # Theta
  zet <- pmax(V - 0.5*(zl+zr), 0.0) # start with midpoint
  zet_mat = zet*tx
  lm0 <- lm(y~x+zet_mat)
  alf.hat <- lm0$coef[alfIndex] # alpha is just the intercept
  bet.hat <- lm0$coef[betIndex] # beta is our alpha
  gam.hat <- lm0$coef[gamIndex] # gamma is beta for our model, coefficients for IC component.
  sig.hat <- summary(lm0)$sigma^2
  # GEL algorithm
  repeat{
    # Updated estimation of omega
    repeat{
      ome.old <- ome.hat
      omat <- matrix(rep(ome.hat, n), byrow=T, nrow=n)
      # line below was added since ymat and relu_mat are matrices
      beta_x_mat = matrix(rep(x%*%bet.hat, m), nrow = n)
      gam_tx = as.vector(tx %*% gam.hat)
      numerator <- alfas*dnorm(ymat, alf.hat+beta_x_mat+gam_tx*relu_mat, sqrt(sig.hat))*omat
      denominator <- matrix(rep(rowSums(numerator), m), nrow=n)
      nuumat <- numerator/denominator
      ome.hat <- colSums(nuumat)/n
      if (sum((ome.hat-ome.old)^2)/sum(ome.old^2)<toler)
        break
    }
    
    # Updated estimation of Theta
    alf.old <- alf.hat
    bet.old <- bet.hat
    gam.old <- gam.hat
    sig.old <- sig.hat
    thet.old <- c(alf.old, bet.old, gam.old, sig.old)
    lnfunc <- function(params){
      alf = params[alfIndex]
      bet = params[betIndex]
      gam = params[gamIndex]
      sig = params[sigIndex]
      beta_x_mat = matrix(rep(x%*%bet, m), nrow = n)
      gam_tx = as.vector(tx %*% gam)
      return(-sum(log(rowSums(alfas*dnorm(ymat, alf+beta_x_mat+gam_tx*relu_mat, sqrt(sig))*omat) +.00000001)))
    }
    name_vector = c("alpha", colnames(x), colnames(tx), "sigma")
    parnames(lnfunc) = name_vector
    m0 <- mle2(lnfunc, start=setNames(c(alf.hat,bet.hat,gam.hat,sig.hat),name_vector), method="L-BFGS-B", vecpar = TRUE)
    sm0 <- summary(m0)@coef
    # index below here had to be updated to account for X being a matrix
    alf.hat <- sm0[alfIndex, 1]
    bet.hat <- sm0[betIndex, 1]
    gam.hat <- sm0[gamIndex, 1]
    sig.hat <- sm0[sigIndex, 1]
    thet.hat <- c(alf.hat, bet.hat, gam.hat, sig.hat)
    if (sum((ome.hat-ome.old)^2)/sum(ome.old^2) + sum((thet.hat-thet.old)^2)/sum(thet.old^2) < toler)
      break
  }
  
  ### modified lnfunc allows for input into hessian() function from numDeriv package.
  ### this is akin the to the numerical approach Langohr and Gomez reports without details.
  modified_lnfunc = function(param){
    alf = param[1]
    bet = param[betIndex]
    gam = param[gamIndex]
    sig = param[sigIndex]
    beta_x_mat = matrix(rep(x%*%bet, m), nrow = n)
    gam_tx = as.vector(tx %*% gam)
    return(-sum(log(rowSums(alfas*dnorm(ymat, alf+beta_x_mat+gam_tx*relu_mat, sqrt(sig))*omat) +.00000001)))
  }
  SigInv = hessian(modified_lnfunc, thet.hat)
  vcov_est = solve(SigInv)
  se_est = sqrt(diag(vcov_est))
  # Distribution function of Z
  cbi <- cbind(sj, omega.hat=round(ome.hat, 4), FZ.hat=round(cumsum(ome.hat), 4))
  ok <- cbi[, 2]>0
  # For theta, we removed rounding to the fourth decimal (which was done by Langohr and Gomez)
  # the vcov and SE return arguments are added, as the original implementation did not
  # provide SE estimation
  return(list(Theta=thet.hat, Omega=cbi[ok, ], VCOV = vcov_est, SE = se_est))
}

