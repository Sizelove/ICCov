
#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
#include <chrono>
using namespace Rcpp;
using namespace std::chrono;

arma::vec normal_calc_indiv_pl_value(arma::vec exp_gamma_X, arma::vec alpha_X, arma::vec beta_X, double sigma_sq, arma::mat function_values_matrix, arma::vec Y, arma::vec lambda_k, arma::vec times, std::vector<double> Left, std::vector<double> Right, arma::vec tau){
  unsigned int n = Y.n_rows;
  unsigned int m = times.n_elem;
  arma::vec pl_val_total(n);
  pl_val_total.zeros();

  arma::mat p_ik_mat(m+2,n);
  p_ik_mat.zeros();

  // parameters
  arma::dvec Big_lambda_k(m+2);
  //

  // update Big_lambda_k
  Big_lambda_k[0] = lambda_k[0];
  for(unsigned int kk = 1; kk < m+1; kk++){
    Big_lambda_k[kk] = Big_lambda_k[kk - 1] + lambda_k[kk];
  }
  Big_lambda_k[m+1] = -999999;

  // initialize
  arma::vec modified_times(m+2);
  modified_times(0) = 0;
  for(unsigned int kk = 0; kk < m; kk++){
    modified_times(kk+1) = times[kk];
  }
  modified_times(m+1) = arma::datum::inf;

  for(unsigned int ii = 0; ii < n ; ii++){
    double pl_val = 0;
    for(unsigned int kk = 0; kk < m+2; kk++){
      if(Left[ii] < modified_times[kk] && modified_times[kk] <= Right[ii]){
        double gaussian_part = pow(2.0*sigma_sq*M_PI, -0.5)*exp( ( -(0.5)*pow(Y.at(ii) - alpha_X[ii] - function_values_matrix(ii, kk)*beta_X[ii],2) /(sigma_sq)) );
        if(kk > 1 && kk < m + 1){
          pl_val += (exp(-Big_lambda_k[kk-1]*exp_gamma_X[ii]) - exp(-Big_lambda_k[kk]*exp_gamma_X[ii]))*gaussian_part;
        }else if(kk == 1){
          pl_val += (1 - exp(-Big_lambda_k[1]*exp_gamma_X[ii]))*gaussian_part;
        }else if(kk == m + 1){
          pl_val += exp(-1*Big_lambda_k[m]*exp_gamma_X[ii])*gaussian_part;
        }
      }
    }
    if(pl_val != 0){
      pl_val_total[ii] = log(pl_val);
    }
  }
  return(pl_val_total);
}

arma::vec fit_pl_score(arma::vec exp_gamma_X, arma::vec alpha_X, arma::vec beta_X, double sigma_sq, arma::mat function_values_matrix, arma::vec Y, arma::vec lambda_tilde_k, double pl_tolerance, unsigned int pl_max_iter, arma::vec times, std::vector<double> Left, std::vector<double> Right, arma::vec tau){
  //

  unsigned int n = Y.n_rows;
  unsigned int m = times.size();
  arma::vec pl_val_total(n);

  // storage
  arma::mat EW_ik_mat(m+2, n);
  EW_ik_mat.zeros();
  arma::mat p_ik_mat(m+2,n);
  p_ik_mat.zeros();

  // containers
  arma::vec denom_sum_iter(m+2);
  arma::vec old_lambda(m+2);
  arma::vec Big_lambda_old(m+2);

  arma::vec sum_pik_vec(m);
  sum_pik_vec.zeros();
  // parameters
  arma::dvec Big_lambda_k(m+2);
  //

  // iterators
  double delta = 0;
  double sum_iter = 0;

  // update Big_lambda_k
  arma::vec lambda_k = lambda_tilde_k;
  Big_lambda_k[0] = lambda_k[0];
  for(unsigned int kk = 1; kk < m + 1; kk++){
    Big_lambda_k[kk] = Big_lambda_k[kk - 1] + lambda_k[kk];
  }
  Big_lambda_k[m + 1] = -999999;

  // initialize
  arma::vec modified_times(m+2);
  modified_times(0) = 0;
  for(unsigned int kk = 0; kk < m; kk++){
    modified_times(kk+1) = times[kk];
  }
  modified_times(m+1) = arma::datum::inf;
  //

  unsigned int em_iter = 0;
  delta = 1000;
  while(delta > pl_tolerance && em_iter < pl_max_iter){
    em_iter ++;
    delta = 0;


     for(unsigned int ii = 0; ii < n; ii++){
       sum_iter = 0;
       for(unsigned int kk = 0; kk < m + 2; kk++){
         if(Left[ii] < modified_times[kk] && modified_times[kk] <= Right[ii]){
           double gaussian_part = pow(2.0*sigma_sq*M_PI, -0.5)*exp( ( -(0.5)*pow(Y.at(ii) - alpha_X[ii] - function_values_matrix(ii, kk)*beta_X[ii],2) /(sigma_sq)) );
           if(kk > 1 && kk < m + 1){
             p_ik_mat(kk, ii) = (exp(-Big_lambda_k[kk-1]*exp_gamma_X[ii]) - exp(-Big_lambda_k[kk]*exp_gamma_X[ii]))*gaussian_part;
           }else if(kk == 1){
             p_ik_mat(kk, ii) = (1 - exp(-Big_lambda_k[1]*exp_gamma_X[ii]))*gaussian_part;
           }else if(kk == m + 1){
             p_ik_mat(kk, ii) = exp(-1*Big_lambda_k[m]*exp_gamma_X[ii])*gaussian_part;
           }
          sum_iter += (p_ik_mat(kk, ii));
         }
       }
       for(unsigned int kk = 0; kk < m+2; kk++){
         if(Left[ii] < modified_times[kk] && modified_times[kk] <= Right[ii]){
           if(sum_iter == 0){sum_iter = 1;}
           p_ik_mat(kk, ii) = p_ik_mat(kk, ii) / sum_iter;
        }
      }
    }


     double inner_val = 0;
     for(unsigned int ii = 0; ii < n; ii++){
       double time_accu = 0;
       for(unsigned int kk = 0; kk < m+1; kk++){
         if(Left[ii] < modified_times[kk] && modified_times[kk] <= Right[ii]){
           inner_val = (1 - exp(-lambda_k[kk]*exp_gamma_X[ii]));
           EW_ik_mat.at(kk, ii) = 0;
           if(inner_val != 0){
             EW_ik_mat.at(kk, ii) = p_ik_mat(kk, ii)*(lambda_k[kk]*exp_gamma_X[ii]/inner_val);
           }
           EW_ik_mat.at(kk, ii) += (time_accu)*(lambda_k[kk]*exp_gamma_X[ii]);
           time_accu += p_ik_mat(kk, ii);
         }
     }
     }

      old_lambda = lambda_k;
      lambda_k.zeros();
      denom_sum_iter.zeros();
      for(unsigned int ii = 0; ii < n; ii++){
        for(unsigned int kk = 0; kk < m + 1; kk++){
          if(Right[ii] >= modified_times[kk]){
            lambda_k[kk] += EW_ik_mat(kk, ii);
            denom_sum_iter[kk] += exp_gamma_X[ii];
          }
        }
      }
      for(unsigned int kk = 0; kk < m + 1; kk++){
        lambda_k[kk] = lambda_k[kk] / denom_sum_iter[kk];
        old_lambda[kk] = std::abs(lambda_k[kk] - old_lambda[kk]);
      }
      R_CheckUserInterrupt();
      Big_lambda_old = Big_lambda_k;
      lambda_k[0] = 0;
      Big_lambda_k[0] = lambda_k[0];
        for(unsigned int kk = 1; kk < m + 1; kk++){
          Big_lambda_k[kk] = Big_lambda_k[kk - 1] + lambda_k[kk];
          if(Big_lambda_k[kk] > 1){
            Big_lambda_old[kk] = std::abs(Big_lambda_k[kk] - Big_lambda_old[kk]) / Big_lambda_k[kk];
          }else{
            Big_lambda_old[kk] = std::abs(Big_lambda_k[kk] - Big_lambda_old[kk]);
          }
        }
      delta = Big_lambda_old.max();
    }

  pl_val_total =normal_calc_indiv_pl_value(exp_gamma_X, alpha_X, beta_X, sigma_sq, function_values_matrix, Y, lambda_k, times, Left, Right, tau);
  return(pl_val_total);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List c_internal_em_normalresponse(arma::vec Y, arma::mat X_HAZ, arma::mat X_MAIN, arma::mat X_IC, arma::mat function_values_matrix, std::vector<double> times, std::vector<double> Left, std::vector<double> Right, double stepsize, double tolerance, double pl_tolerance, arma::vec tau, unsigned int max_iter){
  // changed on 09 / 06 / 2021 to include assumption that the subjects are in order of R_i timel;

  pl_tolerance = tolerance;

  unsigned int n = X_HAZ.n_rows;
  unsigned int m = times.size();
  unsigned int p = X_HAZ.n_cols;

  arma::mat EW_ik_mat(m + 2, n);
  EW_ik_mat.zeros();
  arma::mat p_ik_mat(m+2,n);
  p_ik_mat.zeros();
  arma::dvec exp_gamma_X(n, arma::fill::zeros);

arma::mat X_sums(m+2, p);
  X_sums.zeros();
  arma::vec denom_sum_iter(m + 2);
  denom_sum_iter.zeros();
  arma::rowvec score_gam(p);
  score_gam.zeros();
  arma::cube X_outer(p, p, m+2);
  arma::mat info_gam(p, p);
  arma::vec old_lambda(m+2);
  arma::vec Big_lambda_old(m+2);
  unsigned int alpha_p = X_MAIN.n_cols;
  unsigned int beta_p = X_IC.n_cols;
  unsigned int response_model_p = X_MAIN.n_cols + X_IC.n_cols;
  arma::mat wls_mat(response_model_p, response_model_p);
  arma::vec wls_score(response_model_p);
  wls_mat.zeros();
  wls_score.zeros();
  arma::vec aug_vec(response_model_p);
  aug_vec.zeros();
  arma::vec alpha_beta_res(response_model_p, arma::fill::zeros);
  arma::vec old_alpha_beta_res(response_model_p, arma::fill::zeros);
  bool conv_bool = FALSE;
  arma::vec sum_pik_vec(m + 2);
  sum_pik_vec.zeros();
  // parameters
  arma::dvec lambda_k(m + 2);
  arma::dvec Big_lambda_k(m + 2);
  lambda_k.fill((1/static_cast<double>(m)));
  arma::dvec gamma(p, arma::fill::zeros);
  arma::dvec old_gamma(p, arma::fill::zeros);
  arma::dvec alpha(alpha_p, arma::fill::zeros);
  arma::vec alpha_X(n, arma::fill::zeros);
  arma::dvec beta(beta_p, arma::fill::zeros);
  arma::vec beta_X(n, arma::fill::zeros);
  double sigma_sq = arma::var(Y);
  //

  // iterators
  double delta = 0;
  double sum_iter = 0;

  // update Big_lambda_k
  lambda_k[0] = 0;
  Big_lambda_k[0] = lambda_k[0];
  for(unsigned int kk = 1; kk < m + 1; kk++){
    Big_lambda_k[kk] = Big_lambda_k[kk - 1] + lambda_k[kk];
  }
  Big_lambda_k[m + 1] = -9999; // this should never be accessed, so here is a trip-wire.

  // initialize
  exp_gamma_X = exp(X_HAZ*gamma);
  arma::vec modified_times(m+2);
  modified_times(0) = 0;
  for(unsigned int kk = 0; kk < m; kk++){
      modified_times(kk+1) = times[kk];
  }
  modified_times(m+1) = arma::datum::inf;

  arma::rowvec holder_vec;
  arma::mat holder_mat;


  /////////////// begin loop here //////////////////////
  unsigned int em_iter = 0;
  delta = 1000;

  while(delta > tolerance && em_iter < max_iter && (conv_bool == FALSE)){
   em_iter ++;
   delta = 0;
    // update p_ik_matF
    p_ik_mat.zeros();

    for(unsigned int ii = 0; ii < n; ii++){
      sum_iter = 0;
      for(unsigned int kk = 0; kk < m+2; kk++){
        if(Left[ii] < modified_times[kk] && modified_times[kk] <= Right[ii]){
          double gaussian_part = pow(2.0*sigma_sq*M_PI, -0.5) * exp( ( -(0.5)*pow(Y.at(ii) - alpha_X[ii] - function_values_matrix(ii, kk)*beta_X[ii],2) /(sigma_sq)) );

          if(kk > 1 && kk < m + 1){
            p_ik_mat(kk, ii) = (exp(-Big_lambda_k[kk-1]*exp_gamma_X[ii]) - exp(-Big_lambda_k[kk]*exp_gamma_X[ii]))*gaussian_part;
          }else if(kk == 1){
            p_ik_mat(kk, ii) = (1 - exp(-Big_lambda_k[1]*exp_gamma_X[ii]))*gaussian_part;
          }else if(kk == m + 1){
            p_ik_mat(kk, ii) = exp(-1*Big_lambda_k[m]*exp_gamma_X[ii])*gaussian_part;
          }
       sum_iter += (p_ik_mat(kk, ii));
      }
    }
    for(unsigned int kk = 0; kk < m+2; kk++){
      if(modified_times[kk] > Left[ii] && modified_times[kk] <= Right[ii]){
      p_ik_mat(kk, ii) = p_ik_mat(kk, ii) / sum_iter;
    }
  }
} // end ii

arma::vec p_ik_vec(n);
p_ik_vec.zeros();
for(unsigned int ii = 0; ii < n; ii++){
  for(unsigned int kk = 0; kk < m + 2; kk++){
    p_ik_vec(ii) += p_ik_mat(kk, ii);
  }
}

    double inner_val = 0;
    for(unsigned int ii = 0; ii < n; ii++){
      double time_accu = 0;
      for(unsigned int kk = 0; kk < m+1; kk++){
        if(Left[ii] < modified_times[kk] && Right[ii] >= modified_times[kk]){
          inner_val = (1 - exp(-lambda_k[kk]*exp_gamma_X[ii]));
          EW_ik_mat.at(kk, ii) = 0;
          if(inner_val != 0){
            EW_ik_mat.at(kk, ii) = p_ik_mat(kk, ii)*(lambda_k[kk]*exp_gamma_X[ii]/inner_val);
          }
          EW_ik_mat.at(kk, ii) += (time_accu)*(lambda_k[kk]*exp_gamma_X[ii]);
          time_accu += p_ik_mat(kk, ii);
        }
    }
    }

    /// update gamma
    info_gam.zeros();
    X_outer.zeros();
    X_sums.zeros();
    denom_sum_iter.zeros();
    score_gam.zeros();
    for(unsigned int ii = 0; ii < n; ii++){
      holder_vec = X_HAZ.row(ii)*exp_gamma_X[ii];
      holder_mat = X_HAZ.row(ii).t()*exp_gamma_X[ii]*X_HAZ.row(ii);
      for(unsigned int kk = 1; kk < m + 1; kk++){
        if(Right[ii] >= modified_times[kk]){
          X_sums.row(kk) += holder_vec;
          X_outer.slice(kk) += holder_mat;
          denom_sum_iter[kk] += exp_gamma_X[ii];
        }
      }
    }
    for(unsigned int kk = 0; kk < m + 1; kk++){
      X_sums.row(kk) = X_sums.row(kk) / denom_sum_iter[kk];
      X_outer.slice(kk) = X_outer.slice(kk) / denom_sum_iter[kk];
    }
    for(unsigned int ii = 0; ii < n; ii++){
      for(unsigned int kk = 1; kk < m + 1; kk++){
        if(Right[ii] >= modified_times[kk]){
          score_gam += EW_ik_mat(kk, ii)*(X_HAZ.row(ii) - X_sums.row(kk));
          info_gam += EW_ik_mat(kk, ii)*(X_outer.slice(kk) - (X_sums.row(kk).t()*X_sums.row(kk)));
        }
      }
    }

    old_gamma = gamma;
    gamma = gamma + solve(info_gam, score_gam.t());
    for(unsigned int ii = 0; ii < p; ii++){
      old_gamma[ii] = std::abs(old_gamma[ii] - gamma[ii] );
    }
    exp_gamma_X = exp(X_HAZ*gamma);


    R_CheckUserInterrupt();

    // update alpha and beta
  old_alpha_beta_res = alpha_beta_res;
  aug_vec.zeros();
  wls_mat.zeros();
  wls_score.zeros();
 for(unsigned int ii = 0; ii < n; ii++){
    for(unsigned int pp = 0; pp < alpha_p; pp++){
      aug_vec(pp) = X_MAIN(ii, pp);
    }
   for(unsigned int kk = 0; kk < m+2; kk++){
     if(Left[ii] < modified_times[kk] && modified_times[kk] <= Right[ii]){
      for(unsigned int pp = 0; pp < beta_p; pp++){
        aug_vec(pp + alpha_p) = X_IC(ii, pp)*function_values_matrix(ii, kk);
      }
      wls_mat += p_ik_mat(kk, ii)*aug_vec*aug_vec.t();
      wls_score += aug_vec*p_ik_mat(kk, ii)*Y(ii);
    }
   }
 }
 alpha_beta_res = solve(wls_mat, wls_score);
 alpha = alpha_beta_res.subvec(0, alpha_p - 1);

   alpha_X = X_MAIN*alpha;
   beta = alpha_beta_res.subvec(alpha_p, response_model_p - 1);
   beta_X = X_IC*beta;
   for(unsigned int pp = 0; pp < response_model_p; pp++){
       old_alpha_beta_res[pp] = std::abs(alpha_beta_res[pp] - old_alpha_beta_res[pp]);
  }


  double mse = 0;
  for(unsigned int ii = 0; ii < n; ii++){
    for(unsigned int kk = 0; kk < m+2; kk++){
       if(modified_times[kk] > Left[ii] && modified_times[kk] <= Right[ii]){
      double gaussian_part = pow(Y.at(ii) - alpha_X[ii] - function_values_matrix(ii, kk)*beta_X[ii],2);
       mse += p_ik_mat(kk, ii)*gaussian_part;
     }
    }
  }
  double sigma_sq_delta = sigma_sq;
  sigma_sq = mse / static_cast<double>(n);

  sigma_sq_delta = std::abs(sigma_sq - sigma_sq_delta);

    old_lambda = lambda_k;
    lambda_k.zeros();
    denom_sum_iter.zeros();
    for(unsigned int ii = 0; ii < n; ii++){
      for(unsigned int kk = 0; kk < m + 1; kk++){
        if(Right[ii] >= modified_times[kk]){
          lambda_k[kk] += EW_ik_mat(kk, ii);
          denom_sum_iter[kk] += exp_gamma_X[ii];
        }
      }
    }
    for(unsigned int kk = 0; kk < m + 1; kk++){
      lambda_k[kk] = lambda_k[kk] / denom_sum_iter[kk];
      old_lambda[kk] = std::abs(lambda_k[kk] - old_lambda[kk]);
    }
    R_CheckUserInterrupt();

    Big_lambda_old = Big_lambda_k;
    arma::vec real_Big_lambda_old = Big_lambda_k;
    lambda_k[0] = 0;
    Big_lambda_k[0] = lambda_k[0];
      for(unsigned int kk = 1; kk < m + 1; kk++){
        Big_lambda_k[kk] = Big_lambda_k[kk - 1] + lambda_k[kk];
        if(Big_lambda_k[kk] > 1){
          Big_lambda_old[kk] = std::abs(Big_lambda_k[kk] - Big_lambda_old[kk]) / Big_lambda_k[kk];
        }else{
          Big_lambda_old[kk] = std::abs(Big_lambda_k[kk] - Big_lambda_old[kk]);
        }
      }

    delta = Big_lambda_old.max();
    delta = std::max(delta, old_gamma.max());
    delta = std::max(old_alpha_beta_res.max(), delta);
  }

  unsigned int pl_max_iter = max_iter;
  unsigned int total_d = p + alpha.n_elem + beta.n_elem + 1;

  // pl containers
  arma::mat pl_matrix(total_d, total_d);
  pl_matrix.zeros();
  arma::mat scores(total_d, n);
  scores.zeros();

  // containers for the pl function evaluated at (theta)
  arma::vec f_beta_e(n);
  arma::vec f_beta(n);

  exp_gamma_X = exp(X_HAZ * (gamma));
  alpha_X = X_MAIN * alpha;
  beta_X = X_IC * beta;
  unsigned int counter = 0;
  // pl(theta)
  f_beta =  fit_pl_score(exp_gamma_X, alpha_X, beta_X, sigma_sq,function_values_matrix, Y, lambda_k, pl_tolerance, pl_max_iter, times, Left, Right, tau);
    // pl for gamma, pl(theta + h_n)
    arma::vec col_1_e(p);
    for(unsigned int pp = 0; pp < p; pp++){
      col_1_e.zeros();
      col_1_e(pp) = stepsize;

      exp_gamma_X = exp(X_HAZ * (gamma + col_1_e));
      f_beta_e = fit_pl_score(exp_gamma_X, alpha_X, beta_X, sigma_sq, function_values_matrix, Y, lambda_k, pl_tolerance, pl_max_iter, times, Left, Right, tau);
      for(unsigned int ii = 0; ii < n; ii++){
        scores(counter, ii) = (f_beta_e(ii) - f_beta(ii)) / stepsize;
      }
      counter++;
    }
    exp_gamma_X = exp(X_HAZ * (gamma));

    // pl for alpha
    col_1_e.resize(alpha_p);
    for(unsigned int pp = 0; pp < alpha_p; pp++){
      col_1_e.zeros();
      col_1_e(pp) = stepsize;
      arma::vec temp_alpha = alpha + col_1_e;
      alpha_X = X_MAIN * temp_alpha;
      f_beta_e = fit_pl_score(exp_gamma_X, alpha_X, beta_X, sigma_sq, function_values_matrix, Y, lambda_k, pl_tolerance, pl_max_iter, times, Left, Right, tau);
      for(unsigned int ii = 0; ii < n; ii++){
        scores(counter, ii) = (f_beta_e(ii) - f_beta(ii)) / stepsize;
      }
      counter++;
    }
    alpha_X = X_MAIN * alpha;

    // pl for beta
    col_1_e.resize(beta_p);
    for(unsigned int pp = 0; pp < beta_p; pp++){
      col_1_e.zeros();
      col_1_e(pp) = stepsize;
      arma::vec temp_beta = beta + col_1_e;
      beta_X = X_IC * temp_beta;
      f_beta_e = fit_pl_score(exp_gamma_X, alpha_X, beta_X, sigma_sq, function_values_matrix, Y, lambda_k, pl_tolerance, pl_max_iter, times, Left, Right, tau);
      for(unsigned int ii = 0; ii < n; ii++){
        scores(counter, ii) = (f_beta_e(ii) - f_beta(ii)) / stepsize;
      }
      counter++;
    }
    beta_X = X_IC*beta;

    while(sigma_sq < 2*stepsize){
      stepsize = stepsize / 2;
    }
    f_beta_e = fit_pl_score(exp_gamma_X, alpha_X, beta_X, sigma_sq + stepsize, function_values_matrix, Y, lambda_k, pl_tolerance, pl_max_iter, times, Left, Right, tau);
    f_beta = fit_pl_score(exp_gamma_X, alpha_X, beta_X, sigma_sq - stepsize, function_values_matrix, Y, lambda_k, pl_tolerance, pl_max_iter, times, Left, Right, tau);
    for(unsigned int ii = 0; ii < n; ii++){
      scores(counter, ii) = (f_beta_e(ii) - f_beta(ii)) / (2*stepsize);
    }


    // computing the pl matrix
    for(unsigned int ii = 0; ii < n; ii++){
      pl_matrix += scores.col(ii)*(scores.col(ii).t());
    }
    pl_matrix = pl_matrix.i();


  // calculating the posterior value of Y, E[Y | \theta, \Lambda] = \sum E[Y | \theta, \Lambda, T]*p(T).
  arma::vec post_y(n);
  post_y.zeros();
  for(unsigned int ii = 0; ii < n; ii++){
    for(unsigned int kk = 0; kk < m+2; kk++){
      if(Left[ii] < modified_times[kk] && modified_times[kk] <= Right[ii]){
        post_y(ii)  += p_ik_mat(kk, ii)*(alpha_X(ii) + beta_X(ii)*function_values_matrix(ii, kk));
      }
    }
  }

return(Rcpp::List::create(_["cov"] = pl_matrix, _["gamma"] = gamma, _["alpha"] = alpha, _["beta"] = beta,_["Y"] = Y, _["sigma_sq"] = sigma_sq,
_["Lambda"] = Big_lambda_k, _["times"] = times, _["posterior_Ey"] = post_y));
}
