#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double norm_rs(double a, double b);
double half_norm_rs(double a, double b);
double unif_rs(double a, double b);
double exp_rs(double a, double b);
double rnorm_trunc(double mu, double sigma, double lower, double upper);

// [[Rcpp::export]]
Rcpp::List BaySeqMPeak(IntegerMatrix Y, IntegerMatrix X, NumericVector s, double e, double f, double a, double b, int iter) {
  int n = Y.nrow();
  int W = Y.ncol();
  int KK = X.ncol();
  int K = (KK - 1)/2;
  
  // Hyperparameters
  double phi_s = 10.0;
  double a_pi = 1.0;
  double b_pi = 1.0;
  double a_phi = 0.001;
  double b_phi = 0.001;
  double mu_alpha_ctrl = 0;
  double sigma_alpha_ctrl = 10;
  double mu_alpha_ip = 0;
  double sigma_alpha_ip = 10;
  double a_alpha = a;
  double b_alpha = b;
  double a_beta = a;
  double b_beta = b;
  double a_theta = a;
  double b_theta = b;
  double e_gamma = e;
  double f_gamma = f;
  double e_delta = e;
  double f_delta = f;
  double e_xi = e;
  double f_xi = f;
  
  // Algorithm settings
  int burn = iter/2;
  int total_count_lwr = 10;
  int EE = W*0.05;
  if(EE == 0)
  {
    EE = 2;
  }
  double tau_phi = 1;
  double tau_alpha_0_ctrl = 0.1;
  double tau_alpha_0_ip = 0.1;
  double tau_alpha = 1;
  double tau_beta = 1;
  double tau_theta = 1;
  
  // Temporary variables
  int i, w, k, kk, it, ee;
  int H_sum = 0;
  int count = 0;
  int count_2 = 0;
  double max_temp, sum_temp, count_temp, hastings, phi_temp, alpha_temp, b_temp;
  double pi_mean = 0.0;
  NumericVector prob_temp(2);
  NumericVector loglambda_temp(n);
  IntegerVector temp_3(n, 0);
  NumericVector mh_try(5, 0.0);
  NumericVector mh_accept(5, 0.0);
  NumericVector pi_store(iter);
  IntegerMatrix H(n, W);
  NumericMatrix H_ppi(n, W);
  NumericVector phi(W);
  NumericVector phi_mean(W);
  NumericVector alpha_0_input(W);
  NumericVector alpha_0_ctrl_mean(W);
  NumericVector alpha_0_ip(W);
  NumericVector alpha_0_ip_mean(W);
  NumericMatrix B(KK, W);
  IntegerMatrix Gamma(KK, W);
  NumericMatrix logLambda(n, W);
  NumericVector alpha_mean(W);
  NumericMatrix Beta_mean(K, W);
  NumericMatrix Theta_mean(K, W);
  NumericVector gamma_ppi(W);
  NumericVector gamma_sum(iter);
  NumericMatrix Delta_ppi(K, W);
  NumericVector Delta_sum(iter);
  NumericMatrix Xi_ppi(K, W);
  NumericVector Xi_sum(iter);
  LogicalVector flag(W);
  
  // Preprocessing
  for(w = 0; w < W; w++)
  {
    count = 0;
    count_2 = 0;
    for(i = 0; i < n; i++)
    {
      count = count + Y(i, w);
      if(Y(i, w) != 0)
      {
        count_2++;
      }
    }
    if(count <= total_count_lwr || count_2 <= 2)
    {
      flag(w) = 1; 
    }
    else
    {
      flag(w) = 0;
    }
  }
  
  // Initialization
  double pi = 0.5;
  for(w = 0; w < W; w++)
  {
    if(flag(w) == 0)
    {
      phi(w) = phi_s;
      phi_mean(w) = 0;
      alpha_0_input(w) = runif(1, 0, 2)(0);
      alpha_0_ctrl_mean(w) = 0;
      alpha_0_ip(w) = runif(1, 0, alpha_0_input(w))(0);
      alpha_0_ip_mean(w) = 0;
      alpha_mean(w) = 0;
      for(kk = 0; kk < KK; kk++)
      {
        Gamma(kk, w) = rbinom(1, 1, 0.05)(0);
        if(Gamma(kk, w) == 1)
        {
          if(kk == 0)
          {
            B(kk, w) = alpha_0_ip(w) - alpha_0_input(w) + rnorm(1, 0, 1)(0);
          }
          else
          {
            B(kk, w) = rnorm(1, 0, 1)(0);
          }
        }
        else
        {
          if(kk == 0)
          {
            B(kk, w) = alpha_0_ip(w) - alpha_0_input(w);
          }
          else
          {
            B(kk, w) = 0;
          }
        }
      }
      for(i = 0; i < n; i++)
      {
        if(Y(i, w) == 0)
        {
          H(i, w) = rbinom(1, 1, 0.5)(0);
        }
        H_ppi(i, w) = 0;
        if(H(i, w) == 1)
        {
          H_sum++;
        }
        logLambda(i, w) = alpha_0_input(w);
        for(kk = 0; kk < KK; kk++)
        {
          if(X(i, kk) == 1) 
          {
            if(kk == 0)
            {
              logLambda(i, w) = logLambda(i, w) + B(kk, w);
            }
            else
            {
              if(Gamma(kk, w) == 1)
              {
                logLambda(i, w) = logLambda(i, w) + B(kk, w);
              }
            }
          }
        }
        for(k = 0; k < K; k++)
        {
          Beta_mean(k, w) = 0;
          Theta_mean(k, w) = 0;
          Delta_ppi(k, w) = 0;
          Xi_ppi(k, w) = 0;
        }
      }
    }
    else
    {
      phi(w) = 0;
      phi_mean(w) = 0;
      alpha_0_input(w) = 0;
      alpha_0_ctrl_mean(w) = 0;
      alpha_0_ip(w) = 0;
      alpha_0_ip_mean(w) = 0;
      alpha_mean(w) = 0;
      for(kk = 0; kk < KK; kk++)
      {
        Gamma(kk, w) = 0;
        B(kk, w) = 0;
      }
      for(i = 0; i < n; i++)
      {
        H(i, w) = 0;
        H_ppi(i, w) = 0;
        logLambda(i, w) = 0;
        for(k = 0; k < K; k++)
        {
          Beta_mean(k, w) = 0;
          Theta_mean(k, w) = 0;
          Delta_ppi(k, w) = 0;
          Xi_ppi(k, w) = 0;
        }
      }
    }
  }
  
  // MCMC
  count = 0;
  for(it = 0; it < iter; it++)
  {
    // Update pi
    pi = rbeta(1, a_pi + H_sum, b_pi + n*W - H_sum)(0);
    
    // Update H
    H_sum = 0;
    for(w = 0; w < W; w++)
    {
      if(flag(w) == 0)
      {
        for(i = 0; i < n; i++)
        {
          if(Y(i, w) == 0)
          {
            prob_temp(0) = phi(w)*(log(phi(w)) - log(s(i)*exp(logLambda(i, w)) + phi(w))) + log(1 - pi);
            prob_temp(1) = log(pi);
            max_temp = max(prob_temp);
            prob_temp(0) = prob_temp(0) - max_temp;
            prob_temp(1) = prob_temp(1) - max_temp;
            prob_temp(0) = exp(prob_temp(0));
            prob_temp(1) = exp(prob_temp(1));
            sum_temp = prob_temp(0) + prob_temp(1);
            prob_temp(0) = prob_temp(0)/sum_temp;
            prob_temp(1) = prob_temp(1)/sum_temp;
            H(i, w) = rbinom(1, 1, prob_temp(1))(0);
            if(H(i, w) == 1)
            {
              H_sum++;
            }
          }
        }
      }
    }
    
    // Update phi
    for(w = 0; w < W; w++)
    {
      if(flag(w) == 0)
      {
        count_temp = 0;
        do {
          phi_temp = rgamma(1, phi(w)*phi(w)/tau_phi, tau_phi/phi(w))(0);
          count_temp++;
        } while (phi_temp < 1 && count_temp < 1000);
        if (count_temp == 1000)
        {
          phi_temp = 10;
        }
        hastings = 0;
        for(i = 0; i < n; i++)
        {
          if(H(i, w) == 0) {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + Y(i, w)) - (phi_temp + Y(i, w))*log(phi_temp + s(i)*exp(logLambda(i, w)));
            hastings = hastings - (phi(w)*log(phi(w)) - lgamma(phi(w)) + lgamma(phi(w) + Y(i, w)) - (phi(w) + Y(i, w))*log(phi(w) + s(i)*exp(logLambda(i, w))));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi(w)) - b_phi*phi(w));
        if (it >= burn) {
          mh_try(0)++;
        }
        if(hastings >= log(double(rand()%10001)/10000))
        {
          phi(w) = phi_temp;
          if (it >= burn) {
            mh_accept(0)++;
          }
        }
      }
    }
    
    // Update alpha_0_input
    for(w = 0; w < W; w++)
    {
      if(flag(w) == 0)
      {
        alpha_temp = rnorm_trunc(alpha_0_input(w), tau_alpha_0_ctrl, alpha_0_ip(w), 10);
        hastings = -(alpha_temp - mu_alpha_ctrl)*(alpha_temp - mu_alpha_ctrl)/2/sigma_alpha_ctrl/sigma_alpha_ctrl;
        hastings = hastings - (-(alpha_0_input(w) - mu_alpha_ctrl)*(alpha_0_input(w) - mu_alpha_ctrl)/2/sigma_alpha_ctrl/sigma_alpha_ctrl);
        for(i = 0; i < n; i++)
        {
          if(X(i, 0) == 0 && H(i, w) == 0)
          {
            hastings = hastings + Y(i, w)*log(s(i)*exp(logLambda(i, w) - alpha_0_input(w) + alpha_temp)) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w) - alpha_0_input(w) + alpha_temp));
            hastings = hastings - (Y(i, w)*log(s(i)*exp(logLambda(i, w))) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w))));
          }
        }
        mh_try(1)++;
        if (hastings >= log(double(rand()%10001)/10000))
        {
          for(i = 0; i < n; i++)
          {
            if(X(i, 0) == 0)
            {
              logLambda(i, w) = logLambda(i, w) - alpha_0_input(w) + alpha_temp;
            }
          }
          B(0, w) = B(0, w) + alpha_0_input(w) - alpha_temp;
          alpha_0_input(w) = alpha_temp;
          if(it >= burn) {
            mh_accept(1)++;
          }
        }
      }
    }
    
    // Update alpha_0_ip
    for(w = 0; w < W; w++)
    {
      if(flag(w) == 0)
      {
        alpha_temp = rnorm_trunc(alpha_0_ip(w), tau_alpha_0_ip, -10, alpha_0_input(w));
        hastings = -(alpha_temp - mu_alpha_ip)*(alpha_temp - mu_alpha_ip)/2/sigma_alpha_ip/sigma_alpha_ip;
        hastings = hastings - (-(alpha_0_ip(w) - mu_alpha_ip)*(alpha_0_ip(w) - mu_alpha_ip)/2/sigma_alpha_ip/sigma_alpha_ip);
        for(i = 0; i < n; i++)
        {
          if(X(i, 0) == 1 && H(i, w) == 0)
          {
            hastings = hastings + Y(i, w)*log(s(i)*exp(logLambda(i, w) - alpha_0_ip(w) + alpha_temp)) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w) - alpha_0_ip(w) + alpha_temp));
            hastings = hastings - (Y(i, w)*log(s(i)*exp(logLambda(i, w))) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w))));
          }
        }
        mh_try(2)++;
        if (hastings >= log(double(rand()%10001)/10000))
        {
          for(i = 0; i < n; i++)
          {
            if(X(i, 0) == 1)
            {
              logLambda(i, w) = logLambda(i, w) - alpha_0_ip(w) + alpha_temp;
            }
          }
          B(0, w) = B(0, w) - alpha_0_ip(w) + alpha_temp;
          alpha_0_ip(w) = alpha_temp;
          if(it >= burn) {
            mh_accept(2)++;
          }
        }
      }
    }
    
    // Update B, i.e. alpha, Beta, and Delta
    for(w = 0; w < W; w++)
    {
      if(flag(w) == 0)
      {
        for(kk = 0; kk < KK; kk++)
        {
          if(Gamma(kk, w) == 1)
          {
            if(kk == 0)
            {
              alpha_temp = rnorm_trunc(B(kk, w) - alpha_0_ip(w) + alpha_0_input(w), tau_alpha/10, 0, 10);
              b_temp = alpha_0_ip(w) - alpha_0_input(w) + alpha_temp;
              hastings = (-a_alpha - 1/2)*log(b_alpha + alpha_temp*alpha_temp/2);
              hastings = hastings - (-a_alpha - 1/2)*log(b_alpha + B(kk, w)*B(kk, w)/2);
            }
            else
            {
              if(kk < K + 1)
              {
                b_temp = rnorm_trunc(B(kk, w), tau_beta/10, -10, 10);
                hastings = (-a_beta - 1/2)*log(b_beta + b_temp*b_temp/2);
                hastings = hastings - (-a_beta - 1/2)*log(b_beta + B(kk, w)*B(kk, w)/2);
              }
              else
              {
                b_temp = rnorm_trunc(B(kk, w), tau_theta/10, -10, 10);
                hastings = (-a_theta - 1/2)*log(b_theta + b_temp*b_temp/2);
                hastings = hastings - (-a_theta - 1/2)*log(b_theta + B(kk, w)*B(kk, w)/2);
                hastings =0;
              }
            }
            for(i = 0; i < n; i++)
            {
              if(X(i, kk) == 1 && H(i, w) == 0)
              {
                hastings = hastings + Y(i, w)*log(s(i)*exp(logLambda(i, w) - B(kk, w) + b_temp)) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w) - B(kk, w) + b_temp));
                hastings = hastings - (Y(i, w)*log(s(i)*exp(logLambda(i, w))) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w))));
              }
            }
            mh_try(3)++;
            if(hastings >= log(double(rand()%10001)/10000))
            {
              for(i = 0; i < n; i++)
              {
                if(X(i, kk) == 1)
                {
                  logLambda(i, w) = logLambda(i, w) - B(kk, w) + b_temp;
                }
              }
              B(kk, w) = b_temp;
              if(it >= burn) {
                mh_accept(3)++;
              }
            }
          }
        }
      }
    }
    
    // Update Gamma, i.e. gamma, Delta, and Xi
    for(kk = 0; kk < KK; kk++)
    {
      for(ee = 0; ee < EE; ee++)
      {
        w = rand()%W;
        if(flag(w) == 0)
        {
          hastings = 0;
          if(Gamma(kk, w) == 0) // Add
          {
            if(kk == 0)
            {
              alpha_temp = rnorm_trunc(0, tau_alpha, 0, 10);
              b_temp = alpha_0_ip(w) - alpha_0_input(w) + alpha_temp;
              hastings = (a_alpha*log(b_alpha) - lgamma(a_alpha) + lgamma(a_alpha + 0.5) - (a_alpha + 0.5)*log(b_alpha + alpha_temp*alpha_temp/2));
              hastings = hastings - (-log(tau_alpha) - alpha_temp*alpha_temp/(2*tau_alpha*tau_alpha));
              hastings = hastings + e_gamma;
              if(w < W - 1)
              {
                hastings = hastings + f_gamma*Gamma(kk, w + 1);
              }
              if(w > 0)
              {
                hastings = hastings + f_gamma*Gamma(kk, w - 1);
              }
            }
            else
            {
              if(kk < K + 1)
              {
                b_temp = rnorm_trunc(0, tau_beta, -10, 10);
                hastings = (a_beta*log(b_beta) - lgamma(a_beta) + lgamma(a_beta + 0.5) - (a_beta + 0.5)*log(b_beta + b_temp*b_temp/2));
                hastings = hastings - (-log(tau_beta) - b_temp*b_temp/(2*tau_beta*tau_beta));
                hastings = hastings + e_delta;
                if(w < W - 1)
                {
                  hastings = hastings + f_delta*Gamma(kk, w + 1);
                }
                if(w > 0)
                {
                  hastings = hastings + f_delta*Gamma(kk, w - 1);
                }
              }
              else
              {
                b_temp = rnorm_trunc(0, tau_theta, -10, 10);
                hastings = (a_theta*log(b_theta) - lgamma(a_theta) + lgamma(a_theta + 0.5) - (a_theta + 0.5)*log(b_theta + b_temp*b_temp/2));
                hastings = hastings - (-log(tau_theta) - b_temp*b_temp/(2*tau_theta*tau_theta));
                hastings = hastings + e_xi;
                if(w < W - 1)
                {
                  hastings = hastings + f_xi*Gamma(kk, w + 1);
                }
                if(w > 0)
                {
                  hastings = hastings + f_xi*Gamma(kk, w - 1);
                }
              }
            }
          }
          else // Delete
          {
            if(kk == 0)
            {
              alpha_temp = 0;
              b_temp = alpha_0_ip(w) - alpha_0_input(w) + alpha_temp;
              hastings = -(a_alpha*log(b_alpha) - lgamma(a_alpha) + lgamma(a_alpha + 0.5) - (a_alpha + 0.5)*log(b_alpha + (B(kk, w) - alpha_0_ip(w) + alpha_0_input(w))*(B(kk, w) - alpha_0_ip(w) + alpha_0_input(w))/2));
              hastings = hastings + (-log(tau_alpha) - (B(kk, w) - alpha_0_ip(w) + alpha_0_input(w))*(B(kk, w) - alpha_0_ip(w) + alpha_0_input(w))/(2*tau_alpha*tau_alpha));
              hastings = hastings - e_gamma;
              if(w < W - 1)
              {
                hastings = hastings - f_gamma*Gamma(kk, w + 1);
              }
              if(w > 0)
              {
                hastings = hastings - f_gamma*Gamma(kk, w - 1);
              }
            }
            else
            {
              if(kk < K + 1)
              {
                b_temp = 0;
                hastings = -(a_beta*log(b_beta) - lgamma(a_beta) + lgamma(a_beta + 0.5) - (a_beta + 0.5)*log(b_beta + B(kk, w)*B(kk, w)/2));
                hastings = hastings + (-log(tau_beta) - B(kk, w)*B(kk, w)/(2*tau_beta*tau_beta));
                hastings = hastings - e_delta;
                if(w < W - 1)
                {
                  hastings = hastings - f_delta*Gamma(kk, w + 1);
                }
                if(w > 0)
                {
                  hastings = hastings - f_delta*Gamma(kk, w - 1);
                }
              }
              else
              {
                b_temp = 0;
                hastings = -(a_theta*log(b_theta) - lgamma(a_theta) + lgamma(a_theta + 0.5) - (a_theta + 0.5)*log(b_theta + B(kk, w)*B(kk, w)/2));
                hastings = hastings + (-log(tau_theta) - B(kk, w)*B(kk, w)/(2*tau_theta*tau_theta));
                hastings = hastings - e_xi;
                if(w < W - 1)
                {
                  hastings = hastings - f_xi*Gamma(kk, w + 1);
                }
                if(w > 0)
                {
                  hastings = hastings - f_xi*Gamma(kk, w - 1);
                }
              }
            }        
          }
          for(i = 0; i < n; i++)
          {
            if(X(i, kk) == 1 && H(i, w) == 0)
            {
              hastings = hastings + Y(i, w)*log(s(i)*exp(logLambda(i, w) - B(kk, w) + b_temp)) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w) - B(kk, w) + b_temp));
              hastings = hastings - (Y(i, w)*log(s(i)*exp(logLambda(i, w))) - (Y(i, w) + phi(w))*log(phi(w) + s(i)*exp(logLambda(i, w))));
            }
          }
          mh_try(4)++;
          if (hastings >= log(double(rand()%10001)/10000))
          {
            for(i = 0; i < n; i++)
            {
              if(X(i, kk) == 1)
              {
                logLambda(i, w) = logLambda(i, w) - B(kk, w) + b_temp;
              }
            }
            B(kk, w) = b_temp;
            Gamma(kk, w) = 1 - Gamma(kk, w);
            if(it >= burn) {
              mh_accept(4)++;
            }
          }
        }
      }
    }

    // Monitor the process
    if(it*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    
    pi_store(it) = pi;
    gamma_sum(it) = 0;
    Delta_sum(it) = 0;
    Xi_sum(it) = 0;
    if(it >= burn)
    {
      pi_mean = pi_mean + pi;
    }
    for(w = 0; w < W; w++)
    {
      if(flag(w) == 0)
      {
        for(i = 0; i < n; i++)
        {
          if(it >= burn && H(i, w) == 1)
          {
            H_ppi(i, w)++;
          }
        }
        if(it >= burn)
        {
          phi_mean(w) = phi_mean(w) + phi(w);
          alpha_0_ctrl_mean(w) = alpha_0_ctrl_mean(w) + alpha_0_input(w);
          alpha_0_ip_mean(w) = alpha_0_ip_mean(w) + alpha_0_ip(w);
        }
        if(Gamma(0, w) == 1)
        {
          gamma_sum(it)++;
          if(it >= burn)
          {
            gamma_ppi(w)++;
            alpha_mean(w) = alpha_mean(w) + B(0, w) - alpha_0_ip(w) + alpha_0_input(w);
          }
        }
        for(k = 0; k < K; k++)
        {
          if(Gamma(1 + k, w) == 1)
          {
            Delta_sum(it)++;
            if(it >= burn)
            {
              Beta_mean(k, w) = Beta_mean(k, w) + B(1 + k, w);
              Delta_ppi(k, w)++;
            }
          }
          if(Gamma(1 + K + k, w) == 1)
          {
            Xi_sum(it)++;
            if(it >= burn)
            {
              Theta_mean(k, w) = Theta_mean(k, w) + B(1 + K + k, w);
              Xi_ppi(k, w)++;
            }
          }
        }
      }
    }
  }
  
  for(i = 0; i < 5; i++)
  {
    if(mh_try(i) != 0)
    {
      mh_accept(i) = mh_accept(i)/mh_try(i);
    }
  }
  pi_mean = pi_mean/(iter - burn);
  for(w = 0; w < W; w++)
  {
    if(flag(w) == 0)
    {
      for(i = 0; i < n; i++)
      {
        H_ppi(i, w) = H_ppi(i, w)/(iter - burn);
      }
      phi_mean(w) = phi_mean(w)/(iter - burn);
      alpha_0_ctrl_mean(w) = alpha_0_ctrl_mean(w)/(iter - burn);
      alpha_0_ip_mean(w) = alpha_0_ip_mean(w)/(iter - burn);
      alpha_mean(w) = alpha_mean(w)/(iter - burn);
      gamma_ppi(w) = gamma_ppi(w)/(iter - burn);
      for(k = 0; k < K; k++)
      {
        Beta_mean(k, w) = Beta_mean(k, w)/(iter - burn);
        Delta_ppi(k, w) = Delta_ppi(k, w)/(iter - burn);
        Theta_mean(k, w) = Theta_mean(k, w)/(iter - burn);
        Xi_ppi(k, w) = Xi_ppi(k, w)/(iter - burn);
      }
    }
  }
  
  // return Rcpp::List::create(Rcpp::Named("flag") = flag, Rcpp::Named("logLambda") = logLambda, Rcpp::Named("B") = B, Rcpp::Named("s") = s_mean, Rcpp::Named("pi") = pi_store, Rcpp::Named("H") = H_ppi, Rcpp::Named("phi") = phi_mean, Rcpp::Named("alpha_0_input") = alpha_0_ctrl_mean, Rcpp::Named("alpha_0_ip") = alpha_0_ip_mean, Rcpp::Named("alpha") = alpha_mean, Rcpp::Named("gamma") = gamma_ppi, Rcpp::Named("gamma_sum") = gamma_sum, Rcpp::Named("Theta") = Theta_mean, Rcpp::Named("Xi") = Xi_ppi, Rcpp::Named("Xi_sum") = Xi_sum, Rcpp::Named("Beta") = Beta_mean, Rcpp::Named("Delta") = Delta_ppi, Rcpp::Named("Delta_sum") = Delta_sum, Rcpp::Named("accept") = mh_accept);
  return Rcpp::List::create(Rcpp::Named("flag") = flag, Rcpp::Named("pi") = pi_mean, Rcpp::Named("H") = H_ppi, Rcpp::Named("phi") = phi_mean, Rcpp::Named("alpha_0_input") = alpha_0_ctrl_mean, Rcpp::Named("alpha_0_ip") = alpha_0_ip_mean, Rcpp::Named("alpha") = alpha_mean, Rcpp::Named("gamma") = gamma_ppi, Rcpp::Named("Theta") = Theta_mean, Rcpp::Named("Xi") = Xi_ppi, Rcpp::Named("Beta") = Beta_mean, Rcpp::Named("Delta") = Delta_ppi, Rcpp::Named("accept") = mh_accept);
}

// [[Rcpp::export]]
double rnorm_trunc(double mu, double sigma, double lower, double upper)
{
  int change;
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;
  
  change = 0;
  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;
  
  // First scenario
  if( (a == R_NegInf)||(b == R_PosInf))
  {
    if(a == R_NegInf)
    {
      change = 1;
      a = -b;
      b = R_PosInf;
    }
    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    if(change) z = -z;
  }
  
  // Second scenario
  else if((a*b) <= 0.0)
  {
    // The two possibilities for this scenario
    if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
    {
      z = norm_rs(a, b);
    }
    else z = unif_rs(a,b);
  }
  
  // Third scenario
  else
  {
    if(b < 0)
    {
      tmp = b; b = -a; a = -tmp; change = 1;
    }
    
    lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
    if(lograt <= logt2)
    {
      z = unif_rs(a,b);
    }
    else if((lograt > logt1)&&(a < t3))
    {
      z = half_norm_rs(a,b);
    }
    else
    {
      z = exp_rs(a,b);
    }
    if(change)
    {
      z = -z;
    }
  }
  double output;
  output = sigma*z + mu;
  return (output);
}

// [[Rcpp::export]]
double exp_rs(double a, double b)
{
  double  z, u, rate;
  rate = 1/a;
  
  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a))
  {
    z = R::rexp(rate);
  }
  u = R::runif(0.0, 1.0);
  
  while( log(u) > (-0.5*z*z))
  {
    z = R::rexp(rate);
    while(z > (b-a))
    {
      z = R::rexp(rate);
    }
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}

// [[Rcpp::export]]
double unif_rs(double a, double b)
{
  double xstar, logphixstar, x, logu;
  
  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0) 
  {
    xstar = 0.0;
  }
  else 
  {
    xstar = a;
  }
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
  
  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while(logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
  {
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}

// [[Rcpp::export]]
double half_norm_rs(double a, double b)
{
  double x;
  x = fabs(norm_rand());
  while((x<a)||(x>b)) 
  {
    x = fabs(norm_rand());
  }
  return x;
}

// [[Rcpp::export]]
double norm_rs(double a, double b)
{
  double x;
  x = Rf_rnorm(0.0, 1.0);
  while((x < a)||(x > b)) 
  {
    x = norm_rand();
  }
  return x;
}
