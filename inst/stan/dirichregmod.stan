data { // set up to run a single instance (1 stock) of GSI observations
  int N_samples; // number of samples
  int N_bins; // number of bins, dimensions of X
  matrix[N_samples, N_bins] X; // proportions
  int N_covar; // number of covariates in design matrix X
  matrix[N_samples, N_covar] design_X_mu;
  matrix[N_samples, N_covar] design_X_sigma;
  matrix[N_samples, N_covar] design_X_zero;
  matrix[N_samples, N_covar] design_X_ones;
  int overdisp; // whether or not to include overdispersion term
  int postpred; // whether or not to include posterior predictive samples
  real prior_sd;
}
transformed data {
  int is_zero[N_samples,N_bins]; // indicator for data being 1
  int is_one[N_samples,N_bins]; // indicator for data being 0
  int is_proportion[N_samples,N_bins]; // indicator which elements to estimate
  matrix[N_samples, N_bins] logX; // log proportions
  matrix[N_samples, N_bins] logNX; // log proportions
  vector[N_samples] ESS;
  vector[N_bins] ones;

  // for dirichlet
  for(i in 1:N_bins) {
    ones[i] = 1;
  }
  //int n_est;
  for(i in 1:N_samples) {
    ESS[i] = 0;
    for(j in 1:N_bins) {
      ESS[i] = ESS[i] + X[i,j]; // effective sample size for sample
    }
  }
  for(i in 1:N_samples) {
    for(j in 1:N_bins) {
      if(X[i,j]==0) {
        is_zero[i,j] = 1;
      } else {
        is_zero[i,j] = 0;
      }

      if(X[i,j]==ESS[i]) {
        is_one[i,j] = 1;
      } else {
        is_one[i,j] = 0;
      }

      is_proportion[i,j] = (1-is_zero[i,j]) * (1-is_one[i,j]);
    }

  }

    // do some more book keeping
    for(i in 1:N_samples) {
      for(j in 1:N_bins) {
        if(is_proportion[i,j]==1) {
          logX[i,j] = log(X[i,j]);
          logNX[i,j] = log(ESS[i] - X[i,j]);
        }
      }
    }

}
parameters {
  vector[overdisp] phi_inv; // overdispersion
  matrix[N_bins-1,N_covar] beta_raw_mu;
  matrix[N_bins-1,N_covar] beta_raw_sigma;
  matrix[N_bins-1,N_covar] beta_raw_ones;
  matrix[N_bins-1,N_covar] beta_raw_zero;
}
transformed parameters {
  real phi;
  matrix<lower=0,upper=1>[N_samples, N_bins] p_zero; // probability of 0 for each cell
  matrix<lower=0,upper=1>[N_samples, N_bins] p_one; // probability of 1 for each cell
  matrix[N_bins,N_covar] beta_mu; // coefficients
  matrix[N_bins,N_covar] beta_sigma; // coefficients
  matrix[N_bins,N_covar] beta_ones; // coefficients
  matrix[N_bins,N_covar] beta_zero; // coefficients

  matrix<lower=0,upper=1>[N_samples,N_bins] mu; // estimates, in normal space
  matrix<lower=0,upper=1>[N_samples,N_bins] sigma; // estimates, in normal space
  matrix<lower=0,upper=1>[N_samples,N_bins] sigma2; // estimates, in normal space

  // phi is dynamic
  phi = 1;
  if(overdisp==1) {phi = 1/phi_inv[1];}

  for (l in 1:N_covar) {
    beta_mu[N_bins,l] = 0.0;
    beta_sigma[N_bins,l] = 0.0;
    beta_ones[N_bins,l] = 0.0;
    beta_zero[N_bins,l] = 0.0;
  }
  for (k in 1:(N_bins-1)) {
    for (l in 1:N_covar) {
      beta_mu[k,l] = beta_raw_mu[k,l];
      beta_sigma[k,l] = beta_raw_sigma[k,l];
      beta_ones[k,l] = beta_raw_ones[k,l];
      beta_zero[k,l] = beta_raw_zero[k,l];
    }
  }

  // from betas, we can calculate sample-specific mu
  for (n in 1:N_samples) {
    vector[N_bins] logits_mu;
    vector[N_bins] logits_sigma;
    vector[N_bins] logits_ones;
    vector[N_bins] logits_zero;
    for (m in 1:N_bins){
      logits_mu[m] = design_X_mu[n,] * transpose(beta_mu[m,]);
      logits_zero[m] = design_X_zero[n,] * transpose(beta_zero[m,]);
      logits_ones[m] = design_X_ones[n,] * transpose(beta_ones[m,]);
      logits_sigma[m] = design_X_sigma[n,] * transpose(beta_sigma[m,]);
    }
    logits_mu = softmax(logits_mu);
    logits_ones = softmax(logits_ones);
    logits_zero = softmax(logits_zero);
    //logits_eta = softmax(logits_eta);
    for(m in 1:N_bins) {
      mu[n,m] = logits_mu[m];
      p_zero[n,m] = logits_zero[m];
      p_one[n,m] = logits_ones[m];
      sigma[n,m] = exp(logits_sigma[m]);
      sigma2[n,m] = sigma[n,m]*sigma[n,m];
    }
  }

}
model {
  real alpha_temp; // temp for extended beta
  real beta_temp; // temp for extended beta
  real eta_temp; // temp for extended beta
  if(overdisp==1) {
    phi_inv ~ cauchy(0,5);
  }
  // priors for fixed effects for covariate factors
  for(i in 1:N_covar) {
    for(j in 1:(N_bins-1)) {
      beta_raw_mu[j,i] ~ normal(0,prior_sd);
      beta_raw_sigma[j,i] ~ normal(0,prior_sd);
      beta_raw_ones[j,i] ~ normal(0,prior_sd);
      beta_raw_zero[j,i] ~ normal(0,prior_sd);
    }
  }

  for(i in 1:N_samples) {
    for(j in 1:N_bins) {
      // marginals of the trinomial are independent binomials
      target += bernoulli_lpmf(is_zero[i,j]| p_zero[i,j]);
      target += bernoulli_lpmf(is_one[i,j]| p_one[i,j]);

      if(is_proportion[i,j]==1) {
        eta_temp = ((1 - mu[i,j])/sigma2[i,j] - 1/mu[i,j]);
        alpha_temp = eta_temp * mu[i,j] * mu[i,j];
        beta_temp = eta_temp * mu[i,j] * (1 - mu[i,j]);
        // beta lpmf for 3-parameter model
        target += log(1 - p_zero[i,j] - p_one[i,j]) + (alpha_temp-1)*logX[i,j] + (beta_temp-1)*logNX[i,j] - (alpha_temp+beta_temp-1)*log(ESS[i]) - lbeta(alpha_temp,beta_temp);
      }
    }
  }

}
generated quantities {
  real alpha_temp;
  real beta_temp;
  real eta_temp;
  vector[N_bins] log_lik[N_samples]; // log likelihood
  vector[N_bins*postpred] ynew[N_samples*postpred]; // new data, posterior predictive distribution
  int newy_is_zero[N_samples*postpred,N_bins*postpred];
  int newy_is_one[N_samples*postpred,N_bins*postpred];
  int newy_is_proportion;

  for(i in 1:N_samples) {
    for(j in 1:N_bins) {
      log_lik[i,j] = 0;
      log_lik[i,j] += bernoulli_lpmf(is_zero[i,j]| p_zero[i,j]);
      log_lik[i,j] += bernoulli_lpmf(is_one[i,j]| p_one[i,j]);
      if(is_proportion[i,j]==1) {
        //alpha_temp = mu[i,j]*ESS[i]*phi;
        //beta_temp = (1-mu[i,j])*ESS[i]*phi;
        eta_temp = ((1 - mu[i,j])/sigma2[i,j] - 1/mu[i,j]);
        alpha_temp = eta_temp * mu[i,j] * mu[i,j] * (ESS[i]*phi);
        beta_temp = eta_temp * mu[i,j] * (1 - mu[i,j]) * (ESS[i]*phi);
        log_lik[i,j] += log(1 - p_zero[i,j] - p_one[i,j]) + (alpha_temp-1)*logX[i,j] + (beta_temp-1)*logNX[i,j] - (alpha_temp+beta_temp-1)*log(ESS[i]) - lbeta(alpha_temp,beta_temp);
      }

      // posterior predictive sampling
      if(postpred==1) {
        newy_is_zero[i,j] = bernoulli_rng(p_zero[i,j]);
        if(newy_is_zero[i,j]==1) ynew[i,j] = 0.0;
        newy_is_one[i,j] = bernoulli_rng(p_one[i,j]);
        if(newy_is_one[i,j]==1) ynew[i,j] = 1.0;
        newy_is_proportion = newy_is_zero[i,j] + newy_is_one[i,j];
        if(newy_is_proportion==0) {
          eta_temp = ((1 - mu[i,j])/sigma2[i,j] - 1/mu[i,j]);
          alpha_temp = eta_temp * mu[i,j] * mu[i,j] * (ESS[i]*phi);
          beta_temp = eta_temp * mu[i,j] * (1 - mu[i,j]) * (ESS[i]*phi);
          ynew[i,j] = beta_rng(alpha_temp, beta_temp);
          ynew[i,j] = ynew[i,j]*ESS[i]; // stretch to N
        }
      }
    }
  }
}
