data { // set up to run a single instance (1 stock) of GSI observations
  int N_samples; // number of samples
  int N_bins; // number of bins, dimensions of X
  matrix[N_samples, N_bins] X; // proportions
  int N_covar; // number of covariates in design matrix X
  matrix[N_samples, N_covar] design_X;
  array[N_bins,N_bins-1] int prod_idx;
  int overdisp; // whether or not to include overdispersion term
  int postpred; // whether or not to include posterior predictive samples
  real prior_sd;
  int tot_re; // total number of random intercepts
  int n_groups; // number of re groups
  matrix[N_samples, tot_re] design_Z; // design matrix of random ints
  int re_var_indx[tot_re + 1]; // total index across acll groups
  int n_re_by_group[n_groups + 1]; // indices by group
  int est_re; // 0 or 1 indicator
}
transformed data {
  array[N_samples,N_bins] int is_zero; // indicator for data being 1
  //array[N_samples,N_bins] int is_one; // indicator for data being 0
  array[N_samples,N_bins] int is_proportion; // indicator which elements to estimate
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

      if(X[i,j] < ESS[i] && X[i,j] > 0) {
        is_proportion[i,j] = 1;
      }
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
  vector<lower=0>[overdisp] phi_inv; // overdispersion
  matrix[N_bins-1,N_covar] beta_raw;
  matrix[(N_bins - 1) * est_re, tot_re * est_re] zeta_raw;
  vector<lower=0>[est_re * n_groups] zeta_sds; // sds of random intercepts
}
transformed parameters {
  real<lower=0> phi;
  matrix<lower=0,upper=1>[N_samples, N_bins] p_zero; // probability of 0 for each cell
  matrix<lower=0,upper=1>[N_samples, N_bins] p_one; // probability of 1 for each cell
  matrix[N_bins,N_covar] beta; // coefficients
  matrix[N_bins, tot_re * est_re] zeta; // coefficients
  matrix<lower=0,upper=1>[N_samples,N_bins] mu; // estimates, in normal space
  // phi is dynamic
  phi = 1;
  if(overdisp==1) {phi = 1/phi_inv[1];}
  // identifiability
  for (l in 1:N_covar) {
    beta[N_bins,l] = 0.0;
  }
  for (k in 1:(N_bins-1)) {
    for (l in 1:N_covar) {
      beta[k,l] = beta_raw[k,l];
    }
  }
  if(est_re == 1) {
    for (l in 1:tot_re) {
      zeta[N_bins,l] = 0.0;
    }
    for (k in 1:(N_bins-1)) {
      for (l in 1:tot_re) {
        zeta[k,l] = zeta_raw[k,l];
      }
    }
  }

  // from betas, we can calculate sample-specific mu
  if(est_re == 0) {
    for (n in 1:N_samples) {
      vector[N_bins] logits;
      for (m in 1:N_bins){
        logits[m] = design_X[n,] * transpose(beta[m,]);
      }
      logits = softmax(logits);
      for(m in 1:N_bins) {
        mu[n,m] = logits[m];
      }
    }
  } else {
    for (n in 1:N_samples) {
      vector[N_bins] logits;
      for (m in 1:N_bins){
        logits[m] = design_X[n,] * transpose(beta[m,]) + design_Z[n,] * transpose(zeta[m,]);
      }
      logits = softmax(logits);
      for(m in 1:N_bins) {
        mu[n,m] = logits[m];
      }
    }
  }

  for(i in 1:N_samples) {
    for(j in 1:N_bins) {
      p_zero[i,j] = (1-mu[i,j])^(ESS[i]*phi);
      // calculate probability = 1 as derived
      // = Pr(all other groups = 0 AND group of interest != 0)
    }
    // need 2 loops here because all quantites above used below
    for(j in 1:N_bins) {
      p_one[i,j] = (1 - p_zero[i,j]) * prod(p_zero[i, prod_idx[j,:]]);
    }
  }

}
model {
  real alpha_temp; // temp for extended beta
  real beta_temp; // temp for extended beta

  if(overdisp==1) {
    phi_inv ~ cauchy(0,5);
  }
  // priors for fixed effects for covariate factors
  for(i in 1:N_covar) {
    for(j in 1:(N_bins-1)) {
      beta_raw[j,i] ~ normal(0,prior_sd);
    }
  }
  //
  if(est_re==1) {
    for(i in 1:n_groups) {
      zeta_sds[i] ~ student_t(3,0,2); // priors on standard deviations
    }
    for(i in 1:(N_bins-1)) {
      for(j in 1:tot_re) {
        // normal deviations in mv logit space
        zeta_raw[i,j] ~ normal(0, zeta_sds[re_var_indx[j]]);
      }
    }
  }

  for(i in 1:N_samples) {
    for(j in 1:N_bins) {
      // marginals of the trinomial are independent binomials
      target += bernoulli_lpmf(is_zero[i,j]| p_zero[i,j]);

      if(is_proportion[i,j]==1) {
        alpha_temp = mu[i,j]*ESS[i]*phi;
        beta_temp = (1-mu[i,j])*ESS[i]*phi;
        // beta lpmf for 3-parameter model
        target += log(1 - p_zero[i,j] - p_one[i,j]) + (alpha_temp-1)*logX[i,j] + (beta_temp-1)*logNX[i,j] - (alpha_temp+beta_temp-1)*log(ESS[i]) - lbeta(alpha_temp,beta_temp);
      }
    }
  }

}
generated quantities {
  real alpha_temp;
  real beta_temp;
  array[N_samples] vector[N_bins] log_lik; // log likelihood
  array[N_samples*postpred] vector[N_bins*postpred] ynew; // new data, posterior predictive distribution
  array[N_samples*postpred,N_bins*postpred] int newy_is_zero;
  array[N_samples*postpred,N_bins*postpred] int newy_is_one;
  int newy_is_proportion;

  for(i in 1:N_samples) {
    for(j in 1:N_bins) {
      log_lik[i,j] = 0;
      log_lik[i,j] += bernoulli_lpmf(is_zero[i,j]| p_zero[i,j]);
      if(is_proportion[i,j]==1) {
        alpha_temp = mu[i,j]*ESS[i]*phi;
        beta_temp = (1-mu[i,j])*ESS[i]*phi;
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
          alpha_temp = mu[i,j]*ESS[i];
          beta_temp = (1-mu[i,j])*ESS[i];
          ynew[i,j] = beta_rng(alpha_temp, beta_temp);
          ynew[i,j] = ynew[i,j]*ESS[i]; // stretch to N
        }
      }
    }
  }
}
