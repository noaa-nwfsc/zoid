# Random generation of datasets using the dirichlet broken stick method
#'
#' Random generation of datasets using the dirichlet broken stick method
#'
#' @param n_obs Number of observations (rows of data matrix to simulate). Defaults to 10
#' @param n_groups Number of categories for each observation (columns of data matrix). Defaults to 10
#' @param ess_fraction The effective sample size fraction, defaults to 1
#' @param tot_n The total sample size to simulate for each observation. This is approximate and the actual
#' simulated sample size will be slightly smaller. Defaults to 100
#' @param p The stock proportions to simulate from, as a vector. Optional, and when not included,
#' random draws from the dirichlet are used
#' @return A 2-element list, whose 1st element `X_obs` is the simulated dataset, and whose
#' 2nd element is the underlying vector of proportions `p` used to generate the data
#' @export
#' @importFrom gtools rdirichlet
#' @importFrom stats rbeta rbinom
#'
#' @examples
#' \donttest{
#' y <- broken_stick(n_obs = 3, n_groups = 5, tot_n = 100)
#'
#' # add custom proportions
#' y <- broken_stick(
#'   n_obs = 3, n_groups = 5, tot_n = 100,
#'   p = c(0.1, 0.2, 0.3, 0.2, 0.2)
#' )
#' }
broken_stick <- function(n_obs = 1000,
                         n_groups = 10,
                         ess_fraction = 1,
                         tot_n = 100,
                         p = NULL) {
  if (is.null(p)) {
    p <- gtools::rdirichlet(1, rep(1, n_groups))
  }

  ess <- tot_n * ess_fraction
  # first, determine the presence of zeros using the stick breaking algorithm
  # second, for instances where zeros and tot_n are not observed, rescale the parameters and
  #         use the stick breaking algorithm for the Dirichlet process a second time

  X_mean_prob <- matrix(p, n_obs, n_groups, byrow = T)
  X_var_prob <- matrix((p * (1 - p)) / (ess + 1),
    n_obs, n_groups,
    byrow = T
  )
  X_mean <- X_mean_prob * tot_n
  X_var <- X_var_prob * tot_n^2

  # Simulate occurrence of 0s, 1s
  X_obs_0 <- matrix(NA, n_obs, n_groups) #1 indicates that the observation is a zero
  X_obs_1 <- matrix(NA, n_obs, n_groups)
  for(i in 1:n_obs){
    BREAK = "FALSE"
    repeat{
      for(j in 1:n_groups){
        X_obs_0[i,j] <- rbinom(1,1,(1-p[j])^ess)
      }
      if(sum(X_obs_0[i,])<n_groups){break}
    }
    for(j in 1:n_groups){
      X_obs_1[i,j] <- ifelse(X_obs_0[i,j]==0 & sum(X_obs_0[i,])==(n_groups-1),1,0)
    }
  }

  # Simulate proportions
  X_indicator <- matrix(NA,n_obs,n_groups)
  X_indicator[X_obs_0 == 0] <- 1
  X_indicator[X_obs_0 == 1] <- 0

  X_alpha <- (X_mean_prob) * ((X_mean_prob * (1-X_mean_prob)/ X_var_prob) - 1)
  X_alpha_mod <- X_alpha * 0

  mu_vals <- X_mean * 0 # These will be independent Beta draws, conditioned on being non-zero.
  q_vals <- mu_vals     # These will be equivalent to dirichlet draws (mu_vals modified by stick-breaking algorithm)

  X_obs <- X_indicator*0

  for(i in 1:n_obs){
    if(BREAK == "FALSE"){
      X_alpha_mod[i,] <- X_alpha[i,] * X_indicator[i,]
      for(j in 1:(n_groups-1)){ # Loop over stocks for dirichlet component
        if(j==1){
          if(X_alpha_mod[i,j]>0){
            mu_vals[i,j] <- rbeta(1,X_alpha_mod[i,j],sum(X_alpha_mod[i,(j+1):n_groups]))
            q_vals[i,j] = mu_vals[i,j]
          }
        }else if(j>1 ){
          if(X_alpha_mod[i,j]>0){
            mu_vals[i,j] <- rbeta(1,X_alpha_mod[i,j],sum(X_alpha_mod[i,(j+1):n_groups]))
            q_vals[i,j] <- prod(1 - mu_vals[i,(1:j-1)]) * mu_vals[i,j]
          }
        }
      }
      if(X_alpha_mod[i,n_groups]>0){
        q_vals[i,n_groups] <- 1 - sum(q_vals[i,(1:n_groups-1)])
      }
      X_obs[i,] <- q_vals[i,] * (tot_n)
    }
  }

  return(list(X_obs = X_obs, p = p))
}
