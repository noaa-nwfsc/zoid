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
#'
#' @export
#' @importFrom gtools rdirichlet
#' @importFrom stats median pbeta qbeta quantile rbeta
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
  ESS <- tot_n * ess_fraction
  # This is the mean and variable of the beta component 0 < X < N
  X_mean <- matrix(p, n_obs, n_groups, byrow = T) * tot_n
  X_var <- matrix((p * (1 - p) * tot_n^2) / (ESS + 1),
    n_obs, n_groups,
    byrow = T
  )

  # Scale out the sample size to move to proportion space.
  X_mean_prob <- X_mean / tot_n
  X_var_prob <- X_var / tot_n^2

  # calculate the probability of 0 and probability of 1 first
  p_zero <- (1 - X_mean_prob)^ESS
  p_one <- X_mean_prob^ESS

  # Unconditional mean calculation
  UNCOND.MEAN <- p_one[1, ] * tot_n + (1 - p_zero[1, ] - p_one[1, ]) * X_mean[1, ]
  COND.MEAN <- X_mean[1, ]

  # Method of moments to calculated alpha from the dirichlet
  X_alpha <- (X_mean_prob) * ((X_mean_prob * (1 - X_mean_prob) / X_var_prob) - 1)
  # Calculate the betas for the marginal Beta distribution for later use
  X_beta <- (1 - X_mean_prob) * ((X_mean_prob * (1 - X_mean_prob) / X_var_prob) - 1)

  # use the stick breaking algorithm for the Dirichlet process
  # to generate the values that can be broken into the three categories.

  mu_vals <- X_alpha * 0 # These will be independent Beta draws.
  q_vals <- mu_vals # These will be equivalent to dirichlet draws (mu_vals modified by stick-breaking algorithm)
  for (i in 1:n_obs) {
    for (j in 1:(n_groups - 1)) {
      mu_vals[i, j] <- rbeta(1, X_alpha[i, j], sum(X_alpha[i, (j + 1):n_groups]))
      if (j == 1) {
        q_vals[i, j] <- mu_vals[i, j]
      } else if (j > 1) {
        q_vals[i, j] <- prod(1 - mu_vals[i, (1:j - 1)]) * mu_vals[i, j]
      }
    }
    q_vals[i, n_groups] <- 1 - sum(q_vals[i, (1:n_groups - 1)])
  }

  # Calculate the quantile of each of the q_vals from their respective marginal Betas
  X_cdf <- pbeta(q_vals, X_alpha, X_beta)

  # Define the replicates where all of the observations are in one category
  X_cdf_N <- X_cdf * 0
  X_cdf_N[X_cdf > (1 - p_one)] <- 1
  ROWS <- apply(X_cdf_N, 1, max)
  THESE.rows <- which(ROWS == 1)

  # stretch out the middle section so that the proportion between p(x=zero) and p(x=tot_n)
  # is on the interval 0,1
  X_cdf_mod <- (X_cdf - p_zero) / (1 - p_zero - p_one)

  # Deal with extreme underflow for p(x==0) is very large and X==N cases
  X_cdf_mod[p_zero == 1] <- -99
  X_cdf_mod[THESE.rows, ] <- X_cdf_N[THESE.rows, ]
  X_cdf_mod[X_cdf_mod < 0] <- 0

  # Simulate
  THESE <- which(X_cdf_mod > 0 & X_cdf_mod < 1)
  X_obs <- X_cdf * 0
  X_obs[THESE] <- qbeta(X_cdf_mod[THESE], X_alpha[THESE], X_beta[THESE]) * tot_n
  X_obs[THESE.rows, ] <- X_cdf_mod[THESE.rows, ] * tot_n

  return(list(X_obs = X_obs, p = p))
}
