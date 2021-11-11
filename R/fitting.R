#' Fit a Bayesian Dirichlet regression model, allowing for zero-and-one inflation, covariates, and overdispersion
#'
#' Fit a Bayesian Dirichlet regression model that optionally includes covariates to estimate
#' effects of factor or continuous variables on proportions.
#'
#' @param formula The model formula for the design matrix. Does not need to have a response specified. If =NULL, then
#' the design matrix is ignored and all rows are treated as replicates
#' @param design_matrix A data frame, dimensioned as number of observations, and covariates in columns
#' @param data_matrix A matrix, with observations on rows and number of groups across columns
#' @param chains Number of mcmc chains, defaults to 3
#' @param iter Number of mcmc iterations, defaults to 2000
#' @param warmup Number iterations for mcmc warmup, defaults to 1/2 of the iterations
#' @param overdispersion Whether or not to include overdispersion parameter, defaults to FALSE
#' @param overdispersion_sd Prior standard deviation on 1/overdispersion parameter, Defaults to inv-Cauchy(0,5)
#' @param posterior_predict Whether or not to return draws from posterior predictive distribution (requires more memory)
#' @param moment_match Whether to do moment matching via [loo::loo_moment_match()]. This increases memory by adding all temporary
#' parmaeters to be saved and returned
#' @param prior_sd Parameter to be passed in to use as standard deviation of the normal distribution in transformed space. If
#' covariates are included this defaults to 1, but for models with single replicate, defaults to 1/n_bins.
#' @param ... Any other arguments to pass to [rstan::sampling()].
#'
#' @export
#' @return A list containing the fitted model and arguments and data used
#' to fit the model. These include `model` (the fitted model object of class `stanfit`),
#' `par_names` (the names of monitored parameters), `design_matrix` (the
#' design matrix of covariates), `data_matrix` (the data matrix of responses),
#' `overdispersion` (boolean, whether overdispersion was used),
#' `overdispersion_prior` (the prior used for overdispersion),
#' and `posterior_predict` (boolean, whether posterior prediction was done)
#'
#' @importFrom rstan sampling
#' @importFrom stats model.frame model.matrix rcauchy
#' @import Rcpp
#'
#' @examples
#' \donttest{
#' y <- matrix(c(3.77, 6.63, 2.60, 0.9, 1.44, 0.66, 2.10, 3.57, 1.33),
#'   nrow = 3, byrow = TRUE
#' )
#' # fit a model with no covariates
#' fit <- fit_zoid(data_matrix = y)
#'
#' # fit a model with 1 factor
#' design <- data.frame("y" = c(1, 1, 1), "fac" = c("spring", "spring", "fall"))
#' fit <- fit_zoid(formula = ~fac, design_matrix = design, data_matrix = y)
#' }
#'
fit_zoid <- function(formula = NULL,
                         design_matrix,
                         data_matrix,
                         chains = 3,
                         iter = 2000,
                         warmup = floor(iter / 2),
                         overdispersion = FALSE,
                         overdispersion_sd = 5,
                         posterior_predict = FALSE,
                         moment_match = FALSE,
                         prior_sd = NA,
                         ...) {

  # if a single observation
  if (class(data_matrix)[1] != "matrix") {
    data_matrix <- matrix(data_matrix, nrow = 1)
  }

  if (!is.null(formula)) {
    model_frame <- model.frame(formula, design_matrix)
    model_matrix <- model.matrix(formula, model_frame)
  } else {
    model_matrix <- matrix(1, nrow = nrow(data_matrix))
    colnames(model_matrix) <- "(Intercept)"
  }

  sd_prior <- 1 / ncol(data_matrix) # default if no covariates
  if (ncol(model_matrix) > 1) sd_prior <- 1
  if (!is.na(prior_sd)) sd_prior <- prior_sd

  par_names <- colnames(model_matrix)

  stan_data <- list(
    N_bins = ncol(data_matrix),
    N_samples = nrow(data_matrix),
    X = data_matrix,
    N_covar = ncol(model_matrix),
    design_X = model_matrix,
    overdisp = ifelse(overdispersion == TRUE, 1, 0),
    overdispersion_sd = overdispersion_sd,
    postpred = ifelse(posterior_predict == TRUE, 1, 0),
    prior_sd = sd_prior
  )

  pars <- c("beta", "log_lik", "mu")
  if (overdispersion == TRUE) pars <- c(pars, "phi")
  if (posterior_predict == TRUE) pars <- c(pars, "ynew")
  if (moment_match == TRUE) pars <- c(pars, "phi_inv", "raw_beta", "p_zero", "p_one")
  sampling_args <- list(
    object = stanmodels$dirichregmod,
    chains = chains,
    iter = iter,
    warmup = warmup,
    pars = pars,
    data = stan_data, ...
  )
  fit <- do.call(sampling, sampling_args)

  prior <- NULL
  if (overdispersion) {
    prior <- abs(rcauchy(n = chains * (iter - warmup), location = 0, scale = overdispersion_sd))
  }
  return(list(
    model = fit, par_names = par_names,
    design_matrix = model_matrix,
    data_matrix = data_matrix,
    overdispersion = overdispersion,
    overdispersion_prior = prior,
    posterior_predict = posterior_predict
  ))
}
