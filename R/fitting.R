#' Fit a trinomial mixture model with Stan
#'
#' Fit a trinomial mixture model that optionally includes covariates to estimate
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
#' @param posterior_predict Whether or not to return draws from posterior predictive distribution (requires more memory)
#' @param ... Any other arguments to pass to [rstan::sampling()].
#'
#' @export
#' @importFrom rstan sampling
#' @importFrom stats model.frame model.matrix
#' @import Rcpp
#'
#' @examples
#' \donttest{
#' y <- matrix(c(3.77, 6.63, 2.60, 0.9, 1.44, 0.66, 2.10, 3.57, 1.33),
#'   nrow = 3, byrow = TRUE
#' )
#' # fit a model with no covariates
#' fit <- fit_trinomix(data_matrix = y)
#'
#' # fit a model with 1 factor
#' design <- data.frame("y" = c(1, 1, 1), "fac" = c("spring", "spring", "fall"))
#' fit <- fit_trinomix(formula = ~fac, design_matrix = design, data_matrix = y)
#' }
#'
fit_trinomix <- function(formula = NULL,
                         design_matrix,
                         data_matrix,
                         chains = 3,
                         iter = 2000,
                         warmup = floor(iter / 2),
                         overdispersion = FALSE,
                         posterior_predict = FALSE,
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

  par_names <- colnames(model_matrix)

  stan_data <- list(
    N_stocks = ncol(data_matrix),
    N_samples = nrow(data_matrix),
    X = data_matrix,
    N_covar = ncol(model_matrix),
    design_X = model_matrix,
    overdisp = ifelse(overdispersion == TRUE, 1, 0),
    postpred = ifelse(posterior_predict == TRUE, 1, 0)
  )

  pars <- c("beta", "log_lik", "mu")
  if (overdispersion == TRUE) pars <- c(pars, "theta")
  if (posterior_predict == TRUE) pars <- c(pars, "ynew")

  sampling_args <- list(
    object = stanmodels$dirichreg,
    chains = chains,
    iter = iter,
    warmup = warmup,
    pars = pars,
    data = stan_data, ...
  )
  fit <- do.call(sampling, sampling_args)

  return(list(
    model = fit, par_names = par_names,
    design_matrix = model_matrix,
    data_matrix = data_matrix,
    overdispersion = overdispersion,
    posterior_predict = posterior_predict
  ))
}
