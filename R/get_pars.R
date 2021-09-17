#' Extract parameters from fitted model.
#'
#' Extract estimated parameters from fitted model.
#'
#' @param fitted_model The fitted model returned as an rstan object from the call to trinomix
#' @param conf_int Parameter controlling confidence intervals calculated, defaults to 0.05
#' for 95% intervals
#' @export
#' @importFrom rstan extract
#'
#' @examples
#' \donttest{
#' y = matrix(c(3.77,6.63,2.60,0.9,1.44,0.66,2.10,3.57,1.33),
#' nrow=3, byrow = TRUE)
#' # fit a model with no covariates
#' fit <- fit_trinomix(data_matrix = y)
#' p_hat = get_pars(fit)
#'}
#'
get_pars = function(fitted_model, conf_int = 0.05) {
  if("model" %in% names(fitted_model) == FALSE & class(fitted_model$model)[1] != "stanfit") {
    stop("Error: input isn't an stanfit object")
  }

  p = get_fitted(fitted_model, conf_int = conf_int)

  # add on other parameters
  pars = rstan::extract(fitted_model$model)
  n_group = dim(pars$beta)[2]
  n_cov = dim(pars$beta)[3]
  betas = expand.grid("group" = seq(1,n_group),
                   "cov" = seq(1,n_cov),
                   "par" = NA,
                   "mean" = NA,
                   "median" = NA,
                   "lo" = NA,
                   "hi" = NA)

  for(i in 1:nrow(betas)) {
    betas$mean[i] = mean(pars$beta[,betas$group[i],betas$cov[i]])
    betas$median[i] = median(pars$beta[,betas$group[i],betas$cov[i]])
    betas$lo[i] = quantile(pars$beta[,betas$group[i],betas$cov[i]],conf_int/2.0)
    betas$hi[i] = quantile(pars$beta[,betas$group[i],betas$cov[i]],1 - conf_int/2.0)
    betas$par[i] = fitted_model$par_names[betas$cov[i]]
  }

  par_list = list(p = p, betas = betas)
  if(fitted_model$overdispersion == TRUE) {
    theta = data.frame(
      "mean" = mean(pars$theta),
      "median" = median(pars$theta),
      "lo" = quantile(pars$theta, conf_int/2.0),
      "hi" = quantile(pars$theta, 1 - conf_int/2.0)
    )
    par_list$theta = theta
  }

  return(par_list)
}
