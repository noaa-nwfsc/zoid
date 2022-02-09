#' Extract estimates of predicted latent proportions.
#'
#' Extract point estimates of compositions from fitted model.
#'
#' @param fitted_model The fitted model returned as an rstan object from the call to zoid
#' @param conf_int Parameter controlling confidence intervals calculated, defaults to 0.05
#' for 95% intervals
#' @export
#' @return A list containing the posterior summaries of estimated parameters, with
#' element `mu` (the predicted values in normal space). For predictions
#' in transformed space, or overdispersion, see [get_pars()]
#' @importFrom rstan extract
#' @importFrom stats median quantile
#'
#' @examples
#' \donttest{
#' y <- matrix(c(3.77, 6.63, 2.60, 0.9, 1.44, 0.66, 2.10, 3.57, 1.33),
#'   nrow = 3, byrow = TRUE
#' )
#' # fit a model with no covariates
#' fit <- fit_zoid(data_matrix = y)
#' p_hat <- get_fitted(fit)
#' }
#'
get_fitted <- function(fitted_model, conf_int = 0.05) {
  if ("model" %in% names(fitted_model) == FALSE & class(fitted_model$model)[1] != "stanfit") {
    stop("Error: input isn't an stanfit object")
  }

  pars <- rstan::extract(fitted_model$model)
  n_obs <- dim(pars$mu)[2]
  n_group <- dim(pars$mu)[3]
  df <- expand.grid(
    "obs" = seq(1, n_obs),
    "group" = seq(1, n_group)
  )

  for (i in 1:nrow(df)) {
    df$mean[i] <- mean(pars$mu[, df$obs[i], df$group[i]])
    df$median[i] <- median(pars$mu[, df$obs[i], df$group[i]])
    df$lo[i] <- quantile(pars$mu[, df$obs[i], df$group[i]], conf_int / 2.0)
    df$hi[i] <- quantile(pars$mu[, df$obs[i], df$group[i]], 1 - conf_int / 2.0)
  }

  return(df)
}
