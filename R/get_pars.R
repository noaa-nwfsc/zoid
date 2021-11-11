#' Extract parameters from fitted model.
#'
#' Extract estimated parameters from fitted model.
#'
#' @param fitted_model The fitted model returned as an rstan object from the call to zoid
#' @param conf_int Parameter controlling confidence intervals calculated, defaults to 0.05
#' for 95% intervals
#' @export
#' @return A list containing the posterior summaries of estimated parameters. At minimum,
#' this will include `p` (the estimated proportions) and `betas` (the predicted values in
#' transformed space). For models with overdispersion, an extra
#' element `phi` will also be returned, summarizing overdispersion. For predictions
#' in normal space, see [get_fitted()]
#' @importFrom rstan extract
#'
#' @examples
#' \donttest{
#' y <- matrix(c(3.77, 6.63, 2.60, 0.9, 1.44, 0.66, 2.10, 3.57, 1.33),
#'   nrow = 3, byrow = TRUE
#' )
#' # fit a model with no covariates
#' fit <- fit_zoid(data_matrix = y)
#' p_hat <- get_pars(fit)
#' }
#'
get_pars <- function(fitted_model, conf_int = 0.05) {
  if ("model" %in% names(fitted_model) == FALSE & class(fitted_model$model)[1] != "stanfit") {
    stop("Error: input isn't an stanfit object")
  }

  p <- get_fitted(fitted_model, conf_int = conf_int)

  # add on other parameters
  pars <- rstan::extract(fitted_model$model)
  n_group <- dim(pars$beta)[2]
  n_cov <- dim(pars$beta)[3]
  betas <- expand.grid(
    "group" = seq(1, n_group),
    "cov" = seq(1, n_cov),
    "par" = NA,
    "mean" = NA,
    "median" = NA,
    "lo" = NA,
    "hi" = NA
  )

  for (i in 1:nrow(betas)) {
    betas$mean[i] <- mean(pars$beta[, betas$group[i], betas$cov[i]])
    betas$median[i] <- median(pars$beta[, betas$group[i], betas$cov[i]])
    betas$lo[i] <- quantile(pars$beta[, betas$group[i], betas$cov[i]], conf_int / 2.0)
    betas$hi[i] <- quantile(pars$beta[, betas$group[i], betas$cov[i]], 1 - conf_int / 2.0)
    betas$par[i] <- fitted_model$par_names[betas$cov[i]]
  }

  par_list <- list(p = p, betas = betas)
  if (fitted_model$overdispersion == TRUE) {
    phi <- data.frame(
      "mean" = mean(pars$phi),
      "median" = median(pars$phi),
      "lo" = quantile(pars$phi, conf_int / 2.0),
      "hi" = quantile(pars$phi, 1 - conf_int / 2.0)
    )
    par_list$phi <- phi
  }

  return(par_list)
}
