#' Find appropriate standard deviations for prior
#'
#' @param n_bins Bins for the Dirichlet distribution
#' @param n_draws Numbers of samples to use for doing calculation
#' @param target The goal of the specified prior, e.g. 1 or 1/n_bins
#' @param iterations to try, to ensure robust solution. Defaults to 5
#' @export
#' @importFrom stats optim
#' @return A 3-element list consisting of `sd` (the approximate standard deviation
#' in transformed space that gives a similar prior to that specified), `value` (the
#' value of the root mean squared percent error function being minimized),
#' and `convergence` (0 if convergence occurred, error code from
#' [optim()] otherwise)
#' @examples
#' \donttest{
#' # fit model with 3 components / alpha = 1
#' set.seed(123)
#' f <- fit_prior(n_bins = 3, n_draws = 1000, target = 1)
#' # fit model with 20 components / alpha = 1/20
#' f <- fit_prior(n_bins = 20, n_draws = 1000, target = 1 / 20)
#' }
#'
fit_prior <- function(n_bins, n_draws = 10000, target = 1 / n_bins, iterations = 5) {
  best <- 1.0e10
  best_value <- NA
  for (i in 1:iterations) {
    o <- try(optim(
      par = runif(1), rmspe_calc, n_bins = n_bins,
      target = target,
      n_draws = n_draws,
      method = "BFGS"
    ), silent = TRUE)
    if (class(o) != "try-error") {
      if (o$value < best) {
        best <- o$value
        best_value <- list(
          sd = exp(o$par), value = o$value,
          convergence = o$convergence
        )
      }
    }
  }
  return(best_value)
}


#' Find appropriate prior for a given target distribution.
#'
#' Extract point estimates of compositions from fitted model.
#'
#' @param par The parameter (standard deviation) to be searched over to find a Dirichlet equivalent
#' @param n_bins Bins for the Dirichlet distribution
#' @param n_draws Numbers of samples to use for doing calculation
#' @param target The goal of the specified prior, e.g. 1 or 1/n_bins
#' @importFrom stats rnorm
#' @importFrom compositions fitDirichlet
rmspe_calc <- function(par, n_bins, n_draws, target) {
  x <- matrix(rnorm(n_draws * (n_bins - 1), 0, exp(par)), n_draws, n_bins - 1)
  x <- cbind(x, 0)
  f <- function(x) {
    return(exp(x) / sum(exp(x)))
  }
  p <- t(apply(x, 1, f))

  funct_alpha <- fitDirichlet(p)$alpha
  rmspe <- sqrt(mean(((funct_alpha - target) / target)^2))
  return(rmspe)
}
