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
#' @param overdispersion_sd Prior standard deviation on 1/overdispersion parameter, Defaults to inv-Cauchy(0,5)
#' @param posterior_predict Whether or not to return draws from posterior predictive distribution (requires more memory)
#' @param moment_match Whether to do moment matching via [loo::loo_moment_match()]. This increases memory by adding all temporary
#' parmaeters to be saved and returned
#' @param prior_sd Parameter to be passed in to use as standard deviation of the normal distribution in transformed space. If
#' covariates are included this defaults to 1, but for models with single replicate, defaults to 1/n_bins.
#' @param ... Any other arguments to pass to [rstan::sampling()].
#'
#' @export
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
#' fit <- fit_zoid(data_matrix = y, chains = 1, iter = 100)
#'
#' # fit a model with 1 factor
#' design <- data.frame("fac" = c("spring", "spring", "fall"))
#' fit <- fit_zoid(formula = ~fac, design_matrix = design, data_matrix = y, chains = 1, iter = 100)
#' }
#' # try a model with random effects
#' set.seed(123)
#' y <- matrix(runif(99,1,4), ncol=3)
#' design <- data.frame("fac" = sample(letters[1:5], size=nrow(y), replace=TRUE))
#' design$fac <- as.factor(design$fac)
#' fit <- fit_zoid(formula = ~(1|fac), design_matrix = design, data_matrix = y, chains = 1, iter = 100)
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

  # fill with dummy values
  parsed_res <- list(design_matrix = matrix(0, nrow(data_matrix),ncol=1),
                     var_indx = 1,
                     n_re_by_group = 1,
                     tot_re = 1,
                     n_groups = 1)
  est_re <- FALSE
  re_group_names <- NA
  if (!is.null(formula)) {
    model_frame <- model.frame(formula, design_matrix)
    model_matrix <- model.matrix(formula, model_frame)
    # extract the random effects
    res <- parse_re_formula(formula, design_matrix)
    if(length(res$var_indx) > 0) {
      parsed_res <- res # only update if REs are in formula
      est_re <- TRUE
      model_matrix <- res$fixed_design_matrix
      re_group_names <- res$random_effect_group_names
    }
  } else {
    model_matrix <- matrix(1, nrow = nrow(data_matrix))
    colnames(model_matrix) <- "(Intercept)"
  }

  sd_prior <- 1 / ncol(data_matrix) # default if no covariates
  if (ncol(model_matrix) > 1) sd_prior <- 1
  if (!is.na(prior_sd)) sd_prior <- prior_sd

  par_names <- colnames(model_matrix)

  prod_idx <- matrix(0, ncol(data_matrix), ncol(data_matrix)-1)
  for(j in 1:ncol(data_matrix)){
    prod_idx[j,] <- seq(1,ncol(data_matrix),1)[-j]
  }

  stan_data <- list(
    N_bins = ncol(data_matrix),
    N_samples = nrow(data_matrix),
    X = data_matrix,
    prod_idx = prod_idx,
    N_covar = ncol(model_matrix),
    design_X = model_matrix,
    overdisp = ifelse(overdispersion == TRUE, 1, 0),
    overdispersion_sd = overdispersion_sd,
    postpred = ifelse(posterior_predict == TRUE, 1, 0),
    prior_sd = sd_prior,
    design_Z = parsed_res$design_matrix, # design matrix for Z (random int)
    re_var_indx = c(parsed_res$var_indx, 1), # index of the group for each re
    n_re_by_group = c(parsed_res$n_re_by_group, 1), # number of random ints per group
    tot_re = parsed_res$tot_re, # total number of random ints, across all groups
    n_groups = parsed_res$n_group,
    est_re = as.numeric(est_re)
  )

  pars <- c("beta", "log_lik", "mu")
  if(est_re == TRUE) pars <- c(pars, "zeta", "zeta_sds")
  if (overdispersion == TRUE) pars <- c(pars, "phi")
  if (posterior_predict == TRUE) pars <- c(pars, "ynew")
  if (moment_match == TRUE) pars <- c(pars, "phi_inv", "beta_raw", "p_zero", "p_one")
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
    posterior_predict = posterior_predict,
    stan_data = stan_data,
    re_group_names = re_group_names
  ))
}


#' Fit a trinomial mixture model that optionally includes covariates to estimate
#' effects of factor or continuous variables on proportions.
#'
#' @param formula The model formula for the design matrix.
#' @param data The data matrix used to construct RE design matrix
#' @importFrom stats model.matrix as.formula
parse_re_formula <- function(formula, data) {
  # Convert the formula to a character string
  formula_str <- as.character(formula)
  # Split the formula into parts based on '+' and '-' symbols
  formula_parts <- unlist(strsplit(formula_str, split = "[-+]", perl = TRUE))
  # Trim whitespace from each part
  formula_parts <- trimws(formula_parts)
  # Identify parts containing a bar '|'
  random_effects <- grep("\\|", formula_parts, value = TRUE)
  fixed_effects <- setdiff(formula_parts, random_effects)

  # Create design matrix for fixed effects. Catch the cases where no fixed
  # effects are included, or intercept-only models used
  if (length(fixed_effects) > 1 || (length(fixed_effects) == 1 && fixed_effects != "~")) {
    fixed_formula_str <- paste("~", paste(fixed_effects, collapse = "+"))
  } else {
    fixed_formula_str <- "~ 1" # Only intercept
  }
  fixed_design_matrix <- model.matrix(as.formula(fixed_formula_str), data)

  random_effect_group_names <- sapply(random_effects, function(part) {
    # Extract the part after the '|'
    split_part <- strsplit(part, "\\|", perl = TRUE)[[1]]
    # Remove the closing parenthesis and trim
    group_name <- gsub("\\)", "", split_part[2])
    trimws(group_name)
  })

  # create design matrices by group
  for(i in 1:length(random_effects)) {
    new_formula <- as.formula(paste("~", random_effect_group_names[i], "-1"))
    if(i ==1) {
      design_matrix <- model.matrix(new_formula, data)
      var_indx <- rep(1, ncol(design_matrix))
      n_re <- length(var_indx)
    } else {
      design_matrix <- cbind(design_matrix, model.matrix(new_formula, data))
      var_indx <- c(var_indx, rep(i, ncol(design_matrix)))
      n_re <- c(n_re, length(ncol(design_matrix)))
    }
  }
  n_groups <- 0
  if(length(var_indx) > 0) n_groups <- max(var_indx)
  return(list(design_matrix = design_matrix, var_indx = var_indx, n_re_by_group = n_re,
              tot_re = sum(n_re), n_groups = n_groups,
              fixed_design_matrix = fixed_design_matrix,
              random_effect_group_names = random_effect_group_names))
}
