if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 200
CHAINS <- 1
SEED <- 1234
TOL <- 0.1 # %

# ------------------------------------------------------
# a basic fit


test_that("basic fit", {
  skip_on_cran()
  set.seed(SEED)

  y <- matrix(c(3.77, 6.63, 2.60, 0.9, 1.44, 0.66, 2.10, 3.57, 1.33),
    nrow = 3, byrow = TRUE
  )
  fit <- fit_zoid(data_matrix = y, chains = CHAINS)
  pars <- rstan::extract(fit$model)
  expect_equal(apply(pars$mu, c(2, 3), mean)[1, 1], 0.3052632, tolerance = TOL)
  expect_equal(apply(pars$mu, c(2, 3), mean)[1, 2], 0.4357927, tolerance = TOL)
  expect_equal(apply(pars$mu, c(2, 3), mean)[1, 3], 0.2589441, tolerance = TOL)
})
