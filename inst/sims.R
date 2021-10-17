library(MCMCpack)
library(dplyr)
library(ggplot2)
remotes::install_github("nwfsc-cb/trinomix", "priors")
library(trinomix)

df <- expand.grid(
  n_bins = 10,
  seed = 1:100,
  N = c(10,30,50)
)

for (i in 1:nrow(df)) {
  # simulate data
  set.seed(df$seed[i])
  p <- c(rdirichlet(1, alpha = rep(1, df$n_bins[i])))
  dat <- data.frame("true" = p * df$N[i])

  # fit model with beta prior
  beta_fit <- fit_trinomix(
    data_matrix = t(dat), iter = 4000,
    chains = 1, refresh = 0
  )
  # fit model with improper prior
  imp_fit <- fit_trinomix(
    data_matrix = t(dat), iter = 4000,
    chains = 1, use_beta_priors = FALSE, refresh = 0
  )

  dat$bin <- seq(1, df$n_bins[i])
  dat$p <- p
  dat$N <- df$N[i]
  dat <- rbind(dat, dat)
  dat$est <- c(get_fitted(beta_fit)$mean, get_fitted(beta_fit)$mean)
  dat$prior <- c(rep("beta", df$n_bins[i]), rep("imp", df$n_bins[i]))
  dat$seed <- df$seed[i]
  if (i == 1) {
    out <- dat
  } else {
    out <- rbind(out, dat)
  }
  print(i)
}

# is there consistent bias with either prior
wide <- tidyr::pivot_wider(out, names_from = prior, values_from = est)
ggplot(dplyr::filter(wide,N==10), aes(beta,imp)) + geom_point() +
  xlab("Estimate using beta prior") +
  ylab("Estimate using improper prior")

cor(wide$beta,wide$imp)
