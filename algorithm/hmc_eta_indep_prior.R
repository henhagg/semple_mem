library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

cpp_code = stanc('algorithm/ng_indep_prior_normal_pop.stan')
stanmodel = stan_model(model_code = cpp_code$model_code)

sample_eta_with_nuts = function(c_matrix,
                                prior_alpha,
                                prior_beta,
                                prior_mean_mu,
                                prior_sd_mu,
                                nuts_warmup = 1000) {
  mu = rep(NA, length(prior_mean_mu))
  tau = rep(NA, length(prior_alpha))
  for (j in 1:nrow(c_matrix)) {
    model_data = list(
      M = ncol(c_matrix),
      c = c_matrix[j, ],
      alpha = prior_alpha[j],
      beta = prior_beta[j],
      mu_0 = prior_mean_mu[j],
      sigma_0 = prior_sd_mu[j]
    )
    fit = sampling(
      object = stanmodel,
      data = model_data,
      chains = 1,
      iter = nuts_warmup + 1,
      warmup = nuts_warmup,
      refresh = 0 # do not print progress indicator
    )
    output_df = as.data.frame(fit)
    mu[j] = output_df$mu
    tau[j] = output_df$tau
  }
  return(c(mu, tau))
}