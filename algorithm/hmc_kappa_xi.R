library(cmdstanr)

options(mc.cores = parallel::detectCores())
stan_mod_kappa_xi = cmdstan_model(file.path('algorithm', 'kappa_xi_gllim_normal_prior.stan'))

run_hmc_warmup = function(stan_mod_kappa_xi,
                          gllim_mod,
                          dim_ind_param,
                          dim_kappa_xi,
                          dim_data,
                          ind_param_init,
                          observed_data,
                          prior_means_kappa_xi,
                          prior_sds_kappa_xi,
                          init_kappa_xi,
                          stan_output_dir,
                          iter_warmup = 1000,
                          iter_sampling = 100) {
  
  # the arrays needs to be converted to lists before giving them to Stan
  gllim_Sigma_list = lapply(seq(dim(gllim_mod$Sigma)[3]), function(x)
    gllim_mod$Sigma[, , x])
  gllim_Gamma_list = lapply(seq(dim(gllim_mod$Gamma)[3]), function(x)
    gllim_mod$Gamma[, , x])
  gllim_A_list = lapply(seq(dim(gllim_mod$A)[3]), function(x)
    gllim_mod$A[, , x])
  
  stan_data = list(
    K = gllim_mod$K,
    D = dim_ind_param + dim_kappa_xi,
    L = dim_data,
    M = ncol(observed_data),
    dim_ind_param = dim_ind_param,
    dim_kappa_xi = dim_kappa_xi,
    c_gllim = gllim_mod$c,
    c_model = ind_param_init,
    b = gllim_mod$b,
    pi = gllim_mod$pi,
    A = gllim_A_list,
    Sigma = gllim_Sigma_list,
    Gamma = gllim_Gamma_list,
    y_o = observed_data,
    prior_means_kappa_xi = prior_means_kappa_xi,
    prior_sds_kappa_xi = prior_sds_kappa_xi
  )
  
  if(!dir.exists(stan_output_dir)){
    dir.create(stan_output_dir, recursive = TRUE) # create directory to save Stan output
  }
  
  fit_warmup = stan_mod_kappa_xi$sample(
    data = stan_data,
    init = list(list(kappa_xi = init_kappa_xi)),
    chains = 1,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = 1,
    output_dir = stan_output_dir
  )
  
  return(
    list(
      step_size_adaptation = fit_warmup$metadata()$step_size_adaptation,
      inv_metric = fit_warmup$inv_metric()[[1]],
      metric_type = fit_warmup$metadata()$metric,
      stan_data = stan_data
    )
  )
}

sample_kappa_xi_warmed_up_hmc = function(kappa_xi_init,
                                         ind_param_current,
                                         stan_mod_kappa_xi,
                                         gllim_mod,
                                         num_hmc_samples,
                                         step_size,
                                         inv_metric_vector,
                                         metric_type,
                                         stan_data) {
  stan_data$c_model = ind_param_current # update the current value of indvidual param
  fit_gibbs = stan_mod_kappa_xi$sample(
    data = stan_data,
    init = list(list(kappa_xi = kappa_xi_init)),
    iter_warmup = 0,
    adapt_engaged = FALSE,
    iter_sampling = num_hmc_samples,
    chains = 1,
    show_messages = FALSE,
    show_exceptions = FALSE,
    diagnostics = "", # do not check diagnostics automatically
    refresh = 0,
    step_size = step_size,
    inv_metric = diag(inv_metric_vector),
    metric = metric_type
  )
  
  last_sample = fit_gibbs$draws("kappa_xi")[num_hmc_samples, 1, ]
  return(last_sample)
}
