###################### SETUP AND RUN INFERENCE #################################
run_inference = function(model_name,
                         num_priorpred_samples,
                         num_surrogate_post_samples,
                         num_gibbs_samples,
                         num_rounds,
                         K_start,
                         mixprob_thresh,
                         cov_structure,
                         gllim_maxiter,
                         factor_cvMH,
                         mvn_package,
                         subfolder,
                         alg_version,
                         burnin_ind_param,
                         observation_name,
                         burnin_kappa_xi = 0,
                         gllim_file_path = NULL,
                         parallelize_ind_param_MH = FALSE,
                         random_seed = 1,
                         plot_results = FALSE,
                         sampler_method_eta = "conjugate_prior",
                         sampler_method_kappa_xi = "HMC",
                         iter_hmc_warmup_kappa_xi = 1000,
                         iter_hmc_warmup_eta = 1000) {
  # source files
  source(file.path("models", model_name, paste0(model_name, "_model.R")))
  source(file.path("models", model_name, "prior_definitions.R"))
  source(file.path("algorithm", paste0("semple_functions_", mvn_package, ".R")))
  source(file.path("algorithm", "write_results_to_file.R"))
  
  # read into local environment to save parameter settings
  source(file.path("models", model_name, "prior_definitions.R"), local = TRUE)
  source(file.path("models", model_name, paste0(model_name, "_model.R")), local = TRUE)
  
  # read observed data from file
  observed_data = read_observed_data(model_name, observation_name = observation_name)
  
  # read true parameter values that generated observed data, values are NA if real data is used
  true_param_file_name = file.path(
    "models",
    model_name,
    observation_name,
    "true_param.csv"
  )
  true_param = read.csv(file = true_param_file_name, header = TRUE)
  
  # read model param from file
  model_param_file = file.path(
    "models",
    model_name,
    observation_name,
    "model_param.csv"
  )
  model_param = read.csv(file = model_param_file)
  
  # create output dir
  output_dir = file.path(
    getwd(),
    "results",
    model_name,
    observation_name,
    subfolder
  )
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  assign("output_dir", output_dir, envir = .GlobalEnv)

  # write settings to json
  param_list = mget(setdiff(ls(), lsf.str()))
  settings_json = rjson::toJSON(param_list, indent = 1)
  write(settings_json, file.path(output_dir, "settings.json"))
  
  ########################## RUN ###################################
  if (alg_version == "unperturbed_kappa_xi") {
    source("algorithm/alg_unperturbed_kappa_xi.R")
    output = gibbs_unperturbed_kappa_xi(
      num_priorpred_samples = num_priorpred_samples,
      num_surrogate_post_samples = num_surrogate_post_samples,
      num_gibbs_samples = num_gibbs_samples,
      num_rounds = num_rounds,
      K_start = K_start,
      mixprob_thresh = mixprob_thresh,
      cov_structure = cov_structure,
      gllim_maxiter = gllim_maxiter,
      factor_cvMH = factor_cvMH,
      burnin_ind_param = burnin_ind_param,
      burnin_kappa_xi = burnin_kappa_xi,
      output_dir = output_dir,
      dim_ind_param = dim_ind_param,
      observed_data = observed_data,
      random_seed = random_seed,
      param_names = names(true_param),
      model_param = model_param,
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      prior_mean_mu = prior_mean_mu,
      sampler_method_eta = sampler_method_eta,
      sampler_method_kappa_xi = sampler_method_kappa_xi,
      parallelize_ind_param_MH = parallelize_ind_param_MH,
      iter_hmc_warmup_kappa_xi = iter_hmc_warmup_kappa_xi,
      iter_hmc_warmup_eta = iter_hmc_warmup_eta
    )
    
  } else if (alg_version == "only_random_effects") {
    source("algorithm/alg_only_random_effects.R")
    output = semple_mem_only_random_effects(
      num_priorpred_samples = num_priorpred_samples,
      num_surrogate_post_samples = num_surrogate_post_samples,
      num_gibbs_samples = num_gibbs_samples,
      num_rounds = num_rounds,
      K_start = K_start,
      mixprob_thresh = mixprob_thresh,
      cov_structure = cov_structure,
      gllim_maxiter = gllim_maxiter,
      factor_cvMH = factor_cvMH,
      burnin_ind_param = burnin_ind_param,
      output_dir = output_dir,
      dim_ind_param = dim_ind_param,
      observed_data = observed_data,
      random_seed = random_seed,
      param_names = names(true_param),
      model_param = model_param,
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      prior_mean_mu = prior_mean_mu,
      sampler_method_eta = sampler_method_eta,
      parallelize_ind_param_MH = parallelize_ind_param_MH,
      iter_hmc_warmup_eta = iter_hmc_warmup_eta
    )
    
  }  else {
    stop("Invalid algorithm version (alg_version)")
  }
  
  # optionally plot results
  if(plot_results){
    print("Computing posterior predictive checks")
    source(file.path("plot_scripts", "compute_ppc.R"), local = TRUE)
    compute_ppc_population_and_individual(input_dir = output_dir, round_index = 2)
    
    print("Plotting results")
    source(file.path("plot_scripts", "plot_results.R"), local = TRUE)
    if(alg_version == "only_random_effects"){
      plot_analysis_only_individual_param(output_dir)
    } else {
      plot_full_analysis(output_dir)
    }
    
  }
  
  return(output)
}

read_observed_data = function(model_name, observation_name){
  observed_data_file_name = file.path(
    "models",
    model_name,
    observation_name,
    "observations.csv"
  )
  observed_data = read.csv(file = observed_data_file_name)
  if(!is.null(observed_data$time)){
    observed_data = within(observed_data, rm(time))
  }
  observed_data_matrix = data.matrix(observed_data)
  return(unname(observed_data))
}

simulate_data_parallel = function(ind_param) {
  sim_data = mclapply(
    lapply(seq(ncol(ind_param)), function(i)
      ind_param[, i]),
    model,
    mc.set.seed = FALSE,
    mc.cores = detectCores() - 1
  )
  sim_data = do.call(cbind, sim_data)
  
  return(sim_data)
}

sample_eta_normal_gamma_conjugate = function(c, prior_alpha, prior_beta, prior_mean_mu) {
  tau = sample_tau_normal_gamma_conjugate(c, prior_alpha, prior_beta, prior_mean_mu)
  mu = sample_mu_normal_gamma_conjugate(tau, c, prior_mean_mu)
  eta = c(mu, tau)
  return(eta)
}

sample_tau_normal_gamma_conjugate = function(c, prior_alpha, prior_beta, prior_mean_mu) {
  M = ncol(c)
  tau = rep(NA, length(prior_alpha))
  for (j in 1:nrow(c)){
    c_mean = mean(c[j, ])
    post_alpha = prior_alpha[j] + M * 0.5
    SS = 0.5 * sum((c[j, ] - c_mean) ^ 2)
    post_beta = prior_beta[j] + SS + (M / (1 + M)) * ((c_mean - prior_mean_mu[j]) ^ 2) * 0.5
    # post_beta = prior_beta[j] + SS + (M / (2 * M)) * (c_mean - prior_mean_mu[j])
    tau[j] = rgamma(1, shape = post_alpha, rate = post_beta)
  }
  return(tau)
}

sample_mu_normal_gamma_conjugate = function(tau, c, prior_mean_mu) {
  M = ncol(c)
  mu = rep(NA, length(prior_mean_mu))
  for(j in 1:nrow(c)) {
    c_mean = mean(c[j, ])
    post_mean = (prior_mean_mu[j] + M * c_mean) / (1 + M)
    mu[j] = rnorm(1, mean = post_mean, sd = sqrt(1 / ((1 + M) * tau[j])))
  }
  return(mu)
}

logtarget_kappa_xi = function(kappa_xi,
                              ind_param,
                              gllim_lik_param,
                              pdf_kappa_xi_prior,
                              observed_data) {
  logpdf = log(pdf_kappa_xi_prior(kappa_xi))
  for (i in 1:ncol(observed_data)) {
    likelihood_mix_prob = likelihood_mix_prob(
      c = gllim_lik_param$c,
      Gamma = gllim_lik_param$Gamma,
      pi = gllim_lik_param$prior_mix_prob,
      K = gllim_lik_param$K,
      theta = c(ind_param[, i], kappa_xi)
    )
    
    mean_surrogate_likelihood = mean_surrogate_likelihood(
      A = gllim_lik_param$A,
      b = gllim_lik_param$b,
      observed_data = observed_data[, i],
      theta = c(ind_param[, i], kappa_xi),
      K = gllim_lik_param$K
    )
    
    surrogatelik = mixpdf(
      value = observed_data[, i],
      probs = likelihood_mix_prob,
      mean = mean_surrogate_likelihood,
      covars = gllim_lik_param$Sigma
    )
    
    logpdf = logpdf + log(surrogatelik)
  }
  return(logpdf)
}

sample_kappa_xi_adaptmcmc = function(ind_param,
                                     kappa_xi_old,
                                     observed_data,
                                     gllim_lik_param,
                                     pdf_kappa_xi_prior,
                                     scale_init,
                                     burnin,
                                     target_acc_rate = 0.234,
                                     gamma_adapt_speed = 2 / 3) {
  
  invisible(capture.output( # suppress printing in MCMC function
    list_output <- adaptMCMC::MCMC(
      p = function(kappa_xi) {
        logtarget_kappa_xi(kappa_xi, ind_param, gllim_lik_param, pdf_kappa_xi_prior, observed_data)
      },
      init = kappa_xi_old,
      n = 2 + burnin,
      scale = scale_init,
      acc.rate = target_acc_rate,
      gamma = gamma_adapt_speed,
      showProgressBar = FALSE
    )
  ))
  return(list_output)
}


plot_prior_simulations = function(sim_data) {
  plot(sim_data[, 1], type = "l", ylim = c(min(sim_data), max(sim_data)))
  for (i in 1:dim(sim_data)[2]) {
    lines(sim_data[, i])
  }
}

sample_surrogate_posterior = function(gllim_lik_param,
                                      num_samples,
                                      observed_data,
                                      dim_ind_param,
                                      dim_kappa_xi = 0,
                                      factor_cvMH) {
  M = ncol(observed_data)
  post_param_list = compute_post_param(gllim_lik_param)
  surrogate_posterior_samples = array(NA, dim = c(dim_ind_param + dim_kappa_xi, M, num_samples))
  # ind_param = array(NA, dim = c(dim_ind_param, M, num_samples))
  # kappa_xi_param = array(NA, dim = c(dim_kappa_xi, M, num_samples))
  
  for (i in 1:M) {
    observation_single = observed_data[, i]
    
    # compute surrogate mixture probabilities, means, and covariances
    post_mix_prob = compute_post_mix_prob(gllim_lik_param, observation_single)
    post_mix_mean = compute_mean_proposal(post_param_list, gllim_lik_param$K, observation_single)
    post_mix_cov = post_param_list$Sigma_star
    
    # sample from surrogate posterior mixture and simulate corresponding data
    for (j in 1:num_samples) {
      surrogate_posterior_sample = sample_gmm_inflated(
        probs = post_mix_prob,
        means = post_mix_mean,
        covars = post_mix_cov,
        factor_cvMH = factor_cvMH
      )
      surrogate_posterior_samples[, i, j] = surrogate_posterior_sample
      # ind_param[, i, j] = surrogate_post_sample[1:dim_ind_param]
      # kappa_xi_param[, i, j] = surrogate_post_sample[(dim_ind_param + 1):length(surrogate_post_sample)]
    }
  }
  
  if (dim_kappa_xi > 0) {
    kappa_xi_param = surrogate_posterior_samples[(dim_ind_param + 1):(dim_ind_param + dim_kappa_xi), , ]
    if(dim_kappa_xi == 1){dim(kappa_xi_param) = c(dim_kappa_xi, M, num_samples)}
    
    return(
      list(ind_param = surrogate_posterior_samples[1:dim_ind_param, , ],
           kappa_xi_param = kappa_xi_param)
    )
  }
  else {
    return(list(ind_param = surrogate_posterior_samples))
  }
  
}

compute_surrogate_posterior_mixture_param = function(gllim_lik_param, observed_data) {
  post_mixture_param_list = list()
  post_param_list = compute_post_param(gllim_lik_param)
  for (i in 1:ncol(observed_data)) {
    observation_single = observed_data[, i]
    post_mixture_param_list[[i]] = list(
      post_mix_prob = compute_post_mix_prob(gllim_lik_param, observation_single),
      post_mix_mean = compute_mean_proposal(post_param_list, gllim_lik_param$K, observation_single),
      post_mix_cov = post_param_list$Sigma_star
    )
  }
  return(post_mixture_param_list)
}

remove_kappa_xi_from_surrogate_posterior = function(post_mixture_param_list, dim_ind_param) {
  for (i in 1:length(post_mixture_param_list)) {
    post_mixture_param_list[[i]]$post_mix_mean = post_mixture_param_list[[i]]$post_mix_mean[1:dim_ind_param, ]
    post_mixture_param_list[[i]]$post_mix_cov = post_mixture_param_list[[i]]$post_mix_cov[1:dim_ind_param, 1:dim_ind_param, ]
  }
  return(post_mixture_param_list)
}

sample_ind_param_MH_modified_proposal = function(ind_param_old,
                                                 pdf_ind,
                                                 eta,
                                                 kappa_xi,
                                                 post_mixture_param_list,
                                                 factor_cvMH,
                                                 gllim_lik_param,
                                                 observation_single,
                                                 num_samples_saved,
                                                 burnin) {
  
  K = gllim_lik_param$K
  A = gllim_lik_param$A
  b = gllim_lik_param$b
  Sigma = gllim_lik_param$Sigma
  c = gllim_lik_param$c
  Gamma = gllim_lik_param$Gamma
  prior_mix_prob = gllim_lik_param$prior_mix_prob
  
  post_mix_prob = post_mixture_param_list$post_mix_prob
  post_mix_mean = post_mixture_param_list$post_mix_mean
  post_mix_cov = post_mixture_param_list$post_mix_cov
  
  dim_ind_param = length(ind_param_old)
  
  saved_ind_param = matrix(
    data = NA,
    nrow = dim_ind_param,
    ncol = (burnin + num_samples_saved)
  )
  
  # theta = (c, kappa, xi); ind_param = c
  theta_old = c(ind_param_old, kappa_xi)
  
  likelihood_mix_prob_old = likelihood_mix_prob(c, Gamma, prior_mix_prob, K, theta_old)
  mean_surrogate_likelihood_old = mean_surrogate_likelihood(A, b, observation_single, theta_old, K)
  surrogatelik_old = mixpdf(
    value = observation_single,
    probs = likelihood_mix_prob_old,
    mean = mean_surrogate_likelihood_old,
    covars = Sigma
  )
  
  surrogatepost_old = mixpdf(
    value = ind_param_old,
    probs = post_mix_prob,
    means = post_mix_mean,
    covars = factor_cvMH * post_mix_cov
  )
  
  num_accept = 0
  for (j in 1:(burnin + num_samples_saved)) {
    
    ind_param_proposal = sample_gmm_inflated(
      probs = post_mix_prob,
      means = post_mix_mean,
      covars = post_mix_cov,
      factor_cvMH = factor_cvMH
    )
    
    # add kappa and xi from the previous values of the Gibbs sampler
    theta_proposal = c(ind_param_proposal, kappa_xi)
    
    jac = jacobian(theta_old = ind_param_old, theta = ind_param_proposal) # the ratio of transformation jacobians
    
    # compute surrogate likelihood of proposal
    likelihood_mix_prob = likelihood_mix_prob(c, Gamma, prior_mix_prob, K, theta_proposal)
    mean_surrogate_likelihood = mean_surrogate_likelihood(A, b, observation_single, theta_proposal, K)
    surrogatelik = mixpdf(
      value = observation_single,
      probs = likelihood_mix_prob,
      mean = mean_surrogate_likelihood,
      covars = Sigma
    )
    
    # compute surrogate posterior of proposal
    surrogatepost_proposal = mixpdf(
      value = ind_param_proposal,
      probs = post_mix_prob,
      means = post_mix_mean,
      covars = factor_cvMH * post_mix_cov
    )
    
    MHratio = min(
      1,
      surrogatelik / surrogatelik_old * jac * pdf_ind(ind_param_proposal, eta) / pdf_ind(ind_param_old, eta) *
        surrogatepost_old / surrogatepost_proposal
    )
    
    if (is.nan(MHratio)) {
      saved_ind_param[, j] = ind_param_old
    } else if (runif(1) < MHratio) {
      theta_old = theta_proposal
      surrogatepost_old = surrogatepost_proposal
      surrogatelik_old = surrogatelik
      ind_param_old = ind_param_proposal
      saved_ind_param[, j] = ind_param_proposal
      num_accept = num_accept + 1
    }
    else{
      saved_ind_param[, j] = ind_param_old
    }
  }
  
  accrate = num_accept / (burnin + num_samples_saved)
  
  return(list(param = saved_ind_param, accrate = accrate))
}

sample_ind_param_MH = function(ind_param_old,
                               pdf_ind,
                               eta,
                               kappa_xi,
                               post_mixture_param_list,
                               factor_cvMH,
                               gllim_lik_param,
                               observation_single,
                               num_samples_saved,
                               burnin) {
  
  K = gllim_lik_param$K
  A = gllim_lik_param$A
  b = gllim_lik_param$b
  Sigma = gllim_lik_param$Sigma
  c = gllim_lik_param$c
  Gamma = gllim_lik_param$Gamma
  prior_mix_prob = gllim_lik_param$prior_mix_prob
  
  post_mix_prob_prev = post_mixture_param_list$post_mix_prob
  mean_proposal_prev = post_mixture_param_list$post_mix_mean
  covar_proposal_prev = post_mixture_param_list$post_mix_cov
  
  dim_ind_param = length(ind_param_old)
  
  saved_ind_param = matrix(
    data = NA,
    nrow = dim_ind_param,
    ncol = (burnin + num_samples_saved)
  )
  
  theta_old = c(ind_param_old, kappa_xi)
  
  likelihood_mix_prob_old = likelihood_mix_prob(c, Gamma, prior_mix_prob, K, theta_old)
  mean_surrogate_likelihood_old = mean_surrogate_likelihood(A, b, observation_single, theta_old, K)
  
  surrogatelik_old = mixpdf(
    value = observation_single,
    probs = likelihood_mix_prob_old,
    mean = mean_surrogate_likelihood_old,
    covars = Sigma
  )
  
  surrogatepost_old = mixpdf(
    value = theta_old,
    probs = post_mix_prob_prev,
    means = mean_proposal_prev,
    covars = factor_cvMH * covar_proposal_prev
  )
  
  num_accept = 0
  for (j in 1:(burnin + num_samples_saved)) {
    
    theta_proposal = sample_gmm_inflated(
      probs = post_mix_prob_prev,
      means = mean_proposal_prev,
      covars = covar_proposal_prev,
      factor_cvMH = factor_cvMH
    )
    
    ind_param_proposal = theta_proposal[1:dim_ind_param]
    
    # ignore xi from the proposal and replace it by the previous xi from the Gibbs sampler 
    theta_proposal[(dim_ind_param+1):(dim_ind_param+length(kappa_xi))] = kappa_xi
    
    jac = jacobian(theta_old = theta_old, theta = theta_proposal) # the ratio of transformation jacobians
    
    likelihood_mix_prob = likelihood_mix_prob(c, Gamma, prior_mix_prob, K, theta_proposal)
    mean_surrogate_likelihood = mean_surrogate_likelihood(A, b, observation_single, t(theta_proposal), K)
    surrogatelik = mixpdf(
      value = observation_single,
      probs = likelihood_mix_prob,
      mean = mean_surrogate_likelihood,
      covars = Sigma
    )
    
    surrogatepost_proposal = mixpdf(
      value = theta_proposal,
      probs = post_mix_prob_prev,
      means = mean_proposal_prev,
      covars = factor_cvMH * covar_proposal_prev
    )
    
    MHratio = min(
      1,
      surrogatelik / surrogatelik_old * jac * pdf_ind(ind_param_proposal, eta) / pdf_ind(ind_param_old, eta) *
        surrogatepost_old / surrogatepost_proposal
    )
    
    if (is.nan(MHratio)) {
      saved_ind_param[, j] = ind_param_old
    } else if (runif(1) < MHratio) {
      theta_old = theta_proposal
      surrogatepost_old = surrogatepost_proposal
      surrogatelik_old = surrogatelik
      ind_param_old = ind_param_proposal
      saved_ind_param[, j] = ind_param_proposal
      num_accept = num_accept + 1
    }
    else{
      saved_ind_param[, j] = ind_param_old
    }
  }
  
  accrate = num_accept / (burnin + num_samples_saved)
  
  return(list(param = saved_ind_param, accrate = accrate))
}

sample_MH_surrogatelikelihood = function(proposal_old,
                                         pdf_ind,
                                         eta,
                                         post_mixture_param_list,
                                         factor_cvMH,
                                         gllim_lik_param,
                                         observation_single,
                                         num_samples_saved,
                                         burnin) {
  K = gllim_lik_param$K
  A = gllim_lik_param$A
  b = gllim_lik_param$b
  Sigma = gllim_lik_param$Sigma
  c = gllim_lik_param$c
  Gamma = gllim_lik_param$Gamma
  prior_mix_prob = gllim_lik_param$prior_mix_prob
  
  post_mix_prob_prev = post_mixture_param_list$post_mix_prob
  mean_proposal_prev = post_mixture_param_list$post_mix_mean
  covar_proposal_prev = post_mixture_param_list$post_mix_cov
  
  allparMCMC = matrix(
    data = NA,
    nrow = length(proposal_old),
    ncol = (burnin + num_samples_saved)
  )
  
  likelihood_mix_prob_old = likelihood_mix_prob(c, Gamma, prior_mix_prob, K, proposal_old)
  mean_surrogate_likelihood_old = mean_surrogate_likelihood(A, b, observation_single, proposal_old, K)
  surrogatelik_old = mixpdf(
    value = observation_single,
    probs = likelihood_mix_prob_old,
    mean = mean_surrogate_likelihood_old,
    covars = Sigma
  )
  
  surrogatepost_old = mixpdf(
    value = proposal_old,
    probs = post_mix_prob_prev,
    means = mean_proposal_prev,
    covars = factor_cvMH * covar_proposal_prev
  )
  
  num_accept = 0
  for (j in 1:(burnin + num_samples_saved)) {
    
    proposal = sample_gmm_inflated(
      probs = post_mix_prob_prev,
      means = mean_proposal_prev,
      covars = covar_proposal_prev,
      factor_cvMH = factor_cvMH
    )
    
    jac = jacobian(theta_old = proposal_old, theta = proposal) # the ratio of transformation jacobians
    
    likelihood_mix_prob = likelihood_mix_prob(c, Gamma, prior_mix_prob, K, proposal)
    mean_surrogate_likelihood = mean_surrogate_likelihood(A, b, observation_single, t(proposal), K)
    surrogatelik = mixpdf(
      value = observation_single,
      probs = likelihood_mix_prob,
      mean = mean_surrogate_likelihood,
      covars = Sigma
    )
    
    surrogatepost_proposal = mixpdf(
      value = proposal,
      probs = post_mix_prob_prev,
      means = mean_proposal_prev,
      covars = factor_cvMH * covar_proposal_prev
    )
    
    MHratio = min(
      1,
      surrogatelik / surrogatelik_old * jac * pdf_ind(proposal, eta) / pdf_ind(proposal_old, eta) *
        surrogatepost_old / surrogatepost_proposal
    )
    
    if (is.nan(MHratio)) {
      allparMCMC[, j] = proposal_old
    } else if (runif(1) < MHratio) {
      proposal_old = proposal
      surrogatepost_old = mixpdf(
        value = proposal_old,
        probs = post_mix_prob_prev,
        means = mean_proposal_prev,
        covars = factor_cvMH * covar_proposal_prev
      )
      allparMCMC[, j] = proposal
      
      surrogatelik_old = surrogatelik
      
      num_accept = num_accept + 1
    }
    else{
      allparMCMC[, j] = proposal_old
    }
  }
  
  accrate = num_accept / (burnin + num_samples_saved)
  
  return(list(param = allparMCMC, accrate = accrate))
}