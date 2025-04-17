semple_mem_only_random_effects = function(num_priorpred_samples,
                                      num_surrogate_post_samples,
                                      num_gibbs_samples,
                                      num_rounds,
                                      K_start,
                                      mixprob_thresh,
                                      cov_structure,
                                      gllim_maxiter,
                                      factor_cvMH,
                                      burnin_ind_param,
                                      output_dir,
                                      dim_ind_param,
                                      observed_data,
                                      random_seed,
                                      param_names,
                                      model_param,
                                      prior_alpha,
                                      prior_beta,
                                      prior_mean_mu,
                                      sampler_method_eta,
                                      parallelize_ind_param_MH = FALSE,
                                      iter_hmc_warmup_eta = 1000,
                                      plot_prior_sim = TRUE) {
  
  start_time = Sys.time()
  time_vec = rep(NA, num_rounds + 1)
  K_values = K_start
  dim_data = model_param$dim_data
  
  if (sampler_method_eta == "HMC") {
    source(file.path("algorithm", "hmc_eta.R"), local = TRUE)
  } else if (sampler_method_eta == "HMC_indep_normal_gamma") {
    source(file.path("algorithm", "hmc_eta_indep_prior.R"), local = TRUE)
  }
  
  set.seed(random_seed)
  ################ PRIOR PREDICTIVE STEP #############################
  prior_pop_param = replicate(num_priorpred_samples, sample_eta_prior())
  prior_ind_param = apply(prior_pop_param, MARGIN = 2, sample_ind_param)
  time_vec[1] = difftime(Sys.time(), start_time, units = "secs")

  prior_sim_data = sapply(
    1:num_priorpred_samples,
    FUN = function(i)
      model(prior_ind_param[, i], log_kappa_xi = NULL, model_param = model_param)
  )

  write_prior_pred_to_file(output_dir = output_dir,
                           time_vec = time_vec,
                           pop_param = prior_pop_param,
                           ind_param = prior_ind_param,
                           param_names = param_names)
  
  if(plot_prior_sim){
    plot_prior_simulations(prior_sim_data)
  }
  
  prior_responses = prior_ind_param
  rownames(prior_responses) = NULL
  
  gllim_lik_param = fit_gllim(
    responses = prior_responses,
    covariates = prior_sim_data,
    K = K_start,
    maxiter = gllim_maxiter,
    verbose = 1,
    cov_structure = cov_structure,
    mixprob_thresh = mixprob_thresh,
    init_param = NULL
  )
  saveRDS(object = gllim_lik_param, file = file.path(output_dir, "gllim_model_round0.rds"))
  K_values = c(K_values, gllim_lik_param$K)
  
  ################## SAMPLE FROM SURROGATE POSTERIOR MIXTURE ######################
  M = ncol(observed_data)
  samples_surrogate_post = sample_surrogate_posterior(gllim_lik_param = gllim_lik_param,
                                                      num_samples = num_surrogate_post_samples,
                                                      observed_data = observed_data,
                                                      dim_ind_param = dim_ind_param,
                                                      factor_cvMH = 1)
  ind_param_surrogate_post = samples_surrogate_post$ind_param
  
  time_vec[2] = difftime(Sys.time(), start_time, units = "secs")
  
  sim_data_surrogate_post = array(NA, dim=c(dim_data, M, num_surrogate_post_samples))
  for(i in 1:M){
    for (j in 1:num_surrogate_post_samples) {
      sim_data_surrogate_post[, i, j] = model(ind_param_surrogate_post[, i, j], log_kappa_xi = NULL, model_param = model_param)
    }
  }
  
  write_ind_param_to_csv(
    ind_param = ind_param_surrogate_post,
    dim_ind_param = dim_ind_param,
    M = M,
    num_gibbs_samples = num_surrogate_post_samples,
    output_dir = output_dir,
    round_index = 1
  )
  
  # stack training data from the surrogate posterior samples (D_1) in matrix format
  ind_param_surrogate_post_matrix = matrix(
    ind_param_surrogate_post,
    nrow = dim_ind_param,
    ncol = M * num_surrogate_post_samples
  )
  
  if(num_rounds < 2){return()}
  
  responses = ind_param_surrogate_post_matrix
  covariates = matrix(
    sim_data_surrogate_post,
    nrow = dim_data,
    ncol = M * num_surrogate_post_samples
  )
  
  acceptance_rates_ind_param = array(NA, dim=c(num_gibbs_samples, M, num_rounds - 1))
  
  ind_param_temp = ind_param_surrogate_post[, , dim(ind_param_surrogate_post)[3]] # save last previously accepted values from last round temporarily
  ind_param = array(NA, dim = c(dim_ind_param, M, num_gibbs_samples))
  ind_param[, , 1] = ind_param_temp
  
  # initialize population parameters eta
  eta_matrix = matrix(NA, nrow = dim(prior_pop_param)[1], ncol = num_gibbs_samples)
  if(sampler_method_eta == "conjugate_prior") {
    eta_matrix[, 1] = sample_eta_normal_gamma_conjugate(
      c = ind_param_temp,
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      prior_mean_mu = prior_mean_mu
    )
  }
  else if (sampler_method_eta == "HMC") {
    eta_matrix[, 1] = sample_eta_with_nuts(
      c_matrix = as.matrix(ind_param[, , 1]),
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      prior_mean_mu = prior_mean_mu,
      nuts_warmup = iter_hmc_warmup_eta
    )
  }else if (sampler_method_eta == "HMC_indep_normal_gamma") {
    eta_matrix[, 1] = sample_eta_with_nuts(
      c_matrix = as.matrix(ind_param[, , 1]),
      prior_alpha = prior_alpha,
      prior_beta = prior_beta,
      prior_mean_mu = prior_mean_mu,
      prior_sd_mu = prior_sd_mu,
      nuts_warmup = iter_hmc_warmup_eta
    )
  } else {
    warning("sampler_method_eta is invalid")
  }

  for (r in 2:num_rounds) {
    print(paste("Running round", r))
    
    gllim_lik_param = fit_gllim(
      responses,
      covariates,
      K = gllim_lik_param$K,
      maxiter = gllim_maxiter,
      verbose = 1,
      cov_structure = cov_structure,
      mixprob_thresh = mixprob_thresh,
      init_param = gllim_lik_param$optpar
    )
    saveRDS(object = gllim_lik_param, file = file.path(output_dir, paste0("gllim_model_round", r - 1, ".rds")))
    K_values = c(K_values, gllim_lik_param$K)
    
    post_mixture_param_list = compute_surrogate_posterior_mixture_param(gllim_lik_param, observed_data)
    reduced_post_mixture_param_list = remove_kappa_xi_from_surrogate_posterior(post_mixture_param_list, dim_ind_param)
    
    pb = progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = num_gibbs_samples,
      complete = "=",
      incomplete = "-",
      current = ">",
      clear = FALSE,
      width = 100
    )
    
    print(paste("Running gibbs sampling round", r))
    for (j in 2:num_gibbs_samples) {
      pb$tick()
      
      # SAMPLE NEW C
      if(parallelize_ind_param_MH){
        ind_MH_output = mclapply(1:M, function(i)
          sample_MH_surrogatelikelihood(
            proposal_old = ind_param[, i, j - 1],
            pdf_ind = pdf_ind,
            eta = eta_matrix[, j - 1],
            post_mixture_param_list = post_mixture_param_list[[i]],
            factor_cvMH = factor_cvMH,
            gllim_lik_param = gllim_lik_param,
            observation_single = observed_data[, i],
            num_samples_saved = 1,
            burnin = burnin_ind_param
          )$param[, burnin_ind_param + 1], 
          mc.set.seed = FALSE, # not to mess with random seed from the beginning of the algorithm
          mc.cores = detectCores() - 1)
        ind_param[, , j] = do.call(cbind, ind_MH_output)
      } else {
        for (i in 1:M) {
          ind_MH_output = sample_MH_surrogatelikelihood(
            proposal_old = ind_param[, i, j - 1],
            pdf_ind = pdf_ind,
            eta = eta_matrix[, j - 1],
            post_mixture_param_list = post_mixture_param_list[[i]],
            factor_cvMH = factor_cvMH,
            gllim_lik_param = gllim_lik_param,
            observation_single = observed_data[, i],
            num_samples_saved = 1,
            burnin = burnin_ind_param
          )
          acceptance_rates_ind_param[j, i, r - 1] = ind_MH_output$accrate
          ind_param[, i, j] = ind_MH_output$param[, burnin_ind_param + 1] # remove burnin
        }
      }
      
      # SAMPLE NEW ETA (BY NORMAL-NORMALGAMMA CONJUGACY)
      if(sampler_method_eta == "conjugate_prior"){
        eta_matrix[, j] = sample_eta_normal_gamma_conjugate(
          c = as.matrix(ind_param[, , j]),
          prior_alpha = prior_alpha,
          prior_beta = prior_beta,
          prior_mean_mu = prior_mean_mu
        )
      }else if (sampler_method_eta == "HMC") {
        eta_matrix[, j] = sample_eta_with_nuts(
          c_matrix = as.matrix(ind_param[, , j]),
          prior_alpha = prior_alpha,
          prior_beta = prior_beta,
          prior_mean_mu = prior_mean_mu,
          nuts_warmup = iter_hmc_warmup_eta
        )
      }else if (sampler_method_eta == "HMC_indep_normal_gamma") {
        eta_matrix[, j] = sample_eta_with_nuts(
          c_matrix = as.matrix(ind_param[, , j]),
          prior_alpha = prior_alpha,
          prior_beta = prior_beta,
          prior_mean_mu = prior_mean_mu,
          prior_sd_mu = prior_sd_mu,
          nuts_warmup = iter_hmc_warmup_eta
        )
      } else {
        warning("sampler_method_eta is invalid")
      }
      
    }
    
    time_vec[r + 1] = difftime(Sys.time(), start_time, units = "secs")
    write_gibbs_results_to_csv(
      time_vec = time_vec,
      ind_param = ind_param,
      dim_ind_param = dim_ind_param,
      M = M,
      num_gibbs_samples = num_gibbs_samples,
      output_dir = output_dir,
      eta_matrix = eta_matrix,
      acceptance_rates_ind_param = acceptance_rates_ind_param[, , r - 1],
      round_index = r,
      param_names = param_names
    )
    write.table(
      cbind(K_values, seq(0, r)),
      file = file.path(output_dir, paste0("num_mixture_components.csv")),
      sep = ",",
      row.names = F,
      col.names = c("K", "round"),
      quote = F
    )
    
    if (r < num_rounds) {
      # simulate new data from gibbs samples
      sim_data = array(NA, dim = c(dim_data, M, num_gibbs_samples))
      for (i in 1:M) {
        for (j in 1:num_gibbs_samples) {
          sim_data[, i, j] = model(ind_param[, i, j], model_param)
        }
      }
      
      ind_param_matrix = matrix(ind_param,
                                nrow = dim_ind_param,
                                ncol = M * num_gibbs_samples)
      new_responses = ind_param_matrix
      responses = cbind(responses, new_responses)
      covariates = cbind(covariates,
                         matrix(sim_data, nrow = dim_data, ncol = M * num_gibbs_samples))
      
      # initialize next round to last previously accepted values
      ind_param[, , 1] = ind_param[, , dim(ind_param)[3]]
      eta_matrix[, 1] = eta_matrix[, num_gibbs_samples]
    }
    
  } # end for-loop
}