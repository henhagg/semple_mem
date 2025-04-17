library(xLLiM)
library(rjson)
library(ggplot2)

compute_bic = function(model_name,
                       num_priorpred_samples,
                       K_values,
                       cov_structure,
                       model_param,
                       kappa_xi_exists = TRUE, # there are constant parameters in the model
                       random_seed = 1,
                       gllim_verb = 1,
                       parallelize_simulations = FALSE,
                       output_subfolder = NULL) {
  
  source(file.path("models", model_name, paste0(model_name, "_model.R")), local = T)
  source(file.path("models", model_name, "prior_definitions.R"), local = T)
  
  set.seed(random_seed) # for reproducibility
  
  start_time = Sys.time()
  time_vec = rep(NA, length(K_values))
  
  pop_param = replicate(num_priorpred_samples, sample_eta_prior())
  if(kappa_xi_exists){
    kappa_xi_param = replicate(num_priorpred_samples, sample_kappa_xi_prior())
    if(!is.matrix(kappa_xi_param)) {
      kappa_xi_param = matrix(kappa_xi_param, ncol = length(kappa_xi_param))
    }
  } else {
    kappa_xi_param = NULL
  }

  ind_param = apply(pop_param, MARGIN = 2, sample_ind_param)
  
  if (parallelize_simulations) {
    sim_data = mclapply(
      1:ncol(ind_param),
      FUN = function(i)
        model(ind_param[, i], kappa_xi_param[, i], model_param),
      mc.set.seed = FALSE,
      mc.cores = detectCores() - 1
    )
    sim_data = do.call(cbind, sim_data)
  } else {
    sim_data = sapply(
      1:num_priorpred_samples,
      FUN = function(i)
        model(ind_param[, i], kappa_xi_param[, i], model_param)
    )
  }
  
  bic_values = rep(NA, length(K_values))
  
  for (K in K_values) {
    print(paste("Computing BIC for K =", K))
    mod = gllim(
      ind_param,
      sim_data,
      in_K = K,
      maxiter = 300,
      verb = gllim_verb,
      cstr = list(Sigma = cov_structure)
    )
    bic_values[which(K_values == K)] = mod$nbpar * log(num_priorpred_samples) - 2 *
      mod$LLf
    
    time_vec[which(K_values == K)] = difftime(Sys.time(), start_time, units = "secs")
  }
  
  # save results to file
  if (!is.null(output_subfolder)) {
    output_dir = create_output_dir(model_name, output_subfolder)
    
    # save_settings_to_json_file(cov_structure, num_priorpred_samples, model_param, output_dir)
    # write settings to json
    param_list = mget(setdiff(ls(), lsf.str()))
    param_list[c("pop_param", "kappa_xi_param", "ind_param", "sim_data", "mod")] = NULL
    settings_json = rjson::toJSON(param_list, indent = 1)
    write(settings_json, file.path(output_dir, "settings.json"))
    
    # save bic values to csv
    write.table(
      cbind(K_values, bic_values), file = file.path(output_dir, "bic.csv"), 
      sep = ",",
      row.names = F,
      col.names = c("K", "BIC"),
      quote = F
    )
    
    save_bic_plot_to_pdf(
      K_values = K_values,
      bic_values = bic_values,
      output_dir = output_dir
    )
    
    # save cumulative elapsed time of bic computations to csv
    write.table(
      cbind(K_values, time_vec),
      file = file.path(output_dir, "cumulative_elapsed_time.csv"),
      sep = ",",
      row.names = F,
      col.names = c("K", "seconds"),
      quote = F
    )
  } else {
    plot(
      K_values,
      bic_values,
      type = "b",
      xlab = "K",
      ylab = "BIC"
    )
  }
  
  return(bic_values)
}

compute_bic_perturbed = function(model_name,
                                 num_priorpred_samples,
                                 K_values,
                                 cov_structure,
                                 gllim_verb) {
  source(file.path("models", model_name, paste0(model_name, "_model.R")), local = T)
  source(file.path("models", model_name, "prior_definitions.R"),
         local = T)
  
  set.seed(1) # for reproducibility
  
  pop_param = replicate(num_priorpred_samples, sample_eta_prior())
  ind_param = apply(pop_param, MARGIN = 2, sample_ind_param)
  sim_data = apply(ind_param, MARGIN = 2, model)
  
  bic_values = rep(NA, length(K_values))
  
  for (K in K_values) {
    print(K)
    mod = gllim(
      ind_param,
      sim_data,
      in_K = K,
      maxiter = 300,
      verb = gllim_verb,
      cstr = list(Sigma = cov_structure)
    )
    bic_values[which(K_values == K)] = mod$nbpar * log(num_priorpred_samples) - 2 *
      mod$LLf
  }
  
  return(bic_values)
}

create_output_dir = function(model_name, output_subfolder) {
  output_dir = file.path("results", model_name, "bic", output_subfolder)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  return(output_dir)
}

save_settings_to_json_file = function(gllim_cov_structure,
                                      num_priorpred_samples,
                                      model_param,
                                      output_dir) {
  settings = list(
    gllim_covariance_structure = gllim_cov_structure,
    num_priorpred_samples = num_priorpred_samples,
    model_param = model_param
  )
  settings_json = rjson::toJSON(settings, indent = 1)
  write(settings_json, file.path(output_dir, "settings.json"))
}

save_bic_plot_to_pdf = function(K_values,
                                bic_values,
                                output_dir,
                                pdf_width = 5,
                                pdf_height = 3) {
  bic_data = data.frame(K_values, bic_values)
  ggplot(bic_data, aes(x = K_values, y = bic_values)) + geom_line() + geom_point() + xlab("K") + ylab("BIC")
  ggsave(file.path(output_dir, "bic.pdf"),
         width = pdf_width,
         height = pdf_height)
}

plot_bic_from_file = function(input_file, save_to_file = FALSE) {
  bic_data = read.csv(file = input_file)
  ggplot(bic_data, aes(x = K, y = BIC)) + geom_line() + geom_point()
  
  if(save_to_file){
    ggsave(file.path(dirname(input_file), "bic.pdf"))
  }
}
