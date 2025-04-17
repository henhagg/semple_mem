generate_observations = function(model_name,
                                 mu,
                                 tau,
                                 log_kappa_xi = NULL,
                                 num_individuals,
                                 model_param = NULL,
                                 save_to_file = FALSE,
                                 save_to_file_pepsdi_format = FALSE,
                                 observation_name = NULL,
                                 random_seed = NULL) {
  
  source(file.path("models", model_name, paste0(model_name, "_model.R")), local = TRUE)
  source(file.path("models", model_name, "prior_definitions.R"), local = TRUE)
  
  if(is.null(random_seed)){
    random_seed = num_individuals
  }
  set.seed(random_seed)
  
  eta_true = c(mu, tau)
  
  ind_param = replicate(num_individuals, sample_ind_param(eta_true))
  
  observations = matrix(NA, nrow = model_param$dim_data, ncol = num_individuals)
  for (i in 1:num_individuals) {
    observations[, i] = model(ind_param[, i], log_kappa_xi, model_param)
  }
  time_vec = seq(0, model_param$tt, length.out = model_param$dim_data)
  
  plot(observations[, 1], type = "l", ylim = c(min(observations), max(observations)))
  for (i in 1:num_individuals) {
    lines(observations[, i])
  }
  
  if (save_to_file) {
    if(!is.null(observation_name)){
      output_dir = file.path("models", model_name, observation_name)
    } else{
      output_dir = file.path("models", model_name, paste0("num_observation_", num_individuals))
    }
   
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    pdf(file = file.path(output_dir, "observations.pdf"), width = 7, height = 5)
    plot(observations[, 1], type = "l", ylim = c(min(observations), max(observations)))
    for (i in 1:num_individuals) {
      lines(observations[, i])
    }
    dev.off()
    
    write.table(
      random_seed,
      file = file.path(output_dir, "random_seed.csv"),
      sep = ",",
      row.names = F,
      col.names = "random_seed",
      quote = F
    )
    write.table(
      cbind(time_vec, observations),
      file = file.path(output_dir, "observations.csv"),
      sep = ",",
      row.names = F,
      col.names = c("time", paste0("obs", 1:num_individuals)),
      quote = F
    )
    write.table(
      ind_param,
      file = file.path(output_dir, "ind_param.csv"),
      sep = ",",
      row.names = F,
      col.names = paste0("ind", 1:dim(ind_param)[2]),
      quote = F
    )
    write.table(
      model_param,
      file.path(output_dir, "model_param.csv"),
      sep = ",",
      row.names = F,
      quote = F
    )
    
    true_param = eta_true
    colnames_true_param = c(paste0("mu", 1:length(mu)), paste0("tau", 1:length(tau)))
    
    # append kappa and xi parameters to the vector of true parameter values
    if(!is.null(log_kappa_xi)){
      true_param = c(true_param, log_kappa_xi)
      kappa_exists = exists("dim_kappa")
      xi_exists = exists("dim_xi")
      if(kappa_exists && xi_exists){
        colnames_true_param = c(colnames_true_param, paste0("kappa", 1:dim_kappa), paste0("xi", 1:dim_xi))
      } else if(kappa_exists && !xi_exists){
        colnames_true_param = c(colnames_true_param, paste0("kappa", 1:dim_kappa))
      } else {
        colnames_true_param = c(colnames_true_param, paste0("xi", 1:dim_xi))
      }
    }
    
    write.table(
      t(true_param),
      file.path(output_dir, "true_param.csv"),
      sep = ",",
      row.names = F,
      col.names = colnames_true_param,
      quote = F
    )
    
    if(save_to_file_pepsdi_format) {
      # write to PEPSDI csv file format
      pepsdi_df = data.frame()
      for (i in 1:num_individuals) {
        row_indices = ((i - 1) * model_param$dim_data + 1):(i * model_param$dim_data)
        pepsdi_df[row_indices, 1] = time_vec
        pepsdi_df[row_indices, 2] = observations[, i]
        pepsdi_df[row_indices, 3] = 1
        pepsdi_df[row_indices, 4] = i
      }
      header <- c("time", "obs", "obs_id", "id")
      names(pepsdi_df) <- header
      write.table(
        x = pepsdi_df,
        file = file.path(output_dir, "observations_pepsdi.csv"),
        sep = ",",
        row.names = F,
        col.names = T,
        quote = F
      )
    }
  }
}

generate_observations(model_name = "ornstein_uhlenbeck_unperturbed_noise",
                      mu = c(-0.7, 2.3, -0.9),
                      tau = c(4, 10, 4),
                      log_kappa_xi = -1.2,
                      num_individuals = 40,
                      model_param = list(tt = 10, dim_data = 51, x0 = 0),
                      observation_name = "num_observation_40",
                      save_to_file = T)
generate_observations(model_name = "mrna_fix_tzero",
                      mu = c(-0.694, -3, 0.027),
                      tau = c(10, 10, 10),
                      log_kappa_xi = c(5.704, 0.751, 2.079, -1.6),
                      num_individuals = 40,
                      model_param = list(tt = 30, dt = 0.01, dim_data = 61),
                      save_to_file = T,
                      save_to_file_pepsdi_format = T)
generate_observations(model_name = "mrna_indep_prior",
                      mu = c(-0.694, -3, 0.027, -0.164),
                      tau = c(10, 10, 10, 10),
                      log_kappa_xi = c(5.704, 0.751, 2.079, -3.454),
                      num_individuals = 400,
                      model_param = list(tt = 30, dt = 0.01, dim_data = 60),
                      observation_name = "400ind",
                      save_to_file = T)
generate_observations(model_name = "mrna_indep_prior_only_individual_param",
                      mu = c(-0.694, -3, 0.027, -0.164, 5.704, 0.751, 2.079, -3.454),
                      tau = c(1, 1, 1, 1, 1, 1, 1, 1),
                      num_individuals = 200,
                      model_param = list(tt = 30, dt = 0.01, dim_data = 60),
                      observation_name = "200ind_precision1_all",
                      save_to_file = T)
