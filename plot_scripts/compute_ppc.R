library(Rcpp)
library(rjson)
library(dplyr)
library(ggplot2)

compute_ppc_single_round = function(input_dir, round_index, random_seed = 1) {
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  observed_data = read.csv(
    file.path(
      "models",
      settings$model_name,
      settings$observation_name,
      "observations.csv"
    )
  )
  
  set.seed(random_seed)
  
  source(file.path("models", settings$model_name, paste0(settings$model_name, "_model.R")))
  
  ind_param_samples = read.csv(file.path(input_dir, paste0("ind_param_round", round_index, ".csv")))
  kappa_xi_file_path = file.path(input_dir, paste0("kappa_xi_round", round_index, ".csv"))
  if(file.exists(kappa_xi_file_path)){
    kappa_xi_samples = read.csv(kappa_xi_file_path)
  } else { # there are no constant parameter samples
    kappa_xi_samples = NULL
  }
  
  post_sim_data = matrix(data = NA, nrow = nrow(ind_param_samples), ncol = settings$model_param$dim_data)
  num_ind = length(unique(ind_param_samples$id))
  num_gibbs_cycles = length(unique(ind_param_samples$gibbs_cycle))
  for(i in 1:num_ind){
    ind_param_specific_ind = ind_param_samples %>% filter(id == i) %>% select(!c(id, gibbs_cycle))
    for(j in 1:num_gibbs_cycles) {
      post_sim_data[(i - 1) * num_gibbs_cycles + j, ] = model(
        log_ind_param = unlist(ind_param_specific_ind[j, ]),
        log_kappa_xi = unlist(kappa_xi_samples[j, ]),
        model_param = settings$model_param
      )
    }
  }
  
  probs = c(0, 0.025, 0.05, 0.95, 0.975, 1)
  ppc_qu = t(apply(
    post_sim_data,
    MARGIN = 2,
    FUN = quantile,
    probs = probs
  ))
  
  ppc_qu_df = data.frame(ppc_qu)
  colnames(ppc_qu_df) = paste0("qu", sub("\\.", "", as.character(probs)))
  ppc_qu_df$time = observed_data$time
  
  write.table(
    ppc_qu_df,
    file = file.path(input_dir, paste0("ppc_round", round_index, ".csv")),
    row.names = F,
    sep = ",",
    quote = F
  )
}

compute_ppc_individual_single_round = function(input_dir, round_index, random_seed = 1) {
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  observed_data = read.csv(
    file.path(
      "models",
      settings$model_name,
      settings$observation_name,
      "observations.csv"
    )
  )
  
  set.seed(random_seed)
  
  source(file.path("models", settings$model_name, paste0(settings$model_name, "_model.R")))
  
  ind_param_samples = read.csv(file.path(input_dir, paste0("ind_param_round", round_index, ".csv")))
  kappa_xi_file_path = file.path(input_dir, paste0("kappa_xi_round", round_index, ".csv"))
  if(file.exists(kappa_xi_file_path)){
    kappa_xi_samples = read.csv(kappa_xi_file_path)
  } else { # there are no constant parameter samples
    kappa_xi_samples = NULL
  }
  
  num_ind = length(unique(ind_param_samples$id))
  num_gibbs_cycles = length(unique(ind_param_samples$gibbs_cycle))
  dim_data = settings$model_param$dim_data
  post_sim_data = array(data = NA, dim = c(num_gibbs_cycles, dim_data, num_ind))

  for(i in 1:num_ind){
    ind_param_specific_ind = ind_param_samples %>% filter(id == i) %>% select(!c(id, gibbs_cycle))
    for(j in 1:num_gibbs_cycles) {
      post_sim_data[j, , i] = model(
        log_ind_param = unlist(ind_param_specific_ind[j, ]),
        log_kappa_xi = unlist(kappa_xi_samples[j, ]),
        model_param = settings$model_param
      )
    }
  }

  probs = c(0, 0.025, 0.05, 0.95, 0.975, 1)
  ppc_qu = apply(
    post_sim_data,
    MARGIN = c(2, 3),
    FUN = quantile,
    probs = probs
  )

  ppc_qu_matrix = t(matrix(ppc_qu,
                         nrow = dim(ppc_qu)[1],
                         ncol = dim_data * num_ind))
  ppc_qu_df = data.frame(ppc_qu_matrix)
  colnames(ppc_qu_df) = paste0("qu", sub("\\.", "", as.character(probs)))
  ppc_qu_df$id = rep(1:num_ind, each = dim_data)
  ppc_qu_df$time = rep(observed_data$time, num_ind)

  write.table(
    ppc_qu_df,
    file = file.path(input_dir, paste0("ppc_individual_round", round_index, ".csv")),
    row.names = F,
    sep = ",",
    quote = F
  )
}

compute_ppc_population_and_individual = function(input_dir, round_index) {
  compute_ppc_single_round(input_dir, round_index)
  compute_ppc_individual_single_round(input_dir, round_index)
}

compute_ppc_from_pop_param = function(model_name,
                                      posterior_samples_file,
                                      output_file = NULL,
                                      random_seed = 1) {
  set.seed(random_seed)

  source(file.path("models", model_name, paste(model_name, "_model.R", sep = "")))
  source(file.path("models", model_name, "prior_definitions.R"))

  pop_post_samples = data.matrix(read.csv(posterior_samples_file, header = F))
  dimnames(pop_post_samples) = NULL
  num_post_samples = dim(pop_post_samples)[1]

  ind_param = apply(pop_post_samples, MARGIN = 1, FUN = sample_ind_param)
  post_sim_data = apply(ind_param, MARGIN = 2, FUN = model)

  plot(post_sim_data[, 1], type = "l")
  for (i in 2:num_post_samples) {
    lines(post_sim_data[, i])
  }

  ppc_qu = t(apply(
    post_sim_data,
    MARGIN = 1,
    FUN = quantile,
    probs = c(0, 0.025, 0.05, 0.95, 0.975, 1)
  ))
  lines(ppc_qu[, 1], col = "red", type = "l")
  lines(ppc_qu[, 4], col = "red", type = "l")

  dimnames(ppc_qu) = NULL
  print(ppc_qu)

  if (!is.null(output_file)) {
    write.table(
      ppc_qu,
      file = output_file,
      row.names = FALSE,
      col.names = c("qu0", "qu025", "qu05", "qu95", "qu975", "qu1"),
      sep = ","
    )
  }
}
