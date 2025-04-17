library(readxl)
source("plot_scripts/plot_raw_mrna_data.R", local = TRUE)

generate_observations_from_mrna_data = function(egfp_data,
                                                model_name,
                                                individual_indices,
                                                observation_name,
                                                data_set_name,
                                                only_individual_param = FALSE) {
  
  source(file.path("models", model_name, "prior_definitions.R"), local = TRUE)
  
  # create output directory
  output_dir = file.path("models", model_name, observation_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  colnames(egfp_data) = c("time", 1:ncol(egfp_data)) # name columns
  egfp_data$time = egfp_data$time / 3600 # covert time from seconds to hours
  plot_data(egfp_data, individual_indices) # plot raw data
  
  # log data but not time
  egfp_data[, 2:ncol(egfp_data)] = log(egfp_data[, 2:ncol(egfp_data)])
  
  # plot logged data and save pdf to file
  plot_data(egfp_data, individual_indices)
  ggsave(filename = "plot_log_data.pdf", path = output_dir)
  
  # save logged data to file
  data_to_save = egfp_data[, c(1, ((individual_indices) + 1))]
  colnames(data_to_save) = c("time", paste0("obs", individual_indices))
  write.table(
    data_to_save,
    file = file.path(output_dir, "observations.csv"),
    sep = ",",
    row.names = F,
    quote = F
  )
  
  # save model parameters to file
  model_param = list(tt = 30,
                     dt = 0.01,
                     dim_data = 180)
  write.table(
    model_param,
    file.path(output_dir, "model_param.csv"),
    sep = ",",
    row.names = F,
    quote = F
  )
  
  # write "true parameters"=NA to file, this is needed to recover the param names later
  if (only_individual_param) {
    param_names = c(paste0("mu", 1:dim_ind_param),
                    paste0("tau", 1:dim_ind_param))
  } else{
    param_names = c(
      paste0("mu", 1:dim_ind_param),
      paste0("tau", 1:dim_ind_param),
      paste0("kappa", 1:dim_kappa),
      paste0("xi", 1:dim_xi)
    )
  }

  write.table(
    rbind(param_names, rep(NA, length(param_names))),
    file.path(output_dir, "true_param.csv"),
    sep = ",",
    row.names = F,
    col.names = F,
    quote = F
  )
  
  # write data info to file
  write.table(
    list(data_set = data_set_name, individual_indices = individual_indices),
    file.path(output_dir, "data_info.csv"),
    sep = ",",
    row.names = F,
    quote = F
  )
}

path_data = file.path("models", "mrna_data", "20160427_mean_eGFP.xlsx")
egfp_data = read_xlsx(path = path_data, col_names = FALSE)

generate_observations_from_mrna_data(
  egfp_data,
  model_name = "mrna_indep_prior_only_individual_param",
  individual_indices = 1:40,
  observation_name = "egfp_40ind",
  data_set_name = basename(path_data),
  only_individual_param = TRUE
)

generate_observations_from_mrna_data(
  egfp_data,
  model_name = "mrna_indep_prior",
  individual_indices = 1:40,
  observation_name = "egfp_40ind",
  data_set_name = basename(path_data)
)
