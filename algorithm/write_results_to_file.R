write_prior_pred_to_file = function(output_dir,
                                    time_vec,
                                    pop_param,
                                    ind_param,
                                    param_names,
                                    kappa_xi_param = NULL,
                                    kappa_exists = NULL,
                                    xi_exists = NULL) {
  write.table(
    cbind(time_vec, c(0, rep(NA, length(time_vec)-1))),
    file = file.path(output_dir, "elapsed_time.csv"),
    sep = ",",
    row.names = F,
    col.names = c("seconds", "round"),
    quote = F
  )
  write.table(
    t(pop_param),
    file = file.path(output_dir, "eta_round0.csv"),
    sep = ",",
    row.names = F,
    col.names = param_names[1:nrow(pop_param)],
    quote = F
  )
  write.table(
    t(ind_param),
    file = file.path(output_dir, "ind_param_round0.csv"),
    sep = ",",
    row.names = F,
    col.names = paste0("c", 1:nrow(ind_param)),
    quote = F
  )
  if(!is.null(kappa_xi_param)){
    write_kappa_xi_to_csv(
      output_dir = output_dir,
      kappa_xi_matrix = kappa_xi_param,
      kappa_exists = kappa_exists,
      xi_exists = xi_exists,
      round_index = 0
    )
  }
}

write_gibbs_results_to_csv = function(time_vec,
                                      ind_param,
                                      dim_ind_param,
                                      M,
                                      num_gibbs_samples,
                                      output_dir,
                                      eta_matrix,
                                      acceptance_rates_ind_param,
                                      acceptance_rates_kappa_xi,
                                      round_index,
                                      param_names,
                                      kappa_xi_matrix = NULL,
                                      kappa_exists = NULL,
                                      xi_exists = NULL) {
  write_ind_param_to_csv(ind_param,
                         dim_ind_param,
                         M,
                         num_gibbs_samples,
                         output_dir,
                         round_index)
  
  round_vector = rep(NA, length(time_vec))
  round_vector[1:(round_index+1)] = 0:round_index
  write.table(
    cbind(time_vec, round_vector),
    file = file.path(output_dir, "elapsed_time.csv"),
    sep = ",",
    row.names = F,
    col.names = c("seconds", "round"),
    quote = F
  )
  write.table(
    t(eta_matrix),
    file = file.path(output_dir, paste0("eta_round", round_index, ".csv")),
    sep = ",",
    row.names = F,
    col.names = param_names[1:nrow(eta_matrix)],
    quote = F
  )
  write.table(
    acceptance_rates_ind_param,
    file = file.path(output_dir, paste0("acceptance_rates_ind_param_round", round_index, ".csv")),
    sep = ",",
    row.names = F,
    col.names = paste0("ind", 1:ncol(acceptance_rates_ind_param)),
    quote = F
  )

  if (!is.null(kappa_xi_matrix)) {
    write_kappa_xi_to_csv(
      output_dir = output_dir,
      kappa_xi_matrix = kappa_xi_matrix,
      kappa_exists = kappa_exists,
      xi_exists = xi_exists,
      round_index = round_index
    )
    write.table(
      acceptance_rates_kappa_xi,
      file = file.path(output_dir, paste0("acceptance_rates_kappa_xi.csv")),
      sep = ",",
      row.names = F,
      col.names = paste0("round", 2:(ncol(acceptance_rates_kappa_xi) + 1)),
      quote = F
    )
  }
}

write_ind_param_to_csv = function(ind_param,
                                  dim_ind_param,
                                  M,
                                  num_gibbs_samples,
                                  output_dir,
                                  round_index) {
  dim(ind_param) = c(dim_ind_param, M * num_gibbs_samples)
  df = data.frame(t(ind_param))
  colnames(df) = paste0("c", 1:dim_ind_param)
  df$id = rep(1:M, num_gibbs_samples)
  df$gibbs_cycle = rep(1:num_gibbs_samples, rep(M, num_gibbs_samples))
  df = df[order(df$id),]
  write.csv(
    df,
    file = file.path(output_dir, paste0("ind_param_round", round_index, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
}

write_kappa_xi_to_csv = function(output_dir,
                                 kappa_xi_matrix,
                                 kappa_exists,
                                 xi_exists,
                                 round_index) {
  if(kappa_exists && xi_exists){
    col.names = c(paste0("kappa", 1:dim_kappa), paste0("xi", 1:dim_xi))
  } else if(kappa_exists && !xi_exists){
    col.names = paste0("kappa", 1:dim_kappa)
  } else {
    col.names = paste0("xi", 1:dim_xi)
  }
  write.table(
    t(kappa_xi_matrix),
    file = file.path(output_dir, paste0("kappa_xi_round", round_index, ".csv")),
    sep = ",",
    row.names = F,
    col.names = col.names,
    quote = F
  )
}