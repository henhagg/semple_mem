source(file = "algorithm/bic.R", local = TRUE)

# ::::::::::::: Ornstein-Uhlenbeck ::::::::::::::::::::::::::::::::
bic_values = compute_bic(
  model_name = "ornstein_uhlenbeck_unperturbed_noise",
  num_priorpred_samples = 10000,
  K_values = c(5, 10, 15),
  cov_structure = "",
  model_param = list(tt=10, dim_data = 51, x0 = 0),
  output_subfolder = "10k"
)

#::::::::::::: mRNA fixed t_0 ::::::::::::::::::::::::::::::::
bic_values = compute_bic(
  model_name = "mrna_fix_tzero",
  num_priorpred_samples = 50000,
  K_values = c(5, 10, 15, 20),
  cov_structure = "",
  model_param = list(tt = 30, dt = 0.01, dim_data = 60),
  output_subfolder = "50k"
)

#::::::::::::: mRNA independent population priors ::::::::::::::::::::::::::::::::
bic_values = compute_bic(
  model_name = "mrna_indep_prior",
  num_priorpred_samples = 10000,
  K_values = c(5, 10, 15, 20),
  cov_structure = "",
  model_param = list(tt = 30, dt = 0.01, dim_data = 180),
  output_subfolder = "10k"
)
