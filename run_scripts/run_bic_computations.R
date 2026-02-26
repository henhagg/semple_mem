source(file = file.path("algorithm", "bic.R"), local = TRUE)

# ::::::::::::: Ornstein-Uhlenbeck ::::::::::::::::::::::::::::::::
bic_values = compute_bic(
  model_name = "ornstein_uhlenbeck_unperturbed_noise",
  num_priorpred_samples = 50000,
  K_values = 2:12,
  cov_structure = "",
  model_param = list(tt=10, dim_data = 51, x0 = 0),
  output_subfolder = "50k_2-12"
)

#::::::::::::: mRNA fixed t_0 (simulated data) ::::::::::::::::::::::::::::::::
bic_values = compute_bic(
  model_name = "mrna_fix_tzero",
  num_priorpred_samples = 50000,
  K_values = 2:10,
  cov_structure = "",
  model_param = list(tt = 30, dt = 0.01, dim_data = 60),
  output_subfolder = "50k_2-10"
)

#::::::::::::: mRNA independent population priors (real data) ::::::::::::::::::::::::::::::::
bic_values = compute_bic(
  model_name = "mrna_indep_prior",
  num_priorpred_samples = 50000,
  K_values = 2:10,
  cov_structure = "",
  model_param = list(tt = 30, dt = 0.01, dim_data = 180),
  output_subfolder = "50k"
)
