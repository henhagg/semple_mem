library(mvtnorm)
library(mvnfast)
library(xLLiM)
library(rjson)
library(parallel)
library("Rcpp")
library(progress)
library(adaptMCMC)
library(cmdstanr)

source(file.path("algorithm", "semple_mem.R"))

# Main result in Fig. 2
run_inference(
  model_name = "ornstein_uhlenbeck_unperturbed_noise",
  observation_name = "num_observation_40",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 1250,
  num_gibbs_samples = 50000,
  num_rounds = 4,
  K_start = 10,
  mixprob_thresh = 0,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast",
  subfolder = "50k_samples_R4_seed1_HMC_attempt2",
  alg_version = "unperturbed_kappa_xi",
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC"
)

# Fig. S5 in supplementary material (same as above but with 20 observations instead of 40)
run_inference(
  model_name = "ornstein_uhlenbeck_unperturbed_noise",
  observation_name = "num_observation_20",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 2500,
  num_gibbs_samples = 50000,
  num_rounds = 4,
  K_start = 10,
  mixprob_thresh = 0,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast",
  subfolder = "50k_samples_R4_HMC",
  alg_version = "unperturbed_kappa_xi",
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC"
)

# Fig. S6 in supplementary material (same as above but with 100 observations)
run_inference(
  model_name = "ornstein_uhlenbeck_unperturbed_noise",
  observation_name = "num_observation_100",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 500,
  num_gibbs_samples = 50000,
  num_rounds = 4,
  K_start = 10,
  mixprob_thresh = 0,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast",
  subfolder = "50k_samples_R4_seed1_HMC_21_11",
  alg_version = "unperturbed_kappa_xi",
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC"
)
