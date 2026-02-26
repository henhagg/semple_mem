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

# mRNA model with simulated data and fixed t_0 (Fig. 4 in main paper)
run_inference(
  model_name = "mrna_fix_tzero",
  observation_name = "num_observation_40",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 1250,
  num_gibbs_samples = 50000,
  num_rounds = 4,
  K_start = 7,
  mixprob_thresh = 0,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast", # mvnfast or mvtnorm
  subfolder = "50k_R4",
  alg_version = "unperturbed_kappa_xi", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 6,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC",
  iter_hmc_warmup_eta = 1000
)