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

# Fig. S13 in supplementary material
run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "200ind_precision1_all",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 250,
  num_gibbs_samples = 10000,
  num_rounds = 4,
  K_start = 10,
  mixprob_thresh = 0,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "R4_10k_gibbs_21_11",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)

# Fig. S15 in supplementary material
run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "egfp_40ind",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 1250,
  num_gibbs_samples = 10000,
  num_rounds = 4,
  K_start = 10,
  mixprob_thresh = 0,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "R4_10k_25_11",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
