library(mvtnorm)
library(mvnfast)
library(xLLiM)
library(rjson)
library(parallel)
library("Rcpp")
library(progress)
library(adaptMCMC)
library(cmdstanr)

source("algorithm/semple_mem.R")

run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "200ind_precision1_all",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 1000,
  num_gibbs_samples = 10000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10_burnin0",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)

run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "egfp_40ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 1000,
  num_gibbs_samples = 10000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "burnin10",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
