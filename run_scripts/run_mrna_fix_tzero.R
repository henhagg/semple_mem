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
  model_name = "mrna_fix_tzero",
  observation_name = "num_observation_40",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 1000,
  num_gibbs_samples = 10000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "10k_hmc_randomseed3_K10",
  alg_version = "unperturbed_kappa_xi", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 3,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "conjugate_prior"
)