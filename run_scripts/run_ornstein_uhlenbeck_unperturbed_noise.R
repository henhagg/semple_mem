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
  model_name = "ornstein_uhlenbeck_unperturbed_noise",
  observation_name = "num_observation_40",
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
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast",
  subfolder = "10k_samples_hmc",
  alg_version = "unperturbed_kappa_xi",
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC"
)
