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

# mRNA model with real data
run_inference(
  model_name = "mrna_indep_prior",
  observation_name = "egfp_40ind",
  num_priorpred_samples = 50000,
  num_surrogate_post_samples = 1250,
  num_gibbs_samples = 1000,
  num_rounds = 4,
  K_start = 9,
  mixprob_thresh = 0,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 0,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast",
  subfolder = "R4_1ksamples_21_11",
  alg_version = "unperturbed_kappa_xi",
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
