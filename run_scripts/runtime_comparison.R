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
  model_name = "mrna_indep_prior",
  observation_name = "10ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 4000,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "unperturbed_kappa_xi", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior",
  observation_name = "20ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 2000,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "unperturbed_kappa_xi", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior",
  observation_name = "40ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 1000,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "unperturbed_kappa_xi", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior",
  observation_name = "100ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 400,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  burnin_kappa_xi = 0,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "unperturbed_kappa_xi", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_kappa_xi = "HMC",
  iter_hmc_warmup_kappa_xi = 1000,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)





run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "10ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 4000,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "20ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 2000,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "40ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 1000,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "100ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 400,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "200ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 200,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
run_inference(
  model_name = "mrna_indep_prior_only_individual_param",
  observation_name = "400ind",
  num_priorpred_samples = 10000,
  num_surrogate_post_samples = 100,
  num_gibbs_samples = 1000,
  num_rounds = 2,
  K_start = 10,
  mixprob_thresh = 0.005,
  cov_structure = "",
  gllim_maxiter = 300,
  factor_cvMH = 1,
  burnin_ind_param = 10,
  mvn_package = "mvnfast", # mvnfast, mvtnorm
  subfolder = "K10",
  alg_version = "only_random_effects", #m_models, single_model, single_skip_surrpost, only_gibbs, unperturbed_xi
  random_seed = 1,
  parallelize_ind_param_MH = FALSE,
  plot_results = TRUE,
  sampler_method_eta = "HMC_indep_normal_gamma",
  iter_hmc_warmup_eta = 1000
)
