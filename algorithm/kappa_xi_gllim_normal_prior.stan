data {
  int<lower=0> K;
  int<lower=0> D;
  int<lower=0> L;
  int<lower=0> M;
  int<lower=0> dim_ind_param;
  int<lower=0> dim_kappa_xi;
  matrix[D,K] c_gllim;
  matrix[dim_ind_param, M] c_model;
  matrix[L,K] b;
  simplex[K] pi;
  array[K] matrix[L,D] A;
  array[K] matrix[L,L] Sigma;
  array[K] matrix[D,D] Gamma;
  matrix[L,M] y_o;
  vector[dim_kappa_xi] prior_means_kappa_xi;
  vector[dim_kappa_xi] prior_sds_kappa_xi;
}

parameters {
  vector[dim_kappa_xi] kappa_xi;
}

model {
  // priors
  for(j in 1:dim_kappa_xi) {
    target += normal_lpdf(kappa_xi[j] | prior_means_kappa_xi[j], prior_sds_kappa_xi[j]);
  }
  
  // surrogate likelihoods
  for (i in 1:M) {
    vector[D] theta;
    theta = append_row(c_model[:,i], kappa_xi);
    vector[K] log_eta_non_normalized;
    for(k in 1:K) {
      log_eta_non_normalized[k] = log(pi[k]) + multi_normal_lpdf(theta | c_gllim[:,k], Gamma[k]);
    }
    real log_norm_constant_eta = log_sum_exp(log_eta_non_normalized);
    
    vector[K] lp;
    for (k in 1:K) {
      lp[k] = log_eta_non_normalized[k] - log_norm_constant_eta + multi_normal_lpdf(y_o[:,i] | A[k] * theta + b[:, k], Sigma[k]);
    }
    target += log_sum_exp(lp);
  }
}
