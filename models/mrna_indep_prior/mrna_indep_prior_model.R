library(Rcpp)
sourceCpp("models/mrna.cpp")

model = function(log_ind_param, log_kappa_xi, model_param) {
  # param order: delta, gamma, k, t0, m0, scale, offset, sigma
  cpp_mrna_model(
    logparam = c(log_ind_param, log_kappa_xi),
    tt = model_param$tt,
    dt = model_param$dt,
    dim_data = model_param$dim_data,
    observe_time_zero = FALSE # this is how the real world data works
  )
}

jacobian = function(theta_old, theta) {
  jacobian = 1 # no transformation jacobian needed here
  return(jacobian)
}