#:::: Generative model ::::::::::::::::::::::::::::::::
source(file.path("models", "ornstein_uhlenbeck.R"), local = T)
model = function(log_ind_param, log_kappa_xi, model_param) {
  ind_param = exp(log_ind_param)
  
  x = ornstein_uhlenbeck_model(
    theta = ind_param,
    x0 = model_param$x0,
    dim_data = model_param$dim_data,
    delta = model_param$tt / (model_param$dim_data - 1)
  )
  
  y = x + rnorm(model_param$dim_data,
                mean = 0,
                sd = exp(log_kappa_xi))
  return(y)
}

#::::::::: TRANSFORMATION JACOBIAN ::::::::::::::::::::::
jacobian = function(theta_old, theta) {
  jacobian = 1 # no transformation jacobian needed here
  return(jacobian)
}