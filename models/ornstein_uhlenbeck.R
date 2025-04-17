ornstein_uhlenbeck_model = function(theta, x0, dim_data, delta){
  x = rep(NA, dim_data)
  x[1] = x0
  diffusion_coeff = sqrt((theta[3] ^ 2 / (2 * theta[1])) * (1 - exp(-2 * theta[1] * delta)))
  for (i in 2:(dim_data)) {
    x[i] =  theta[2] + (x[i - 1] - theta[2]) * exp(-theta[1] * delta) + diffusion_coeff * rnorm(1)
  }
  return(x)
}
