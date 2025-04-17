data {
  int<lower=1> M;
  vector[M] c;
  real<lower=0> alpha;
  real<lower=0> beta;
  real mu_0;
  real sigma_0;
}

parameters {
  real<lower=0> tau;
  real mu;
}

model {
  target += gamma_lpdf(tau | alpha, beta); // prior log-gamma
  target += normal_lpdf(mu | mu_0, sigma_0); // prior log-normal
  
  for (i in 1:M) {
    target += normal_lpdf(c[i] | mu, 1/sqrt(tau));
  }
}
