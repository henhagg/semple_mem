data {
  int<lower=1> M;
  vector[M] c;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> lambda;
  real mu_0;
}

parameters {
  real<lower=0> tau;
  real mu;
}

model {
  target += gamma_lpdf(tau | alpha, beta); // prior log-gamma
  target += normal_lpdf(mu | mu_0, 1/sqrt(lambda * tau)); // prior log-normal
  
  for (i in 1:M) {
    target += normal_lpdf(c[i] | mu, 1/sqrt(tau));
  }
}
