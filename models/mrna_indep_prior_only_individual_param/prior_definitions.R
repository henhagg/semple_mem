############### PRIOR DEFINITIONS ############################
# true values mu c(-0.694, -3, 0.027, -0.164, 5.704, 0.751, 2.079, -1.6)
prior_mean_mu = c(-1, -5, 0.5, 0, 5, 1, 3, -1)
prior_sd_mu = c(1, 2, 1, 1, 1, 1, 1, 1)

prior_alpha = c(2, 2, 2, 2, 2, 2, 2, 2)
prior_beta =  c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5) #rate

dim_ind_param = 8
dim_eta = 16

sample_mu_prior = function() {
  c(
    rnorm(1, mean = prior_mean_mu[1], sd = prior_sd_mu[1]),
    rnorm(1, mean = prior_mean_mu[2], sd = prior_sd_mu[2]),
    rnorm(1, mean = prior_mean_mu[3], sd = prior_sd_mu[3]),
    rnorm(1, mean = prior_mean_mu[4], sd = prior_sd_mu[4]),
    rnorm(1, mean = prior_mean_mu[5], sd = prior_sd_mu[5]),
    rnorm(1, mean = prior_mean_mu[6], sd = prior_sd_mu[6]),
    rnorm(1, mean = prior_mean_mu[7], sd = prior_sd_mu[7]),
    rnorm(1, mean = prior_mean_mu[8], sd = prior_sd_mu[8])
  )
}

sample_tau_prior = function() {
  c(
    rgamma(1, shape = prior_alpha[1], rate = prior_beta[1]),
    rgamma(1, shape = prior_alpha[2], rate = prior_beta[2]),
    rgamma(1, shape = prior_alpha[3], rate = prior_beta[3]),
    rgamma(1, shape = prior_alpha[4], rate = prior_beta[4]),
    rgamma(1, shape = prior_alpha[5], rate = prior_beta[5]),
    rgamma(1, shape = prior_alpha[6], rate = prior_beta[6]),
    rgamma(1, shape = prior_alpha[7], rate = prior_beta[7]),
    rgamma(1, shape = prior_alpha[8], rate = prior_beta[8])
  )
}

sample_eta_prior = function() {
  tau = sample_tau_prior()
  mu = sample_mu_prior()
  return(c(mu,tau))
}

pdf_ind = function(c, eta) {
  dnorm(c[1], mean = eta[1], sd = 1 / sqrt(eta[9])) *
    dnorm(c[2], mean = eta[2], sd = 1 / sqrt(eta[10])) *
    dnorm(c[3], mean = eta[3], sd = 1 / sqrt(eta[11])) *
    dnorm(c[4], mean = eta[4], sd = 1 / sqrt(eta[12])) *
    dnorm(c[5], mean = eta[5], sd = 1 / sqrt(eta[13])) *
    dnorm(c[6], mean = eta[6], sd = 1 / sqrt(eta[14])) *
    dnorm(c[7], mean = eta[7], sd = 1 / sqrt(eta[15])) *
    dnorm(c[8], mean = eta[8], sd = 1 / sqrt(eta[16]))
}

sample_ind_param = function(eta) {
  c(
    rnorm(1, mean = eta[1], sd = 1 / sqrt(eta[9])),
    rnorm(1, mean = eta[2], sd = 1 / sqrt(eta[10])),
    rnorm(1, mean = eta[3], sd = 1 / sqrt(eta[11])),
    rnorm(1, mean = eta[4], sd = 1 / sqrt(eta[12])),
    rnorm(1, mean = eta[5], sd = 1 / sqrt(eta[13])),
    rnorm(1, mean = eta[6], sd = 1 / sqrt(eta[14])),
    rnorm(1, mean = eta[7], sd = 1 / sqrt(eta[15])),
    rnorm(1, mean = eta[8], sd = 1 / sqrt(eta[16]))
  )
}