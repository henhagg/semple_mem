############### PRIOR DEFINITIONS ############################
prior_mean_mu = c(-0.694, -3, 0.027)

prior_alpha = c(2, 2, 2)
prior_beta =  c(0.5, 0.5, 0.5) #rate

prior_mean_kappa = c(5, 1, 3)
prior_sd_kappa = c(1, 1, 1)

prior_mean_xi = -1.5
prior_sd_xi = 1

dim_ind_param = 3
dim_eta = 6
dim_kappa = 3
dim_xi = 1
dim_kappa_xi = 4

sample_mu_prior = function(tau) {
  c(
    rnorm(1, mean = prior_mean_mu[1], sd = sqrt(1 / tau[1])),
    rnorm(1, mean = prior_mean_mu[2], sd = sqrt(1 / tau[2])),
    rnorm(1, mean = prior_mean_mu[3], sd = sqrt(1 / tau[3]))
  )
}

sample_tau_prior = function() {
  c(
    rgamma(1, shape = prior_alpha[1], rate = prior_beta[1]),
    rgamma(1, shape = prior_alpha[2], rate = prior_beta[2]),
    rgamma(1, shape = prior_alpha[3], rate = prior_beta[3])
  )
}

sample_kappa_prior = function(){
  c(
    rnorm(1, mean=prior_mean_kappa[1], sd=prior_sd_kappa[1]),
    rnorm(1, mean=prior_mean_kappa[2], sd=prior_sd_kappa[2]),
    rnorm(1, mean=prior_mean_kappa[3], sd=prior_sd_kappa[3])
  )
}

sample_xi_prior = function() {
  rnorm(1, mean = prior_mean_xi, sd = prior_sd_xi)
}

sample_kappa_xi_prior = function() {
  c(sample_kappa_prior(), sample_xi_prior())
}

pdf_kappa_xi_prior = function(kappa_xi) {
  dnorm(kappa_xi[1], mean = prior_mean_kappa[1], sd = prior_sd_kappa[1]) *
    dnorm(kappa_xi[2], mean = prior_mean_kappa[2], sd = prior_sd_kappa[2]) *
    dnorm(kappa_xi[3], mean = prior_mean_kappa[3], sd = prior_sd_kappa[3]) *
    dnorm(kappa_xi[4], mean = prior_mean_xi, sd = prior_sd_xi)
}

sample_eta_prior = function() {
  tau = sample_tau_prior()
  mu = sample_mu_prior(tau)
  return(c(mu,tau))
}

pdf_ind = function(c, eta) {
  dnorm(c[1], mean = eta[1], sd = 1 / sqrt(eta[4])) *
    dnorm(c[2], mean = eta[2], sd = 1 / sqrt(eta[5])) *
    dnorm(c[3], mean = eta[3], sd = 1 / sqrt(eta[6]))
}

sample_ind_param = function(eta) {
  c(
    rnorm(1, mean = eta[1], sd = 1 / sqrt(eta[4])),
    rnorm(1, mean = eta[2], sd = 1 / sqrt(eta[5])),
    rnorm(1, mean = eta[3], sd = 1 / sqrt(eta[6])) 
  )
}