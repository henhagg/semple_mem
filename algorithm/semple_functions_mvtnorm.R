########################### HELP FUNCTIONS #####################################
fit_gllim = function(responses,covariates,K,maxiter,verbose,cov_structure,mixprob_thresh,init_param){
  #remove potential NA values before fitting
  col_index_na = which(colSums(is.na(covariates))!=0)
  if(length(col_index_na) > 0){
    print("col_index_na")
    print(col_index_na)
    responses = responses[,-col_index_na]
    covariates = covariates[,-col_index_na]
  }
  
  #remove potential Inf values before fitting
  col_index_inf = which(colSums(is.infinite(covariates))!=0)
  if(length(col_index_inf) > 0){
    print("col_index_inf")
    print(col_index_inf)
    responses = responses[,-col_index_inf]
    covariates = covariates[,-col_index_inf]
  }
  
  #remove potential NaN values before fitting
  col_index_nan = which(colSums(is.nan(covariates))!=0)
  if(length(col_index_nan) > 0){
    print("col_index_nan")
    print(col_index_nan)
    responses = responses[,-col_index_nan]
    covariates = covariates[,-col_index_nan]
  }
  
  # fit prior predictive training data via GLLIM
  mod = gllim(responses,covariates,in_K=K,maxiter=maxiter,in_theta=init_param,verb=verbose,cstr=list(Sigma=cov_structure))
  
  # retrieve MLE parameters for the mixture of f(sims|allpar)  (easier to fit than f(allpar|sims))
  A = mod$A
  b = mod$b
  Sigma = mod$Sigma
  c = mod$c
  Gamma = mod$Gamma
  prior_mix_prob = mod$pi  # these are the prior mixture probabilities denote with \pi_k in Forbes et al
  K = length(prior_mix_prob)
  
  # check if exist components having prior_mix_prob < mixprob_thresh for numerical stability
  id_eliminate = which(prior_mix_prob<mixprob_thresh)
  if(length(id_eliminate)>0){ # if >0 there exist components having small mixture probabilities
    K = K - length(id_eliminate)  # reduce number of components accordingly
    warning('mixprob_thresh reduced the mixture number of components from ', K + length(id_eliminate),  ' to ',K)
    A = A[,,-id_eliminate]
    b = b[,-id_eliminate]
    Sigma = Sigma[,,-id_eliminate]
    c = c[,-id_eliminate]
    Gamma = Gamma[,,-id_eliminate]
    prior_mix_prob = prior_mix_prob[-id_eliminate]
    prior_mix_prob = prior_mix_prob/sum(prior_mix_prob) # normalize after removing components
  }
  
  # store the MLEs so we can reuse those as starting values at next optimization round
  optpar=NULL
  optpar$A = A
  optpar$b = b
  optpar$Sigma = Sigma
  optpar$c = c
  optpar$Gamma = Gamma
  optpar$pi = prior_mix_prob
  
  return(list(K=K, A=A, b=b, Sigma=Sigma, c=c, Gamma=Gamma, prior_mix_prob=prior_mix_prob, pi = prior_mix_prob, optpar=optpar))
}

# compute mixture posterior probabilities
compute_post_mix_prob = function(gllim_param, observed_data){
  K = gllim_param$K
  A = gllim_param$A
  b = gllim_param$b
  Sigma = gllim_param$Sigma
  c = gllim_param$c
  Gamma = gllim_param$Gamma
  pi = gllim_param$prior_mix_prob
  
  post_mix_prob = vector(length=K)
  for (k in 1:K){
    if(K==1){
      c_star = A%*%c + b
      Gamma_star = Sigma + A%*%Gamma%*%t(A)
      post_mix_prob = pi*mvtnorm::dmvnorm(observed_data, mean=c_star, sigma=Gamma_star)
    }
    else{
      c_star_k = A[,,k]%*%c[,k] + b[,k]
      Gamma_star_k = Sigma[,,k] + A[,,k]%*%Gamma[,,k]%*%t(A[,,k])
      post_mix_prob[k] = pi[k]*mvtnorm::dmvnorm(observed_data, mean=c_star_k, sigma=Gamma_star_k)
    }
  }
  post_mix_prob = post_mix_prob/sum(post_mix_prob)
  return(post_mix_prob)
}

# compute mean of mixture proposal distribution
compute_mean_proposal = function(post_param_list, K, observed_data){
  A_star = post_param_list$A_star
  b_star = post_param_list$b_star
  
  mean_proposal = array(NA, dim=c(dim(b_star)[1], K))
  for (k in 1:K){
    if(K==1){
      mean_proposal = A_star%*%observed_data + b_star
    }
    else{
      mean_proposal[,k] = A_star[,,k]%*%observed_data + b_star[,k]
    }
  }
  return(mean_proposal)
}

# compute mixture posterior parameters mean and covariance parameters A_star, b_star and Sigma_star
compute_post_param = function(gllim_param){
  K = gllim_param$K
  A = gllim_param$A 
  b = gllim_param$b
  Sigma = gllim_param$Sigma
  c = gllim_param$c
  Gamma = gllim_param$Gamma
  
  dim_data = dim(A)[1]
  dim_param = dim(A)[2]
  
  A_star = array(NA,dim=c(dim_param,dim_data,K)) # L x d x K
  b_star = array(NA,dim=c(dim_param,K))  # L x K
  Sigma_star = array(NA,dim=c(dim_param,dim_param,K)) # L x L x K
  
  for (k in 1:K){
    if(K==1){
      Sigma_star = solve(solve(Gamma) + t(A)%*%solve(Sigma)%*%A)
      A_star = Sigma_star %*% t(A) %*% solve(Sigma)
      b_star = Sigma_star %*% (solve(Gamma)%*%c - t(A)%*%solve(Sigma)%*%b)
    }
    else
    {
      Sigma_star_k = solve(solve(Gamma[,,k]) + t(A[,,k])%*%solve(Sigma[,,k])%*%A[,,k])
      Sigma_star[,,k] = Sigma_star_k
      A_star[,,k] = Sigma_star_k %*% t(A[,,k]) %*% solve(Sigma[,,k])
      b_star[,k] = Sigma_star_k %*% (solve(Gamma[,,k])%*%c[,k] - t(A[,,k])%*%solve(Sigma[,,k])%*%b[,k])
    }
    
  }
  return(list(A_star=A_star, b_star=b_star, Sigma_star=Sigma_star))
}

# function returning a sample from a GMM with inflated covariance. 
sample_gmm_inflated = function(probs, means, covars, factor_cvMH=1){
  numcomp = length(probs)
  
  k = which(rmultinom(1,1,probs)==1)  # samples an index from the multinomial distribution with associated probabilities probs
  if(numcomp>1)
  {
    mvtnorm::rmvnorm(1, mean=means[,k], sigma=factor_cvMH*covars[,,k])
  }else{
    mvtnorm::rmvnorm(1, mean=means[,k], sigma=factor_cvMH*covars)
  }
}

# compute mean of surrogate likelihood 
mean_surrogate_likelihood = function(A, b, observed_data, theta, K){
  mean_surrogate_likelihood = array(NA, dim=c(length(observed_data), K))
  for (k in 1:K){
    if(K==1){
      mean_surrogate_likelihood = A%*%theta + b
    }
    else{
      mean_surrogate_likelihood[,k] = A[,,k]%*%theta + b[,k]
    }
  }
  return(mean_surrogate_likelihood)
}

# compute surrogate likelihood mixture probabilities
likelihood_mix_prob = function(c, Gamma, pi, K, theta){
  likelihood_mix_prob = vector(length=K)
  for (k in 1:K){
    if(K==1){
      likelihood_mix_prob = pi*mvtnorm::dmvnorm(theta, mean=c, sigma=Gamma)
    }
    else{
      likelihood_mix_prob[k] = pi[k]*mvtnorm::dmvnorm(theta, mean=c[,k], sigma=Gamma[,,k])
    }
  }
  likelihood_mix_prob = likelihood_mix_prob/sum(likelihood_mix_prob)
  return(likelihood_mix_prob)
}

mixpdf = function(value,probs,means,covars){
  pdf = 0
  numcomp = length(probs)
  if(numcomp>1){
    for(k in 1:numcomp){
      pdf = pdf + probs[k]*mvtnorm::dmvnorm(value,mean=means[,k],sigma=covars[,,k])
    }
  }else{
    pdf = mvtnorm::dmvnorm(value,mu=means,sigma=covars)
  }
  return(pdf)
}