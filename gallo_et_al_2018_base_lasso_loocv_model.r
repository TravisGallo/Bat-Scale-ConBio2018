model {
  # process model of occupancy status
  # we just have a random year/session effect here
  for(i in 1:n){
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- beta_0[yrs[i]] 
  }

  
  # random effect for intercept of occupancy
  for(i in 1:nsession){
    beta_0[i] ~ dnorm(int_yr_mu[yr_mu_vec[i]], int_tau_session)
  }
  # the above random effect is nested within the year of sampling
  for(i in 1:nyr){
  	int_yr_mu[i] ~ dnorm(int_mu, int_tau_year)
  }
  
 # hyperparameters for intercept of process model
  int_mu ~ dnorm(0, 3) # mean
  int_tau_session ~ dgamma(1,1) # inverse gamma 
  int_sd_session <- 1 / sqrt(int_tau_session) # sd
  int_tau_year ~ dgamma(1,1) # inverse gamma
  int_sd_year <- 1 / sqrt(int_tau_year) # sd

  # create beta parameters for covariates
  for(i in 1:ncov){
    laplace[i] ~ ddexp(0, lambda) # lasso regression
    #pi[i] ~ dbern(pp) # indicator variable
    beta[i] <- laplace[i] * pi[i] # mixture of priors
  }
  pp~dbeta(4,8)
  # scaling parameter for double exponential priors
  # smaller mean more mass at zero
  lambda ~ dgamma(1,1)
  ####################################
  
  # observation model

  # random session within year effect for detection
  for(i in 1:nsession){
    lp_0[i] ~ dnorm(lp_mu_yr[yr_mu_vec[i]], lp_tau_session)
  }
  # the year effect
  for(i in 1:nyr){
    lp_mu_yr[i] ~ dnorm(lp_mu, lp_tau_year)
  }
  # hyperparameters for detection
  lp_mu ~ dnorm(0,3) # mean logit detection prob across study
  lp_tau_session ~ dgamma(1,1) # inverse gamma
  lp_sd_session <- 1 / sqrt(lp_tau_session) # sd
  lp_tau_year ~ dgamma(1,1) # inverse gamma
  lp_sd_year <- 1 / sqrt(lp_tau_year) # sd
  
  # observation model
  for(i in 1:n){
    logit(detect_prob[i]) <- lp_0[yrs[i]] + inprod(x[sites[i],], beta[])
    y[i] ~ dbin(z[i] * detect_prob[i], jmat[i])
  }
  
  # code for a binomial coefficient
  BinCo <- exp(logfact(jmat[loo]) - (logfact(pred_y) + 
      logfact(jmat[loo] - pred_y)))
  # the likelihood of our held out data given the model
  lik <- ifelse(equals(pred_y,0), 
    psi[loo]*((1-detect_prob[loo])^jmat[loo]) + (1-psi[loo]),
    BinCo*psi[loo]*(detect_prob[loo]^pred_y) * 
      (1-detect_prob[loo])^(jmat[loo]-pred_y))
}