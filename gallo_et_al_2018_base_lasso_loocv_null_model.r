model {
  # process model of occupancy status
  # we just have a random year/session effect here
  for(i in 1:n){
    z[i] ~ dbern(psi)
  }
psi ~ dbeta(1,1)
  ####################################
  
  # observation model

  detect_prob ~ dbeta(1,1)
  # observation model
  for(i in 1:n){
    y[i] ~ dbin(z[i] * detect_prob, jmat[i])
  }
  
  # code for a binomial coefficient
  BinCo <- exp(logfact(jmat[loo]) - (logfact(pred_y) + 
      logfact(jmat[loo] - pred_y)))
  # the likelihood of our held out data given the model
  lik <- ifelse(equals(pred_y,0), 
    psi*((1-detect_prob)^jmat[loo]) + (1-psi),
    BinCo*psi*(detect_prob^pred_y) * 
      (1-detect_prob)^(jmat[loo]-pred_y))
}