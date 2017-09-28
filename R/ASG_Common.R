ASGCommon <- function(y, x, count, group, priors, niter, nchains=3, burnin=niter/10, thin=10){
  # Load Library
  require(rjags)
  
  # Setup data for model
  dat = list()
  dat$y     <- y
  dat$x     <- x
  dat$num   <- count
  dat$n     <- length(y)
  dat$nG    <- length(unique(group))
  dat$group <- as.numeric(group)
  # Set priors
  dat$b1 <- priors$b1
  dat$b2 <- priors$b2
  dat$e  <- priors$e
  dat$m  <- priors$m
  dat$s1 <- priors$s1
  dat$s2 <- priors$s2
  
  
  # Set up the model in Jags
  ASGCommon = "
  model{
  
  for (i in 1:n) {
    y[i] ~ dbinom(theta[i], num[i])
    logit(theta[i]) <- muTheta[i]
    u[i] = ifelse(x[i] < mu, 1, 0)
    muTheta[i] = u[i]*(beta1 + (eta - beta1)*exp(-(x[i] - mu)^2 / (2*sigma1^2))) + (1-u[i])*(beta2 + (eta-beta2)*exp(-(x[i] - mu)^2 / (2*sigma2^2)))
  }
  
  beta1  ~ dnorm(0, 1/b1)
  beta2  ~ dnorm(0, 1/b2)
  eta    ~ dnorm(0, 1/e)
  mu     ~ dnorm(15, 1/m)
  sigma1 ~ dt(0, 1/s1, 1) T(0,)
  sigma2 ~ dt(0, 1/s1, 1) T(0,)
  
  }"
  m = jags.model(textConnection(ASGCommon), data=dat, n.chains=nchains, n.adapt=burnin)
  res = coda.samples(m, c("eta","mu","muTheta","beta1","beta2","theta","sigma1","sigma2"), niter, thin=thin)
  return(res)
}