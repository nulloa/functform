ASGIndep <- function(y, x, count, group, priors, niter, nchains=3, burnin=niter/10, thin=10){
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
  dat$tb1 <- priors$b1
  dat$tb2 <- priors$b2
  dat$te  <- priors$e
  dat$tm  <- priors$m
  dat$ts1 <- priors$s1
  dat$ts2 <- priors$s2
  
  
  # Set up the model in Jags
  ASGHier = "
  model{
  
  for (i in 1:n) {
  y[i] ~ dbinom(theta[i], num[i])
  logit(theta[i]) <- muTheta[i]
  u[i] = ifelse(x[i] < mu[group[i]], 1, 0)
  muTheta[i] = u[i]*(beta1[group[i]] + (eta[group[i]] - beta1[group[i]])*exp(-(x[i] - mu[group[i]])^2 / (2*sigma1[group[i]]^2))) + (1-u[i])*(beta2[group[i]] + (eta[group[i]]-beta2[group[i]])*exp(-(x[i] - mu[group[i]])^2 / (2*sigma2[group[i]]^2)))
  }
  
  for (g in 1:nG) {
    beta1[g]  ~ dnorm(0, t_b1)
    beta2[g]  ~ dnorm(0, t_b2)
    eta[g]    ~ dnorm(0, t_e)
    mu[g]     ~ dnorm(15, t_m)
    sigma1[g] ~ dt(0, t_s1, 1) T(0,)
    sigma2[g] ~ dt(0, t_s1, 1) T(0,)
  }

  t_b1 <- 1/sqrt(tau_b1)
  tau_b1 ~ dt(0, 1/tb1, 1) T(0,)
  
  t_b2 <- 1/sqrt(tau_b2)
  tau_b2 ~ dt(0, 1/tb2, 1) T(0,)
  
  t_e <- 1/sqrt(tau_e)
  tau_pp ~ dt(0, 1/te, 1) T(0,)
  
  t_m <- 1/sqrt(tau_m)
  tau_m ~ dt(0, 1/tm, 1) T(0,)
  
  t_s1 <- 1/sqrt(tau_s1)
  tau_s1 ~ dt(0, 1/ts1, 1) T(0,)
  
  t_s2 <- 1/sqrt(tau_s2)
  tau_s2 ~ dt(0, 1/ts2, 1) T(0,)
  
  }"
  m = jags.model(textConnection(ASGHier), data=dat, n.chains=nchains, n.adapt=burnin)
  res = coda.samples(m, c("pp","mu","muTheta","beta1","beta2","theta","sigma1","sigma2","tau_b1","tau_b2","tau_pp","tau_m","tau_s1","tau_s2","t_b1","t_b2","t_pp","t_m","t_s1","t_s2"), niter, thin=thin)
  return(res)
}