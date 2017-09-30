#' lognormal_hier
#' 
#' Runs Asymmetric Gaussian MCMC with a hierarchical mean structure accross the groups
#'
#' @param y response variable which follows binomial dist
#' @param x explanatory variable
#' @param count n in binomial dist
#' @param group groups of response
#' @param priors list of priors
#' @param niter number of interations to be run
#' @param nchains number of chains to be run (default=3)
#' @param burnin number of samples to be used as burnin (technically adaption, see link below)
#' @param thin when you want to thin (default=10)
#' 
#' @seealso \url{http://www.mikemeredith.net/blog/2016/Adapt_or_burn.htm}
#'
#' @return A MCMC object
#'
#' @examples
#' priors <- list()
#' priors$vm  <- 10
#' priors$mx  <- 15
#' priors$vmx <- 10
#' priors$vs  <- 10
#' 
#'
#' @export


lognormal_hier <- function(y, x, count, group, priors, niter, nchains=3, burnin=niter/10, thin=10){
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
  dat$vm  <- priors$vm
  dat$mx  <- priors$mx
  dat$vmx <- priors$vmx
  dat$vs  <- priors$vs
  
  
  # Set up the model in Jags
  Hier = "
  model{
  
  for (i in 1:n) {
    y[i] ~ dbinom(theta[i], num[i])
    logit(theta[i]) <- ltheta[i]
    ltheta[i] ~ dlnorm(mu[group[i]], sig[group[i]])
  }
  
  for (g in 1:nG) {
    mu[g]    ~ dnorm(m_x, t_m)
    sig[g] <- 1/sigma[g]
    sigma[g] ~ dt(0, t_s, 1) T(0,)
  }
  
  m_x ~ dnorm(mx, 1/vmx)

  t_m <- 1/sqrt(tau_m)
  tau_m ~ dt(0, 1/vm, 1) T(0,)
  
  t_s <- 1/tau_s
  tau_s ~ dt(0, 1/vs, 1) T(0,)

  }"
  m = jags.model(textConnection(Hier), data=dat, n.chains=nchains, n.adapt=burnin)
  res = coda.samples(m, c("theta","ltheta","mu","sigma","tau_m","tau_s","t_m","t_s"), niter, thin=thin)
  return(res)
}