#' asg_hierln
#' 
#' Runs Asymmetric Gaussian MCMC with a hierarchical mean structure accross the groups using log-normal hierarchical models
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
#' priors = list()
#' priors$vtb1
#' priors$vtb2
#' priors$vte
#' priors$mx
#' priors$vmx
#' priors$vtm
#' priors$vts1
#' priors$vts2
#' 
#'
#' @export


asg_hierln <- function(y, x, count, group, priors, niter, nchains=3, burnin=niter/10, thin=10){
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
  #betas
  dat$vtb1 <- priors$vtb1
  dat$vtb2 <- priors$vtb2
  # etas
  dat$mn   <- priors$mn
  dat$vmn  <- priors$vmn
  dat$vtn  <- priors$vtn
  # mus
  dat$mx   <- priors$mx
  dat$vmx  <- priors$vmx
  dat$vtm  <- priors$vtm
  # sigmas
  dat$ms1 <- priors$ms1
  dat$vms1 <- priors$vms1
  dat$ms2 <- priors$ms2
  dat$vms2 <- priors$vms2
  dat$vts1 <- priors$vts1
  dat$vts2 <- priors$vts2
  
  
  # Set up the model in Jags
  ASGHier = "
  model{
  
  for (i in 1:n) {
  y[i] ~ dbinom(theta[i], num[i])
  logit(theta[i]) <- ltheta[i]
  u[i] = ifelse(x[i] < mu[group[i]], 1, 0)
  ltheta[i] = u[i]*(beta1[group[i]] + (nu[group[i]] - beta1[group[i]])*exp(-(x[i] - mu[group[i]])^2 / (2*sigma1[group[i]]^2))) + (1-u[i])*(beta2[group[i]] + (nu[group[i]]-beta2[group[i]])*exp(-(x[i] - mu[group[i]])^2 / (2*sigma2[group[i]]^2)))
  }
  
  for (g in 1:nG) {
  beta1[g] ~ dnorm(0, t_b1)
  beta2[g] ~ dnorm(0, t_b2)
  nu[g]   ~ dnorm(m_n, t_n)
  mu[g]    ~ dnorm(m_x, t_m)
  sigma1[g] ~ dlnorm(m_s1, t_s1)
  sigma2[g] ~ dlnorm(m_s2, t_s2)
  }
  
  m_n ~ dnorm(mn, 1/vmn)
  m_x ~ dnorm(mx, 1/vmx)
  m_s1 ~ dnorm(ms1, 1/vms1)
  m_s2 ~ dnorm(ms2, 1/vms2)
  
  t_b1 <- 1/sqrt(tau_b1)
  tau_b1 ~ dt(0, 1/vtb1, 1) T(0,)
  
  t_b2 <- 1/sqrt(tau_b2)
  tau_b2 ~ dt(0, 1/vtb2, 1) T(0,)
  
  t_n <- 1/sqrt(tau_n)
  tau_n ~ dt(0, 1/vtn, 1) T(0,)
  
  t_m <- 1/sqrt(tau_m)
  tau_m ~ dt(0, 1/vtm, 1) T(0,)
  
  t_s1 <- 1/sqrt(tau_s1)
  tau_s1 ~ dt(0, 1/vts1, 1) T(0,)
  
  t_s2 <- 1/sqrt(tau_s2)
  tau_s2 ~ dt(0, 1/vts2, 1) T(0,)
  
  }"
  m = jags.model(textConnection(ASGHier), data=dat, n.chains=nchains, n.adapt=burnin)
  res = coda.samples(m, c("nu","mu","ltheta","beta1","beta2","theta","sigma1","sigma2","tau_b1","tau_b2","tau_n","tau_m","tau_s1","tau_s2","t_b1","t_b2","t_n","t_m","t_s1","t_s2","m_x","m_n","m_s1","m_s2"), niter, thin=thin)
  return(res)
}