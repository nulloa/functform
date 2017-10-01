#' asg_multi
#' 
#' Runs Asymmetric Gaussian MCMC with an multivariate, independent mean structure accross the groups
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
#' priors$vb1
#' priors$vb2
#' priors$ve
#' priors$mx
#' priors$vm
#' priors$vs1
#' priors$vs2
#' 
#'
#' @export

asg_multi <- function(y, x, count, group, priors, niter, nchains=3, burnin=niter/10, thin=10){
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
  dat$mu0  <- priors$mu0
  dat$phi0 <- priors$phi0

  
  
  # Set up the model in Jags
  Multi = "
  model{
  
  for (i in 1:n) {
    y[i] ~ dbinom(theta[i], num[i])
    logit(theta[i]) <- ltheta[i]
    u[i] = ifelse(x[i] < mu[group[i]], 1, 0)
    ltheta[i] = u[i]*(beta1[group[i]] + (eta[group[i]] - beta1[group[i]])*exp(-(x[i] - mu[group[i]])^2 / (2*sigma1[group[i]]^2))) + (1-u[i])*(beta2[group[i]] + (eta[group[i]]-beta2[group[i]])*exp(-(x[i] - mu[group[i]])^2 / (2*sigma2[group[i]]^2)))
  }
  
  for (g in 1:nG) {
    log(eta[g]) <- leta[g]
    log(sigma1[g]) <- lsigma1[g]
    log(sigma2[g]) <- lsigma2[g]
    ctheta[g,] <- c(beta1[g], beta1[g], leta[g], mu[g], lsigma1[g], lsigma2[g])
    ctheta[g, 1:6] ~ dmnorm(mu0, phi0)
  }
  
  }"
  m = jags.model(textConnection(Multi), data=dat, n.chains=nchains, n.adapt=burnin)
  res = coda.samples(m, c("eta","mu","ltheta","beta1","beta2","theta","sigma1","sigma2"), niter, thin=thin)
  return(res)
}