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
#' mean.vector = c(beta1, beta2, leta, mu, lsigma1, lsigma2)
#' 
#' priors <- list()
#' priors$mu0  # A 6x1 vector
#' priors$phi0 # A 6x6 inverse var-cov matrix
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
    u[i] = ifelse(x[i] < ctheta[3], 1, 0)
    ltheta[i] = u[i]*(ctheta[1] + (exp(ctheta[3]) - ctheta[1])*exp(-(x[i] - ctheta[4])^2 / (2*exp(ctheta[5]^2)))) + (1-u[i])*(ctheta[2] + (exp(ctheta[3])-ctheta[2])*exp(-(x[i] - ctheta[4])^2 / (2*exp(ctheta[6]^2))))
  }
  
  ctheta[1:6] ~ dmnorm(mu0, phi0) # c(beta1[g], beta2[g], leta[g], mu[g], lsigma1[g], lsigma2[g])
  
  }"
  m = jags.model(textConnection(Multi), data=dat, n.chains=nchains, n.adapt=burnin)
  res = coda.samples(m, c("ctheta","theta","ltheta"), niter, thin=thin)
  return(res)
}