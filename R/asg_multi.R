#' asg_multi
#' 
#' Runs Asymmetric Gaussian MCMC with an multivariate, independent mean structure accross the groups
#'
#' @param y response variable which follows binomial dist
#' @param x explanatory variable
#' @param count n in binomial dist
#' @param group groups of response
#' @param priors list of priors
#' @param niter number of interations to be run (default=2000)
#' @param nchains number of chains to be run (default=3)
#' @param nclusters number of clusters to be used (default=nchains)
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
#' priors$C    # A 6x6 precision matrix for mean MVN
#' priors$phi0 # A 6x6 precision matrix for Wish Dist
#' 
#'
#' @export

asg_multi <- function(y, x, count, group, priors, niter=2000, nchains=3, nclusters=nchains, burnin=niter/2, thin=10){
  # Load Library
  require(R2jags)
  
  # Setup data for model
  dat = list()
  dat$y     <- y
  dat$x     <- x
  dat$num   <- count
  dat$n     <- length(y)
  dat$group <- as.numeric(group)
  dat$nG    <- length(unique(group))
  
  # Set priors
  dat$mu0  <- priors$mu0
  dat$C    <- priors$C
  dat$phi0 <- priors$phi0
  
  list2env(dat, envir=globalenv() )

  # Set up the model in Jags
  m = jags.parallel(data=dat, 
                    inits=NULL,
                    parameters.to.save=c("ctheta","theta","mu.theta","phi.theta"), 
                    model.file = system.file("model", "multi_hier.txt", package = "functform"),
                    n.chains = nchains, 
                    n.iter = niter,
                    n.burnin=burnin,
                    n.thin=thin,
                    n.cluster= nclusters
  )
  return(coda::as.mcmc(m))
}