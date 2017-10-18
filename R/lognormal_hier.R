#' lognormal_hier
#' 
#' Runs Asymmetric Gaussian MCMC with a hierarchical mean structure accross the groups
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
#' priors <- list()
#' priors$mx
#' priors$vmx
#' priors$vm
#' priors$vs
#' 
#'
#' @export


lognormal_hier <- function(y, x, count, group, priors, niter=2000, nchains=3, nclusters=nchains, burnin=niter/2, thin=10){
  # Load Library
  require(R2jags)
  
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
  
  list2env(dat, envir=globalenv() )
  
  # Set up the model in Jags
  m = jags.parallel(data=dat, 
                    inits=NULL,
                    parameters.to.save=c("theta","mu","sigma","t_m","t_s"),
                    model.file = system.file("model", "ln_hier.txt", package = "functform"),
                    n.chains = nchains, 
                    n.iter = niter,
                    n.burnin=burnin,
                    n.thin=thin,
                    n.cluster= nclusters
  )
  return(coda::as.mcmc(m))
}