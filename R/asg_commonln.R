#' asg_commonln
#' 
#' Runs Asymmetric Gaussian MCMC with a common mean structure accross the groups with log-normal priors
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
#'   priors$vb1
#'   priors$vb2
#'   priors$mn
#'   priors$vn
#'   priors$mx
#'   priors$vm
#'   priors$ms1
#'   priors$vs1
#'   priors$ms2
#'   priors$vs2
#'
#' @export


asg_commonln <- function(y, x, count, group, priors, niter=2000, nchains=3, nclusters=nchains, burnin=niter/2, thin=10){
  # Load Library
  require(R2jags)
  
  # Setup data for model
  dat <- list(y=y, x=x, num=count, n=length(y), nG=length(unique(group)), group=as.numeric(group))
  # Set priors
  dat <- c(dat, priors)
  
  list2env(dat, envir=globalenv() )
  
  # Set up the model in Jags
  m = jags.parallel(data=dat, 
                    inits=NULL,
                    parameters.to.save=c("beta1","beta2","nu","mu","sigma1","sigma2","theta"), 
                    model.file = system.file("model", "asg_commonln.txt", package = "functform"),
                    n.chains = nchains, 
                    n.iter = niter,
                    n.burnin=burnin,
                    n.thin=thin,
                    n.cluster= nclusters
  )
  return(coda::as.mcmc(m))
}