#' asg_indep
#' 
#' Runs Asymmetric Gaussian MCMC with an independent mean structure accross the groups
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
#' priors$vb1
#' priors$vb2
#' priors$mn
#' priors$vn
#' priors$mx
#' priors$vm
#' priors$vs1
#' priors$vs2
#' 
#'
#' @export

asg_indep <- function(y, x, count, group, priors, niter=2000, nchains=3, nclusters=nchains, burnin=niter/2, thin=10){
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
  dat$vb1 <- priors$vb1
  dat$vb2 <- priors$vb2
  dat$mn  <- priors$mn
  dat$vn  <- priors$vn
  dat$mx  <- priors$mx
  dat$vm  <- priors$vm
  dat$vs1 <- priors$vs1
  dat$vs2 <- priors$vs2
  
  list2env(dat, envir=globalenv() )
  
  # Set up the model in Jags
  m = jags.parallel(data=dat, 
                    inits=NULL,
                    parameters.to.save=c("beta1","beta2","nu","mu","sigma1","sigma2","theta"), 
                    model.file = system.file("model", "asg_indep.txt", package = "functform"),
                    n.chains = nchains, 
                    n.iter = niter,
                    n.burnin=burnin,
                    n.thin=thin,
                    n.cluster= nclusters
  )
  return(coda::as.mcmc(m))
}