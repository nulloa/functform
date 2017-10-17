#' asg_hierln
#' 
#' Runs Asymmetric Gaussian MCMC with a hierarchical mean structure accross the groups using log-normal hierarchical models
#'
#' @param y response variable which follows binomial dist
#' @param x explanatory variable
#' @param count n in binomial dist
#' @param group groups of response
#' @param priors list of priors
#' @param niter number of interations to be run (default=2000)
#' @param nchains number of chains to be run (default=3)
#' @param nclusters number of clusters to be used (default=nchains)
#' @param burnin number of samples to be used as burnin (technically adaption, see link below) (default=niter/2)
#' @param thin when you want to thin (default=10)
#' 
#' @seealso \url{http://www.mikemeredith.net/blog/2016/Adapt_or_burn.htm}
#'
#' @return A MCMC object
#'
#' @examples
#' priors = list()
#' #betas
#' priors$vtb1
#' priors$vtb2
#' # etas
#' priors$mn
#' priors$vmn
#' priors$vtn
#' # mus
#' priors$mx
#' priors$vmx
#' priors$vtm
#' # sigmas
#' priors$ms1
#' priors$vms1
#' priors$ms2
#' priors$vms2
#' priors$vts1
#' priors$vts2
#' 
#'
#' @export


asg_hierln <- function(y, x, count, group, priors, niter=2000, nchains=3, ncluster=nchains, burnin=niter/2, thin=10){
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
  
  list2env(dat, envir=globalenv() )
  
  # Set up the model in Jags
  m = jags.parallel(data=dat, 
                    inits=NULL,
                    parameters.to.save=c("beta1","beta2","nu","mu","sigma1","sigma2","theta","t_b1","t_b2","t_n","t_m","t_s1","t_s2","m_x","m_n","m_s1","m_s2"), 
                    model.file = "inst/model/asg_hierln.txt",
                    n.chains = nchains, 
                    n.iter = niter,
                    n.burnin=burnin,
                    n.thin=thin,
                    n.cluster= nclusters
  )
  return(coda:as.mcmc(m))
}