#' drake_mcmc
#' 
#' A wrapper function to make running mcmc fns on drake easier
#'
#' @param dat response variable which follows binomial dist
#' @param iter number of interations to be run (default=2000)
#' @param chains number of chains to be run (default=3)
#' @param clusters number of clusters to be used (default=nchains)
#' @param burnin number of samples to be used as burnin (technically adaption, see link below)
#' @param thin when you want to thin (default=10)
#' @param whichmodel characeter of which model you want to run
#' 
#' @seealso \url{http://www.mikemeredith.net/blog/2016/Adapt_or_burn.htm}
#'
#' @return A MCMC object
#'
#' @examples
#' drake_mcmc(d, iter = 1000, chains = 3, clusters=3, burnin = iter/2, thin=10, whichmodel=="asg_common")
#'
#' @export


drake_mcmc <- function(d, iter = 1000, chains = 3, clusters=chains, burnin = iter/2, thin=10, whichmodel=NULL){
    if(is.null(whichmodel)) stop("No model was chosen.")
    require(functform)
    dat <- d[c(1,2,3,4)]
    priors <- d[-c(1,2,3,4)]
    
    if(whichmodel=="asg_common"){
      s <- asg_common(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="asg_common_ln"){
      s <- asg_commonln(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="asg_indep"){
      s <- asg_indep(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="asg_indep_ln"){
      s <- asg_indepln(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="asg_hier"){
      s <- asg_hier(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="asg_hier_ln"){
      s <- asg_hierln(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="ln_common"){
      s <- lognormal_common(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="ln_indep"){
      s <- lognormal_indep(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="ln_hier"){
      s <- lognormal_hier(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="n_common"){
      s <- normal_common(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="n_indep"){
      s <- normal_indep(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    if(whichmodel=="n_hier"){
      s <- normal_hier(dat$y, dat$week, dat$num, dat$group, priors, nchains=chains, nclusters=clusters, niter=iter, burnin=burnin, thin=thin)
    }
    
    return(s)
  }