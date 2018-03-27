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
#' @param inits expects TRUE/FALSE, if TRUE will use maximum likelihood to get starting values. Needs 3 chains.
#' @param whichmodel characeter of which model you want to run
#' 
#' @seealso \url{http://www.mikemeredith.net/blog/2016/Adapt_or_burn.htm}
#'
#' @return A MCMC object
#'
#' @examples
#' drake_mcmc(d, iter = 1000, chains = 3, clusters=3, warmup = iter/2, thin=1, whichmodel=="asg_common")
#'
#' @export


drake_mcmc <- function(d, iter = 4000, warmup=iter/2, chains = 3, thin=1, whichmodel=NULL){
    if(is.null(whichmodel)) stop("No model was chosen.")
    require(functform)
    dat <- d$dat
    priors <- d$priors
    inits <- d$inits

    
    if(whichmodel=="asg_common"){
      s <- asg_common(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
    }
    if(whichmodel=="asg_indep"){
      s <- asg_indep(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
    }
    if(whichmodel=="asg_simple_cs"){
      s <- asg_simple_cs(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
    }
    if(whichmodel=="asg_group_cs"){
      s <- asg_group_cs(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
    }
    if(whichmodel=="asg_simple_sc"){
      s <- asg_simple_sc(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
    }
    if(whichmodel=="asg_group_sc"){
      s <- asg_group_sc(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
    }
    if(whichmodel=="asg_mix"){
      s <- asg_mix(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
    }
    if(whichmodel=="asg_mix2"){
      s <- asg_mix2(dat$y, dat$x, dat$num, dat$group, dat$seas, priors, nchains=chains, nwarmup=warmup, niter=iter, thin=thin, inits=inits)
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