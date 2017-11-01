#' asg_indepln5
#' 
#' Runs Asymmetric Gaussian MCMC with an independent mean structure accross the groups using log-normal priors with only one intercept
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
#' @param mle expects TRUE/FALSE, if TRUE will use maximum likelihood to get starting values. Needs 3 chains.
#' @param inits Add specific initial values
#' 
#' @seealso \url{http://www.mikemeredith.net/blog/2016/Adapt_or_burn.htm}
#'
#' @return A MCMC object
#'
#' @examples
#'   priors$vb1
#'   priors$mn
#'   priors$vn
#'   priors$mx
#'   priors$vm
#'   priors$ms1
#'   priors$vs1
#'   priors$ms2
#'   priors$vs2
#' 
#'
#' @export

asg_indepln5 <- function(y, x, count, group, priors, niter=2000, nchains=3, nclusters=nchains, burnin=niter/2, thin=10, inits=NULL, mle=FALSE){
  # Load Library
  require(R2jags)
  
  # Setup data for model
  dat <- list(y=y, x=x, num=count, n=length(y), nG=length(unique(group)), group=as.numeric(group))
  # Set priors
  dat <- c(dat, priors)
  
  if(isTRUE(mle)){
    
    asg <- Vectorize(function(x, beta1, mu, h, sigma1, sigma2){
      top <- beta1 + (h - beta1)*exp(-((x - mu)^2)/(2*sigma1^2))
      bot <- beta1 + (h - beta1)*exp(-((x - mu)^2)/(2*sigma2^2))
      ifelse(x < mu,return(top),return(bot))
    })
    
    log.lik <- function(y, num, x, par){
      theta <- boot::inv.logit(asg(x, par[1], par[2], par[3], par[4], par[5]))
      ll <- sum(lchoose(num, y)) + sum(y*log(theta)) + sum((num-y)*log(1-theta))
      if(par[2] <= .Machine$double.xmin | par[4] <= .Machine$double.xmin | par[5] <= .Machine$double.xmin){ll <- -Inf}
      return(-ll)
    }
    
    c1 <- c2 <- c3 <- matrix(data=NA, ncol=5, nrow=length(unique(dat$group)))
    for(g in 1:length(unique(dat$group))){
      suby   <- dat$y[dat$group==g]
      subx   <- dat$x[dat$group==g]
      subnum <- dat$num[dat$group==g]
      t <- try(optim(par=c(-5, 12, -2.5, 10, 10), log.lik, y=suby, x=subx, num=subnum, hessian=TRUE))
      if("try-error" %in% class(t)){
        tt <- try(optim(par=c(-5, 12, -2.5, 1, 10), log.lik, y=suby, x=subx, num=subnum, hessian=TRUE))
        if("try-error" %in% class(tt)){
          fit <- optim(par=c(-5, 12, -2.5, 1, 2), log.lik, y=suby, x=subx, num=subnum, hessian=TRUE)
        }
        else{
          fit <- optim(par=c(-5, 12, -2.5, 1, 10), log.lik, y=suby, x=subx, num=subnum, hessian=TRUE)
        }
      }else{
        fit <- optim(par=c(-5, 12, -2.5, 10, 10), log.lik, y=suby, x=subx, num=subnum, hessian=TRUE)
      }
      fisher_info <- MASS::ginv(fit$hessian)
      prop_sigma  <- sqrt(abs(diag(fisher_info)))
      prop_sigma[prop_sigma==0] <- 1
      upper <- fit$par+1.96*prop_sigma
      lower <- fit$par-1.96*prop_sigma
      c1[g,] <- fit$par
      c2[g,] <- c(lower[1:3], abs(lower[4:5]))
      c3[g,] <- c(upper[1:3], abs(upper[4:5]))
    }
    
    init <- list(list("beta1"=c1[,1],"mu"=c1[,2],"nu"=c1[,3],"sigma1"=c1[,4],"sigma2"=c1[,5]),
                 list("beta1"=c2[,1],"mu"=c2[,2],"nu"=c2[,3],"sigma1"=c2[,4],"sigma2"=c2[,5]),
                 list("beta1"=c3[,1],"mu"=c3[,2],"nu"=c3[,3],"sigma1"=c3[,4],"sigma2"=c3[,5])
    )
    rm(c1,c2,c3, upper, lower, suby, subx, subnum, fit, fisher_info, prop_sigma, asg, log.lik)
  }else{init=inits}
  
  list2env(dat, envir=globalenv() )
  
  # Set up the model in Jags
  m = jags(data=dat, 
                    inits=init,
                    parameters.to.save=c("beta1","nu","mu","sigma1","sigma2","theta"), 
                    model.file = system.file("model", "asg_indepln5.txt", package = "functform"),
                    n.chains = nchains, 
                    n.iter = niter,
                    n.burnin=burnin,
                    n.thin=thin
                    #n.cluster= nclusters
  )
  return(coda::as.mcmc(m))
}