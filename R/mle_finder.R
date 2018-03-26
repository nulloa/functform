#' mle_finder
#' 
#' A function which finds and formats the mle of ASG ff
#' 
#' @param y response variable which follows binomial dist
#' @param x explanatory variable
#' @param count n in binomial dist
#' @param group groups of response
#' @param season seasons in response
#' @param inits initial values for optim (defualt=NULL)
#' @param ll likelihood to be optimized (defualt=ASG)
#'
#' @return A function
#'
#' @examples
#' asg_mle(y, x, num, group, seas)
#'
#' @export


asg_mle <- function(y, x, num, group, season, inits=NULL, ll=NULL){
  
  asg <- Vectorize(function(x, beta1, beta2, mu, h, sigma1, sigma2){
    top <- beta1 + (h - beta1)*exp(-((x - mu)^2)/(2*sigma1^2))
    bot <- beta2 + (h - beta2)*exp(-((x - mu)^2)/(2*sigma2^2))
    ifelse(x < mu,return(top),return(bot))
  })
  
  log.lik <- function(y, num, x, par){
    theta <- boot::inv.logit(asg(x, par[1], par[2], par[3], par[4], par[5], par[6]))
    ll <- sum(lchoose(num, y)) + sum(y*log(theta)) + sum((num-y)*log(1-theta))
    if(par[3] <= 0 | par[5] <= 0 | par[6] <= 0){ll <- -Inf}
    return(-ll)
  }
  
  if(is.null(ll)){ll <- log.lik}
  if(is.null(inits)){inits <- c(-4.5, -4.5, 15, -2.5, 4, 5)}
  
  df <- data.frame(y=y, x=x, num=num, group=group, season=season)
  opt <- array(NA, dim=c(length(unique(df$group)), length(unique(df$season)), 6))
  
  for(s in 1:length(unique(df$season))){
    for(g in 1:length(unique(df$group))){
      dat <- subset(df, season==paste(unique(df$season)[s], sep="") & group==paste(unique(df$group)[g], sep=""))
      opts <- optim(par=inits, ll, y=dat$y, x=dat$x, num=dat$num)
      opt[unique(df$group)[g], unique(df$season)[s], ] <- c(opts$par[1], opts$par[2], opts$par[4], opts$par[3], opts$par[5], opts$par[6])
    }
  }
  
  return(opt)
}
