#' get_res_df
#' 
#' Use mcmc objects from package to get posterior estimates of theta
#'
#' @param dat data in same format as used in functions
#' @param x mcmc object
#' @param niter number of interations to be ran
#' @param model.name name of the model associated with df
#'
#' @return A data.frame
#'
#' @examples
#' 
#' 
#' 
#'
#' @export

get_res_df <- function(dat, x, niter, model.name){
  # Generate True mean theta 95 Credible intervals and est mean theta lines for each model
  PThetaEst  <- UpQ <- LwQ <- NULL
  PlThetaEst <- UplTQ <- LwlTQ <- NULL
  
  # Get predictions and 95% Credible prediction intervals
  nsim <- (niter/10)/10
  ypred <- UPpred <- LWpred <- PostPredSim <- NULL
  
  for (i in 1:dat$n){
    ####### Prediction and 95% PI for Y ######
    for(sim in 1:nsim){
      PostPredSim[sim] <- rbinom(1, dat$num[i],  x[,c(paste("theta[", i, "]", sep=""))][[1]][(niter/10)-nsim+sim])
    }
    # Model 1 - ASG Indep
    ypred[i]  <- mean(PostPredSim)
    LWpred[i] <- quantile(PostPredSim, probs=c(.025))
    UPpred[i] <- quantile(PostPredSim, probs=c(.975))
    
    ####### Est and 95% CI for Theta ######
    # For Model 1 - ASG Indep
    LwQ[i]      <- quantile(x[,c(paste("theta[", i, "]", sep=""))][[1]],probs=c(.025))
    UpQ[i]      <- quantile(x[,c(paste("theta[", i, "]", sep=""))][[1]],probs=c(.975))
    PThetaEst[i] <- mean(x[,c(paste("theta[", i, "]", sep=""))][[1]])
    
    ####### Est and 95% CI for muTheta ######
    LwlTQ[i]      <- quantile(x[,c(paste("ltheta[", i, "]", sep=""))][[1]],probs=c(.025))
    UplTQ[i]      <- quantile(x[,c(paste("ltheta[", i, "]", sep=""))][[1]],probs=c(.975))
    PlThetaEst[i] <- mean(x[,c(paste("ltheta[", i, "]", sep=""))][[1]])
  }
  
  
  # Create simplified data frame
  df <- data.frame(Prop = dat$y/dat$num,
                   nTot = dat$num,
                   nILI = dat$y,
                   Region = dat$region,
                   Week = dat$week,
                   Model = as.factor(c(rep(paste(model.name, sep=""),dat$n))),
                   PredY = c(ypred),
                   PYLB =  c(LWpred),
                   PYUB =  c(UPpred),
                   EstTheta = c(PThetaEst),
                   ThetaLB =  c(LwQ),
                   ThetaUB =  c(UpQ),
                   EstmuTheta = c(PlThetaEst),
                   muThetaLB =  c(LwlTQ),
                   muThetaUB =  c(UplTQ))
  return(df)
}