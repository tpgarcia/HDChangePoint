#' Logistic model for the parametric nonlinear mixed effects model (nlme) procedure
#'
#' @param theta1 a scale parameter, which determines steepness of the logistic model.
#' @param theta2 a parameter for an inflection point.
#' @param theta3 a parameter, which determines the maximum  (asymptote) of the logistic model.
#' @param logage a vector of time points (log-scaled ages), which is generated around \code{theta2}. It can be an any length of numeric vector.
#'
#'
#'
#' @return a vector of response variables under the deterministic logistic model, whose length is the same as the length of \code{logage}.
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ## Specify parameters
#' theta1=6;
#' theta2=0;
#' theta3=1;
#'
#'
#' ## Specify the time points (log-scaled ages)
#' logage<-seq(-0.5, 0.5, by=0.1)
#'
#' ## Obtain a vector of response variables when parametric NLME assumes the logistic model
#' logit.ft<-logistft(6, 0, 1, logage)
#'
#' ## Plot the S-shaped curve under the logistic model
#' plot(logage, logit.ft, type='l')
#'
#'
logistft<- function(theta1,theta2, theta3, logage){

  ###########################################
  ## nonlinear parametrized omega function ##
  ###########################################

  denom <- 1+exp(-theta1*(logage-theta2))
  denom2 <- denom*denom

  logist <- theta3/denom

  return(logist)

  #  analytical dervivatives
  logistgrad <- array(0,c(length(logage),3),list(NULL,c("theta1","theta2", "theta3")))
  logistgrad[,"theta1"] <- ((logage-theta2)*(theta3)*(denom-1)/denom2)
  logistgrad[,"theta2"] <- -(theta1*theta3*(denom-1)/denom2)
  logistgrad[,"theta3"] <- (1/denom)

  ## attributes of function logist
  attr(logist,"gradient") <- logistgrad


}



#' First derivative of the logistic model function for the parametric nonlinear mixed effects model (nlme) procedure
#'
#'
#' @param theta1 a scale parameter, which determines steepness of the logistic model.
#' @param theta2 a parameter for an inflection point.
#' @param theta3 a parameter, which determines the maximum  (asymptote) of the logistic model.
#' @param logage a vector of time points (log-scaled ages), which is generated around \code{theta2}. It can be an any length of numeric vector.
#'
#' @return a vector of the first derivative of the logistic function with respect to (logage-theta2).
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ## Specify parameters
#' theta1=6;
#' theta2=0;
#' theta3=1;
#'
#'
#' ## Specify the time points (log-scaled ages)
#' logage<-seq(-0.5, 0.5, by=0.1)
#'
#' ## Obtain a vector of the first derivative of the logistic model
#' logit_first_derv<-logist_first_deriv_ft(6, 0, 1,  logage);
#'
#' ## Plot the first derivative of the arctangent model
#' plot(logage, logit_first_derv, type='l')
#'
#'
#'
#'
logist_first_deriv_ft<- function(theta1,theta2, theta3, logage){

  ###########################################
  ## nonlinear parametrized omega function ##
  ###########################################

  numer<-theta3*theta1*exp(-theta1*(logage-theta2))
  denom <- 1+exp(-theta1*(logage-theta2))
  denom2 <- denom*denom

  logist_first_deriv <- numer/denom2

  return(logist_first_deriv)


}




#' Second derivative of the logistic model function for the parametric nonlinear mixed effects model (nlme) procedure
#'
#'
#'
#' @param theta1 a scale parameter, which determines steepness of the logistic model.
#' @param theta2 a parameter for an inflection point.
#' @param theta3 a parameter, which determines the maximum  (asymptote) of the logistic model.
#' @param logage a vector of time points (log-scaled ages), which is generated around \code{theta2}. It can be an any length of numeric vector.
#'
#'
#' @return   a vector of the second derivative of the logistic function with respect to (logage-theta2).
#' @export
#'
#' @examples
#'
#'
#' library(HDChangePoint)
#'
#' ## Specify parameters
#' theta1=6;
#' theta2=0;
#' theta3=1;
#'
#' ## Specify the time points (log-scaled ages)
#' logage<-seq(-0.5, 0.5, by=0.1)
#'
#' ## Obtain a vector of the second derivative of the logistic model
#' logit_sec_derv<-logist_second_deriv_ft(6, 0, 1, logage);
#'
#' ## Plot the second derivative of the logistic model
#' plot(logage, logit_sec_derv, type='l')
#'
#'
#'
#'
logist_second_deriv_ft<- function(theta1,theta2, theta3, logage){

  ###########################################
  ## nonlinear parametrized omega function ##
  ###########################################
  denom <- 1+exp(-theta1*(logage-theta2))
  numer1<--theta3*(theta1)^2*exp(-theta1*(logage-theta2))
  numer2<-1-exp(-theta1*(logage-theta2))

  denom3 <- denom*denom*denom
  numer<-numer1*numer2

  logist_second_deriv <- numer/denom3


  return(logist_second_deriv)
}


