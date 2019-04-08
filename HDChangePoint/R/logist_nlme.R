
###############################
# logistic function for nlme ##
###############################

#' Title
#'
#' @param theta1
#' @param theta2
#' @param theta3
#' @param logage
#'
#' @return a vector of variables from the logistic function.
#' @export
#'
#' @examples  logage<-seq(-0.5, 0.5, by=0.1);
#'            logit.ft<-logistft(6, 0, 1, logage);
#'            plot(logage, logit.ft, type='l')
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



###################################################################
# first derivative of the logistic function for parametric NLME  ##
# with respect to  logage-theta2                                 ##
###################################################################

#' Title
#'
#' @param theta1
#' @param theta2
#' @param theta3
#' @param logage
#'
#' @return a vector of the first derivatives of the logistic function.
#' @export
#'
#' @examples  logage<-seq(-0.5, 0.5, by=0.1);
#'            logit_first_derv<-logist_first_deriv_ft(2.45/pi, 0, pi/1.1, 0.8, logage);
#'            plot(logage, logit_first_derv, type='l')
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




###################################################################
# second derivative of the logistic function for parametric NLME ##
# with respect to logage-theta2                                  ##
###################################################################

#' Title
#'
#' @param theta1
#' @param theta2
#' @param theta3
#' @param logage
#'
#' @return  a vector of the second derivatives of the logistic function
#' @export
#'
#' @examples  logage<-seq(-0.5, 0.5, by=0.1);
#'            logit_sec_derv<-logist_second_deriv_ft(2.45/pi, 0, pi/1.1, 0.8, logage);
#'            plot(logage, logit_sec_derv, type='l')
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


