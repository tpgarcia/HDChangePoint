#' Arctanget model for the parametric nonlinear mixed effects model (nlme) procedure
#'
#' @param gam1 a scale parameter.
#' @param gam2 a parameter for an inflection point.
#' @param gam3 a parameter, which determines steepness of the arctangent model.
#' @param gam4 a vertical shift parameter that can shift the maximum and the minimum of the function values.
#' @param logage a vector of time points (log-scaled ages), which is generated around \code{gam2}. It can be an any length of vector.
#'
#' @return a vector of response variables under the deterministic arctangent model, whose length is the same as the length of \code{logage}.
#'
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ## Specify parameters
#' gam1=2.45/pi
#' gam2=0
#' gam3=pi/1.1
#' gam4=0.8
#'
#' ## Specify the time points (log-scaled ages)
#' logage<-seq(-0.5, 0.5, by=0.1)
#'
#' ## Obtain a vector of response variables when parametric NLME assumes the arctangent model
#' arctanft<-arctanf(2.45/pi, 0, pi/1.1, 0.8, logage)
#'
#' ## Plot the S-shaped curve under the arctangent model
#' plot(logage, arctanft, type='l')
#'
#'
arctanf<- function(gam1,gam2,gam3, gam4, logage){


  ## objective function ##
  arctan<-(gam1)*atan(gam3*(logage-gam2))+gam4

  return(arctan=arctan)

  #  analytical dervivatives
  denom<-{1+(gam3*(logage-gam2))^{2}}

  arctangrad <- array(0,c(length(logage),4),list(NULL,c("gam1","gam2", "gam3", "gam4")))
  arctangrad[,"gam1"] <-atan(gam3*(x-gam2))
  arctangrad[,"gam2"] <-(gam1)*(1/denom)*(-gam3)
  arctangrad[,"gam3"] <-(gam1)*(1/denom)*(logage-gam2)
  arctangrad[,"gam4"] <-1

  attr(arctan,"gradient") <- arctangrad

}


#' First derivative of the arctangent model function for the parametric nonlinear mixed effects model (nlme) procedure
#'
#' @param gam1 a scale parameter.
#' @param gam2 a parameter for an inflection point.
#' @param gam3 a parameter, which determines steepness of the arctangent model.
#' @param logage a vector of time points (log-scaled ages), which is generated around \code{gam2}. It can be an any length of vector.
#'
#' @return a vector of the first derivative of the arctangent function with respect to (logage-gam2).
#'
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ## Specify parameters
#' gam1=2.45/pi
#' gam2=0
#' gam3=pi/1.1
#'
#' ## Specify the time points (log-scaled ages)
#' logage<-seq(-0.5, 0.5, by=0.1)
#'
#' ## Obtain a vector of the first derivative of the arctangent model
#' arctan_first_derv<-arctan_first_deriv_ft(2.45/pi, 0, pi/1.1, 0.8, logage)
#'
#' ## Plot the first derivative of the arctangent model
#' plot(logage, arctan_first_derv, type='l')
#'
arctan_first_deriv_ft<- function(gam1,gam2,gam3,logage){

  ## objective function ##
  #arctan<-(theta1)*atan(theta3*(logage-theta2))+theta4

  denom<-1+((gam3*(logage-gam2))^{2})
  arctan_first_deriv<-(gam1*gam3)/denom


  return(arctan_first_deriv)

}


#' Second derivative of the arctangent model function for the parametric nonlinear mixed effects model (nlme) procedure
#'
#' @param gam1 a scale parameter.
#' @param gam2 a parameter for an inflection point.
#' @param gam3 a parameter, which determines steepness of the arctangent model.
#' @param logage a vector of time points (log-scaled ages), which is generated around \code{gam2}. It can be an any length of numeric vector.
#'
#' @return  a vector of the second derivative of the arctangent function with respect to (logage-gam2).
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ## Specify parameters
#' gam1=2.45/pi
#' gam2=0
#' gam3=pi/1.1
#'
#' ## Specify the time points (log-scaled ages)
#' logage<-seq(-0.5, 0.5, by=0.1)
#'
#' ## Obtain a vector of the second derivative of the arctangent model
#' arct_sec_derv<-arctan_second_deriv_ft(2.45/pi, 0, pi/1.1, logage)
#'
#' ## Plot the second derivative of the arctangent model
#' plot(logage, arct_sec_derv, type='l')
#'
#'
#'
arctan_second_deriv_ft<- function(gam1,gam2,gam3,logage){

  ## objective function ##
  #arctan<-(theta1)*atan(theta3*(logage-theta2))+theta4
  numer<--2*gam1*(gam3^{3})*(logage-gam2)
  denom<-1+(gam3*(logage-gam2))^{2}

  arctan_second_deriv<-numer/denom


  return(arctan_second_deriv)

}
