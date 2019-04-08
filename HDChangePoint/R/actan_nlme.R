
######################################################
# s-shaped model for nlme:     arctangent function   #
######################################################

#' Title
#'
#' @param gam1
#' @param gam2
#' @param gam3
#' @param gam4
#' @param logage
#'
#' @return a vector of variables from the arctangent function.
#' @export
#'
#' @examples  logage<-seq(-0.5, 0.5, by=0.1);
#'            arctanft<-arctanf(2.45/pi, 0, pi/1.1, 0.8, logage);
#'            plot(logage, arctanft, type='l')
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



##########################################################
# first derivative of arctangent model with respect to  ##
# (logage-gam2)                                         ##
##########################################################


#' Title
#'
#' @param gam1
#' @param gam2
#' @param gam3
#' @param logage
#'
#' @return a vector of the first derivatives of the arctangent function.
#' @export
#'
#' @examples  logage<-seq(-0.5, 0.5, by=0.1);
#'            arctan_first_derv<-arctan_first_deriv_ft(2.45/pi, 0, pi/1.1, 0.8, logage);
#'            plot(logage, arctan_first_derv, type='l')
#'
arctan_first_deriv_ft<- function(gam1,gam2,gam3,logage){

  ## objective function ##
  #arctan<-(theta1)*atan(theta3*(logage-theta2))+theta4

  denom<-1+((gam3*(logage-gam2))^{2})
  arctan_first_deriv<-(gam1*gam3)/denom


  return(arctan_first_deriv)

}




#####################################################################
# second derivative of the arctangent function for parametric NLME ##
# with respect to (logage-gam2)                                    ##
#####################################################################


#' Title
#'
#' @param gam1
#' @param gam2
#' @param gam3
#' @param logage
#'
#' @return  a vector of the second derivatives of the arctangent function
#' @export
#'
#' @examples  logage<-seq(-0.5, 0.5, by=0.1);
#'            arct_sec_derv<-arctan_second_deriv_ft(2.45/pi, 0, pi/1.1, 0.8, logage);
#'            plot(logage, arct_sec_derv, type='l')
#'
arctan_second_deriv_ft<- function(gam1,gam2,gam3,logage){

  ## objective function ##
  #arctan<-(theta1)*atan(theta3*(logage-theta2))+theta4
  numer<--2*gam1*(gam3^{3})*(logage-gam2)
  denom<-1+(gam3*(logage-gam2))^{2}

  arctan_second_deriv<-numer/denom


  return(arctan_second_deriv)

}
