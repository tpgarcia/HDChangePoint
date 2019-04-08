############################################################
## This is the true omega function: logsitic or arctangent #
## logistic: true fixed effects theta1=1, theta2=6         #
## arctangent: alpha1=2.45/pi, alpha2=pi/1.1, alpha3=0.8   #
## random inflection point logT should be generated        #
############################################################


#' Title
#'
#' @param x
#' @param logT
#' @param model
#'
#' @return a vector of response values corresponding time points x, whose trajecotry follows an S-shape.
#'         In practice, the locations of the inflection points (logT) are latent values.
#' @export
#'
#' @examples  x<-seq(-0.5, 0.5, by=0.1); logT<-0;
#'            out<-w(x, logT, model="logist"); plot(out)
#'
#'
#'
w<-function(x, logT, model) {

  if (model=="logist"){

    w<-1/(1+exp(-6*(x-logT)))

  } else {

    w<-(2.45/pi)*atan((pi/1.1)*(x-logT))+0.8

  }

  return(w)
}


