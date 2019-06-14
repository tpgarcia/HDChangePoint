#' Nonlinear S-shaped Model
#'
#' @param x a vector of time points. It can be any length of a vector.
#' @param logT a parameter for an inflection point.
#' @param model a character string for a nonlinear model: \code{"logist"} or \code{"arctan"}.
#'
#' @return a vector of response values corresponding a vector of time points \eqn{x}, whose trajecotry follows an S-shape.
#'         Two deterministic \code{"S"}\code{-}shaped functions are the logstic and the aractangent functions and user can choose either one.
#'
#' @export
#'
#'
#'
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#'
#' ## Generate a vector of time points and an inflection point in the domain of the time points
#' x<-seq(-0.5, 0.5, by=0.1); logT<-0;
#'
#' ## Obtain an S-shaped function under the logistic model
#' out<-w(x, logT, model="logist");
#' plot(out)
#'
#' ## Obtain an S-shaped function under the arctangent model
#' out1<-w(x, logT, model="arctan");
#' plot(out1)
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


