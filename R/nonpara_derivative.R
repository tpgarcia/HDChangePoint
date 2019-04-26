#' First Derivative Using A Finite Differences Method
#'
#' @param x a m-length of numeric vector.
#' @param y a m-length of numeric vector.
#'
#' @return Given vectors x and y, returns a m-length of the first derivatives of vector y corresponding a vector x. The reference is \url{https://www.r-bloggers.com/numerical-differentiation-with-finite-differences-in-r/}.
#'
#'
#' @export
#'
#' @examples
#'
#'
#' library(HDChangePoint)
#'
#' ##Generate the same length of vectors x and y.
#' x<-seq(-0.5, 0.5, by=0.01);
#' y<-w(x, 0, model="logist");
#'
#' ## Obtain the first derivatives of y using the forward difference method.
#' out<-finite.differences(x, y); plot(out);
#'
#'
#'
finite.differences <- function(x, y) {


  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }

  n <- length(x)

  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = n)

  # Iterate through the values using the forward differencing method
  for (i in 2:n) {
    fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])
  }

  # Use the backward differencing approach for the last value
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])

  return(fdx)
}



#' Second Derivatives Using the Midpoint Formula
#'
#' @param x a m-length of numeric vector.
#' @param y a m-length of numeric vector.
#'
#' @return Given m-length of vectors x and y, returns a m-length of the second derivatives of y
#'         corresponding a vector x.
#'
#' @export
#'
#' @examples
#'
#' library(HDChangePoint)
#'
#' ##Generate the same length of vectors x and y.
#'
#' x<-seq(-0.5, 0.5, by=0.01);
#' y<-w(x, 0, model="logist");
#'
#'
#' ## Obtain the second derivatives of y using the midpoint formula.
#' out<-sec.der.mid(x, y); plot(out);
#'
#'

sec.der.mid <- function(x, y) {



  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }

  n <- length(x)

  # Initialize a vector of length n to enter the derivative approximations
  fddx<-rep(100, n)

  # Iterate through the values using the Midpoint Formula
  for (i in 2:n-1) {

    fddx[i-1] <- (y[i-1] -2*y[i]+ y[i+1]) / (x[i-1] - x[i])^{2}
  }


  fddx[n] <- (y[n-2] -2*y[n-1]+y[n]) / (x[n-2] - x[n-1])^{2}
  fddx[n-1] <- (y[n-3] -2*y[n-2]+y[n-1]) / (x[n-3] - x[n-2])^{2}

  return(fddx)
}
