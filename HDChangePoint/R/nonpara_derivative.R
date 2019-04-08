#' Title
#'
#' @param x
#' @param y
#'
#' @return Given vectors x and y, returns first derivatives of vector y corresponding a vector x
#' @export
#'
#' @examples
#'           x<-seq(-0.5, 0.5, by=0.1); y<-w(x, 0, model="logist");
#'           out<-finite.differences(x, y); plot(out);
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

  # For the last value, since we are unable to perform the forward differencing method
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration
  # in the forward differencing method, but it is used as an approximation as we
  # don't have any more information

  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])

  return(fdx)
}


#' Title
#'
#' @param x
#' @param y
#'
#' @return Given vectors x and y, returns second derivatives of vector y corresponding a vector x
#' @export
#'
#' @examples x<-seq(-0.5, 0.5, by=0.01); y<-w(x, 0, model="logist");
#'           out<-sec.der.mid(x, y); plot(out);
#'
#'

sec.der.mid <- function(x, y) {



  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }

  n <- length(x)

  # Iterate through the values using the forward differencing method
  fddx<-rep(100, n)

  for (i in 2:n-1) {

    fddx[i-1] <- (y[i-1] -2*y[i]+ y[i+1]) / (x[i-1] - x[i])^{2}
  }

  # For the last value, since we are unable to perform the forward differencing method
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration
  # in the forward differencing method, but it is used as an approximation as we
  # don't have any more information
  fddx[n] <- (y[n-2] -2*y[n-1]+y[n]) / (x[n-2] - x[n-1])^{2}
  fddx[n-1] <- (y[n-3] -2*y[n-2]+y[n-1]) / (x[n-3] - x[n-2])^{2}

  return(fddx)
}
