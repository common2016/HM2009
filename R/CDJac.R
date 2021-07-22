#' CDJac
#'
#' @description computes a central difference approximation of the Jacobian
#'           matrix of a system of n non-linear functions y_i=f^i(x), where
#'           x is a column vector of dimension m.
#'
#' @param f pointer to the routine that returns the m-vector f(x)
#' @param x0 m vector x0, the point at which the derivatives are to be evaluated
#' @param ... any additional arguments passed to \code{f}.
#'
#' @details algorithm: based on (A.2.8) in Heer and Maussner, see also Dennis and Schnabel (1983),
#'           Algorithm A5.6.4.
#' @return Jac n by m matrix of parital derivatives
#' @export

CDJac <- function(f,x0,...){
  m <- length(x0)
  n <- length(f(x0,...))
  df <- matlab::zeros(n,m)
  eps <- .Machine$double.eps^(1/3)
  x1 <- x2 <- x0

  for (i in 1:m){
    if (x0[i] < 0){
      h <- -eps*max(c(abs(x0[i]),1.0))
    } else {
      h <- eps*max(c(abs(x0[i]),1.0))
    }
    temp <- x0[i]
    x1[i] <- temp + h
    x2[i] <- temp - h
    h <- x1[i]-temp # Trick to increase precision slightly, see Dennis and Schnabel (1983), p. 99
    f1 <- f(x1,...)
    f2 <- f(x2,...)
    for (j in 1:n){
      df[j,i]=(f1[j]-f2[j])/(2*h)
    }
    x1[i] <- x0[i]
    x2[i] <- x0[i]
  }
  return(df)
}

